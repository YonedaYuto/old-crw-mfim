# =============================================================================
# 05_m2_clm.R
#   解析計画書 §8.3「M2：累積ロジット（比例オッズ）モデル」（plan.summary §4-2 の
#   ①②③、E-3）の実装。plan.detail §11 のスクリプト12 に相当。
#
#   入力 : output/04_ridit_spline.rda（dat_m2、spline_terms、M2_SETTINGS）
#          01_labels.R（表示用の文字列）
#   出力 : output/table_m2_odds_ratio.csv    Table 2 に併記する共通オッズ比（比例オッズ、E-3）
#          output/table_m2_ppo_odds_ratio.csv Table 2 の閾値別オッズ比（部分比例オッズ、§7-2）
#          output/table_m2_ppo_profile_check.csv 自前の profile 計算の照合
#          output/table_m2_coefficients.csv  補足表4（部分比例オッズモデルの係数。
#                                            年齢群以外は解釈しない）
#          output/table_m2_po_coefficients.csv 比例オッズモデルの係数（参考）
#          output/table_m2_fit.csv           当てはめの諸元と収束診断（4本）
#          output/table_checks_05.csv        点検結果
#          output/05_m2_clm.rda              fit_G1 / fit_G2 / fit_ppo_G1 / fit_ppo_G2 ほか
#                                            （06 が読む）
#          output/log_05_m2_clm.txt          実行ログ
#
#   改訂 : 2026-09-22 plan.summary.txt §7-1・§7-2 により、年齢群について比例オッズ
#          仮定を外した部分比例オッズモデル（M2-PPO）を加えた。
#          ・アウトカムは 70 点未満・70〜76 点・77 点以上の3区分とし、年齢群を
#            nominal 項に置く。79 段階のまま nominal にすると、年齢群内で観測の
#            ない水準（G1 の 69・74 点、G3 の 71 点など）をはさむ閾値が発散し、
#            最尤推定値が存在しない（著者の判断、2026-09-22）。
#          ・70 点・77 点での年齢群のオッズ比の 95% 区間は profile likelihood
#            とする（§4-2）。ordinal の profile は nominal 項を扱わないため、
#            対数尤度を自前で書いて計算する（著者の判断、2026-09-22）。06 で同じ
#            再標本からブートストラップ区間も求め、両者を突き合わせる。
#          ・比較として併記する比例オッズモデル（79 段階）の共通オッズ比は
#            従来どおり残す。M2_FORMULA_TXT はこのモデルの式のまま変えない
#            （10 の照合を保つため）。部分比例オッズモデルの式は
#            M2_PPO_FORMULA_TXT / M2_PPO_NOMINAL_TXT に置く。
#          ・m2_coefficients（補足表4）は部分比例オッズモデルの係数に替えた。
#            列の構成は変えていない（11 を変えずに済ませるため）。
#
#   計画との対応
#     ・モデル式 : logit P(Y ≥ y | X) = α_y + β_G2・G2 + β_G3・G3 + γ'Z（§8.3）
#       M2 = 年齢群 + 性別 + 疾患区分 + 病前の要介護状態
#            + 入院時運動FIM（ridit スプライン） + 入院時認知FIM（ridit スプライン）
#            + 入院時期区分
#     ・区間 : profile likelihood（未調整95%）。G3対G2 は G2 を参照とした
#       再当てはめから得る（§8.3）。線形対比の Wald 区間は用いない
#     ・大域尤度比検定（H0: β_G2 = β_G3 = 0、自由度2）は計算してログにのみ残す（§8.5）
#     ・対比較の尤度比検定は行わない（§8.3）
#     ・p値は出力表に入れない（§8.5-1：Table 2 にはp値を載せない）
#     ・部分比例オッズモデル（§7-2）:
#       logit P(Y ≥ t | X) = α_t + β_G2,t・G2 + β_G3,t・G3 + γ'Z、t = 70, 77。
#       年齢群の係数だけを閾値ごとに自由にし、共変量の係数は2つの閾値で共通。
#       Y は 3 区分（<70、70–76、≥77）。区間は profile likelihood（未調整95%）で、
#       G3対G2 は G2 を参照とした再当てはめから得る（比例オッズモデルと同じ手順）
#
#   `ordinal::clm` の符号について
#     clm は logit P(Y ≤ j) = θ_j − x'β と書く。したがって coef の符号は
#     「Y が大きくなる向き」に揃っており、exp(β) はそのまま
#     「Y ≥ y となるオッズ」の比（§8.3 の β）になる。この対応は当てはめ後に
#     予測確率を使って数値で点検する（下記 §7）。
#     nominal 項は閾値の側に入り、logit P(Y ≤ j) = θ_j + τ_j,g − x'β となる。
#     したがって閾値 j での「Y ≥ y となるオッズ」の比は exp(−τ_j,g) であり、
#     本スクリプトは符号を反転して報告する。この対応も §7 で数値で点検する。
#
#   文字コード : UTF-8（Windows の R では source(..., encoding = "UTF-8")）
# =============================================================================


# -----------------------------------------------------------------------------
# 0. 設定
# -----------------------------------------------------------------------------

OUT_DIR <- "output"
IN_RDA  <- file.path(OUT_DIR, "04_ridit_spline.rda")

CONF_LEVEL <- 0.95              # 未調整95%（§8.5）

## --- 対比（§8.3、§8.4）。高齢側を被減数・分子に置く -------------------------
CONTRASTS <- list(
  list(key = "G2_vs_G1", ref = "G1", cmp = "G2", role = "secondary"),
  list(key = "G3_vs_G1", ref = "G1", cmp = "G3", role = "primary"),
  list(key = "G3_vs_G2", ref = "G2", cmp = "G3", role = "secondary")
)
PRIMARY_KEY <- "G3_vs_G1"       # 主対比（§1.2）

## --- 部分比例オッズモデル（plan.summary §7-1・§7-2）--------------------------
## 70 点：入浴と階段昇降以外が見守りレベル（旧設定 65 点を修正）。
## 77 点：入浴と階段昇降以外が修正自立レベル。昇順に置く。
PPO_THRESHOLDS <- c(70, 77)

#' 退院時運動FIM を閾値で区切った順序因子（3区分）を作る。06 も同じ規則を使う。
#' 水準名は閾値だけから決まり、値域（13〜91、感度分析Sでは 12〜91）に依存しない。
make_ppo_levels <- function(thresholds) {
  th <- sort(thresholds)
  mid <- if (length(th) > 1L) sprintf("%g-%g", head(th, -1L), tail(th, -1L) - 1) else character(0)
  c(sprintf("<%g", th[1]), mid, sprintf(">=%g", th[length(th)]))
}
make_ppo_outcome <- function(y, thresholds) {
  lev <- make_ppo_levels(thresholds)
  factor(lev[findInterval(y, sort(thresholds)) + 1L], levels = lev, ordered = TRUE)
}
PPO_LEVELS <- make_ppo_levels(PPO_THRESHOLDS)
PPO_CUTS   <- paste(head(PPO_LEVELS, -1L), tail(PPO_LEVELS, -1L), sep = "|")
names(PPO_CUTS) <- as.character(PPO_THRESHOLDS)   # 閾値 → clm の閾値名

## --- 自前の profile likelihood（§4-2、著者の判断 2026-09-22）------------------
PROFILE_UNIROOT_TOL <- 1e-9     # 根の位置の許容誤差（対数オッズ比の尺度）
PROFILE_MAX_SE      <- 50       # 根を探す範囲の上限（Wald 標準誤差の倍数）
PROFILE_AGREE_TOL   <- 1e-3     # ordinal の profile 区間との照合の許容差（位置の係数）
LOGLIK_AGREE_TOL    <- 1e-6     # clm の対数尤度との照合の許容差

## --- 収束と発散の判定（§8.3「本体の当てはめが収束しない場合の事前規則」）----
## 許容する対処は反復回数上限の引き上げと初期値の変更だけである。
CLM_MAX_ITER     <- 200L        # clm.control(maxIter)
CLM_MAX_LINE     <- 50L         # clm.control(maxLineIter)
CLM_MAX_MOD_ITER <- 10L         # clm.control(maxModIter)
## |β| の判定はスプライン基底以外の係数に当てる。ns() の基底は値域が狭いため、
## 対応する係数は分離がなくても大きくなる。分離の兆候は、カテゴリ共変量の
## 係数の大きさと、すべての係数の標準誤差の大きさで見る。
DIVERGE_COEF <- 10              # スプライン以外の |β| がこれを超えたら REVIEW
DIVERGE_SE   <- 10              # se がこれを超えたら同上
IS_SPLINE_TERM <- function(nm) grepl("^ns\\(", nm)


# -----------------------------------------------------------------------------
# 1. 表示用ラベル
# -----------------------------------------------------------------------------

if (!file.exists("01_labels.R")) {
  stop("01_labels.R が見つからない。05_m2_clm.R は 01_labels.R と同じ",
       "フォルダを作業ディレクトリにして実行すること。")
}
source("01_labels.R", encoding = "UTF-8")

contrast_label <- function(ref, cmp) {
  sprintf("%s vs %s (%s vs %s)", cmp, ref,
          relabel_levels("age_group", cmp), relabel_levels("age_group", ref))
}


# -----------------------------------------------------------------------------
# 2. 実行の準備（出力先とログ）
# -----------------------------------------------------------------------------

if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

.log_con <- file(file.path(OUT_DIR, "log_05_m2_clm.txt"), open = "wt",
                 encoding = "UTF-8")
sink(.log_con, split = TRUE)

.old_error <- getOption("error")
.log_close <- function() {
  while (sink.number() > 0) sink()
  try(close(.log_con), silent = TRUE)
  options(error = .old_error)
}
options(error = function() .log_close())

say  <- function(...) cat(..., "\n", sep = "")
rule <- function(title) say("\n", strrep("-", 74), "\n", title, "\n", strrep("-", 74))

.checks <- list()
check <- function(name, ok, detail = "") {
  .checks[[length(.checks) + 1L]] <<-
    data.frame(check = name, result = if (isTRUE(ok)) "OK" else "REVIEW",
               detail = detail, stringsAsFactors = FALSE)
  say(sprintf("  [%-6s] %s%s", if (isTRUE(ok)) "OK" else "REVIEW", name,
              if (nzchar(detail)) paste0(" : ", detail) else ""))
  invisible(ok)
}

#' 関数に実在する引数だけを渡して呼ぶ（パッケージの版で引数名が変わっても落ちない）
call_known_args <- function(fun, args) {
  fml  <- setdiff(names(formals(fun)), "...")
  drop <- setdiff(names(args), fml)
  if (length(drop)) {
    say("    （この版にない引数は渡さない: ", paste(drop, collapse = ", "), "）")
  }
  do.call(fun, args[names(args) %in% fml])
}

rule(paste0("05_m2_clm.R  ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
say("R : ", R.version.string)

has_ordinal <- requireNamespace("ordinal", quietly = TRUE)
check("ordinal が利用できる", has_ordinal,
      if (has_ordinal) paste0("version ", as.character(utils::packageVersion("ordinal")))
      else "install.packages(\"ordinal\") が必要")
if (!has_ordinal) stop("ordinal が無ければ M2 は当てはめられない。")
library(ordinal)
library(splines)


# -----------------------------------------------------------------------------
# 3. 入力の取り込み
# -----------------------------------------------------------------------------

rule("3. Input")

if (!exists("dat_m2", inherits = TRUE) || !exists("spline_terms", inherits = TRUE)) {
  if (!file.exists(IN_RDA)) {
    stop(IN_RDA, " が見つからない。先に 04_ridit_spline.R を実行すること。")
  }
  load(IN_RDA)
  say("loaded : ", IN_RDA)
} else {
  say("dat_m2 / spline_terms はワークスペース上のものを使う")
}

OUTCOME    <- M2_SETTINGS$outcome
AGE_LEVELS <- M2_SETTINGS$age_levels
COVARS_CAT <- M2_SETTINGS$covars_cat

dat <- as.data.frame(dat_m2, stringsAsFactors = FALSE)
say("dat_m2 (M2 complete-case population) : ", nrow(dat), " rows")

n_by_g <- table(dat$age_group)
for (g in AGE_LEVELS) {
  say(sprintf("  %-3s %-14s n = %5d", g, relabel_levels("age_group", g),
              as.integer(n_by_g[[g]])))
}

check("順序因子のアウトカムがある", is.ordered(dat$y_ord),
      sprintf("levels = %d", nlevels(dat$y_ord)))
check("年齢群の参照が G1 である", identical(levels(dat$age_group)[1], "G1"),
      paste(levels(dat$age_group), collapse = ", "))

## --- 部分比例オッズモデルのアウトカム（3区分、§7-1・§7-2）-------------------
dat$y_ppo <- make_ppo_outcome(dat[[OUTCOME]], PPO_THRESHOLDS)
ppo_cells <- table(factor(dat$age_group, levels = AGE_LEVELS), dat$y_ppo)
say(sprintf("\n  PPO outcome : %s (thresholds %s)", paste(PPO_LEVELS, collapse = " / "),
            paste(PPO_THRESHOLDS, collapse = ", ")))
print(ppo_cells)
check("3区分のアウトカムに欠測がない", !anyNA(dat$y_ppo), "")
## nominal 項の推定値が存在するには、年齢群 × 区分 のどのセルにも観測が要る
check("年齢群 × 3区分 のすべてのセルに観測がある（nominal 項の推定に必要）",
      all(ppo_cells > 0L),
      sprintf("最小セル = %d", min(ppo_cells)))


# -----------------------------------------------------------------------------
# 4. モデル式（§8.3）
#    スプラインの項はノットの数値が埋め込まれた文字列である（04 で作成）。
#    ブートストラップの各回でも同じ文字列を使い、ノットを再算出しない（§6.4-5）。
# -----------------------------------------------------------------------------

rule("4. Model formula")

M2_RHS <- paste(c("age_group", "sex", "class", "support_in",
                  unname(spline_terms["mFIM_in_r"]),
                  unname(spline_terms["cFIM_in_r"]),
                  "period"), collapse = " + ")
M2_FORMULA_TXT <- paste("y_ord ~", M2_RHS)
M2_FORMULA     <- stats::as.formula(M2_FORMULA_TXT)

## --- 部分比例オッズモデル（§7-2）-------------------------------------------
## 位置の項は M2 から年齢群を除いたもの。年齢群は nominal 項に置く。
## M2_FORMULA_TXT は比例オッズモデルの式のまま残す（10 が照合に使う）。
M2_PPO_RHS         <- sub("^age_group \\+ ", "", M2_RHS)
M2_PPO_FORMULA_TXT <- paste("y_ppo ~", M2_PPO_RHS)
M2_PPO_NOMINAL_TXT <- "~ age_group"
M2_PPO_FORMULA     <- stats::as.formula(M2_PPO_FORMULA_TXT)
M2_PPO_NOMINAL     <- stats::as.formula(M2_PPO_NOMINAL_TXT)

say("  proportional odds (M2)            : ", M2_FORMULA_TXT)
say("  partial proportional odds (M2-PPO): ", M2_PPO_FORMULA_TXT,
    " ; nominal = ", M2_PPO_NOMINAL_TXT)
check("M2-PPO の位置の項が M2 から年齢群だけを除いたものである",
      !identical(M2_PPO_RHS, M2_RHS) &&
        identical(paste("age_group +", M2_PPO_RHS), M2_RHS), "")
say("\n  変数の対応:")
for (v in c("age_group", "sex", "class", "support_in",
            "mFIM_in_r", "cFIM_in_r", "period")) {
  say(sprintf("    %-12s %s", v, relabel_vars(v)))
}

## 年齢群以外の係数は解釈しない（§8.3、Table 2 fallacy［8］）
say("\n  解釈するのは年齢群の係数だけである（§8.3）。",
    "他の係数は補足表4に載せ、解釈しない。")


# -----------------------------------------------------------------------------
# 5. 当てはめ
# -----------------------------------------------------------------------------

rule("5. Fitting")

ctrl <- call_known_args(ordinal::clm.control,
                        list(maxIter = CLM_MAX_ITER, maxLineIter = CLM_MAX_LINE,
                             maxModIter = CLM_MAX_MOD_ITER))

#' 年齢群の参照水準を入れ替えたデータで M2 を当てはめる
fit_m2 <- function(data, ref) {
  d <- data
  d$age_group <- stats::relevel(factor(as.character(d$age_group),
                                       levels = AGE_LEVELS), ref = ref)
  t0 <- proc.time()[["elapsed"]]
  fit <- ordinal::clm(M2_FORMULA, data = d, link = "logit", control = ctrl)
  attr(fit, "seconds") <- proc.time()[["elapsed"]] - t0
  attr(fit, "ref")     <- ref
  fit
}

say("  fitting with age_group reference = G1 ...")
fit_G1 <- fit_m2(dat, "G1")
say(sprintf("    done in %.1f s", attr(fit_G1, "seconds")))

say("  fitting with age_group reference = G2 （G3対G2 の profile 区間のため。§8.3）...")
fit_G2 <- fit_m2(dat, "G2")
say(sprintf("    done in %.1f s", attr(fit_G2, "seconds")))

#' 部分比例オッズモデル（§7-2）。年齢群を nominal 項に置く。参照の入れ替えは M2 と同じ。
fit_m2_ppo <- function(data, ref) {
  d <- data
  d$age_group <- stats::relevel(factor(as.character(d$age_group),
                                       levels = AGE_LEVELS), ref = ref)
  t0 <- proc.time()[["elapsed"]]
  fit <- ordinal::clm(M2_PPO_FORMULA, nominal = M2_PPO_NOMINAL, data = d,
                      link = "logit", control = ctrl)
  attr(fit, "seconds") <- proc.time()[["elapsed"]] - t0
  attr(fit, "ref")     <- ref
  fit
}

say("  fitting M2-PPO with age_group reference = G1 ...")
fit_ppo_G1 <- fit_m2_ppo(dat, "G1")
say(sprintf("    done in %.1f s", attr(fit_ppo_G1, "seconds")))
say("  fitting M2-PPO with age_group reference = G2 （G3対G2 の profile 区間のため）...")
fit_ppo_G2 <- fit_m2_ppo(dat, "G2")
say(sprintf("    done in %.1f s", attr(fit_ppo_G2, "seconds")))

## --- 収束診断（§8.3）--------------------------------------------------------
n_of <- function(fit) {
  n <- tryCatch(as.integer(stats::nobs(fit))[1], error = function(e) NA_integer_)
  if (is.na(n)) {
    nn <- fit[["n"]]
    if (!is.null(nn) && length(nn) == 1L) n <- suppressWarnings(as.integer(nn))
  }
  n
}

## clm の $convergence は長さ3のリスト（コード・要約・詳細）である。
conv_code <- function(fit) {
  cv <- fit[["convergence"]]
  if (is.null(cv)) return(NA_integer_)
  if (is.list(cv)) cv <- cv[[1]]
  suppressWarnings(as.integer(cv[1]))
}
conv_msg <- function(fit) {
  cv <- fit[["convergence"]]
  if (is.null(cv)) return("")
  paste(trimws(unlist(cv)), collapse = " / ")
}

#' nominal 項に置いた年齢群の母数の名前（比例オッズモデルでは空）
age_nominal_names <- function(fit) grep("\\.age_group", names(fit[["alpha"]]), value = TRUE)

#' 発散の判定に使う係数：位置の係数（β）＋ nominal 項の年齢群の母数（M2-PPO のみ）。
#' 比例オッズモデルでは β だけになり、従来の判定と同じである。
fit_info <- function(fit, tag) {
  an <- age_nominal_names(fit)
  b  <- c(fit[["beta"]], fit[["alpha"]][an])
  se <- tryCatch(sqrt(diag(vcov(fit)))[names(b)], error = function(e) rep(NA_real_, length(b)))
  mg <- fit[["maxGradient"]]
  yl <- fit[["y.levels"]]
  data.frame(
    model        = tag,
    reference    = attr(fit, "ref"),
    n            = n_of(fit),
    n_thresholds = if (length(yl)) length(yl) - 1L else length(fit[["alpha"]]),
    n_beta       = length(b),
    logLik       = as.numeric(stats::logLik(fit)),
    AIC          = as.numeric(stats::AIC(fit)),
    convergence  = conv_code(fit),
    convergence_message = conv_msg(fit),
    max_abs_gradient = if (is.null(mg) || length(mg) != 1L) NA_real_ else as.numeric(mg),
    max_abs_beta = max(abs(b)),
    max_abs_beta_nonspline = if (any(!IS_SPLINE_TERM(names(b))))
                               max(abs(b[!IS_SPLINE_TERM(names(b))])) else NA_real_,
    max_se       = if (all(is.na(se))) NA_real_ else max(se, na.rm = TRUE),
    seconds      = attr(fit, "seconds"),
    stringsAsFactors = FALSE
  )
}
m2_fit_info <- rbind(fit_info(fit_G1, "M2 (ref G1)"), fit_info(fit_G2, "M2 (ref G2)"),
                     fit_info(fit_ppo_G1, "M2-PPO (ref G1)"),
                     fit_info(fit_ppo_G2, "M2-PPO (ref G2)"))
say("")
print(m2_fit_info, row.names = FALSE, digits = 6)

for (i in seq_len(nrow(m2_fit_info))) {
  r <- m2_fit_info[i, ]
  check(sprintf("[%s] 収束した", r$model),
        is.na(r$convergence) || r$convergence == 0L,
        sprintf("convergence = %s (%s), max|gradient| = %s",
                as.character(r$convergence), r$convergence_message,
                formatC(r$max_abs_gradient, format = "g", digits = 3)))
  check(sprintf("[%s] 係数・標準誤差に発散の兆候がない", r$model),
        all(is.finite(c(r$max_abs_beta, r$max_abs_beta_nonspline))) &&
          r$max_abs_beta_nonspline < DIVERGE_COEF &&
          (is.na(r$max_se) || r$max_se < DIVERGE_SE),
        sprintf("max|beta| = %.3f（うちスプライン以外 %.3f）, max se = %s（閾値 %g / %g）",
                r$max_abs_beta, r$max_abs_beta_nonspline,
                formatC(r$max_se, format = "g", digits = 3),
                DIVERGE_COEF, DIVERGE_SE))
}

## ordinal::convergence() の診断（あれば）をログに残す
.conv <- tryCatch(ordinal::convergence(fit_G1), error = function(e) NULL)
if (!is.null(.conv)) {
  say("\n  --- ordinal::convergence(fit_G1) ---")
  print(.conv)
}
.conv_ppo <- tryCatch(ordinal::convergence(fit_ppo_G1), error = function(e) NULL)
if (!is.null(.conv_ppo)) {
  say("\n  --- ordinal::convergence(fit_ppo_G1) ---")
  print(.conv_ppo)
}

check("2つの当てはめの対数尤度が一致する（参照水準の入れ替えのみであるため）",
      abs(m2_fit_info$logLik[1] - m2_fit_info$logLik[2]) < 1e-6,
      sprintf("%.8f vs %.8f", m2_fit_info$logLik[1], m2_fit_info$logLik[2]))
check("M2-PPO の2つの当てはめの対数尤度が一致する（参照水準の入れ替えのみであるため）",
      abs(m2_fit_info$logLik[3] - m2_fit_info$logLik[4]) < 1e-6,
      sprintf("%.8f vs %.8f", m2_fit_info$logLik[3], m2_fit_info$logLik[4]))
check("当てはめに使われた例数が完全ケース集団と一致する",
      all(m2_fit_info$n == nrow(dat)),
      sprintf("%s vs %d", paste(m2_fit_info$n, collapse = ", "), nrow(dat)))
check("M2-PPO の閾値が2つ（70 点・77 点）で、年齢群の nominal 母数が 2 × 2 個ある",
      m2_fit_info$n_thresholds[3] == length(PPO_THRESHOLDS) &&
        length(age_nominal_names(fit_ppo_G1)) ==
          length(PPO_THRESHOLDS) * (length(AGE_LEVELS) - 1L),
      sprintf("thresholds = %d, age nominal = %s", m2_fit_info$n_thresholds[3],
              paste(age_nominal_names(fit_ppo_G1), collapse = ", ")))

## nominal 項は閾値の順序を保証しない。年齢群を置き換えた予測で、
## 累積確率 P(Y ≤ <70) ≤ P(Y ≤ 70-76) がすべての患者で成り立つことを確かめる。
.nd_ppo <- dat
.nd_ppo$y_ord <- NULL; .nd_ppo$y_ppo <- NULL
.mono_ppo <- vapply(AGE_LEVELS, function(g) {
  d <- .nd_ppo
  d$age_group <- stats::relevel(factor(rep(g, nrow(d)), levels = AGE_LEVELS), ref = "G1")
  p  <- stats::predict(fit_ppo_G1, newdata = d, type = "cum.prob")
  cp <- if (is.list(p) && !is.null(p$cprob1)) p$cprob1 else p
  all(apply(cp, 1L, function(r) all(diff(r) >= -1e-12)))
}, logical(1))
check("M2-PPO の累積確率が全患者・全年齢群で単調である（閾値の交差がない）",
      all(.mono_ppo),
      paste(sprintf("%s=%s", AGE_LEVELS, ifelse(.mono_ppo, "ok", "NG")), collapse = ", "))


# -----------------------------------------------------------------------------
# 6. profile likelihood 95%区間（§8.3）
# -----------------------------------------------------------------------------

rule("6. Profile likelihood 95% confidence intervals")

#' profile likelihood 区間を取り出す。版によって confint.clm の引数が異なるため、
#' parm 指定 → 全 beta の profile → 失敗、の順に降りる。Wald への差し替えはしない。
get_profile_ci <- function(fit, parm, level = CONF_LEVEL) {
  as_mat <- function(x, nm) {
    if (is.null(x)) return(NULL)
    if (is.null(dim(x))) x <- matrix(x, nrow = 1L,
                                     dimnames = list(nm, names(x)))
    as.matrix(x)
  }
  ci <- tryCatch(as_mat(confint(fit, parm = parm, level = level, type = "profile"),
                        parm),
                 error = function(e) {
                   say("    confint(parm=...) でエラー: ", conditionMessage(e))
                   NULL
                 })
  src <- "confint(type = \"profile\", parm)"
  if (is.null(ci) || !all(parm %in% rownames(ci))) {
    ci2 <- tryCatch(as_mat(confint(fit, level = level, type = "profile"), parm),
                    error = function(e) {
                      say("    confint(type=\"profile\") でエラー: ",
                          conditionMessage(e)); NULL
                    })
    if (!is.null(ci2) && all(parm %in% rownames(ci2))) {
      ci  <- ci2[parm, , drop = FALSE]
      src <- "confint(type = \"profile\", all beta)"
    } else {
      return(list(ci = NULL, source = "failed"))
    }
  }
  list(ci = ci[parm, , drop = FALSE], source = src)
}

#' Wald 区間（ログでの対照のみ。本文・図表には用いない。§8.3）
get_wald <- function(fit, parm, ref, level = CONF_LEVEL) {
  b  <- fit$beta[parm]
  se <- tryCatch(sqrt(diag(vcov(fit)))[parm], error = function(e) rep(NA_real_, length(parm)))
  z  <- stats::qnorm(1 - (1 - level) / 2)
  data.frame(group_reference = ref, term = parm,
             beta = as.numeric(b), se = as.numeric(se),
             lower = as.numeric(b - z * se), upper = as.numeric(b + z * se),
             stringsAsFactors = FALSE)
}

## 参照 G1 の当てはめから G2対G1・G3対G1、参照 G2 の当てはめから G3対G2 を取る
parm_G1 <- c("age_groupG2", "age_groupG3")
parm_G2 <- c("age_groupG3")

say("  profile likelihood（参照 G1）...")
pl_G1 <- get_profile_ci(fit_G1, parm_G1)
say("    source : ", pl_G1$source)
say("  profile likelihood（参照 G2）...")
pl_G2 <- get_profile_ci(fit_G2, parm_G2)
say("    source : ", pl_G2$source)

check("参照 G1 の profile 区間が得られた", !is.null(pl_G1$ci), pl_G1$source)
check("参照 G2 の profile 区間が得られた", !is.null(pl_G2$ci), pl_G2$source)

wald_all <- rbind(get_wald(fit_G1, parm_G1, "G1"),
                  get_wald(fit_G2, parm_G2, "G2"))

## --- 対比ごとに1行を組み立てる ----------------------------------------------
or_rows <- list()
for (ct in CONTRASTS) {
  key <- ct$key
  if (identical(ct$ref, "G1")) {
    fit <- fit_G1; pl <- pl_G1; term <- paste0("age_group", ct$cmp)
  } else {
    fit <- fit_G2; pl <- pl_G2; term <- paste0("age_group", ct$cmp)
  }
  b  <- unname(fit$beta[term])
  se <- tryCatch(unname(sqrt(diag(vcov(fit)))[term]), error = function(e) NA_real_)
  ci <- if (!is.null(pl$ci) && term %in% rownames(pl$ci)) pl$ci[term, ] else c(NA, NA)

  or_rows[[length(or_rows) + 1L]] <- data.frame(
    contrast         = key,
    label            = contrast_label(ct$ref, ct$cmp),
    group_reference  = ct$ref,
    group_comparison = ct$cmp,
    reference_label  = relabel_levels("age_group", ct$ref),
    comparison_label = relabel_levels("age_group", ct$cmp),
    n_reference      = as.integer(n_by_g[[ct$ref]]),
    n_comparison     = as.integer(n_by_g[[ct$cmp]]),
    role             = ct$role,
    term             = term,
    term_label       = relabel_vars(term),
    beta             = b,
    se               = se,
    odds_ratio       = exp(b),
    ci_lower         = exp(as.numeric(ci[1])),
    ci_upper         = exp(as.numeric(ci[2])),
    conf_level       = CONF_LEVEL,
    ci_method        = "profile likelihood, unadjusted",
    fitted_with      = sprintf("age_group reference = %s", ct$ref),
    p_value          = NA_real_,   # §8.5-1：E-3 に p 値は付けない
    stringsAsFactors = FALSE
  )
}
m2_odds_ratio <- do.call(rbind, or_rows)

say("\n  E-3 : common odds ratio with profile likelihood 95% CI (Table 2)")
print(m2_odds_ratio[, c("contrast", "label", "n_reference", "n_comparison", "role",
                        "odds_ratio", "ci_lower", "ci_upper")],
      row.names = FALSE, digits = 4)

## --- 点検 -------------------------------------------------------------------
check("対比が3つある", nrow(m2_odds_ratio) == 3L,
      sprintf("rows = %d", nrow(m2_odds_ratio)))
check("すべての対比で profile 区間が得られている",
      all(is.finite(m2_odds_ratio$ci_lower) & is.finite(m2_odds_ratio$ci_upper)),
      paste(m2_odds_ratio$contrast[!is.finite(m2_odds_ratio$ci_lower) |
                                   !is.finite(m2_odds_ratio$ci_upper)], collapse = ", "))
check("区間が点推定値を挟む",
      all(is.na(m2_odds_ratio$ci_lower) |
          (m2_odds_ratio$ci_lower <= m2_odds_ratio$odds_ratio &
           m2_odds_ratio$odds_ratio <= m2_odds_ratio$ci_upper)), "")
check("Table 2 に載せる E-3 に p 値を置いていない（§8.5-1）",
      all(is.na(m2_odds_ratio$p_value)), "")

## 推移性の点検：β(G3|ref G1) − β(G2|ref G1) = β(G3|ref G2)（点推定値は一致する）
.b31 <- unname(fit_G1$beta["age_groupG3"])
.b21 <- unname(fit_G1$beta["age_groupG2"])
.b32 <- unname(fit_G2$beta["age_groupG3"])
check("β(G3 vs G1) − β(G2 vs G1) = β(G3 vs G2)（参照の入れ替えの整合）",
      abs((.b31 - .b21) - .b32) < 1e-6,
      sprintf("%.8f − %.8f = %.8f vs %.8f", .b31, .b21, .b31 - .b21, .b32))

## Wald との対照（ログのみ。本文・図表には用いない）
say("\n  --- profile と Wald の対照（ログのみ、§8.3）---")
cmpw <- merge(m2_odds_ratio[, c("contrast", "group_reference", "term",
                                "odds_ratio", "ci_lower", "ci_upper")],
              wald_all[, c("group_reference", "term", "lower", "upper")],
              by = c("group_reference", "term"), all.x = TRUE)
cmpw$wald_or_lower <- exp(cmpw$lower)
cmpw$wald_or_upper <- exp(cmpw$upper)
cmpw$lower <- NULL; cmpw$upper <- NULL
cmpw <- cmpw[match(m2_odds_ratio$contrast, cmpw$contrast), , drop = FALSE]
print(cmpw, row.names = FALSE, digits = 4)
say("  本文・Table 2 に載せるのは profile likelihood 区間である（§8.3）。")


# -----------------------------------------------------------------------------
# 6b. M2-PPO：70 点・77 点での年齢群のオッズ比と profile likelihood 95%区間
#     （plan.summary §7-2、§4-2。著者の判断 2026-09-22）
#     ordinal の profile.clm が扱うのは位置の係数（β）とスケールの係数だけで、
#     nominal 項（閾値の側の母数）は対象外である。そこで M2-PPO の対数尤度を
#     自前で書き、対象の母数を固定して残りを最大化する profile を計算する。
#     実装の正しさは次の3点で確かめる。
#       (i)  clm の推定値での対数尤度が logLik(fit) と一致する
#       (ii) clm の推定値から全母数を最大化し直しても対数尤度が増えない
#       (iii) 位置の係数について、自前の profile 区間が ordinal の
#             confint(type = "profile") と一致する
# -----------------------------------------------------------------------------

rule("6b. M2-PPO: threshold-specific odds ratios with profile likelihood 95% CI")

#' M2-PPO の対数尤度の計算に要る行列をまとめる。
#' 母数の並びは clm と同じ（alpha：nominal の列ごとに閾値を並べたもの、続いて beta）。
ppo_parts <- function(fit, data, resp = "y_ppo") {
  X <- stats::model.matrix(stats::delete.response(stats::terms(fit)), data = data)
  bnm <- names(fit[["beta"]])
  if (!all(bnm %in% colnames(X))) stop("位置の計画行列の列が係数名と対応しない")
  X <- X[, bnm, drop = FALSE]
  N <- stats::model.matrix(fit[["nom.terms"]], data = data)
  lev  <- levels(data[[resp]])
  cuts <- paste(head(lev, -1L), tail(lev, -1L), sep = "|")
  anm  <- as.vector(outer(cuts, colnames(N), paste, sep = "."))
  if (!identical(anm, names(fit[["alpha"]]))) {
    stop("閾値側の母数の並びが想定と異なる: ", paste(names(fit[["alpha"]]), collapse = ", "))
  }
  list(X = X, N = N, y = as.integer(data[[resp]]), K = length(lev),
       C = length(cuts), q = ncol(N), p = ncol(X), cuts = cuts,
       names = c(anm, bnm))
}

#' 対数尤度（grad = TRUE なら勾配）。logit P(Y ≤ c) = θ_c + τ_c'n − x'β。
#' 観測された区分の確率が正でない母数では -Inf（勾配は NA）を返す。
ppo_loglik <- function(par, P, grad = FALSE) {
  C <- P$C; q <- P$q
  A <- matrix(par[seq_len(C * q)], nrow = C, ncol = q)
  b <- par[C * q + seq_len(P$p)]
  H <- P$N %*% t(A) - as.vector(P$X %*% b)          # n × C の logit P(Y ≤ c)
  n <- length(P$y); k <- P$y
  up <- rep(1, n); lo <- rep(0, n)
  iu <- which(k < P$K); il <- which(k > 1L)
  up[iu] <- stats::plogis(H[cbind(iu, k[iu])])
  lo[il] <- stats::plogis(H[cbind(il, k[il] - 1L)])
  pr <- up - lo
  if (any(!is.finite(pr)) || any(pr <= 0)) {
    return(if (grad) rep(NA_real_, length(par)) else -Inf)
  }
  if (!grad) return(sum(log(pr)))
  G <- matrix(0, n, C)                               # dℓ/dH
  G[cbind(iu, k[iu])] <- stats::dlogis(H[cbind(iu, k[iu])]) / pr[iu]
  G[cbind(il, k[il] - 1L)] <- G[cbind(il, k[il] - 1L)] -
    stats::dlogis(H[cbind(il, k[il] - 1L)]) / pr[il]
  c(as.vector(crossprod(G, P$N)), -as.vector(crossprod(P$X, rowSums(G))))
}

#' 母数 idx を val に固定し、残りを最大化する（BFGS、解析的な勾配）。
#' 初期値はつねに clm の推定値とする（評価の順序によらず同じ結果になる）。
ppo_max_fixed <- function(P, start, idx = integer(0), val = numeric(0)) {
  full <- function(r) {
    z <- start
    if (length(idx)) { z[idx] <- val; z[-idx] <- r } else z <- r
    z
  }
  r0 <- if (length(idx)) start[-idx] else start
  fn <- function(r) { v <- ppo_loglik(full(r), P); if (is.finite(v)) -v else Inf }
  gr <- function(r) {
    g <- ppo_loglik(full(r), P, grad = TRUE)
    if (length(idx)) -g[-idx] else -g
  }
  o <- stats::optim(r0, fn, gr, method = "BFGS",
                    control = list(reltol = 1e-14, maxit = 5000))
  list(loglik = -o$value, par = full(o$par), convergence = o$convergence)
}

#' profile likelihood 区間（母数 idx の尺度）。
#' 2{ℓ(θ̂) − ℓ_p(ψ)} = χ²_1(level) となる ψ を推定値の両側で求める。
ppo_profile_ci <- function(P, mle, idx, se, level = CONF_LEVEL) {
  ## se は根を探す最初の刻みにだけ使う。得られない場合は 0.5 から探し始める。
  if (length(se) != 1L || !is.finite(se) || se <= 0) se <- 0.5
  l_hat  <- ppo_loglik(mle, P)
  crit   <- stats::qchisq(level, df = 1)
  n_eval <- 0L; n_bad <- 0L
  dev <- function(v) {
    m <- ppo_max_fixed(P, mle, idx, v)
    n_eval <<- n_eval + 1L
    if (m$convergence != 0L) n_bad <<- n_bad + 1L
    2 * (l_hat - m$loglik) - crit
  }
  one_side <- function(dir) {
    step <- 2 * se
    v  <- mle[idx] + dir * step
    dv <- dev(v)
    while (is.finite(dv) && dv < 0 && step < PROFILE_MAX_SE * se) {
      step <- 2 * step
      v  <- mle[idx] + dir * step
      dv <- dev(v)
    }
    if (!is.finite(dv) || dv < 0) return(NA_real_)
    stats::uniroot(dev, lower = min(mle[idx], v), upper = max(mle[idx], v),
                   tol = PROFILE_UNIROOT_TOL)$root
  }
  ci <- c(one_side(-1), one_side(+1))
  list(ci = ci, n_eval = n_eval, n_not_converged = n_bad, loglik_max = l_hat)
}

## --- 参照 G1・G2 の当てはめごとに、対数尤度の実装を確かめる -------------------
.dat_ref <- function(ref) {
  d <- dat
  d$age_group <- stats::relevel(factor(as.character(d$age_group), levels = AGE_LEVELS),
                                ref = ref)
  d
}
ppo_P   <- list(G1 = ppo_parts(fit_ppo_G1, .dat_ref("G1")),
                G2 = ppo_parts(fit_ppo_G2, .dat_ref("G2")))
ppo_mle <- list(G1 = c(fit_ppo_G1$alpha, fit_ppo_G1$beta),
                G2 = c(fit_ppo_G2$alpha, fit_ppo_G2$beta))
ppo_fit <- list(G1 = fit_ppo_G1, G2 = fit_ppo_G2)

ppo_impl_check <- do.call(rbind, lapply(c("G1", "G2"), function(r) {
  own   <- ppo_loglik(ppo_mle[[r]], ppo_P[[r]])
  refit <- ppo_max_fixed(ppo_P[[r]], ppo_mle[[r]])
  data.frame(reference = r,
             loglik_clm = as.numeric(stats::logLik(ppo_fit[[r]])),
             loglik_own = own,
             loglik_own_remaximised = refit$loglik,
             remax_convergence = refit$convergence,
             stringsAsFactors = FALSE)
}))
print(ppo_impl_check, row.names = FALSE, digits = 10)
check("(i) 自前の対数尤度が clm の推定値で logLik(fit) と一致する",
      all(abs(ppo_impl_check$loglik_own - ppo_impl_check$loglik_clm) < LOGLIK_AGREE_TOL),
      sprintf("max |diff| = %.3e (tol %.0e)",
              max(abs(ppo_impl_check$loglik_own - ppo_impl_check$loglik_clm)),
              LOGLIK_AGREE_TOL))
check("(ii) clm の推定値から最大化し直しても対数尤度が増えない（clm の推定値が最大）",
      all(ppo_impl_check$loglik_own_remaximised - ppo_impl_check$loglik_own <
            LOGLIK_AGREE_TOL),
      sprintf("max gain = %.3e",
              max(ppo_impl_check$loglik_own_remaximised - ppo_impl_check$loglik_own)))

## --- (iii) 位置の係数で ordinal の profile 区間と照合する -------------------
## スプライン以外の位置の係数（性別・疾患区分・病前の要介護状態・入院時期区分）
.loc_parm <- names(fit_ppo_G1$beta)[!IS_SPLINE_TERM(names(fit_ppo_G1$beta))]
say("\n  (iii) own profile vs ordinal::confint(type = \"profile\") for location coefficients ...")
.pl_loc  <- get_profile_ci(fit_ppo_G1, .loc_parm)
.se_ppo1 <- tryCatch(sqrt(diag(vcov(fit_ppo_G1))),
                     error = function(e) stats::setNames(rep(NA_real_, length(ppo_mle$G1)),
                                                         ppo_P$G1$names))
ppo_profile_check <- do.call(rbind, lapply(.loc_parm, function(nm) {
  i  <- match(nm, ppo_P$G1$names)
  pr <- ppo_profile_ci(ppo_P$G1, ppo_mle$G1, i, unname(.se_ppo1[nm]))
  oc <- if (!is.null(.pl_loc$ci) && nm %in% rownames(.pl_loc$ci)) .pl_loc$ci[nm, ] else c(NA, NA)
  data.frame(term = nm, term_label = relabel_vars(nm, strict = FALSE),
             estimate = unname(ppo_mle$G1[i]),
             own_lower = pr$ci[1], own_upper = pr$ci[2],
             ordinal_lower = as.numeric(oc[1]), ordinal_upper = as.numeric(oc[2]),
             stringsAsFactors = FALSE)
}))
ppo_profile_check$max_abs_diff <- pmax(abs(ppo_profile_check$own_lower -
                                             ppo_profile_check$ordinal_lower),
                                       abs(ppo_profile_check$own_upper -
                                             ppo_profile_check$ordinal_upper))
print(ppo_profile_check[, c("term_label", "estimate", "own_lower", "ordinal_lower",
                            "own_upper", "ordinal_upper", "max_abs_diff")],
      row.names = FALSE, digits = 6)
check("(iii) 位置の係数で、自前の profile 区間が ordinal の profile 区間と一致する",
      !is.null(.pl_loc$ci) && all(is.finite(ppo_profile_check$max_abs_diff)) &&
        max(ppo_profile_check$max_abs_diff) < PROFILE_AGREE_TOL,
      if (is.null(.pl_loc$ci)) "ordinal の profile 区間が得られなかった"
      else sprintf("max |diff| = %.2e（係数の尺度、許容 %.0e。ordinal は補間による近似）",
                   max(ppo_profile_check$max_abs_diff), PROFILE_AGREE_TOL))

## --- 年齢群の閾値別オッズ比 --------------------------------------------------
## clm の nominal 母数 τ は logit P(Y ≤ j) に足される。オッズ比は exp(−τ) であり、
## τ の区間 [τ_L, τ_U] はオッズ比の区間 [exp(−τ_U), exp(−τ_L)] に対応する。
ppo_or_rows <- list()
ppo_wald    <- list()
z_wald <- stats::qnorm(1 - (1 - CONF_LEVEL) / 2)
for (thr in PPO_THRESHOLDS) {
  cut <- PPO_CUTS[[as.character(thr)]]
  for (ct in CONTRASTS) {
    r    <- ct$ref
    fit  <- ppo_fit[[r]]
    term <- paste0(cut, ".age_group", ct$cmp)
    idx  <- match(term, ppo_P[[r]]$names)
    if (is.na(idx)) stop("M2-PPO の母数が見つからない: ", term)
    tau  <- unname(ppo_mle[[r]][idx])
    se   <- tryCatch(unname(sqrt(diag(vcov(fit)))[term]), error = function(e) NA_real_)
    t0   <- proc.time()[["elapsed"]]
    pr   <- ppo_profile_ci(ppo_P[[r]], ppo_mle[[r]], idx, se)
    sec  <- proc.time()[["elapsed"]] - t0
    say(sprintf("    %-9s y >= %g : %d evaluations, %.1f s", ct$key, thr, pr$n_eval, sec))
    ppo_or_rows[[length(ppo_or_rows) + 1L]] <- data.frame(
      contrast         = ct$key,
      label            = contrast_label(ct$ref, ct$cmp),
      group_reference  = ct$ref,
      group_comparison = ct$cmp,
      reference_label  = relabel_levels("age_group", ct$ref),
      comparison_label = relabel_levels("age_group", ct$cmp),
      n_reference      = as.integer(n_by_g[[ct$ref]]),
      n_comparison     = as.integer(n_by_g[[ct$cmp]]),
      role             = ct$role,
      threshold        = thr,
      cut              = cut,
      term             = term,
      term_label       = sprintf("%s, P(%s >= %g)",
                                 relabel_vars(paste0("age_group", ct$cmp)),
                                 relabel_vars(OUTCOME), thr),
      beta             = -tau,               # log OR for Y >= threshold
      se               = se,
      odds_ratio       = exp(-tau),
      ci_lower         = exp(-pr$ci[2]),
      ci_upper         = exp(-pr$ci[1]),
      conf_level       = CONF_LEVEL,
      ci_method        = "profile likelihood (own implementation), unadjusted",
      fitted_with      = sprintf("age_group reference = %s", ct$ref),
      profile_evaluations   = pr$n_eval,
      profile_not_converged = pr$n_not_converged,
      p_value          = NA_real_,           # §8.5-1：p 値は付けない
      stringsAsFactors = FALSE
    )
    ppo_wald[[length(ppo_wald) + 1L]] <- data.frame(
      contrast = ct$key, threshold = thr,
      wald_or_lower = exp(-tau - z_wald * se), wald_or_upper = exp(-tau + z_wald * se),
      stringsAsFactors = FALSE)
  }
}
m2_ppo_odds_ratio <- do.call(rbind, ppo_or_rows)
ppo_wald <- do.call(rbind, ppo_wald)

say("\n  M2-PPO : threshold-specific odds ratios with profile likelihood 95% CI (Table 2)")
print(m2_ppo_odds_ratio[, c("threshold", "contrast", "role", "odds_ratio",
                            "ci_lower", "ci_upper")],
      row.names = FALSE, digits = 4)

## --- 点検 -------------------------------------------------------------------
check("M2-PPO のオッズ比が 閾値 × 対比 = 6 行ある",
      nrow(m2_ppo_odds_ratio) == length(PPO_THRESHOLDS) * length(CONTRASTS),
      sprintf("rows = %d", nrow(m2_ppo_odds_ratio)))
check("M2-PPO のすべてのオッズ比で profile 区間が得られている",
      all(is.finite(m2_ppo_odds_ratio$ci_lower) & is.finite(m2_ppo_odds_ratio$ci_upper)),
      paste(with(m2_ppo_odds_ratio,
                 paste0(contrast, "@", threshold)[!is.finite(ci_lower) | !is.finite(ci_upper)]),
            collapse = ", "))
check("M2-PPO の profile の最大化がすべて収束した",
      all(m2_ppo_odds_ratio$profile_not_converged == 0L),
      sprintf("未収束 %d / %d 回", sum(m2_ppo_odds_ratio$profile_not_converged),
              sum(m2_ppo_odds_ratio$profile_evaluations)))
check("M2-PPO の区間が点推定値を挟む",
      all(is.na(m2_ppo_odds_ratio$ci_lower) |
          (m2_ppo_odds_ratio$ci_lower <= m2_ppo_odds_ratio$odds_ratio &
           m2_ppo_odds_ratio$odds_ratio <= m2_ppo_odds_ratio$ci_upper)), "")
for (thr in PPO_THRESHOLDS) {
  .o <- m2_ppo_odds_ratio[m2_ppo_odds_ratio$threshold == thr, , drop = FALSE]
  .l <- stats::setNames(.o$beta, .o$contrast)
  check(sprintf("M2-PPO（y >= %g）：β(G3 vs G1) − β(G2 vs G1) = β(G3 vs G2)（参照の入れ替えの整合）",
                thr),
        abs((.l[["G3_vs_G1"]] - .l[["G2_vs_G1"]]) - .l[["G3_vs_G2"]]) < 1e-6,
        sprintf("%.8f − %.8f = %.8f vs %.8f", .l[["G3_vs_G1"]], .l[["G2_vs_G1"]],
                .l[["G3_vs_G1"]] - .l[["G2_vs_G1"]], .l[["G3_vs_G2"]]))
}

## Wald との対照（ログのみ。本文・図表には用いない）
say("\n  --- M2-PPO：profile と Wald の対照（ログのみ）---")
print(merge(m2_ppo_odds_ratio[, c("threshold", "contrast", "odds_ratio",
                                  "ci_lower", "ci_upper")],
            ppo_wald, by = c("contrast", "threshold"))[, c("threshold", "contrast",
            "odds_ratio", "ci_lower", "ci_upper", "wald_or_lower", "wald_or_upper")],
      row.names = FALSE, digits = 4)
say("  Table 2 に載せるのは profile likelihood 区間である（§4-2）。",
    "ブートストラップ区間との突き合わせは 06 で行う。")


# -----------------------------------------------------------------------------
# 7. 符号の向きの点検（clm の母数化と §8.3 のモデル式の対応）
#    logit P(Y ≤ j) = θ_j − x'β であれば、β が大きいほど P(Y ≥ y) は大きい。
#    予測確率で直接確かめる。
# -----------------------------------------------------------------------------

rule("7. Direction of the coefficients (numerical check)")

## predict.clm は、応答を含まない newdata に対して cprob1（n × 水準数の
## P(Y ≤ 水準) の行列）を返す。列を指定して取り出す（06 も同じ扱いにする）。
.k  <- max(1L, floor(nlevels(dat$y_ord) / 2))
.j  <- levels(dat$y_ord)[.k]
.nd <- dat
.nd$y_ord <- NULL
.pr <- function(g) {
  d <- .nd
  d$age_group <- stats::relevel(factor(rep(g, nrow(d)), levels = AGE_LEVELS), ref = "G1")
  p  <- stats::predict(fit_G1, newdata = d, type = "cum.prob")
  cp <- if (is.list(p) && !is.null(p$cprob1)) p$cprob1 else p
  v  <- if (is.matrix(cp)) cp[, .k] else as.numeric(cp)
  mean(1 - v)
}
.p1 <- .pr("G1"); .p3 <- .pr("G3")
check("β の符号と予測確率の向きが一致する（exp(β) が Y ≥ y のオッズ比である）",
      is.finite(.p1) && is.finite(.p3) &&
        ((.b31 > 0) == (.p3 > .p1) || abs(.b31) < 1e-8),
      sprintf("β(G3 vs G1) = %.4f ; 全員 G1 とした P(Y ≥ %s) = %.4f ; 全員 G3 = %.4f",
              .b31, .j, .p1, .p3))

## --- M2-PPO：nominal 母数の符号（オッズ比 = exp(−τ)）を予測確率で確かめる -----
## 各患者について logit P(Y ≥ t | 比較群) − logit P(Y ≥ t | G1) を予測確率から求め、
## 報告する対数オッズ比（−τ）と全員で一致することを確かめる。
.cp_ppo <- function(g) {
  d <- .nd_ppo
  d$age_group <- stats::relevel(factor(rep(g, nrow(d)), levels = AGE_LEVELS), ref = "G1")
  p <- stats::predict(fit_ppo_G1, newdata = d, type = "cum.prob")
  if (is.list(p) && !is.null(p$cprob1)) p$cprob1 else p
}
.cp_ref <- .cp_ppo("G1")
.sign_diff <- c()
for (thr in PPO_THRESHOLDS) {
  k <- match(thr, PPO_THRESHOLDS)               # 1 − P(Y ≤ 区分 k) = P(Y ≥ thr)
  for (g in setdiff(AGE_LEVELS, "G1")) {
    lo  <- stats::qlogis(1 - .cp_ppo(g)[, k]) - stats::qlogis(1 - .cp_ref[, k])
    rep_b <- m2_ppo_odds_ratio$beta[m2_ppo_odds_ratio$threshold == thr &
                                      m2_ppo_odds_ratio$contrast == paste0(g, "_vs_G1")]
    .sign_diff <- c(.sign_diff, max(abs(lo - rep_b)))
  }
}
check("M2-PPO：報告する対数オッズ比（−τ）が予測確率の logit 差と全患者で一致する",
      length(.sign_diff) == length(PPO_THRESHOLDS) * (length(AGE_LEVELS) - 1L) &&
        all(.sign_diff < 1e-8),
      sprintf("max |diff| = %.3e", max(.sign_diff)))


# -----------------------------------------------------------------------------
# 8. 大域尤度比検定（H0: β_G2 = β_G3 = 0、自由度2）  [LOG ONLY]
#    §8.5 により本文にも補足資料にも載せない。ログにのみ残す。
#    M2-PPO についても、年齢群の母数すべてが 0 という検定（自由度4）をログに残す。
# -----------------------------------------------------------------------------

rule("8. Global likelihood ratio tests (M2: df = 2 ; M2-PPO: df = 4)  [LOG ONLY]")
say("  §8.5：大域検定は本文にも補足資料にも報告しない。記録としてログにのみ残す。")

M2_NULL_TXT <- sub("^y_ord ~ age_group \\+ ", "y_ord ~ ", M2_FORMULA_TXT)
m2_lrt <- NULL
if (identical(M2_NULL_TXT, M2_FORMULA_TXT)) {
  check("帰無モデルの式を作れた", FALSE, "式から age_group を外せなかった")
} else {
  fit_null <- tryCatch(ordinal::clm(stats::as.formula(M2_NULL_TXT), data = dat,
                                    link = "logit", control = ctrl),
                       error = function(e) { say("  帰無モデルでエラー: ",
                                                 conditionMessage(e)); NULL })
  if (!is.null(fit_null)) {
    m2_lrt <- tryCatch(anova(fit_null, fit_G1), error = function(e) NULL)
    if (!is.null(m2_lrt)) { say(""); print(m2_lrt) }
    .st <- 2 * (as.numeric(stats::logLik(fit_G1)) - as.numeric(stats::logLik(fit_null)))
    say(sprintf("\n  LR statistic = %.4f, df = 2, p = %.6g（ログのみ）",
                .st, stats::pchisq(.st, df = 2, lower.tail = FALSE)))
  }
}
## M2-PPO：年齢群の母数すべて（閾値 2 × 年齢群 2 = 自由度4）が 0 という帰無仮説（ログのみ）
m2_ppo_lrt <- NULL
fit_ppo_null <- tryCatch(ordinal::clm(M2_PPO_FORMULA, data = dat, link = "logit",
                                      control = ctrl),
                         error = function(e) { say("  M2-PPO の帰無モデルでエラー: ",
                                                   conditionMessage(e)); NULL })
if (!is.null(fit_ppo_null)) {
  m2_ppo_lrt <- tryCatch(anova(fit_ppo_null, fit_ppo_G1), error = function(e) NULL)
  if (!is.null(m2_ppo_lrt)) { say("\n  M2-PPO:"); print(m2_ppo_lrt) }
  .df_ppo <- length(age_nominal_names(fit_ppo_G1))
  .st_ppo <- 2 * (as.numeric(stats::logLik(fit_ppo_G1)) -
                    as.numeric(stats::logLik(fit_ppo_null)))
  say(sprintf("\n  M2-PPO LR statistic = %.4f, df = %d, p = %.6g（ログのみ）",
              .st_ppo, .df_ppo, stats::pchisq(.st_ppo, df = .df_ppo, lower.tail = FALSE)))
}
check("大域尤度比検定を出力表・図表に入れていない", TRUE,
      "ログと rda にのみ保存する（§8.5）")


# -----------------------------------------------------------------------------
# 9. 補足表4：共変量の係数（解釈しない。§8.3、Table 2 fallacy［8］）
#    §7-2 により、補足表4 は M2-PPO の係数とする（m2_coefficients）。
#    比例オッズモデルの係数は参考として m2_po_coefficients に残す。
#    列の構成（term, term_label, type, estimate, std_error, odds_ratio,
#    is_age_group, interpreted, note）は従来のまま変えない（11 が転記するため）。
#    M2-PPO の行の type は次の3種類である。
#      "threshold"                      閾値（切片）
#      "threshold-specific coefficient" 年齢群の nominal 母数（閾値別）。
#                                       符号を反転し、exp(estimate) を Y ≥ t の
#                                       オッズ比として報告する
#      "regression coefficient"         位置の係数（共変量。2つの閾値で共通）
# -----------------------------------------------------------------------------

rule("9. Supplementary Table 4 : coefficients of M2-PPO (covariates not interpreted)")

NOTE_NOT_INTERPRETED <- "Not interpreted (Table 2 fallacy). Reported for transparency only."

#' 比例オッズモデルの係数表（従来の m2_coefficients と同じ作り方）
coef_table_po <- function(fit) {
  co  <- summary(fit)$coefficients
  nm  <- rownames(co)
  is_beta <- nm %in% names(fit$beta)
  data.frame(
    term       = nm,
    term_label = ifelse(is_beta, relabel_vars(nm, strict = FALSE), nm),
    type       = ifelse(is_beta, "regression coefficient", "threshold"),
    estimate   = as.numeric(co[, 1]),
    std_error  = as.numeric(co[, 2]),
    odds_ratio = ifelse(is_beta, exp(as.numeric(co[, 1])), NA_real_),
    is_age_group = grepl("^age_group", nm),
    interpreted  = grepl("^age_group", nm),
    note = ifelse(grepl("^age_group", nm), "", NOTE_NOT_INTERPRETED),
    stringsAsFactors = FALSE
  )
}

#' M2-PPO の係数表。年齢群の nominal 母数は符号を反転する（exp(−τ) がオッズ比）。
coef_table_ppo <- function(fit) {
  co  <- summary(fit)$coefficients
  nm  <- rownames(co)
  is_beta <- nm %in% names(fit$beta)
  is_age  <- nm %in% age_nominal_names(fit)
  cut_of  <- sub("\\.age_group.*$", "", nm)
  grp_of  <- sub("^.*\\.age_group", "", nm)
  thr_of  <- names(PPO_CUTS)[match(cut_of, PPO_CUTS)]
  est <- as.numeric(co[, 1])
  est[is_age] <- -est[is_age]
  lab <- nm
  lab[is_beta] <- relabel_vars(nm[is_beta], strict = FALSE)
  lab[is_age]  <- sprintf("%s, P(%s >= %s)",
                          relabel_vars(paste0("age_group", grp_of[is_age]), strict = FALSE),
                          relabel_vars(OUTCOME), thr_of[is_age])
  data.frame(
    term       = nm,
    term_label = lab,
    type       = ifelse(is_beta, "regression coefficient",
                        ifelse(is_age, "threshold-specific coefficient", "threshold")),
    estimate   = est,
    std_error  = as.numeric(co[, 2]),
    odds_ratio = ifelse(is_beta | is_age, exp(est), NA_real_),
    is_age_group = is_age,
    interpreted  = is_age,
    note = ifelse(is_age,
                  paste0("Threshold-specific (nominal) effect of age group; sign reversed ",
                         "from the clm parameter so that exp(estimate) is the odds ratio ",
                         "for mFIM at discharge >= ", thr_of, "."),
                  ifelse(is_beta, NOTE_NOT_INTERPRETED, "")),
    stringsAsFactors = FALSE
  )
}

m2_po_coefficients <- coef_table_po(fit_G1)
m2_coefficients    <- coef_table_ppo(fit_ppo_G1)
is_beta <- m2_coefficients$type == "regression coefficient"

say("  M2-PPO coefficients (thresholds are omitted from this printout):")
print(m2_coefficients[m2_coefficients$type != "threshold",
                      c("type", "term_label", "estimate", "std_error", "odds_ratio")],
      row.names = FALSE, digits = 4)
say(sprintf("\n  thresholds : %d （表には含めるが印字は省略）",
            sum(m2_coefficients$type == "threshold")))

check("係数表に年齢群の閾値別の係数が 2 × 2 個含まれる",
      sum(m2_coefficients$is_age_group) ==
        length(PPO_THRESHOLDS) * (length(AGE_LEVELS) - 1L),
      sprintf("n = %d", sum(m2_coefficients$is_age_group)))
check("補足表4 の年齢群の係数が M2-PPO のオッズ比（参照 G1）と一致する",
      {
        .a <- m2_coefficients[m2_coefficients$is_age_group, , drop = FALSE]
        .o <- m2_ppo_odds_ratio[m2_ppo_odds_ratio$group_reference == "G1", , drop = FALSE]
        nrow(.a) == nrow(.o) && all(abs(sort(.a$estimate) - sort(.o$beta)) < 1e-10)
      }, "")
check("補足表4に「解釈しない」旨の注記がある",
      any(m2_coefficients$note == NOTE_NOT_INTERPRETED), "")
check("係数名がすべて英語ラベルに変換できた",
      !any(is_beta & m2_coefficients$term_label == m2_coefficients$term),
      paste(m2_coefficients$term[is_beta &
        m2_coefficients$term_label == m2_coefficients$term], collapse = ", "))
check("補足表4 の列の構成が従来と同じである（11 がそのまま転記できる）",
      identical(names(m2_coefficients), names(m2_po_coefficients)) &&
        identical(names(m2_coefficients),
                  c("term", "term_label", "type", "estimate", "std_error", "odds_ratio",
                    "is_age_group", "interpreted", "note")),
      paste(names(m2_coefficients), collapse = ", "))


# -----------------------------------------------------------------------------
# 10. 出来上がりの点検
# -----------------------------------------------------------------------------

rule("10. Verification")

check("主対比が G3 対 G1 である",
      PRIMARY_KEY %in% m2_odds_ratio$contrast &&
        m2_odds_ratio$role[m2_odds_ratio$contrast == PRIMARY_KEY] == "primary", "")
check("3対比とも同じ構成法・同じ水準の区間である",
      length(unique(m2_odds_ratio$ci_method)) == 1L &&
        length(unique(m2_odds_ratio$conf_level)) == 1L,
      sprintf("method = %s, level = %s",
              paste(unique(m2_odds_ratio$ci_method), collapse = " / "),
              paste(unique(m2_odds_ratio$conf_level), collapse = " / ")))
check("各対比の例数が完全ケース集団の年齢群の例数と一致する",
      all(m2_odds_ratio$n_reference == as.integer(n_by_g[m2_odds_ratio$group_reference])) &&
      all(m2_odds_ratio$n_comparison == as.integer(n_by_g[m2_odds_ratio$group_comparison])), "")
check("M2-PPO：主対比が G3 対 G1 である（両閾値）",
      all(vapply(PPO_THRESHOLDS, function(t)
        any(m2_ppo_odds_ratio$threshold == t & m2_ppo_odds_ratio$contrast == PRIMARY_KEY &
              m2_ppo_odds_ratio$role == "primary"), logical(1))), "")
check("M2-PPO：6つのオッズ比が同じ構成法・同じ水準の区間である",
      length(unique(m2_ppo_odds_ratio$ci_method)) == 1L &&
        length(unique(m2_ppo_odds_ratio$conf_level)) == 1L,
      sprintf("method = %s", paste(unique(m2_ppo_odds_ratio$ci_method), collapse = " / ")))
check("M2-PPO：p 値を置いていない（§8.5-1）", all(is.na(m2_ppo_odds_ratio$p_value)), "")

checks <- do.call(rbind, .checks)
say("\n  REVIEW 項目: ", sum(checks$result == "REVIEW"), " / ", nrow(checks))
if (any(checks$result == "REVIEW")) {
  print(checks[checks$result == "REVIEW", ], row.names = FALSE)
}


# -----------------------------------------------------------------------------
# 11. 保存
# -----------------------------------------------------------------------------

rule("11. Output")

w <- function(x, f) {
  p <- file.path(OUT_DIR, f)
  write.csv(x, p, row.names = FALSE, fileEncoding = "UTF-8")
  say("  written: ", p)
}
w(m2_odds_ratio,      "table_m2_odds_ratio.csv")
w(m2_ppo_odds_ratio,  "table_m2_ppo_odds_ratio.csv")
w(ppo_profile_check,  "table_m2_ppo_profile_check.csv")
w(m2_coefficients,    "table_m2_coefficients.csv")
w(m2_po_coefficients, "table_m2_po_coefficients.csv")
w(m2_fit_info,        "table_m2_fit.csv")
w(checks,             "table_checks_05.csv")

## 06 が使うもの：fit_G1（比例オッズ）・fit_ppo_G1（部分比例オッズ）、両モデルの式、
## PPO_THRESHOLDS / PPO_LEVELS / PPO_CUTS / make_ppo_outcome、m2_odds_ratio、
## m2_ppo_odds_ratio。
save(fit_G1, fit_G2, m2_odds_ratio, m2_coefficients, m2_fit_info, m2_lrt,
     M2_FORMULA, M2_FORMULA_TXT, ctrl, CONTRASTS, PRIMARY_KEY, CONF_LEVEL,
     fit_ppo_G1, fit_ppo_G2, m2_ppo_odds_ratio, m2_po_coefficients, m2_ppo_lrt,
     M2_PPO_FORMULA, M2_PPO_FORMULA_TXT, M2_PPO_NOMINAL, M2_PPO_NOMINAL_TXT,
     PPO_THRESHOLDS, PPO_LEVELS, PPO_CUTS, make_ppo_levels, make_ppo_outcome,
     ppo_cells, ppo_impl_check, ppo_profile_check, ppo_wald,
     checks,
     file = file.path(OUT_DIR, "05_m2_clm.rda"))
say("  written: ", file.path(OUT_DIR, "05_m2_clm.rda"))

say("\n--- sessionInfo ---")
print(utils::sessionInfo())

say("\nDone. M2 N = ", nrow(dat), " ; contrasts = ", nrow(m2_odds_ratio),
    " ; M2-PPO odds ratios = ", nrow(m2_ppo_odds_ratio),
    " (thresholds ", paste(PPO_THRESHOLDS, collapse = ", "), ")",
    " ; primary contrast = ", PRIMARY_KEY)

.log_close()
