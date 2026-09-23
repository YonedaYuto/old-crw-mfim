# =============================================================================
# 06_standardize_boot.R
#   解析計画書 §8.4「周辺標準化」（plan.summary §4-2 の④、E-4）の実装。
#   plan.detail §11 のスクリプト14 に相当。Table 2 の数値をすべて出力する。
#
#   入力 : output/04_ridit_spline.rda（dat_m2、spline_terms、M2_SETTINGS）
#          output/05_m2_clm.rda      （fit_G1、fit_ppo_G1、M2_FORMULA_TXT、
#                                     M2_PPO_FORMULA_TXT、M2_PPO_NOMINAL_TXT、
#                                     PPO_THRESHOLDS、m2_odds_ratio、m2_ppo_odds_ratio）
#          01_labels.R（表示用の文字列）
#   出力 : output/table2.csv                     Table 2（閾値別OR・標準化確率・RD、共通OR）
#          output/table2_standardized_risk.csv   年齢群別の標準化 P(Y ≥ 70)・P(Y ≥ 77)
#                                                （部分比例オッズモデル）
#          output/table2_risk_difference.csv     3対比のリスク差RD（同上、2閾値）
#          output/table2_ppo_odds_ratio.csv      閾値別OR（profile 区間とブートストラップ区間）
#          output/table_po_standardized.csv      比例オッズモデル由来の標準化確率とRD（比較用）
#          output/table_ppo_vs_po.csv            2つのモデルの標準化量の食い違い
#          output/table_or_ci_profile_vs_boot.csv OR の profile 区間とブートストラップ区間の突き合わせ
#          output/table_boot_summary.csv         ブートストラップの諸元と失敗回
#          output/table2_meta.csv                Table 2 の脚注に必要な諸元
#          output/table2_footnotes.csv           Table 2 の脚注（下書き）
#          output/table_checks_06.csv            点検結果
#          output/06_standardize_boot.rda        上記＋再標本の分布
#          output/log_06_standardize_boot.txt    実行ログ
#
#   計画との対応
#     ・標準化の対象は M2 の完全ケース集団。全員の年齢群を G1・G2・G3 に置き換え、
#       他の共変量は観測値のまま動かさない（§8.4-1）
#     ・`predict.clm(type = "cum.prob")` は P(Y ≤ j) を返すため、
#       P(Y ≥ t) = 1 − P(Y ≤ j)。比例オッズモデル（79 段階）では j は
#       「t 点未満の最大の観測水準」（§8.4）、部分比例オッズモデル（3区分）では
#       j は t 点のすぐ下の区分である。用いた j はログに出力する
#     ・RD の向きは「A対B」で RD = リスク(A) − リスク(B)、A は高齢側（§8.4）
#     ・区間はノンパラメトリック・ブートストラップ（2,000回、パーセンタイル法）。
#       再標本化は年齢群で層化して患者単位で行う（§8.4）
#     ・ridit とノットは基準集団で固定し、各回で再算出しない（§6.4-5）。
#       04 が作った項の文字列をそのまま使うことで、これを構造的に保証する
#     ・失敗した回は除外し、補充の再標本化は行わない。失敗回数と割合を
#       モデル別に報告する（§8.4「ブートストラップが失敗した回の扱い」）
#     ・Table 2 に p 値は載せない（§8.5-1）。RR は算出しない（§8.4-3）
#
#   改訂 : 2026-09-22 plan.summary.txt §7-1・§7-2 により、次のとおり変更した。
#          ・閾値を 65 点から 70 点・77 点の2つにした（THRESHOLDS）。
#          ・部分比例オッズモデル（M2-PPO：3区分、年齢群を nominal 項）と
#            比例オッズモデル（M2：79 段階）の両方で、同じ再標本を使って標準化と
#            ブートストラップを行う。Table 2 の標準化確率・RD は M2-PPO 由来とし、
#            M2 由来の値は比較用の表に置く。
#          ・閾値 65 の二値ロジスティックモデル（§8.6(2)）を削除した。
#          ・"fast" 経路と plogis(θ_j − x'β) による手計算の点検は、年齢群の効果が
#            全閾値で共通であることを前提にしていた。nominal 項（閾値ごとの年齢群の
#            母数 τ_j,g）を含む形 plogis(θ_j + τ_j,g − x'β) に一般化した。
#          ・M2-PPO の閾値別オッズ比の区間は 05 の profile likelihood 区間を Table 2 に
#            載せる。同じ再標本からブートストラップ区間も求め、両者に大きな差が
#            ないことを確かめる（著者の判断、2026-09-22）。
#          ・Table 2 の組み立てを、閾値別OR・2閾値の標準化確率とRD・共通OR（併記）
#            に直した。
#
#   文字コード : UTF-8（Windows の R では source(..., encoding = "UTF-8")）
# =============================================================================


# -----------------------------------------------------------------------------
# 0. 設定
# -----------------------------------------------------------------------------

OUT_DIR   <- "output"
IN_RDA_04 <- file.path(OUT_DIR, "04_ridit_spline.rda")
IN_RDA_05 <- file.path(OUT_DIR, "05_m2_clm.rda")

CONF_LEVEL <- 0.95          # 未調整95%（§8.5）
## 標準化する到達確率 P(運動FIM ≥ t) の閾値（plan.summary §7-1）。昇順に置く。
## 70 点：入浴と階段昇降以外が見守りレベル（旧設定 65 点を修正）。
## 77 点：入浴と階段昇降以外が修正自立レベル。05 の PPO_THRESHOLDS と一致させる。
THRESHOLDS <- c(70, 77)

## --- ブートストラップ（§8.4）------------------------------------------------
NBOOT      <- 2000L         # 本番。§8.10・§13 のパイロット計測では 50 にする
BOOT_SEED  <- 20260920L     # 03_m0_ranktest.R の SEED と揃える（著者の決定、2026-09-20）
                            # ※ plan.detail §8.4 の記載は 20260918。計画書側を改訂すること
CI_TYPE    <- 7L            # stats::quantile の type。パーセンタイル法

## 標準化の対象集団をブートストラップの各回でどう取るか。
##   "resample" : その回の再標本を対象集団とする（推定手続き全体を再標本化する）
##   "fixed"    : 元の完全ケース集団に固定する（共変量分布の不確かさを区間に含めない）
## 計画（§8.4）は「各回でモデルを当てはめ直す」「層化により年齢群構成を固定する」と
## 述べており、推定手続き全体の再標本化にあたる "resample" を既定とする。
STD_POP_IN_BOOT <- "resample"

## --- 並列実行（既定は逐次）--------------------------------------------------
## 再標本のインデックスは先に一括生成するため、PAR_CORES の値によらず結果は同一。
PAR_CORES      <- 3L        # 2以上で parallel::makePSOCKcluster を使う
PROGRESS_EVERY <- 50L       # 進捗をログに出す間隔

## --- 標準化リスクの計算方法 --------------------------------------------------
##   "predict.clm" : §8.4 の記載どおり `predict.clm(type = "cum.prob")` を使う（既定）
##   "fast"        : 同じ量を plogis(θ_j(g) − x'β) で直接計算する。θ_j(g) は
##                   比例オッズモデルでは θ_j、部分比例オッズモデルでは
##                   θ_j + τ_j,g（nominal 項の年齢群の母数。参照群は 0）である。
##                   x'β は年齢群の列を 0 にした線形予測子に、比例オッズモデルでは
##                   年齢群の係数を足したもの。predict.clm の 1/1000 以下の時間で済む。
## "fast" は推定量を変えない実装上の変更である（§8.10）。誤りを黙って持ち込まない
## ため、本体の当てはめで両者が一致することを両モデルで毎回確かめ、一致しなければ
## "predict.clm" に戻して REVIEW を立てる。所要時間が問題になる場合にだけ使うこと。
PREDICT_METHOD <- "predict.clm"
PREDICT_AGREE_TOL <- 1e-10

## --- 収束と発散の判定（§8.4）-----------------------------------------------
CLM_MAX_ITER     <- 200L
CLM_MAX_LINE     <- 50L
CLM_MAX_MOD_ITER <- 10L
## |β| の判定はスプライン基底以外の係数に当てる（05_m2_clm.R と同じ扱い）。
## 部分比例オッズモデルでは、nominal 項の年齢群の母数も判定に含める。
DIVERGE_COEF <- 10
DIVERGE_SE   <- 10
IS_SPLINE_TERM <- function(nm) grepl("^ns\\(", nm)
FAIL_WARN_RATE <- 0.05      # 失敗が5%を超えたら解釈を保留する旨を記す（§8.4）

## --- profile 区間とブートストラップ区間の突き合わせ（著者の判断、2026-09-22）--
## 対数オッズ比の尺度で、区間の端のずれの大きい方を profile 区間の幅で割った値。
## これがこの値を超えたら REVIEW とする（差の大きさは表に数値で残す）。
BOOT_PROFILE_TOL <- 0.10

## --- 対比（§8.4「対比の向き」）----------------------------------------------
CONTRASTS <- list(
  list(key = "G2_vs_G1", ref = "G1", cmp = "G2", role = "secondary"),
  list(key = "G3_vs_G1", ref = "G1", cmp = "G3", role = "primary"),
  list(key = "G3_vs_G2", ref = "G2", cmp = "G3", role = "secondary")
)
PRIMARY_KEY <- "G3_vs_G1"

## --- モデルの呼び名 ----------------------------------------------------------
MODEL_PPO <- "partial proportional odds"
MODEL_PO  <- "proportional odds"


# -----------------------------------------------------------------------------
# 1. 表示用ラベル
# -----------------------------------------------------------------------------

if (!file.exists("01_labels.R")) {
  stop("01_labels.R が見つからない。06_standardize_boot.R は 01_labels.R と同じ",
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

.log_con <- file(file.path(OUT_DIR, "log_06_standardize_boot.txt"), open = "wt",
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

call_known_args <- function(fun, args) {
  fml  <- setdiff(names(formals(fun)), "...")
  drop <- setdiff(names(args), fml)
  if (length(drop)) {
    say("    （この版にない引数は渡さない: ", paste(drop, collapse = ", "), "）")
  }
  do.call(fun, args[names(args) %in% fml])
}

rule(paste0("06_standardize_boot.R  ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
say("R : ", R.version.string)
say("nboot : ", NBOOT, " ; seed : ", BOOT_SEED, " ; conf.level : ", CONF_LEVEL)
say("thresholds : P(mFIM at discharge >= t), t = ", paste(THRESHOLDS, collapse = ", "))
say("standardisation population inside the bootstrap : ", STD_POP_IN_BOOT)
say("parallel cores : ", PAR_CORES)

has_ordinal <- requireNamespace("ordinal", quietly = TRUE)
check("ordinal が利用できる", has_ordinal, "")
if (!has_ordinal) stop("ordinal が無ければ E-4 は算出できない。")
library(ordinal)
library(splines)


# -----------------------------------------------------------------------------
# 3. 入力の取り込み
# -----------------------------------------------------------------------------

rule("3. Input")

if (!exists("dat_m2", inherits = TRUE) || !exists("spline_terms", inherits = TRUE)) {
  if (!file.exists(IN_RDA_04)) {
    stop(IN_RDA_04, " が見つからない。先に 04_ridit_spline.R を実行すること。")
  }
  load(IN_RDA_04); say("loaded : ", IN_RDA_04)
}
if (!exists("fit_G1", inherits = TRUE) || !exists("M2_FORMULA_TXT", inherits = TRUE) ||
    !exists("fit_ppo_G1", inherits = TRUE) || !exists("M2_PPO_FORMULA_TXT", inherits = TRUE)) {
  if (!file.exists(IN_RDA_05)) {
    stop(IN_RDA_05, " が見つからない。先に 05_m2_clm.R を実行すること。")
  }
  load(IN_RDA_05); say("loaded : ", IN_RDA_05)
}
if (!exists("fit_ppo_G1", inherits = TRUE) || !exists("m2_ppo_odds_ratio", inherits = TRUE)) {
  stop("05_m2_clm.rda に部分比例オッズモデル（fit_ppo_G1 ほか）が無い。",
       "§7 版の 05_m2_clm.R を先に実行すること。")
}

OUTCOME    <- M2_SETTINGS$outcome
AGE_LEVELS <- M2_SETTINGS$age_levels

## 部分比例オッズモデルのアウトカム（3区分）。05 と同じ規則で作る。
make_ppo_levels <- function(thresholds) {
  th <- sort(thresholds)
  mid <- if (length(th) > 1L) sprintf("%g-%g", head(th, -1L), tail(th, -1L) - 1) else character(0)
  c(sprintf("<%g", th[1]), mid, sprintf(">=%g", th[length(th)]))
}
make_ppo_outcome <- function(y, thresholds) {
  lev <- make_ppo_levels(thresholds)
  factor(lev[findInterval(y, sort(thresholds)) + 1L], levels = lev, ordered = TRUE)
}
PPO_LEVELS_06 <- make_ppo_levels(THRESHOLDS)

## 年齢群の参照を G1 に固定したデータ（05 の fit_G1・fit_ppo_G1 と同じ並び）
DAT <- as.data.frame(dat_m2, stringsAsFactors = FALSE)
DAT$age_group <- factor(as.character(DAT$age_group), levels = AGE_LEVELS)
DAT$y_ppo <- make_ppo_outcome(DAT[[OUTCOME]], THRESHOLDS)
N_M2 <- nrow(DAT)
n_by_g <- table(DAT$age_group)

say("dat_m2 (M2 complete-case population) : ", N_M2, " rows")
for (g in AGE_LEVELS) {
  say(sprintf("  %-3s %-14s n = %5d", g, relabel_levels("age_group", g),
              as.integer(n_by_g[[g]])))
}
say("\nmodel formula (M2-PPO) : ", M2_PPO_FORMULA_TXT, " ; nominal = ", M2_PPO_NOMINAL_TXT)
say("model formula (M2)     : ", M2_FORMULA_TXT)

check("06 の閾値が 05 の部分比例オッズモデルの閾値と一致する",
      identical(as.numeric(THRESHOLDS), as.numeric(PPO_THRESHOLDS)) &&
        identical(PPO_LEVELS_06, PPO_LEVELS),
      sprintf("06: %s / 05: %s", paste(THRESHOLDS, collapse = ", "),
              paste(PPO_THRESHOLDS, collapse = ", ")))
check("3区分のアウトカムが 05 の当てはめと同じ度数である",
      identical(as.integer(table(DAT$y_ppo)), as.integer(table(fit_ppo_G1$y))),
      paste(sprintf("%s=%d", PPO_LEVELS_06, as.integer(table(DAT$y_ppo))), collapse = ", "))
check("標準化の対象集団は M2 の完全ケース集団である（§8.4-4）", TRUE,
      sprintf("N = %d", N_M2))
check("ridit とノットを再算出していない（04 の項の文字列をそのまま使う。§6.4-5）",
      all(grepl("Boundary.knots = c\\(", unname(spline_terms))) &&
        all(vapply(unname(spline_terms), function(s)
          grepl(s, M2_FORMULA_TXT, fixed = TRUE) && grepl(s, M2_PPO_FORMULA_TXT, fixed = TRUE),
          logical(1))),
      "2つのモデル式にノットの数値が literal で埋め込まれている")


# -----------------------------------------------------------------------------
# 4. 閾値の水準 j の決定（§8.4「P(Y ≥ 65) の算出」を 70 点・77 点に適用）
#    P(Y ≥ t) = 1 − P(Y ≤ j)。
#    比例オッズモデル（79 段階）：j は「t 点未満の最大の観測水準」。
#    部分比例オッズモデル（3区分）：j は t 点のすぐ下の区分（<70、70-76）。
# -----------------------------------------------------------------------------

rule("4. The level j used for P(Y >= t) = 1 - P(Y <= j)")

#' 観測水準（数値）から、閾値未満の最大の水準を返す。無ければ NA。
pick_j <- function(levels_num, threshold) {
  lt <- levels_num[levels_num < threshold]
  if (!length(lt)) return(NA_real_)
  max(lt)
}

THR_KEYS <- paste0("t", THRESHOLDS)            # 列名などに使う閾値の鍵（"t70"、"t77"）

Y_LEVELS_MAIN <- as.numeric(levels(DAT$y_ord))
J_PO <- stats::setNames(vapply(THRESHOLDS, function(t) pick_j(Y_LEVELS_MAIN, t), numeric(1)),
                        THR_KEYS)
J_PO_CHR  <- stats::setNames(as.character(J_PO), THR_KEYS)
J_PPO_CHR <- stats::setNames(PPO_LEVELS_06[seq_along(THRESHOLDS)], THR_KEYS)

say(sprintf("  observed outcome levels : %d (range %g - %g)",
            length(Y_LEVELS_MAIN), min(Y_LEVELS_MAIN), max(Y_LEVELS_MAIN)))
for (i in seq_along(THRESHOLDS)) {
  t <- THRESHOLDS[i]
  say(sprintf("  threshold %d :", t))
  say(sprintf("    M2     (79 levels)  : P(Y >= %d) = 1 - P(Y <= %s)  （%d 点未満の最大の観測水準）",
              t, J_PO_CHR[[i]], t))
  say(sprintf("    M2-PPO (3 categories): P(Y >= %d) = 1 - P(Y <= \"%s\")", t, J_PPO_CHR[[i]]))
  if (!(t %in% Y_LEVELS_MAIN)) {
    say(sprintf("    注意：%d 点そのものは観測されていない。§8.4 の規則により j = %s を使う。",
                t, J_PO_CHR[[i]]))
  }
}

check("M2 の j が閾値ごとに決まった", all(is.finite(J_PO)),
      paste(sprintf("t = %g: j = %s", THRESHOLDS, J_PO_CHR), collapse = ", "))
check("M2 の j が閾値未満の最大の観測水準である",
      all(vapply(seq_along(THRESHOLDS), function(i) {
        j <- J_PO[[i]]; t <- THRESHOLDS[i]
        is.finite(j) && j < t && !any(Y_LEVELS_MAIN > j & Y_LEVELS_MAIN < t)
      }, logical(1))), "")
check("M2 の j が最大水準ではない（P(Y ≤ j) = 1 にならない）",
      all(is.finite(J_PO) & J_PO < max(Y_LEVELS_MAIN)), "")
check("M2-PPO の j が閾値のすぐ下の区分である",
      all(vapply(seq_along(THRESHOLDS), function(i) {
        lv <- as.character(make_ppo_outcome(c(THRESHOLDS[i] - 1, THRESHOLDS[i]), THRESHOLDS))
        identical(lv[1], J_PPO_CHR[[i]]) && !identical(lv[2], J_PPO_CHR[[i]])
      }, logical(1))),
      paste(sprintf("t = %g: \"%s\"", THRESHOLDS, J_PPO_CHR), collapse = ", "))

## 素の（調整しない）到達割合。記述であり、標準化リスクとは別物である（Fig. 2(a) と同じ量）
crude <- vapply(THRESHOLDS, function(t) vapply(AGE_LEVELS, function(g)
  mean(DAT[[OUTCOME]][DAT$age_group == g] >= t), numeric(1)),
  numeric(length(AGE_LEVELS)))
crude <- matrix(crude, nrow = length(AGE_LEVELS), dimnames = list(AGE_LEVELS, THR_KEYS))
say("\n  crude (unadjusted) proportion reaching each threshold, for reference only:")
for (g in AGE_LEVELS) {
  say(sprintf("    %-3s %-14s %s", g, relabel_levels("age_group", g),
              paste(sprintf("P(Y >= %d) = %.4f", THRESHOLDS, crude[g, ]), collapse = " ; ")))
}


# -----------------------------------------------------------------------------
# 5. 周辺標準化の関数（§8.4）
# -----------------------------------------------------------------------------

rule("5. Marginal standardisation")

ctrl <- call_known_args(ordinal::clm.control,
                        list(maxIter = CLM_MAX_ITER, maxLineIter = CLM_MAX_LINE,
                             maxModIter = CLM_MAX_MOD_ITER))

M2_PPO_NOMINAL_F <- stats::as.formula(M2_PPO_NOMINAL_TXT)

## clm の $convergence は長さ3のリスト（コード・要約・詳細）である。
conv_code <- function(fit) {
  cv <- fit[["convergence"]]
  if (is.null(cv)) return(NA_integer_)
  if (is.list(cv)) cv <- cv[[1]]
  suppressWarnings(as.integer(cv[1]))
}

#' nominal 項に置いた年齢群の母数の名前（比例オッズモデルでは空）
age_nominal_names <- function(fit) grep("\\.age_group", names(fit[["alpha"]]), value = TRUE)
is_nominal_fit    <- function(fit) !is.null(fit[["nom.terms"]]) && length(age_nominal_names(fit)) > 0L

#' 応答の水準から clm の閾値名を作る（"69|70"、"<70|70-76" など）
cut_names <- function(ylev) paste(head(ylev, -1L), tail(ylev, -1L), sep = "|")

#' 閾値 j・年齢群 g の θ_j(g)。比例オッズモデルでは θ_j、部分比例オッズモデルでは
#' θ_j + τ_j,g（参照群は τ = 0）。clm の母数の名前で取り出す。
theta_jg <- function(fit, j_chr, g, ylev) {
  k <- match(j_chr, ylev)
  if (is.na(k) || k >= length(ylev)) stop("j が応答の水準に見つからない: ", j_chr)
  cut <- cut_names(ylev)[k]
  a <- fit[["alpha"]]
  if (is_nominal_fit(fit)) {
    int <- paste0(cut, ".(Intercept)")
    if (!(int %in% names(a))) stop("閾値の母数が見つからない: ", int)
    tn <- paste0(cut, ".age_group", g)
    unname(a[[int]]) + (if (tn %in% names(a)) unname(a[[tn]]) else 0)
  } else {
    if (cut %in% names(a)) unname(a[[cut]]) else unname(a[[k]])
  }
}

#' 累積ロジットモデルからの標準化リスク P(Y ≥ t)（年齢群 × 閾値 の行列）
#' @param fit   clm の当てはめ（比例オッズ・部分比例オッズのどちらでもよい）
#' @param data  標準化の対象集団（モデルの共変量をもつ data.frame）
#' @param j_chr 1 − P(Y ≤ j) の j（文字列のベクトル。名前は閾値の鍵）
#' @param ylev  当てはめに使った応答の水準（昇順の文字列）
#'
#' `predict.clm(type = "cum.prob")` は、応答を含まない newdata に対して
#' cprob1 = P(Y ≤ 水準) の n × 水準数の行列を返す。j に対応する列を取り、
#' P(Y ≥ t) = 1 − P(Y ≤ j) とする（§8.4）。nominal 項の年齢群も newdata の
#' 年齢群から計算される。
std_risk_clm <- function(fit, data, j_chr, groups, ylev) {
  k <- match(j_chr, ylev)
  if (anyNA(k)) stop("j が応答の水準に見つからない: ", paste(j_chr[is.na(k)], collapse = ", "))
  nd0 <- data
  nd0$y_ord <- NULL                       # 応答は予測に渡さない
  nd0$y_ppo <- NULL
  out <- matrix(NA_real_, nrow = length(groups), ncol = length(j_chr),
                dimnames = list(groups, names(j_chr)))
  for (g in groups) {
    nd <- nd0
    nd$age_group <- factor(rep(g, nrow(nd)), levels = levels(data$age_group))
    p  <- stats::predict(fit, newdata = nd, type = "cum.prob")
    cp <- if (is.list(p) && !is.null(p$cprob1)) p$cprob1 else p
    if (!is.matrix(cp)) stop("predict.clm が行列を返さなかった")
    out[g, ] <- colMeans(1 - cp[, k, drop = FALSE])
  }
  out
}

#' 同じ量を plogis(θ_j(g) − x'β) で直接計算する（PREDICT_METHOD = "fast"）。
#' 位置の項は年齢群の列を 0 にした線形予測子 eta0 に、比例オッズモデルでは
#' 年齢群の係数を足す（参照群は 0）。閾値の側は θ_j(g)（theta_jg）を使う。
#' 部分比例オッズモデルでは位置の項に年齢群がないため、年齢群は θ_j(g) にだけ入る。
std_risk_clm_fast <- function(fit, data, j_chr, groups, ylev) {
  b  <- fit[["beta"]]
  mm <- stats::model.matrix(stats::delete.response(stats::terms(fit)), data = data)
  if (!all(names(b) %in% colnames(mm))) {
    stop("model.matrix の列が係数名と対応しない")
  }
  X <- mm[, names(b), drop = FALSE]
  ageco <- grep("^age_group", names(b), value = TRUE)
  if (length(ageco)) X[, ageco] <- 0
  eta0 <- as.numeric(X %*% b)
  out <- matrix(NA_real_, nrow = length(groups), ncol = length(j_chr),
                dimnames = list(groups, names(j_chr)))
  for (g in groups) {
    tm  <- paste0("age_group", g)
    add <- if (tm %in% names(b)) unname(b[[tm]]) else 0
    for (i in seq_along(j_chr)) {
      th <- theta_jg(fit, j_chr[[i]], g, ylev)
      out[g, i] <- mean(1 - stats::plogis(th - (eta0 + add)))
    }
  }
  out
}

#' 設定に従って標準化リスクを返す
std_risk <- function(fit, data, j_chr, groups, ylev) {
  if (identical(PREDICT_METHOD, "fast")) {
    std_risk_clm_fast(fit, data, j_chr, groups, ylev)
  } else {
    std_risk_clm(fit, data, j_chr, groups, ylev)
  }
}

#' 群別リスク（年齢群 × 閾値）から3対比のRDを作る（RD = リスク(高齢側) − リスク(若年側)）
risk_to_rd <- function(risk) {
  out <- t(vapply(CONTRASTS, function(ct) risk[ct$cmp, ] - risk[ct$ref, ],
                  numeric(ncol(risk))))
  if (ncol(risk) == 1L) out <- t(out)
  dimnames(out) <- list(RD_KEYS, colnames(risk))
  out
}
RD_KEYS <- vapply(CONTRASTS, function(ct) ct$key, character(1))

#' 行列（行 = 年齢群または対比、列 = 閾値）を名前つきベクトルにする（"G1_t70" など）
flat <- function(m) stats::setNames(as.vector(m),
                                    as.vector(outer(rownames(m), colnames(m), paste, sep = "_")))
RISK_KEYS <- as.vector(outer(AGE_LEVELS, THR_KEYS, paste, sep = "_"))
RDT_KEYS  <- as.vector(outer(RD_KEYS,    THR_KEYS, paste, sep = "_"))

#' 部分比例オッズモデルの閾値別の対数オッズ比（Y ≥ t）。3対比 × 閾値。
#' clm の nominal 母数 τ は logit P(Y ≤ j) に足されるため、log OR = −(τ_cmp − τ_ref)。
logor_ppo <- function(fit, ylev) {
  m <- matrix(NA_real_, nrow = length(RD_KEYS), ncol = length(THRESHOLDS),
              dimnames = list(RD_KEYS, THR_KEYS))
  for (i in seq_along(THRESHOLDS)) {
    for (ct in CONTRASTS) {
      m[ct$key, i] <- -(theta_jg(fit, J_PPO_CHR[[i]], ct$cmp, ylev) -
                          theta_jg(fit, J_PPO_CHR[[i]], ct$ref, ylev))
    }
  }
  m
}
#' 比例オッズモデルの共通対数オッズ比（3対比）
logor_po <- function(fit) {
  b <- fit[["beta"]]
  bg <- function(g) { tm <- paste0("age_group", g); if (tm %in% names(b)) unname(b[[tm]]) else 0 }
  stats::setNames(vapply(CONTRASTS, function(ct) bg(ct$cmp) - bg(ct$ref), numeric(1)), RD_KEYS)
}

## --- 本体の推定値 -----------------------------------------------------------
Y_LEVELS_CHR <- levels(DAT$y_ord)
PPO_LEVELS_CHR <- levels(DAT$y_ppo)

FITS  <- list(ppo = fit_ppo_G1, po = fit_G1)
JCHR  <- list(ppo = J_PPO_CHR,  po = J_PO_CHR)
YLEVS <- list(ppo = PPO_LEVELS_CHR, po = Y_LEVELS_CHR)
MODEL_NAME <- c(ppo = MODEL_PPO, po = MODEL_PO)

check("fit_ppo_G1 が nominal 項をもち、fit_G1 がもたない",
      is_nominal_fit(fit_ppo_G1) && !is_nominal_fit(fit_G1),
      paste(age_nominal_names(fit_ppo_G1), collapse = ", "))

## --- 2つの計算方法が一致することを本体の当てはめで確かめる（§8.10）----------
risk_pkg <- list(); agree <- logical(0)
for (m in names(FITS)) {
  t0 <- proc.time()[["elapsed"]]
  risk_pkg[[m]] <- std_risk_clm(FITS[[m]], DAT, JCHR[[m]], AGE_LEVELS, YLEVS[[m]])
  t_pkg <- proc.time()[["elapsed"]] - t0
  t0 <- proc.time()[["elapsed"]]
  rf <- tryCatch(std_risk_clm_fast(FITS[[m]], DAT, JCHR[[m]], AGE_LEVELS, YLEVS[[m]]),
                 error = function(e) { say("  fast path error: ", conditionMessage(e)); NULL })
  t_fast <- proc.time()[["elapsed"]] - t0
  agree[m] <- !is.null(rf) && max(abs(risk_pkg[[m]] - rf)) < PREDICT_AGREE_TOL
  say(sprintf("  [%s] predict.clm : %.3f s / fast : %.3f s （1回の標準化あたり、3群 × %d 閾値）",
              m, t_pkg, t_fast, length(THRESHOLDS)))
  check(sprintf("[%s] 2つの計算方法（predict.clm と fast）が一致する", MODEL_NAME[[m]]),
        agree[[m]],
        if (is.null(rf)) "fast の計算に失敗した"
        else sprintf("max |diff| = %.3e (tol %.0e)", max(abs(risk_pkg[[m]] - rf)),
                     PREDICT_AGREE_TOL))
}
if (identical(PREDICT_METHOD, "fast") && !all(agree)) {
  PREDICT_METHOD <- "predict.clm"
  check("PREDICT_METHOD を predict.clm に戻した", FALSE,
        "fast が一致しなかったため、§8.4 の記載どおりの実装に戻す")
}
say("  PREDICT_METHOD in the bootstrap : ", PREDICT_METHOD)

risk_point    <- risk_pkg$ppo                  # Table 2（M2-PPO）
rd_point      <- risk_to_rd(risk_point)
risk_point_po <- risk_pkg$po                   # 比較用（M2）
rd_point_po   <- risk_to_rd(risk_point_po)
logor_point_ppo <- logor_ppo(fit_ppo_G1, PPO_LEVELS_CHR)
logor_point_po  <- logor_po(fit_G1)

for (m in names(FITS)) {
  rk <- risk_pkg[[m]]; rd <- risk_to_rd(rk)
  say(sprintf("\n  standardised risk, %s (point estimates):", MODEL_NAME[[m]]))
  for (g in AGE_LEVELS) {
    say(sprintf("    %-3s %-14s %s", g, relabel_levels("age_group", g),
                paste(sprintf("P(Y >= %d) = %.6f", THRESHOLDS, rk[g, ]), collapse = " ; ")))
  }
  say("  risk differences:")
  for (k in RD_KEYS) {
    say(sprintf("    %-10s %s", k,
                paste(sprintf("t = %d: %+.6f", THRESHOLDS, rd[k, ]), collapse = " ; ")))
  }
  check(sprintf("[%s] 標準化リスクが (0, 1) に収まる", MODEL_NAME[[m]]),
        all(is.finite(rk) & rk > 0 & rk < 1),
        paste(sprintf("%s=%.4f", names(flat(rk)), flat(rk)), collapse = ", "))
  check(sprintf("[%s] 標準化リスクが閾値の高い方で小さい（P(Y ≥ 77) ≤ P(Y ≥ 70)）", MODEL_NAME[[m]]),
        ncol(rk) < 2L || all(apply(rk, 1L, function(r) all(diff(r) <= 1e-12))), "")
  check(sprintf("[%s] RD の推移性（RD(G3,G1) = RD(G3,G2) + RD(G2,G1)、閾値ごと）", MODEL_NAME[[m]]),
        all(abs(rd["G3_vs_G1", ] - (rd["G3_vs_G2", ] + rd["G2_vs_G1", ])) < 1e-10), "")
}

## --- predict.clm の返り値を手計算で裏づける（母数化の取り違えを防ぐ）--------
## clm は logit P(Y ≤ j) = θ_j + τ_j'n − x'β。比例オッズモデルでは nominal 項がない。
## 観測された年齢群のまま、全患者の P(Y ≤ j) を両モデル・両閾値で照合する。
manual_cum <- function(fit, data, j_chr, ylev) {
  b   <- fit[["beta"]]
  mm  <- stats::model.matrix(stats::delete.response(stats::terms(fit)), data = data)
  if (!all(names(b) %in% colnames(mm))) stop("model.matrix を作れない")
  eta <- as.numeric(mm[, names(b), drop = FALSE] %*% b)
  a   <- fit[["alpha"]]
  vapply(j_chr, function(j) {
    k   <- match(j, ylev)
    cut <- cut_names(ylev)[k]
    th  <- if (is_nominal_fit(fit)) {
      Nm <- stats::model.matrix(fit[["nom.terms"]], data = data)
      as.numeric(Nm %*% a[paste0(cut, ".", colnames(Nm))])
    } else rep(unname(a[[k]]), nrow(data))
    stats::plogis(th - eta)
  }, numeric(nrow(data)))
}
manual_max_diff <- numeric(0)
for (m in names(FITS)) {
  manual_max_diff[m] <- tryCatch({
    pman <- manual_cum(FITS[[m]], DAT, JCHR[[m]], YLEVS[[m]])
    nd   <- DAT; nd$y_ord <- NULL; nd$y_ppo <- NULL
    p    <- stats::predict(FITS[[m]], newdata = nd, type = "cum.prob")
    cp   <- if (is.list(p) && !is.null(p$cprob1)) p$cprob1 else p
    max(abs(pman - cp[, match(JCHR[[m]], YLEVS[[m]]), drop = FALSE]))
  }, error = function(e) NA_real_)
  check(sprintf("[%s] predict.clm(type = \"cum.prob\") が手計算 plogis(θ_j(g) − x'β) と一致する",
                MODEL_NAME[[m]]),
        is.finite(manual_max_diff[[m]]) && manual_max_diff[[m]] < 1e-8,
        if (!is.finite(manual_max_diff[[m]])) "手計算の照合を行えなかった（model.matrix を作れない）"
        else sprintf("max |diff| = %.3e", manual_max_diff[[m]]))
}

## --- 05 の閾値別オッズ比と、ここで θ_j(g) から求めた値が一致すること ----------
.o5 <- m2_ppo_odds_ratio
.d5 <- vapply(seq_len(nrow(.o5)), function(i)
  abs(.o5$beta[i] - logor_point_ppo[.o5$contrast[i], paste0("t", .o5$threshold[i])]),
  numeric(1))
check("M2-PPO の閾値別対数オッズ比が 05 の m2_ppo_odds_ratio と一致する",
      nrow(.o5) == length(RDT_KEYS) && all(.d5 < 1e-6),
      sprintf("max |diff| = %.3e", max(.d5)))
.d5po <- abs(m2_odds_ratio$beta - logor_point_po[m2_odds_ratio$contrast])
check("M2 の共通対数オッズ比が 05 の m2_odds_ratio と一致する",
      all(.d5po < 1e-6), sprintf("max |diff| = %.3e", max(.d5po)))


# -----------------------------------------------------------------------------
# 6. ブートストラップ（§8.4）
#    年齢群で層化した患者単位の復元抽出。インデックスは先に一括生成するため、
#    逐次でも並列でも結果は同一である。
#    各回で M2-PPO と M2 の2本を当てはめ、両閾値の標準化リスクと、
#    オッズ比（M2-PPO の閾値別、M2 の共通）を記録する。
# -----------------------------------------------------------------------------

rule("6. Stratified nonparametric bootstrap")

idx_by_g <- lapply(AGE_LEVELS, function(g) which(DAT$age_group == g))
names(idx_by_g) <- AGE_LEVELS

set.seed(BOOT_SEED)
BOOT_IDX <- vector("list", NBOOT)
for (b in seq_len(NBOOT)) {
  BOOT_IDX[[b]] <- unlist(lapply(AGE_LEVELS, function(g)
    sample(idx_by_g[[g]], length(idx_by_g[[g]]), replace = TRUE)),
    use.names = FALSE)
}
say(sprintf("  %d resamples generated (stratified by age group; sizes fixed at %s)",
            NBOOT, paste(sprintf("%s=%d", AGE_LEVELS,
                                 vapply(idx_by_g, length, integer(1))), collapse = ", ")))
check("再標本の大きさが各回で観測値に固定されている",
      all(vapply(BOOT_IDX, length, integer(1)) == N_M2),
      sprintf("N = %d", N_M2))

#' 当てはめの成否（収束と発散）を判定する。発散の判定に使う係数は、
#' スプライン以外の位置の係数と、nominal 項の年齢群の母数（M2-PPO のみ）。
fit_status <- function(f, tag) {
  cc   <- conv_code(f)
  conv <- is.na(cc) || cc == 0L
  an   <- age_nominal_names(f)
  bb   <- c(f[["beta"]], f[["alpha"]][an])
  bns  <- bb[!IS_SPLINE_TERM(names(bb))]
  se   <- tryCatch(sqrt(diag(vcov(f)))[names(bb)],
                   error = function(e) rep(NA_real_, length(bb)))
  div  <- !all(is.finite(bb)) ||
          (length(bns) && max(abs(bns)) > DIVERGE_COEF) ||
          (all(is.na(se)) || max(se, na.rm = TRUE) > DIVERGE_SE)
  if (!conv) return(paste(tag, "did not converge"))
  if (div)   return(paste(tag, "coefficients or standard errors diverged"))
  ""
}

## --- 1回分の処理 ------------------------------------------------------------
## 返り値：M2-PPO・M2 それぞれの成否、標準化リスク（年齢群 × 閾値）、対数オッズ比。
## 失敗の条件（§8.4）：収束しない／係数・標準誤差の発散／共変量の水準の欠落／
## いずれかの閾値の上または下に観測が存在しない。M2-PPO では、3区分のいずれかに
## 観測がない回も失敗とする（nominal 項が推定できない）。
boot_one <- function(idx) {
  out <- list(ok_ppo = FALSE, ok_po = FALSE,
              risk_ppo = stats::setNames(rep(NA_real_, length(RISK_KEYS)), RISK_KEYS),
              risk_po  = stats::setNames(rep(NA_real_, length(RISK_KEYS)), RISK_KEYS),
              logor_ppo = stats::setNames(rep(NA_real_, length(RDT_KEYS)), RDT_KEYS),
              logor_po  = stats::setNames(rep(NA_real_, length(RD_KEYS)), RD_KEYS),
              j_used = stats::setNames(rep(NA_real_, length(THRESHOLDS)), THR_KEYS),
              msg_ppo = "", msg_po = "")

  d <- DAT[idx, , drop = FALSE]

  ## (i) 各閾値の両側に観測があるか
  y <- d[[OUTCOME]]
  if (!all(vapply(THRESHOLDS, function(t) any(y >= t) && any(y < t), logical(1)))) {
    out$msg_ppo <- out$msg_po <- "no observation on one side of a threshold"
    return(out)
  }

  ## (ii) 共変量の水準の欠落
  for (v in c("age_group", "sex", "class", "support_in", "period")) {
    if (is.factor(DAT[[v]])) {
      d[[v]] <- factor(as.character(d[[v]]), levels = levels(DAT[[v]]))
      if (any(table(d[[v]]) == 0L)) {
        out$msg_ppo <- out$msg_po <- paste0("empty level in ", v)
        return(out)
      }
    }
  }

  ## 標準化の対象集団（§8.4）。応答は予測に使わないため水準を揃える必要はない。
  pop <- if (identical(STD_POP_IN_BOOT, "fixed")) DAT else d

  ## --- M2-PPO（3区分、年齢群を nominal 項）-----------------------------------
  d$y_ppo <- make_ppo_outcome(y, THRESHOLDS)
  if (any(table(d$y_ppo) == 0L)) {
    out$msg_ppo <- "empty category of the PPO outcome"
  } else {
    f <- tryCatch(ordinal::clm(stats::as.formula(M2_PPO_FORMULA_TXT),
                               nominal = M2_PPO_NOMINAL_F, data = d,
                               link = "logit", control = ctrl),
                  error = function(e) e)
    if (inherits(f, "error")) {
      out$msg_ppo <- paste("clm error:", conditionMessage(f))
    } else {
      st <- fit_status(f, "clm (PPO)")
      if (nzchar(st)) {
        out$msg_ppo <- st
      } else {
        r <- tryCatch(std_risk(f, pop, J_PPO_CHR, AGE_LEVELS, PPO_LEVELS_CHR),
                      error = function(e) e)
        if (inherits(r, "error")) {
          out$msg_ppo <- paste("predict error:", conditionMessage(r))
        } else if (!all(is.finite(r))) {
          out$msg_ppo <- "non-finite standardised risk"
        } else {
          out$risk_ppo  <- flat(r)[RISK_KEYS]
          out$logor_ppo <- flat(logor_ppo(f, PPO_LEVELS_CHR))[RDT_KEYS]
          out$ok_ppo    <- TRUE
        }
      }
    }
  }

  ## --- M2（79 段階、比例オッズ）----------------------------------------------
  ## 応答の順序因子を再標本の観測水準で作り直す
  ylev <- sort(unique(y))
  d$y_ord <- factor(as.character(y), levels = as.character(ylev), ordered = TRUE)
  j <- vapply(THRESHOLDS, function(t) pick_j(ylev, t), numeric(1))
  if (!all(is.finite(j)) || any(j >= max(ylev))) {
    out$msg_po <- "j is not defined in this resample"
    return(out)
  }
  out$j_used <- stats::setNames(j, THR_KEYS)
  j_chr    <- stats::setNames(as.character(j), THR_KEYS)
  ylev_chr <- as.character(ylev)

  f <- tryCatch(ordinal::clm(stats::as.formula(M2_FORMULA_TXT), data = d,
                             link = "logit", control = ctrl),
                error = function(e) e)
  if (inherits(f, "error")) {
    out$msg_po <- paste("clm error:", conditionMessage(f))
  } else {
    st <- fit_status(f, "clm (PO)")
    if (nzchar(st)) {
      out$msg_po <- st
    } else {
      r <- tryCatch(std_risk(f, pop, j_chr, AGE_LEVELS, ylev_chr),
                    error = function(e) e)
      if (inherits(r, "error")) {
        out$msg_po <- paste("predict error:", conditionMessage(r))
      } else if (!all(is.finite(r))) {
        out$msg_po <- "non-finite standardised risk"
      } else {
        out$risk_po  <- flat(r)[RISK_KEYS]
        out$logor_po <- logor_po(f)[RD_KEYS]
        out$ok_po    <- TRUE
      }
    }
  }

  out
}

## --- 実行 -------------------------------------------------------------------
t_start <- proc.time()[["elapsed"]]

if (PAR_CORES > 1L && requireNamespace("parallel", quietly = TRUE)) {
  say(sprintf("  running in parallel on %d cores ...", PAR_CORES))
  cl <- parallel::makePSOCKcluster(PAR_CORES)
  boot_res <- tryCatch({
    parallel::clusterEvalQ(cl, { library(ordinal); library(splines); NULL })
    parallel::clusterExport(cl, varlist = c(
      "DAT", "AGE_LEVELS", "OUTCOME", "THRESHOLDS", "THR_KEYS", "M2_FORMULA_TXT",
      "M2_PPO_FORMULA_TXT", "M2_PPO_NOMINAL_F", "PPO_LEVELS_CHR", "J_PPO_CHR",
      "ctrl", "CONTRASTS", "RD_KEYS", "RISK_KEYS", "RDT_KEYS", "STD_POP_IN_BOOT",
      "DIVERGE_COEF", "DIVERGE_SE", "IS_SPLINE_TERM", "pick_j", "conv_code",
      "age_nominal_names", "is_nominal_fit", "cut_names", "theta_jg", "fit_status",
      "make_ppo_levels", "make_ppo_outcome", "flat", "logor_ppo", "logor_po",
      "PREDICT_METHOD", "std_risk", "std_risk_clm", "std_risk_clm_fast",
      "risk_to_rd", "boot_one"),
      envir = environment())
    parallel::parLapply(cl, BOOT_IDX, boot_one)
  }, error = function(e) {
    say("  parallel error: ", conditionMessage(e)); NULL
  })
  try(parallel::stopCluster(cl), silent = TRUE)
  if (is.null(boot_res)) {
    stop("並列実行に失敗した。PAR_CORES <- 1L にして逐次で再実行すること。")
  }
} else {
  if (PAR_CORES > 1L) say("  parallel が使えないため逐次で実行する。")
  say("  running sequentially ...")
  boot_res <- vector("list", NBOOT)
  for (b in seq_len(NBOOT)) {
    boot_res[[b]] <- boot_one(BOOT_IDX[[b]])
    if (b %% PROGRESS_EVERY == 0L || b == NBOOT) {
      el  <- proc.time()[["elapsed"]] - t_start
      eta <- el / b * (NBOOT - b)
      say(sprintf("    %5d / %d  (%5.1f s elapsed, ETA %5.1f s)", b, NBOOT, el, eta))
    }
  }
}
t_boot <- proc.time()[["elapsed"]] - t_start
say(sprintf("  bootstrap finished in %.1f s (%.3f s per resample)",
            t_boot, t_boot / max(NBOOT, 1L)))


# -----------------------------------------------------------------------------
# 7. 区間の算出（パーセンタイル法、§8.4）
# -----------------------------------------------------------------------------

rule("7. Percentile confidence intervals")

ok_ppo <- vapply(boot_res, function(x) isTRUE(x$ok_ppo), logical(1))
ok_po  <- vapply(boot_res, function(x) isTRUE(x$ok_po),  logical(1))
j_used <- do.call(rbind, lapply(boot_res, function(x) x$j_used))

risk_boot_ppo  <- do.call(rbind, lapply(boot_res, function(x) x$risk_ppo))
risk_boot_po   <- do.call(rbind, lapply(boot_res, function(x) x$risk_po))
logor_boot_ppo <- do.call(rbind, lapply(boot_res, function(x) x$logor_ppo))
logor_boot_po  <- do.call(rbind, lapply(boot_res, function(x) x$logor_po))

#' 再標本ごとの群別リスク（列 "G1_t70" ...）から RD（列 "G2_vs_G1_t70" ...）を作る
rd_from_boot <- function(rb) {
  do.call(cbind, lapply(THR_KEYS, function(tk) {
    m <- vapply(CONTRASTS, function(ct)
      rb[, paste(ct$cmp, tk, sep = "_")] - rb[, paste(ct$ref, tk, sep = "_")],
      numeric(nrow(rb)))
    m <- matrix(m, nrow = nrow(rb))
    colnames(m) <- paste(RD_KEYS, tk, sep = "_")
    m
  }))[, RDT_KEYS, drop = FALSE]
}
rd_boot_ppo <- rd_from_boot(risk_boot_ppo)
rd_boot_po  <- rd_from_boot(risk_boot_po)

pct_ci <- function(x, level = CONF_LEVEL) {
  x <- x[is.finite(x)]
  if (!length(x)) return(c(lower = NA_real_, upper = NA_real_))
  a <- (1 - level) / 2
  q <- stats::quantile(x, c(a, 1 - a), names = FALSE, type = CI_TYPE)
  c(lower = q[1], upper = q[2])
}

n_ok_ppo <- sum(ok_ppo); n_ok_po <- sum(ok_po)
say(sprintf("  successful resamples : M2-PPO %d / %d (failed %d, %.2f%%)",
            n_ok_ppo, NBOOT, NBOOT - n_ok_ppo, 100 * (NBOOT - n_ok_ppo) / NBOOT))
say(sprintf("                         M2     %d / %d (failed %d, %.2f%%)",
            n_ok_po, NBOOT, NBOOT - n_ok_po, 100 * (NBOOT - n_ok_po) / NBOOT))

## 失敗の理由の内訳
fail_reason <- function(msgs, ok) {
  m <- msgs[!ok]
  m <- m[nzchar(m)]
  if (!length(m)) return(data.frame(reason = character(0), n = integer(0),
                                    stringsAsFactors = FALSE))
  ## エラー文のうち可変部分を落として集約する
  m <- sub(":.*$", "", m)
  tb <- sort(table(m), decreasing = TRUE)
  data.frame(reason = names(tb), n = as.integer(tb), stringsAsFactors = FALSE)
}
fr_ppo <- fail_reason(vapply(boot_res, function(x) x$msg_ppo, character(1)), ok_ppo)
fr_po  <- fail_reason(vapply(boot_res, function(x) x$msg_po,  character(1)), ok_po)
if (nrow(fr_ppo)) { say("\n  failure reasons (M2-PPO):"); print(fr_ppo, row.names = FALSE) }
if (nrow(fr_po))  { say("\n  failure reasons (M2):");     print(fr_po,  row.names = FALSE) }

## 再標本で用いられた j の分布（§8.4：用いた j をログに出力する。M2 のみ）
say("\n  j used in the resamples (M2):")
for (tk in THR_KEYS) {
  say("  ", tk, ":")
  print(table(j_used[, tk], useNA = "ifany"))
}
check("再標本の j が主解析の j と一致している（M2、閾値ごと）",
      all(vapply(THR_KEYS, function(tk)
        all(is.na(j_used[, tk]) | j_used[, tk] == J_PO[[tk]]), logical(1))),
      paste(vapply(THR_KEYS, function(tk)
        sprintf("%s: main j = %s ; different in %d resamples", tk, J_PO_CHR[[tk]],
                sum(!is.na(j_used[, tk]) & j_used[, tk] != J_PO[[tk]])), character(1)),
        collapse = " / "))

check("M2-PPO の失敗が 5% 以下（§8.4）",
      (NBOOT - n_ok_ppo) / NBOOT <= FAIL_WARN_RATE,
      sprintf("failed %d / %d = %.2f%%。超えた場合は区間の解釈を保留する旨を結果と限界に記す",
              NBOOT - n_ok_ppo, NBOOT, 100 * (NBOOT - n_ok_ppo) / NBOOT))
check("M2 の失敗が 5% 以下（§8.4）",
      (NBOOT - n_ok_po) / NBOOT <= FAIL_WARN_RATE,
      sprintf("failed %d / %d = %.2f%%",
              NBOOT - n_ok_po, NBOOT, 100 * (NBOOT - n_ok_po) / NBOOT))
check("失敗した回を補充していない（§8.4-2）", TRUE,
      "追加の再標本化は行っていない")


# -----------------------------------------------------------------------------
# 8. 出力表の組み立て
# -----------------------------------------------------------------------------

rule("8. Assembling the output tables")

## --- (E-4a) 年齢群別の標準化リスク（M2-PPO、閾値 × 年齢群）------------------
std_risk_tab <- do.call(rbind, lapply(seq_along(THRESHOLDS), function(i) {
  t <- THRESHOLDS[i]; tk <- THR_KEYS[i]
  do.call(rbind, lapply(AGE_LEVELS, function(g) {
    key <- paste(g, tk, sep = "_")
    ci  <- pct_ci(risk_boot_ppo[ok_ppo, key])
    data.frame(
      threshold       = t,
      age_group       = g,
      age_group_label = relabel_levels("age_group", g),
      n               = as.integer(n_by_g[[g]]),
      estimand        = sprintf("Standardised P(mFIM at discharge >= %d)", t),
      estimate        = risk_point[g, tk],
      ci_lower        = unname(ci[["lower"]]),
      ci_upper        = unname(ci[["upper"]]),
      percent         = 100 * risk_point[g, tk],
      percent_lower   = 100 * unname(ci[["lower"]]),
      percent_upper   = 100 * unname(ci[["upper"]]),
      conf_level      = CONF_LEVEL,
      ci_method       = "bootstrap percentile, unadjusted",
      n_boot_used     = n_ok_ppo,
      model           = MODEL_PPO,
      crude_proportion = unname(crude[g, tk]),
      stringsAsFactors = FALSE
    )
  }))
}))

## --- (E-4b) 3対比のリスク差（M2-PPO、閾値 × 対比）---------------------------
rd_tab <- do.call(rbind, lapply(seq_along(THRESHOLDS), function(i) {
  t <- THRESHOLDS[i]; tk <- THR_KEYS[i]
  do.call(rbind, lapply(CONTRASTS, function(ct) {
    k  <- ct$key
    ci <- pct_ci(rd_boot_ppo[ok_ppo, paste(k, tk, sep = "_")])
    data.frame(
      threshold        = t,
      contrast         = k,
      label            = contrast_label(ct$ref, ct$cmp),
      group_reference  = ct$ref,
      group_comparison = ct$cmp,
      reference_label  = relabel_levels("age_group", ct$ref),
      comparison_label = relabel_levels("age_group", ct$cmp),
      role             = ct$role,
      estimand         = sprintf("Risk difference in standardised P(mFIM >= %d), older minus younger",
                                 t),
      estimate         = rd_point[k, tk],
      ci_lower         = unname(ci[["lower"]]),
      ci_upper         = unname(ci[["upper"]]),
      percentage_point       = 100 * rd_point[k, tk],
      percentage_point_lower = 100 * unname(ci[["lower"]]),
      percentage_point_upper = 100 * unname(ci[["upper"]]),
      conf_level       = CONF_LEVEL,
      ci_method        = "bootstrap percentile, unadjusted",
      n_boot_used      = n_ok_ppo,
      model            = MODEL_PPO,
      p_value          = NA_real_,    # §8.5-1：Table 2 に p 値は載せない
      stringsAsFactors = FALSE
    )
  }))
}))

## --- (E-3) M2-PPO の閾値別オッズ比：profile 区間（05）とブートストラップ区間 ---
ppo_or_tab <- m2_ppo_odds_ratio
ppo_or_tab$boot_ci_lower <- NA_real_
ppo_or_tab$boot_ci_upper <- NA_real_
for (i in seq_len(nrow(ppo_or_tab))) {
  key <- paste(ppo_or_tab$contrast[i], paste0("t", ppo_or_tab$threshold[i]), sep = "_")
  ci  <- pct_ci(logor_boot_ppo[ok_ppo, key])
  ppo_or_tab$boot_ci_lower[i] <- exp(unname(ci[["lower"]]))
  ppo_or_tab$boot_ci_upper[i] <- exp(unname(ci[["upper"]]))
}
ppo_or_tab$boot_ci_method  <- "bootstrap percentile, unadjusted"
ppo_or_tab$n_boot_used     <- n_ok_ppo
ppo_or_tab$model           <- MODEL_PPO
ppo_or_tab <- ppo_or_tab[order(ppo_or_tab$threshold,
                               match(ppo_or_tab$contrast, RD_KEYS)), , drop = FALSE]
rownames(ppo_or_tab) <- NULL

## --- profile 区間とブートストラップ区間の突き合わせ（著者の判断、2026-09-22）---
## 対数オッズ比の尺度で、端のずれの大きい方を profile 区間の幅で割る。
## M2 の共通オッズ比（ordinal の profile 区間）も同じ再標本で求め、参考として並べる。
rel_shift <- function(pl, pu, bl, bu) {
  pmax(abs(log(bl) - log(pl)), abs(log(bu) - log(pu))) / (log(pu) - log(pl))
}
.po_boot <- t(vapply(m2_odds_ratio$contrast, function(k) {
  ci <- pct_ci(logor_boot_po[ok_po, k]); exp(c(ci[["lower"]], ci[["upper"]]))
}, numeric(2)))
or_ci_compare <- rbind(
  data.frame(model = MODEL_PPO, threshold = ppo_or_tab$threshold,
             contrast = ppo_or_tab$contrast, role = ppo_or_tab$role,
             odds_ratio = ppo_or_tab$odds_ratio,
             profile_lower = ppo_or_tab$ci_lower, profile_upper = ppo_or_tab$ci_upper,
             profile_source = "own implementation (05)",
             boot_lower = ppo_or_tab$boot_ci_lower, boot_upper = ppo_or_tab$boot_ci_upper,
             n_boot_used = n_ok_ppo, used_for_check = TRUE,
             stringsAsFactors = FALSE),
  data.frame(model = MODEL_PO, threshold = NA_real_,
             contrast = m2_odds_ratio$contrast, role = m2_odds_ratio$role,
             odds_ratio = m2_odds_ratio$odds_ratio,
             profile_lower = m2_odds_ratio$ci_lower, profile_upper = m2_odds_ratio$ci_upper,
             profile_source = "ordinal::confint (05)",
             boot_lower = unname(.po_boot[, 1]), boot_upper = unname(.po_boot[, 2]),
             n_boot_used = n_ok_po, used_for_check = FALSE,
             stringsAsFactors = FALSE)
)
or_ci_compare$ratio_lower <- or_ci_compare$boot_lower / or_ci_compare$profile_lower
or_ci_compare$ratio_upper <- or_ci_compare$boot_upper / or_ci_compare$profile_upper
or_ci_compare$relative_shift <- rel_shift(or_ci_compare$profile_lower,
                                          or_ci_compare$profile_upper,
                                          or_ci_compare$boot_lower,
                                          or_ci_compare$boot_upper)
or_ci_compare$within_tolerance <- or_ci_compare$relative_shift <= BOOT_PROFILE_TOL
say("  odds ratios: profile likelihood vs bootstrap percentile intervals")
print(or_ci_compare[, c("model", "threshold", "contrast", "odds_ratio", "profile_lower",
                        "boot_lower", "profile_upper", "boot_upper", "relative_shift")],
      row.names = FALSE, digits = 4)
.cmp_ppo <- or_ci_compare[or_ci_compare$used_for_check, , drop = FALSE]
.cmp_po  <- or_ci_compare[!or_ci_compare$used_for_check, , drop = FALSE]
check(sprintf(paste0("M2-PPO の閾値別オッズ比で、profile 区間とブートストラップ区間に",
                     "大きな差がない（端のずれ ≤ profile 区間の幅の %.0f%%、対数尺度）"),
              100 * BOOT_PROFILE_TOL),
      all(is.finite(.cmp_ppo$relative_shift)) && all(.cmp_ppo$within_tolerance),
      sprintf("max = %.1f%%（%s）", 100 * max(.cmp_ppo$relative_shift),
              with(.cmp_ppo[which.max(.cmp_ppo$relative_shift), ],
                   sprintf("%s, y >= %g", contrast, threshold))))
check("（参考）M2 の共通オッズ比での同じ比較（ordinal の profile 区間を基準）", TRUE,
      sprintf("max = %.1f%%。判定には使わない。差の大きさの目安として並べる",
              100 * max(.cmp_po$relative_shift)))

## --- 比例オッズモデル（M2）由来の標準化量（比較用。Table 2 には載せない）-----
po_tab <- rbind(
  do.call(rbind, lapply(seq_along(THRESHOLDS), function(i) {
    t <- THRESHOLDS[i]; tk <- THR_KEYS[i]
    do.call(rbind, lapply(AGE_LEVELS, function(g) {
      ci <- pct_ci(risk_boot_po[ok_po, paste(g, tk, sep = "_")])
      data.frame(quantity = "standardised risk", threshold = t, key = g,
                 label = relabel_levels("age_group", g),
                 estimate = risk_point_po[g, tk],
                 ci_lower = unname(ci[["lower"]]), ci_upper = unname(ci[["upper"]]),
                 j_used = J_PO_CHR[[tk]],
                 stringsAsFactors = FALSE)
    }))
  })),
  do.call(rbind, lapply(seq_along(THRESHOLDS), function(i) {
    t <- THRESHOLDS[i]; tk <- THR_KEYS[i]
    do.call(rbind, lapply(CONTRASTS, function(ct) {
      ci <- pct_ci(rd_boot_po[ok_po, paste(ct$key, tk, sep = "_")])
      data.frame(quantity = "risk difference", threshold = t, key = ct$key,
                 label = contrast_label(ct$ref, ct$cmp),
                 estimate = rd_point_po[ct$key, tk],
                 ci_lower = unname(ci[["lower"]]), ci_upper = unname(ci[["upper"]]),
                 j_used = J_PO_CHR[[tk]],
                 stringsAsFactors = FALSE)
    }))
  }))
)
po_tab$conf_level  <- CONF_LEVEL
po_tab$ci_method   <- "bootstrap percentile, unadjusted"
po_tab$n_boot_used <- n_ok_po
po_tab$model       <- MODEL_PO

## --- M2-PPO と M2 の食い違い（§8.6「事前の規則」に準じて数値で示す）----------
diff_tab <- do.call(rbind, lapply(seq_along(THRESHOLDS), function(i) {
  tk <- THR_KEYS[i]
  data.frame(
    threshold = THRESHOLDS[i],
    key       = c(AGE_LEVELS, RD_KEYS),
    quantity  = c(rep("standardised risk", length(AGE_LEVELS)),
                  rep("risk difference", length(RD_KEYS))),
    ppo_model = c(unname(risk_point[AGE_LEVELS, tk]), unname(rd_point[RD_KEYS, tk])),
    po_model  = c(unname(risk_point_po[AGE_LEVELS, tk]), unname(rd_point_po[RD_KEYS, tk])),
    stringsAsFactors = FALSE
  )
}))
diff_tab$difference <- diff_tab$po_model - diff_tab$ppo_model
say("\n  M2-PPO vs M2 (proportional odds) standardised quantities:")
print(diff_tab, row.names = FALSE, digits = 4)
say("  Table 2 に載せるのは M2-PPO 由来の値である（plan.summary §7-2）。")

## --- Table 2（閾値別OR・標準化確率・RD を M2-PPO から、共通OR を M2 から）------
## 共通オッズ比・閾値別オッズ比は 05 の出力をそのまま使う（profile likelihood 95%区間）。
if (!exists("m2_odds_ratio", inherits = TRUE)) {
  stop("m2_odds_ratio が見つからない。先に 05_m2_clm.R を実行すること。")
}

t2_risk <- data.frame(
  block      = "Standardised risk",
  threshold  = std_risk_tab$threshold,
  row_key    = std_risk_tab$age_group,
  row_label  = std_risk_tab$age_group_label,
  role       = "",
  n          = std_risk_tab$n,
  quantity   = std_risk_tab$estimand,
  estimate   = std_risk_tab$estimate,
  ci_lower   = std_risk_tab$ci_lower,
  ci_upper   = std_risk_tab$ci_upper,
  scale      = "probability",
  ci_method  = std_risk_tab$ci_method,
  model      = MODEL_PPO,
  stringsAsFactors = FALSE
)
t2_rd <- data.frame(
  block      = "Risk difference",
  threshold  = rd_tab$threshold,
  row_key    = rd_tab$contrast,
  row_label  = rd_tab$label,
  role       = rd_tab$role,
  n          = NA_integer_,
  quantity   = sprintf("Difference in standardised P(mFIM >= %d)", rd_tab$threshold),
  estimate   = rd_tab$estimate,
  ci_lower   = rd_tab$ci_lower,
  ci_upper   = rd_tab$ci_upper,
  scale      = "probability difference",
  ci_method  = rd_tab$ci_method,
  model      = MODEL_PPO,
  stringsAsFactors = FALSE
)
t2_or_ppo <- data.frame(
  block      = "Odds ratio",
  threshold  = ppo_or_tab$threshold,
  row_key    = ppo_or_tab$contrast,
  row_label  = ppo_or_tab$label,
  role       = ppo_or_tab$role,
  n          = NA_integer_,
  quantity   = sprintf("Odds ratio for P(Y >= %d)", ppo_or_tab$threshold),
  estimate   = ppo_or_tab$odds_ratio,
  ci_lower   = ppo_or_tab$ci_lower,
  ci_upper   = ppo_or_tab$ci_upper,
  scale      = "odds ratio",
  ci_method  = ppo_or_tab$ci_method,
  model      = MODEL_PPO,
  stringsAsFactors = FALSE
)
t2_or <- data.frame(
  block      = "Common odds ratio",
  threshold  = NA_real_,
  row_key    = m2_odds_ratio$contrast,
  row_label  = m2_odds_ratio$label,
  role       = m2_odds_ratio$role,
  n          = NA_integer_,
  quantity   = "Common odds ratio for P(Y >= y)",
  estimate   = m2_odds_ratio$odds_ratio,
  ci_lower   = m2_odds_ratio$ci_lower,
  ci_upper   = m2_odds_ratio$ci_upper,
  scale      = "odds ratio",
  ci_method  = m2_odds_ratio$ci_method,
  model      = MODEL_PO,
  stringsAsFactors = FALSE
)
table2 <- rbind(t2_risk, t2_rd, t2_or_ppo, t2_or)
table2$conf_level   <- CONF_LEVEL
table2$analysis_set <- "M2 complete-case population"
table2$analysis_n   <- N_M2
table2$p_value      <- NA_real_     # §8.5-1：Table 2 にはp値を載せない

say("\n  Table 2 :")
print(table2[, c("block", "threshold", "row_label", "role", "estimate", "ci_lower", "ci_upper")],
      row.names = FALSE, digits = 4)


# -----------------------------------------------------------------------------
# 9. 諸元と脚注
# -----------------------------------------------------------------------------

rule("9. Table 2 metadata and footnotes")

## 複数の値をもつ項目（threshold など）は ";" で区切った1つの値にする。
## 項目を1行に保つのは、12 の meta_get() が同名の項目の先頭しか返さないため。
table2_meta <- data.frame(
  item = c("analysis_set", "analysis_n", paste0("n_", AGE_LEVELS),
           "outcome", "outcome_label", "threshold", "ppo_categories",
           "j_used_ppo", "j_used_po",
           "estimand_risk", "estimand_rd", "estimand_or", "estimand_or_common",
           "conf_level", "ci_method_risk_rd", "ci_method_or", "ci_method_or_common",
           "ci_method_or_check",
           "nboot", "nboot_used_ppo", "nboot_used_po",
           "boot_seed", "boot_stratified_by", "standardisation_population",
           "standardisation_population_in_bootstrap",
           "primary_contrast", "model", "model_nominal", "model_common_or",
           "p_values", "risk_ratio",
           "predict_method", "boot_profile_tolerance",
           "seconds_bootstrap", "seconds_per_resample"),
  value = c("M2 complete-case population", as.character(N_M2),
            as.character(as.integer(n_by_g[AGE_LEVELS])),
            OUTCOME, relabel_vars(OUTCOME),
            paste(THRESHOLDS, collapse = ";"), paste(PPO_LEVELS_CHR, collapse = ";"),
            paste(J_PPO_CHR, collapse = ";"), paste(J_PO_CHR, collapse = ";"),
            sprintf("Standardised P(%s >= t) by age group, t = %s",
                    relabel_vars(OUTCOME), paste(THRESHOLDS, collapse = ", ")),
            "Risk difference, older minus younger",
            sprintf("Odds ratio for P(Y >= t), age group, t = %s (partial proportional odds)",
                    paste(THRESHOLDS, collapse = ", ")),
            "Common odds ratio for P(Y >= y), age group (proportional odds, for comparison)",
            as.character(CONF_LEVEL),
            "bootstrap percentile, unadjusted",
            unique(m2_ppo_odds_ratio$ci_method)[1],
            unique(m2_odds_ratio$ci_method)[1],
            "bootstrap percentile (same resamples), compared with the profile likelihood interval",
            as.character(NBOOT), as.character(n_ok_ppo), as.character(n_ok_po),
            as.character(BOOT_SEED), "age group",
            "M2 complete-case population (all covariates held at observed values)",
            STD_POP_IN_BOOT,
            PRIMARY_KEY,
            M2_PPO_FORMULA_TXT, M2_PPO_NOMINAL_TXT, M2_FORMULA_TXT,
            "not reported (plan 8.5)", "not computed (plan 8.4)",
            PREDICT_METHOD, as.character(BOOT_PROFILE_TOL),
            sprintf("%.1f", t_boot),
            sprintf("%.3f", t_boot / max(NBOOT, 1L))),
  stringsAsFactors = FALSE
)
print(table2_meta, row.names = FALSE)

table2_footnotes <- data.frame(
  tag = c("i", "ii", "iii", "iv", "v"),
  text = c(
    paste0("The primary contrast is ", contrast_label("G1", "G3"),
           "; the other two contrasts are secondary."),
    paste("The models are conditional on sex, disease category, pre-admission care need,",
          "motor and cognitive FIM at admission (ridit, restricted cubic spline)",
          "and admission period; estimates are not to be compared in magnitude",
          "with estimates made under a different conditioning set."),
    paste0("Odds ratios at ", paste(THRESHOLDS, collapse = " and "), " points and ",
           "standardised probabilities are from a partial proportional odds model in ",
           "which the effect of age group was allowed to differ between the thresholds ",
           "(outcome grouped as ", paste(PPO_LEVELS_CHR, collapse = ", "), " points) and ",
           "covariate effects were common to both thresholds. Odds ratios have profile ",
           "likelihood 95% intervals."),
    paste0("Standardised risks and risk differences are averaged over the M2 ",
           "complete-case population (N = ", N_M2, ") with every patient assigned ",
           "in turn to each age group; 95% intervals are bootstrap percentile ",
           "intervals from ", NBOOT, " resamples stratified by age group ",
           "(", n_ok_ppo, " usable)."),
    paste0("The common odds ratio is from the cumulative logit model with proportional ",
           "odds for all variables (motor FIM at discharge in ",
           length(Y_LEVELS_CHR), " levels), shown for comparison, with profile ",
           "likelihood 95% intervals. All intervals are unadjusted for multiplicity ",
           "and are to be read as compatibility intervals.")
  ),
  stringsAsFactors = FALSE
)
say("")
for (i in seq_len(nrow(table2_footnotes))) {
  say(sprintf("  (%s) %s", table2_footnotes$tag[i], table2_footnotes$text[i]))
}

boot_summary <- data.frame(
  model = c(MODEL_PPO, MODEL_PO),
  n_boot = c(NBOOT, NBOOT),
  n_used = c(n_ok_ppo, n_ok_po),
  n_failed = c(NBOOT - n_ok_ppo, NBOOT - n_ok_po),
  percent_failed = c(100 * (NBOOT - n_ok_ppo) / NBOOT,
                     100 * (NBOOT - n_ok_po) / NBOOT),
  exceeds_5_percent = c((NBOOT - n_ok_ppo) / NBOOT > FAIL_WARN_RATE,
                        (NBOOT - n_ok_po) / NBOOT > FAIL_WARN_RATE),
  seed = c(BOOT_SEED, BOOT_SEED),
  stratified_by = c("age group", "age group"),
  seconds = c(t_boot, NA_real_),   # 2つのモデルを同じ回で当てはめた合計
  stringsAsFactors = FALSE
)


# -----------------------------------------------------------------------------
# 10. 出来上がりの点検
# -----------------------------------------------------------------------------

rule("10. Verification")

.n_t2 <- 3L * length(THRESHOLDS) * 3L + length(CONTRASTS)
check(sprintf("Table 2 の行数が 標準化リスク 6 + RD 6 + 閾値別OR 6 + 共通OR 3 = %d である", .n_t2),
      nrow(table2) == .n_t2, sprintf("rows = %d", nrow(table2)))
check("Table 2 に p 値を置いていない（§8.5-1）", all(is.na(table2$p_value)), "")
check("Table 2 にリスク比（RR）を置いていない（§8.4-3）",
      !any(grepl("risk ratio", c(table2$quantity, table2$scale), ignore.case = TRUE)),
      "RR は算出しない。群別標準化リスクから読者が計算できる")
check("Table 2 に閾値65の値を置いていない（§7-1）",
      !any(grepl("65", table2$quantity)) && all(table2$threshold %in% c(THRESHOLDS, NA)), "")
check("すべての区間が点推定値を挟む",
      all(is.na(table2$ci_lower) |
          (table2$ci_lower <= table2$estimate + 1e-12 &
           table2$estimate <= table2$ci_upper + 1e-12)),
      paste(table2$row_key[!(is.na(table2$ci_lower) |
        (table2$ci_lower <= table2$estimate + 1e-12 &
         table2$estimate <= table2$ci_upper + 1e-12))], collapse = ", "))
check("M2-PPO のオッズ比のブートストラップ区間も点推定値を挟む",
      all(ppo_or_tab$boot_ci_lower <= ppo_or_tab$odds_ratio + 1e-12 &
            ppo_or_tab$odds_ratio <= ppo_or_tab$boot_ci_upper + 1e-12), "")
check("標準化リスクの区間が (0, 1) に収まる",
      all(is.na(std_risk_tab$ci_lower) |
          (std_risk_tab$ci_lower > 0 & std_risk_tab$ci_upper < 1)), "")
check("RD の区間が (−1, 1) に収まる",
      all(is.na(rd_tab$ci_lower) | (rd_tab$ci_lower > -1 & rd_tab$ci_upper < 1)), "")
check("主対比が G3 対 G1 である（両閾値）",
      all(vapply(THRESHOLDS, function(t)
        any(rd_tab$threshold == t & rd_tab$contrast == PRIMARY_KEY & rd_tab$role == "primary"),
        logical(1))), "")
check("Table 2 の閾値別OR が 05 の m2_ppo_odds_ratio（profile 区間）そのものである",
      {
        .a <- table2[table2$block == "Odds ratio", , drop = FALSE]
        .b <- m2_ppo_odds_ratio[match(paste(.a$row_key, .a$threshold),
                                      paste(m2_ppo_odds_ratio$contrast,
                                            m2_ppo_odds_ratio$threshold)), , drop = FALSE]
        nrow(.a) == nrow(m2_ppo_odds_ratio) &&
          isTRUE(all.equal(.a$estimate, .b$odds_ratio)) &&
          isTRUE(all.equal(.a$ci_lower, .b$ci_lower)) &&
          isTRUE(all.equal(.a$ci_upper, .b$ci_upper))
      }, "")
check("解析集団のNを Table 2 の諸元に載せている（§6.5「報告の規則」）",
      all(table2$analysis_n == N_M2), sprintf("N = %d", N_M2))
check("E-2（Fig. 2）と E-4（Table 2）で N が異なりうることを記録した", TRUE,
      sprintf("M2 complete-case N = %d。E-1・E-2 は解析対象の全例を使う（§6.5）", N_M2))
check("ブートストラップの再標本が年齢群で層化されている（§8.4）", TRUE,
      paste(sprintf("%s=%d", AGE_LEVELS, vapply(idx_by_g, length, integer(1))),
            collapse = ", "))

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
w(table2,           "table2.csv")
w(std_risk_tab,     "table2_standardized_risk.csv")
w(rd_tab,           "table2_risk_difference.csv")
w(ppo_or_tab,       "table2_ppo_odds_ratio.csv")
w(po_tab,           "table_po_standardized.csv")
w(diff_tab,         "table_ppo_vs_po.csv")
w(or_ci_compare,    "table_or_ci_profile_vs_boot.csv")
w(boot_summary,     "table_boot_summary.csv")
w(table2_meta,      "table2_meta.csv")
w(table2_footnotes, "table2_footnotes.csv")
w(checks,           "table_checks_06.csv")

## risk_point / rd_point は M2-PPO（Table 2）、*_po は M2（比較用）。
## いずれも 行 = 年齢群または対比、列 = 閾値（"t70"、"t77"）の行列。
## 再標本の分布は列名 "G1_t70"、"G2_vs_G1_t70" などの行列。
save(table2, std_risk_tab, rd_tab, ppo_or_tab, po_tab, diff_tab, or_ci_compare,
     boot_summary, table2_meta, table2_footnotes,
     risk_point, rd_point, risk_point_po, rd_point_po, logor_point_ppo, logor_point_po,
     risk_boot_ppo, rd_boot_ppo, logor_boot_ppo,
     risk_boot_po, rd_boot_po, logor_boot_po,
     ok_ppo, ok_po, j_used, J_PO, J_PO_CHR, J_PPO_CHR, THRESHOLDS, THR_KEYS,
     NBOOT, BOOT_SEED, fr_ppo, fr_po, checks,
     file = file.path(OUT_DIR, "06_standardize_boot.rda"))
say("  written: ", file.path(OUT_DIR, "06_standardize_boot.rda"))

say("\n--- sessionInfo ---")
print(utils::sessionInfo())

say("\nDone. M2 N = ", N_M2, " ; bootstrap M2-PPO ", n_ok_ppo, " / ", NBOOT,
    " usable, M2 ", n_ok_po, " / ", NBOOT, " usable ; thresholds = ",
    paste(THRESHOLDS, collapse = ", "), " ; M2 j = ", paste(J_PO_CHR, collapse = ", "))

.log_close()
