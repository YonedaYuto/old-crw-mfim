# =============================================================================
# 10_age_continuous.R
#   解析計画書 §8.9「補足解析：連続年齢のモデル（層C）」の実装。
#   補足図2 を描くための数値をすべて出力する（作図は別スクリプト）。
#
#   入力 : output/04_ridit_spline.rda  dat_m2 / main / spline_terms / M2_SETTINGS
#          output/05_m2_clm.rda        M2_FORMULA_TXT（共変量の構成の照合に使う）
#          01_labels.R（表示用の文字列）
#   出力 : output/table_supp_fig2_age_or.csv     補足図2 の本体（profile 95%区間）
#          output/table_supp_fig2_curve_wald.csv 細かい格子の点推定値と Wald 区間
#          output/table_supp_fig2_knots.csv      年齢スプラインのノット位置
#          output/table_supp_fig2_fit.csv        当てはめの諸元と収束診断
#          output/table_supp_fig2_meta.csv       補足図2 の脚注に必要な諸元
#          output/table_supp_fig2_footnotes.csv  補足図2 の脚注（下書き）
#          output/table_checks_10.csv            点検結果
#          output/10_age_continuous.rda
#          output/log_10_age_continuous.txt      実行ログ
#
#   計画との対応（§8.9）
#     ・曝露 : 年齢を連続量とし、制限付き3次スプライン（ノット4個、**解析対象**に
#       おける年齢の第5・35・65・95百分位。`ns()` への渡し方は §6.4 と同じ）で表す
#     ・モデル : 累積ロジット（最尤）。共変量の構成は M2 と同じ1本に限る
#     ・報告 : 70歳を参照とした年齢の共通オッズ比を、年齢を横軸にした曲線
#       （95%信頼区間つき）で補足図2 に示す
#     ・参照点を70歳とするのは、65歳が年齢の最小値であり、分布の境界を参照点に
#       すると推定が不安定になるためである
#     ・層Cであり、本文の結論は3区分にもとづく（§7.5）
#
#   区間の構成法（著者の決定、2026-09-20）
#     §8.9 は「95%信頼区間つき」とだけ記し、方法を明記していない。**profile
#     likelihood** を用いる。§8.3（年齢群の共通オッズ比）と構成法が揃う。
#
#     profile 区間は1つの母数についてしか取れないため、年齢 a ごとに
#     「log OR(a vs 70) そのものが第1係数になる」向きへモデルを再母数化して
#     当てはめ直す。
#       c = B(a) − B(70)                       （B は年齢スプラインの基底、長さ K）
#       M = rbind(c', Q')                      （Q は c に直交する K−1 本。M は正則）
#       A = M^(-1),  Z = B A                   （Z は B と同じ列空間を張る）
#     Z を説明変数に使うと、係数は γ = M β となり γ₁ = c'β = log OR(a vs 70)。
#     当てはめは元のモデルと同値（対数尤度が一致することを毎回点検する）。
#     Wald（デルタ法）区間も細かい格子で算出し、比較としてログと別表に残す。
#     本文・補足図2 に載せるのは profile 区間である。
#
#   文字コード : UTF-8（Windows の R では source(..., encoding = "UTF-8")）
# =============================================================================


# -----------------------------------------------------------------------------
# 0. 設定
# -----------------------------------------------------------------------------

OUT_DIR   <- "output"
IN_RDA_04 <- file.path(OUT_DIR, "04_ridit_spline.rda")
IN_RDA_05 <- file.path(OUT_DIR, "05_m2_clm.rda")

CONF_LEVEL <- 0.95

## --- 年齢スプライン（§8.9、§6.4 と同じ渡し方）-------------------------------
AGE_VAR        <- "age"
AGE_REF        <- 70            # 参照点（§8.9）
KNOT_PROBS_4   <- c(0.05, 0.35, 0.65, 0.95)   # 4ノット
KNOT_PROBS_3   <- c(0.10, 0.50, 0.90)         # 退化した場合の代替（§6.4-4）
KNOT_MIN_GAP   <- 1e-8
QUANTILE_TYPE  <- 1L            # 04 の設定があればそちらに合わせる
## ノットを決める母集団。§8.9 は「解析対象における年齢の…百分位」と記す。
##   "analysis_set"  : 解析対象（main。E1〜E3 を除いた集団）  ← 計画どおり
##   "complete_case" : M2 の完全ケース集団（当てはめに使う集団）
KNOT_POPULATION <- "analysis_set"

## --- 格子（§8.9：年齢を横軸にした曲線）--------------------------------------
## profile 区間は年齢の点ごとにモデルを当てはめ直したうえで尤度を profile する。
## 1点あたり本体の当てはめの数十回ぶんに相当する（実測：本体 0.4 秒のとき1点 30 秒
## 前後）。65〜100歳を1歳刻みにすると 20〜40 分かかる見込みである。時間が問題に
## なる場合は PROFILE_BY を 2 や 5 にして格子を粗くする。点推定値の曲線は
## CURVE_BY の細かい格子のまま別表に出力されるため、図の線は粗くならない。
PROFILE_BY   <- 1               # profile 区間を算出する年齢の間隔（歳）
CURVE_BY     <- 0.5             # 点推定値と Wald 区間だけを算出する細かい格子
AGE_MAX_GRID <- NA_real_        # 上限を切る場合に指定（NA なら観測された最大年齢）
## 代表年齢。ログに数値で出す（本文に1行書くときの材料）
REPORT_AGES  <- c(65, 75, 80, 85, 90, 95)

## --- 収束と発散の判定（05・06 と同じ）---------------------------------------
CLM_MAX_ITER     <- 200L
CLM_MAX_LINE     <- 50L
CLM_MAX_MOD_ITER <- 10L
DIVERGE_COEF   <- 10
DIVERGE_SE     <- 10
IS_SPLINE_TERM <- function(nm) grepl("^ns\\(|^az[0-9]+$", nm)


# -----------------------------------------------------------------------------
# 1. 表示用ラベル
# -----------------------------------------------------------------------------

if (!file.exists("01_labels.R")) {
  stop("01_labels.R が見つからない。10_age_continuous.R は 01_labels.R と同じ",
       "フォルダを作業ディレクトリにして実行すること。")
}
source("01_labels.R", encoding = "UTF-8")


# -----------------------------------------------------------------------------
# 2. 実行の準備（出力先とログ）
# -----------------------------------------------------------------------------

if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

.log_con <- file(file.path(OUT_DIR, "log_10_age_continuous.txt"), open = "wt",
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
  fml <- setdiff(names(formals(fun)), "...")
  do.call(fun, args[names(args) %in% fml])
}

rule(paste0("10_age_continuous.R  ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
say("R : ", R.version.string)
say("reference age : ", AGE_REF, " ; CI : profile likelihood (", CONF_LEVEL, ")")

has_ordinal <- requireNamespace("ordinal", quietly = TRUE)
check("ordinal が利用できる", has_ordinal, "")
if (!has_ordinal) stop("ordinal が無ければ補足解析は当てはめられない。")
library(ordinal)
library(splines)


# -----------------------------------------------------------------------------
# 3. 入力の取り込み
# -----------------------------------------------------------------------------

rule("3. Input")

load_env <- function(path, required = TRUE) {
  if (!file.exists(path)) {
    if (required) stop(path, " が見つからない。先に該当のスクリプトを実行すること。")
    say("注意：", path, " が無い。")
    return(NULL)
  }
  e <- new.env(parent = emptyenv())
  load(path, envir = e)
  say("loaded : ", path)
  e
}
take <- function(e, nm, path) {
  miss <- nm[!vapply(nm, exists, logical(1), envir = e, inherits = FALSE)]
  if (length(miss)) stop(path, " に必要なオブジェクトがない: ",
                         paste(miss, collapse = ", "))
  for (n in nm) assign(n, get(n, envir = e), envir = globalenv())
  invisible(NULL)
}

e04 <- load_env(IN_RDA_04)
take(e04, c("dat_m2", "main", "spline_terms", "M2_SETTINGS"), IN_RDA_04)
e05 <- load_env(IN_RDA_05)
take(e05, c("M2_FORMULA_TXT"), IN_RDA_05)

OUTCOME    <- M2_SETTINGS$outcome
AGE_LEVELS <- M2_SETTINGS$age_levels
COVARS_CAT <- M2_SETTINGS$covars_cat
RIDIT_VARS <- M2_SETTINGS$ridit_vars
if (!is.null(M2_SETTINGS$quantile_type)) QUANTILE_TYPE <- M2_SETTINGS$quantile_type

DAT <- as.data.frame(dat_m2, stringsAsFactors = FALSE)
DAT$age_group <- factor(as.character(DAT$age_group), levels = AGE_LEVELS)
MAIN <- as.data.frame(main, stringsAsFactors = FALSE)
N_M2 <- nrow(DAT)

say("dat_m2 (fitted population, M2 complete cases) : ", N_M2, " rows")
say("main   (analysis set; knots are taken here)   : ", nrow(MAIN), " rows")
say(sprintf("age in dat_m2 : min %g, max %g, median %g",
            min(DAT[[AGE_VAR]]), max(DAT[[AGE_VAR]]),
            stats::median(DAT[[AGE_VAR]])))

check("年齢に欠測がない（当てはめ集団）", !anyNA(DAT[[AGE_VAR]]),
      sprintf("missing = %d", sum(is.na(DAT[[AGE_VAR]]))))
check("参照年齢が観測範囲の内側にある",
      AGE_REF > min(DAT[[AGE_VAR]]) && AGE_REF < max(DAT[[AGE_VAR]]),
      sprintf("ref = %g, range = [%g, %g]", AGE_REF,
              min(DAT[[AGE_VAR]]), max(DAT[[AGE_VAR]])))


# -----------------------------------------------------------------------------
# 4. 年齢スプラインのノット（§8.9、§6.4 と同じ規則）
# -----------------------------------------------------------------------------

rule("4. Knots of the age spline")

knot_source <- if (identical(KNOT_POPULATION, "analysis_set")) MAIN else DAT
knot_label  <- if (identical(KNOT_POPULATION, "analysis_set"))
  "analysis set (E1-E3 excluded)" else "M2 complete-case population"
say("  knot population : ", knot_label, " (n = ", nrow(knot_source), ")")

knots_ok <- function(k) all(is.finite(k)) && all(diff(k) > KNOT_MIN_GAP)

x_knot <- knot_source[[AGE_VAR]]
x_knot <- x_knot[!is.na(x_knot)]
k4 <- stats::quantile(x_knot, KNOT_PROBS_4, names = FALSE, type = QUANTILE_TYPE)
if (knots_ok(k4)) {
  AGE_KNOTS <- k4; AGE_KNOT_PROBS <- KNOT_PROBS_4; AGE_NKNOT <- 4L; AGE_FALLBACK <- FALSE
} else {
  k3 <- stats::quantile(x_knot, KNOT_PROBS_3, names = FALSE, type = QUANTILE_TYPE)
  if (!knots_ok(k3)) {
    stop("年齢のノットが重なる。4ノットでも3ノットでも狭義単調にならない。")
  }
  AGE_KNOTS <- k3; AGE_KNOT_PROBS <- KNOT_PROBS_3; AGE_NKNOT <- 3L; AGE_FALLBACK <- TRUE
}
AGE_INNER <- AGE_KNOTS[-c(1, length(AGE_KNOTS))]
AGE_BOUND <- AGE_KNOTS[c(1, length(AGE_KNOTS))]

say(sprintf("  %d knots at percentiles %s", AGE_NKNOT,
            paste(100 * AGE_KNOT_PROBS, collapse = "/")))
say(sprintf("    knots          = %s", paste(sprintf("%.6g", AGE_KNOTS), collapse = ", ")))
say(sprintf("    inner (knots)  = %s", paste(sprintf("%.6g", AGE_INNER), collapse = ", ")))
say(sprintf("    Boundary.knots = %s", paste(sprintf("%.6g", AGE_BOUND), collapse = ", ")))
if (AGE_FALLBACK) say("    （4ノットが重なったため §6.4-4 により3ノットに減らした）")

## 参考：もう一方の母集団で取ったノット（差が大きければログで気づける）
x_alt <- (if (identical(KNOT_POPULATION, "analysis_set")) DAT else MAIN)[[AGE_VAR]]
k_alt <- stats::quantile(x_alt[!is.na(x_alt)], AGE_KNOT_PROBS, names = FALSE,
                         type = QUANTILE_TYPE)
say(sprintf("  (reference) the same percentiles in the other population : %s",
            paste(sprintf("%.6g", k_alt), collapse = ", ")))

check("ノットが狭義単調増加である", knots_ok(AGE_KNOTS),
      paste(sprintf("%.6g", AGE_KNOTS), collapse = ", "))
check("ノットを解析対象から取った（§8.9）",
      identical(KNOT_POPULATION, "analysis_set"),
      sprintf("population = %s", knot_label))
.out_b <- sum(DAT[[AGE_VAR]] < AGE_BOUND[1] | DAT[[AGE_VAR]] > AGE_BOUND[2])
check("境界ノットの外側にある当てはめ集団の例数を記録した", TRUE,
      sprintf("n = %d (%.1f%%)（ns() は境界の外では線形に外挿する）",
              .out_b, 100 * .out_b / N_M2))

num_txt <- function(x) paste0("c(", paste(vapply(x, function(z)
  sprintf("%.17g", z), character(1)), collapse = ", "), ")")
AGE_TERM_TXT <- sprintf("ns(%s, knots = %s, Boundary.knots = %s)",
                        AGE_VAR, num_txt(AGE_INNER), num_txt(AGE_BOUND))
say("\n  age term in the model formula:")
say("    ", AGE_TERM_TXT)

#' 年齢の値から基底行列を作る（ノットが明示されているため x に依存しない）
#' 使うノットは AGE_INNER_USE・AGE_BOUND_USE。6節で、当てはめに実際に使われた
#' ノット（有効数字15桁に丸められた値）に置き換える。
AGE_INNER_USE <- AGE_INNER
AGE_BOUND_USE <- AGE_BOUND
age_basis <- function(x) {
  B <- splines::ns(x, knots = AGE_INNER_USE, Boundary.knots = AGE_BOUND_USE)
  matrix(as.numeric(B), nrow = length(x),
         dimnames = list(NULL, paste0("b", seq_len(ncol(B)))))
}
K <- ncol(age_basis(DAT[[AGE_VAR]][1:2]))
say(sprintf("  spline degrees of freedom (columns of the basis) : %d", K))
check("自由度がノットの数 − 1 である", K == AGE_NKNOT - 1L,
      sprintf("df = %d, knots = %d", K, AGE_NKNOT))

knot_tab <- data.frame(
  variable       = AGE_VAR,
  variable_label = relabel_vars(AGE_VAR),
  knot_population = knot_label,
  n_knots        = AGE_NKNOT,
  knot_index     = seq_along(AGE_KNOTS),
  percentile     = 100 * AGE_KNOT_PROBS,
  age_value      = AGE_KNOTS,
  passed_as      = ifelse(seq_along(AGE_KNOTS) %in% c(1L, length(AGE_KNOTS)),
                          "Boundary.knots", "knots"),
  fallback_3knot = AGE_FALLBACK,
  reference_age  = AGE_REF,
  stringsAsFactors = FALSE
)


# -----------------------------------------------------------------------------
# 5. モデル式（共変量の構成は M2 と同じ1本。年齢群は入れない）
# -----------------------------------------------------------------------------

rule("5. Model formula")

COV_RHS <- paste(c("sex", "class", "support_in",
                   unname(spline_terms["mFIM_in_r"]),
                   unname(spline_terms["cFIM_in_r"]),
                   "period"), collapse = " + ")

## 05 のモデル式から年齢群を外したものと一致するかを点検する
.m2_rhs   <- sub("^y_ord ~ ", "", M2_FORMULA_TXT)
.m2_nogrp <- trimws(sub("^age_group \\+ ", "", .m2_rhs))
check("共変量の構成が M2 と同じである（年齢群を年齢スプラインに置き換えただけ）",
      identical(trimws(COV_RHS), .m2_nogrp),
      sprintf("10: %s\n                 05: %s", COV_RHS, .m2_nogrp))
check("ridit とノットを再算出していない（04 の項の文字列をそのまま使う。§6.4-5）",
      all(grepl("Boundary.knots = c\\(", unname(spline_terms))), "")

AGE_FORMULA_TXT <- paste("y_ord ~", AGE_TERM_TXT, "+", COV_RHS)
say("  ", AGE_FORMULA_TXT)
say("\n  年齢群は入れない。年齢の総合的な効果として読まない（§1.4、§8.9）。")

ctrl <- call_known_args(ordinal::clm.control,
                        list(maxIter = CLM_MAX_ITER, maxLineIter = CLM_MAX_LINE,
                             maxModIter = CLM_MAX_MOD_ITER))


# -----------------------------------------------------------------------------
# 6. 本体の当てはめ
# -----------------------------------------------------------------------------

rule("6. Fitting the continuous-age model")

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

t0 <- proc.time()[["elapsed"]]
fit_age <- ordinal::clm(stats::as.formula(AGE_FORMULA_TXT), data = DAT,
                        link = "logit", control = ctrl)
t_fit <- proc.time()[["elapsed"]] - t0
say(sprintf("  fitted in %.1f s", t_fit))

beta  <- fit_age[["beta"]]
V     <- tryCatch(stats::vcov(fit_age), error = function(e) NULL)

## ordinal::clm は受け取ったモデル式を内部で deparse して作り直す
## （ordinal:::get_clmFormulas）。deparse は数値を有効数字15桁で書くため、
##   (1) 係数名の中のノットは丸めた表示になり、AGE_TERM_TXT（%.17g で書いた文字列）
##       とは一致しない（例：85.549999999999955 → 85.55）。
##   (2) 当てはめに実際に使われるノットも有効数字15桁に丸めた値になり、4節の
##       AGE_INNER・AGE_BOUND と最下位の桁（相対 1e-15 程度）で食い違うことがある。
## (1) には、R と同じ規則で項の表示名を作り直して照合することで対処する。
AGE_TERM_LABEL <- paste(deparse(parse(text = AGE_TERM_TXT, keep.source = FALSE)[[1]],
                                width.cutoff = 500L, backtick = TRUE),
                        collapse = " ")
AGE_COEFS <- paste0(AGE_TERM_LABEL, seq_len(K))
check("年齢スプラインの係数名がモデル式の項と対応する",
      all(AGE_COEFS %in% names(beta)),
      sprintf("expected = %s", paste(AGE_COEFS, collapse = ", ")))
if (!all(AGE_COEFS %in% names(beta))) {
  say("  係数名の一覧: ", paste(names(beta), collapse = ", "))
  stop("年齢スプラインの係数を特定できない。")
}
if (!identical(AGE_TERM_LABEL, AGE_TERM_TXT)) {
  say("  係数名での年齢の項（ノットは有効数字15桁に丸めて表示される）:")
  say("    ", AGE_TERM_LABEL)
}

## (2) には、当てはめに実際に使われたノットをモデルフレームから取り出し、7節以降の
## 対比 B(a) − B(70) を作る age_basis() に使うことで対処する。こうすると対比が β と
## 厳密に対応する。丸めの差そのものは相対 1e-15 程度であり、結果には影響しない。
.mf <- fit_age[["model"]]
.B_fit <- if (is.null(.mf)) NULL else .mf[[AGE_TERM_LABEL]]
if (is.null(.B_fit)) {
  ## モデルフレームが保存されていない場合は、係数名と同じ（丸めた）項を評価して作る
  say("  注意：モデルフレームから年齢の基底を取り出せない。係数名の項を評価して代用する。")
  .B_fit <- eval(parse(text = AGE_TERM_LABEL, keep.source = FALSE)[[1]],
                 envir = DAT, enclos = environment())
}
AGE_INNER_FIT <- as.numeric(attr(.B_fit, "knots"))
AGE_BOUND_FIT <- as.numeric(attr(.B_fit, "Boundary.knots"))
.len_ok <- length(AGE_INNER_FIT) == length(AGE_INNER) && length(AGE_BOUND_FIT) == 2L
.dk <- if (.len_ok) max(abs(c(AGE_INNER_FIT - AGE_INNER, AGE_BOUND_FIT - AGE_BOUND)) /
                         abs(c(AGE_INNER, AGE_BOUND))) else NA_real_
.knots_ok <- .len_ok && is.finite(.dk) && .dk < 1e-12
check("当てはめに使われたノットが4節で算出したノットと一致する（15桁の丸めの範囲）",
      .knots_ok,
      sprintf("max relative diff = %s ; fitted knots = %s",
              formatC(.dk, format = "e", digits = 3),
              paste(sprintf("%.17g", c(AGE_BOUND_FIT[1], AGE_INNER_FIT,
                                       AGE_BOUND_FIT[2])), collapse = ", ")))
if (!.knots_ok) {
  stop("当てはめに使われた年齢スプラインのノットが4節で算出したものと一致しない。")
}
AGE_INNER_USE <- AGE_INNER_FIT
AGE_BOUND_USE <- AGE_BOUND_FIT

.dB <- if (NROW(.B_fit) != N_M2) NA_real_ else
  max(abs(matrix(as.numeric(.B_fit), nrow = N_M2) - age_basis(DAT[[AGE_VAR]])))
check("age_basis() が当てはめに使われた年齢の基底を再現する",
      is.finite(.dB) && .dB < 1e-12,
      sprintf("max |basis diff| = %s", formatC(.dB, format = "e", digits = 3)))
if (!(is.finite(.dB) && .dB < 1e-12)) {
  stop("age_basis() が当てはめに使われた年齢の基底を再現しない。対比が無効になる。")
}

.se  <- if (is.null(V)) rep(NA_real_, length(beta)) else sqrt(diag(V))[names(beta)]
.bns <- beta[!IS_SPLINE_TERM(names(beta))]
cc   <- conv_code(fit_age)
check("収束した", is.na(cc) || cc == 0L,
      sprintf("convergence = %s (%s)", as.character(cc), conv_msg(fit_age)))
check("係数・標準誤差に発散の兆候がない",
      all(is.finite(beta)) && (!length(.bns) || max(abs(.bns)) < DIVERGE_COEF) &&
        (all(is.na(.se)) || max(.se, na.rm = TRUE) < DIVERGE_SE),
      sprintf("max|beta(non-spline)| = %.3f, max se = %s",
              if (length(.bns)) max(abs(.bns)) else NA_real_,
              formatC(max(.se, na.rm = TRUE), format = "g", digits = 3)))

fit_info <- data.frame(
  model        = "cumulative logit, continuous age (restricted cubic spline)",
  n            = as.integer(stats::nobs(fit_age)),
  n_thresholds = length(fit_age[["alpha"]]),
  n_beta       = length(beta),
  logLik       = as.numeric(stats::logLik(fit_age)),
  AIC          = as.numeric(stats::AIC(fit_age)),
  convergence  = cc,
  convergence_message = conv_msg(fit_age),
  max_abs_beta = max(abs(beta)),
  max_se       = if (all(is.na(.se))) NA_real_ else max(.se, na.rm = TRUE),
  seconds      = t_fit,
  stringsAsFactors = FALSE
)
print(fit_info, row.names = FALSE, digits = 6)


# -----------------------------------------------------------------------------
# 7. 点推定値と Wald（デルタ法）区間の曲線（細かい格子）
#    log OR(a vs 70) = (B(a) − B(70))' β_age
# -----------------------------------------------------------------------------

rule("7. Point estimates and Wald intervals on the fine grid")

b_age  <- beta[AGE_COEFS]
V_age  <- if (is.null(V)) NULL else V[AGE_COEFS, AGE_COEFS, drop = FALSE]
B_ref  <- as.numeric(age_basis(AGE_REF))

age_min <- min(DAT[[AGE_VAR]])
age_max <- if (is.finite(AGE_MAX_GRID)) AGE_MAX_GRID else max(DAT[[AGE_VAR]])
grid_fine <- sort(unique(c(seq(age_min, age_max, by = CURVE_BY), AGE_REF, REPORT_AGES)))
grid_fine <- grid_fine[grid_fine >= age_min & grid_fine <= age_max]

contrast_vec <- function(a) as.numeric(age_basis(a)) - B_ref

z <- stats::qnorm(1 - (1 - CONF_LEVEL) / 2)
curve_rows <- lapply(grid_fine, function(a) {
  cvec <- contrast_vec(a)
  lo   <- sum(cvec * b_age)
  se   <- if (is.null(V_age)) NA_real_ else
    sqrt(as.numeric(t(cvec) %*% V_age %*% cvec))
  data.frame(
    age            = a,
    n_at_age       = sum(DAT[[AGE_VAR]] == a),
    log_odds_ratio = lo,
    odds_ratio     = exp(lo),
    se_log_odds_ratio = se,
    wald_lower     = exp(lo - z * se),
    wald_upper     = exp(lo + z * se),
    stringsAsFactors = FALSE
  )
})
age_curve_wald <- do.call(rbind, curve_rows)
age_curve_wald$reference_age <- AGE_REF
age_curve_wald$conf_level    <- CONF_LEVEL
age_curve_wald$ci_method     <- "Wald (delta method), unadjusted"

check("参照年齢でオッズ比が 1 になる",
      isTRUE(abs(age_curve_wald$odds_ratio[age_curve_wald$age == AGE_REF] - 1) < 1e-10),
      sprintf("OR(%g) = %.12f", AGE_REF,
              age_curve_wald$odds_ratio[age_curve_wald$age == AGE_REF][1]))
check("曲線に非有限の値がない", all(is.finite(age_curve_wald$odds_ratio)), "")

say("  odds ratios at the reporting ages (Wald intervals, for reference):")
print(age_curve_wald[age_curve_wald$age %in% REPORT_AGES,
                     c("age", "n_at_age", "odds_ratio", "wald_lower", "wald_upper")],
      row.names = FALSE, digits = 4)


# -----------------------------------------------------------------------------
# 8. profile likelihood 区間（補足図2 に載せる区間）
#    年齢 a ごとに、log OR(a vs 70) が第1係数になるようモデルを再母数化して
#    当てはめ直し、その係数の profile 区間を取る。
# -----------------------------------------------------------------------------

rule("8. Profile likelihood intervals")

grid_prof <- seq(ceiling(age_min), floor(age_max), by = PROFILE_BY)
grid_prof <- sort(unique(c(grid_prof, AGE_REF, REPORT_AGES)))
grid_prof <- grid_prof[grid_prof >= age_min & grid_prof <= age_max]
say(sprintf("  profile grid : %d ages from %g to %g (step %g)",
            length(grid_prof), min(grid_prof), max(grid_prof), PROFILE_BY))
say(sprintf("  本体の当てはめは %.1f s であった。", t_fit))
say("  profile 区間は1点につき当てはめを数十回くり返すため、1点あたりの所要は",
    "本体の当てはめの数十倍になる。")
say("  最初の数点の実測から下に ETA を出す。時間が問題になる場合は PROFILE_BY を",
    "大きくして格子を粗くすること（点推定値の曲線は細かい格子のまま出力される）。")

AZ_NAMES <- paste0("az", seq_len(K))
REPAR_FORMULA_TXT <- paste("y_ord ~", paste(AZ_NAMES, collapse = " + "), "+", COV_RHS)
say("  reparameterised formula : ", REPAR_FORMULA_TXT)

B_all <- age_basis(DAT[[AGE_VAR]])      # n x K

#' c に直交する K−1 本を足して正則な M を作り、A = M^{-1} を返す
repar_matrix <- function(cvec) {
  Q <- qr.Q(qr(matrix(cvec, ncol = 1L)), complete = TRUE)   # K x K
  M <- rbind(as.numeric(cvec), t(Q[, -1, drop = FALSE]))
  if (abs(det(M)) < 1e-12) return(NULL)
  list(M = M, A = solve(M))
}

prof_rows <- list()
t_start <- proc.time()[["elapsed"]]
for (i in seq_along(grid_prof)) {
  a    <- grid_prof[i]
  cvec <- contrast_vec(a)
  lo_pt <- sum(cvec * b_age)

  if (isTRUE(all.equal(a, AGE_REF)) || max(abs(cvec)) < 1e-12) {
    ## 参照点では対比が 0 であり、再母数化できない。OR = 1（区間は定義しない）。
    prof_rows[[length(prof_rows) + 1L]] <- data.frame(
      age = a, log_odds_ratio = 0, odds_ratio = 1,
      ci_lower = 1, ci_upper = 1, ok = TRUE,
      logLik_repar = as.numeric(stats::logLik(fit_age)),
      seconds = 0, note = "reference age", stringsAsFactors = FALSE)
    next
  }

  rp <- repar_matrix(cvec)
  if (is.null(rp)) {
    prof_rows[[length(prof_rows) + 1L]] <- data.frame(
      age = a, log_odds_ratio = lo_pt, odds_ratio = exp(lo_pt),
      ci_lower = NA_real_, ci_upper = NA_real_, ok = FALSE,
      logLik_repar = NA_real_, seconds = 0,
      note = "reparameterisation matrix is singular", stringsAsFactors = FALSE)
    next
  }

  d <- DAT
  Z <- B_all %*% rp$A
  for (k in seq_len(K)) d[[AZ_NAMES[k]]] <- Z[, k]

  ta <- proc.time()[["elapsed"]]
  f  <- tryCatch(ordinal::clm(stats::as.formula(REPAR_FORMULA_TXT), data = d,
                              link = "logit", control = ctrl),
                 error = function(e) e)
  if (inherits(f, "error")) {
    prof_rows[[length(prof_rows) + 1L]] <- data.frame(
      age = a, log_odds_ratio = lo_pt, odds_ratio = exp(lo_pt),
      ci_lower = NA_real_, ci_upper = NA_real_, ok = FALSE,
      logLik_repar = NA_real_, seconds = proc.time()[["elapsed"]] - ta,
      note = paste("clm error:", conditionMessage(f)), stringsAsFactors = FALSE)
    next
  }

  ci <- tryCatch({
    x <- suppressMessages(confint(f, parm = "az1", level = CONF_LEVEL,
                                  type = "profile"))
    if (is.null(dim(x))) x <- matrix(x, nrow = 1L, dimnames = list("az1", names(x)))
    as.matrix(x)["az1", ]
  }, error = function(e) {
    say("    age ", a, " : confint でエラー: ", conditionMessage(e)); c(NA, NA)
  })

  prof_rows[[length(prof_rows) + 1L]] <- data.frame(
    age = a,
    log_odds_ratio = unname(f[["beta"]]["az1"]),
    odds_ratio     = exp(unname(f[["beta"]]["az1"])),
    ci_lower       = exp(as.numeric(ci[1])),
    ci_upper       = exp(as.numeric(ci[2])),
    ok             = all(is.finite(as.numeric(ci))),
    logLik_repar   = as.numeric(stats::logLik(f)),
    seconds        = proc.time()[["elapsed"]] - ta,
    note           = "", stringsAsFactors = FALSE)

  if (i <= 3L || i %% 5L == 0L || i == length(grid_prof)) {
    el <- proc.time()[["elapsed"]] - t_start
    say(sprintf("    %3d / %d  (%5.1f s elapsed, ETA %5.1f s)",
                i, length(grid_prof), el, el / i * (length(grid_prof) - i)))
  }
}
age_profile <- do.call(rbind, prof_rows)
say(sprintf("  profile finished in %.1f s", proc.time()[["elapsed"]] - t_start))

## --- 再母数化が同値であることの点検 -----------------------------------------
ll0  <- as.numeric(stats::logLik(fit_age))
dll  <- max(abs(age_profile$logLik_repar - ll0), na.rm = TRUE)
check("再母数化した当てはめの対数尤度が本体と一致する（同値な母数化である）",
      is.finite(dll) && dll < 1e-6, sprintf("max |diff| = %.3e", dll))

pt_ref <- vapply(age_profile$age, function(a) sum(contrast_vec(a) * b_age), numeric(1))
dpt <- max(abs(age_profile$log_odds_ratio - pt_ref))
check("再母数化した第1係数が (B(a) − B(70))'β と一致する",
      is.finite(dpt) && dpt < 1e-6, sprintf("max |diff| = %.3e", dpt))

check("すべての格子点で profile 区間が得られた", all(age_profile$ok),
      sprintf("failed at %d ages: %s", sum(!age_profile$ok),
              paste(age_profile$age[!age_profile$ok], collapse = ", ")))
check("区間が点推定値を挟む",
      all(is.na(age_profile$ci_lower) |
            (age_profile$ci_lower <= age_profile$odds_ratio &
               age_profile$odds_ratio <= age_profile$ci_upper)), "")


# -----------------------------------------------------------------------------
# 9. 補足図2 の出力表
# -----------------------------------------------------------------------------

rule("9. Supplementary Figure 2 table")

## Wald（デルタ法）の値は格子の一致に頼らず、その場で計算する
wald_at <- do.call(rbind, lapply(age_profile$age, function(a) {
  cvec <- contrast_vec(a)
  lo   <- sum(cvec * b_age)
  se   <- if (is.null(V_age)) NA_real_ else
    sqrt(as.numeric(t(cvec) %*% V_age %*% cvec))
  data.frame(se_log_odds_ratio = se,
             wald_lower = exp(lo - z * se), wald_upper = exp(lo + z * se),
             stringsAsFactors = FALSE)
}))

age_or <- data.frame(
  age            = age_profile$age,
  age_label      = sprintf("%g years", age_profile$age),
  reference_age  = AGE_REF,
  n_at_age       = vapply(age_profile$age, function(a) sum(DAT[[AGE_VAR]] == a), integer(1)),
  n_at_or_above  = vapply(age_profile$age, function(a) sum(DAT[[AGE_VAR]] >= a), integer(1)),
  log_odds_ratio = age_profile$log_odds_ratio,
  odds_ratio     = age_profile$odds_ratio,
  ci_lower       = age_profile$ci_lower,
  ci_upper       = age_profile$ci_upper,
  conf_level     = CONF_LEVEL,
  ci_method      = "profile likelihood, unadjusted",
  ci_ok          = age_profile$ok,
  se_log_odds_ratio = wald_at$se_log_odds_ratio,
  wald_lower     = wald_at$wald_lower,
  wald_upper     = wald_at$wald_upper,
  note           = age_profile$note,
  p_value        = NA_real_,      # §8.5：層Cに p 値は付けない
  stringsAsFactors = FALSE
)

say("  odds ratios relative to age ", AGE_REF, " (profile likelihood 95% CI):")
print(age_or[age_or$age %in% REPORT_AGES,
             c("age", "n_at_age", "odds_ratio", "ci_lower", "ci_upper",
               "wald_lower", "wald_upper")],
      row.names = FALSE, digits = 4)

## profile と Wald の食い違い（ログのみ）
.d <- with(age_or[is.finite(age_or$ci_lower) & is.finite(age_or$wald_lower), ],
           max(c(abs(ci_lower - wald_lower), abs(ci_upper - wald_upper))))
say(sprintf("\n  max |profile - Wald| on the odds ratio scale : %.4f", .d))
say("  補足図2 に載せるのは profile 区間である（著者の決定、2026-09-20）。")

check("補足図2 に p 値を置いていない（§8.5）", all(is.na(age_or$p_value)), "")
check("層Cであり本文の結論には用いない（§7.5、§8.9）", TRUE,
      "3区分にもとづく解析（E-2・E-3・E-4）を本文の結論とする")

## G2（75〜89歳）の内部で年齢差が残っているか（§8.9 の目的）
g2_lo <- 75; g2_hi <- 89
.i <- age_or$age %in% c(g2_lo, g2_hi)
if (sum(.i) == 2L) {
  .or <- age_or$odds_ratio[.i]
  say(sprintf("\n  within G2 (%g to %g years): OR(%g vs %g) = %.3f / OR(%g vs %g) = %.3f ; ratio = %.3f",
              g2_lo, g2_hi, g2_lo, AGE_REF, .or[1], g2_hi, AGE_REF, .or[2],
              .or[2] / .or[1]))
  say("  この比は、G2 の15年幅の内部で年齢差がどれだけ残っているかの目安である（§8.9）。")
}


# -----------------------------------------------------------------------------
# 10. 諸元と脚注
# -----------------------------------------------------------------------------

meta <- data.frame(
  item = c("population", "n", "outcome", "exposure", "reference_age",
           "knot_population", "knot_percentiles", "knots", "spline_df",
           "model_formula", "ci_method", "conf_level", "profile_grid_step",
           "curve_grid_step", "layer"),
  value = c("M2 complete-case population", as.character(N_M2),
            relabel_vars(OUTCOME), "Age, continuous (restricted cubic spline)",
            as.character(AGE_REF), knot_label,
            paste(100 * AGE_KNOT_PROBS, collapse = ", "),
            paste(sprintf("%.6g", AGE_KNOTS), collapse = ", "),
            as.character(K), AGE_FORMULA_TXT,
            "profile likelihood, unadjusted", sprintf("%.2f", CONF_LEVEL),
            as.character(PROFILE_BY), as.character(CURVE_BY),
            "C (supplementary; the main conclusions rest on the three age groups)"),
  stringsAsFactors = FALSE
)

footnotes <- data.frame(
  n = seq_len(6),
  footnote = c(
    sprintf(paste("Odds ratios are for P(discharge motor FIM >= y) at any y, relative to",
                  "age %g, from a cumulative logit model with age as a restricted cubic",
                  "spline (%d knots at the %s percentiles of the analysis set)."),
            AGE_REF, AGE_NKNOT, paste(100 * AGE_KNOT_PROBS, collapse = "/")),
    paste("Covariates are the same as in the main model: sex, disease category,",
          "pre-admission care need, admission motor and cognitive FIM (ridit, spline),",
          "and admission period."),
    sprintf(paste("Age %g is the reference because %g years is the lower bound of the",
                  "eligible age range and an estimate anchored at the boundary is unstable."),
            AGE_REF, min(DAT[[AGE_VAR]])),
    paste("Intervals are unadjusted 95% profile likelihood intervals, obtained by",
          "refitting the model in a parameterisation in which the log odds ratio at",
          "each age is a single coefficient."),
    paste("This is a supplementary analysis. The conclusions of the paper rest on the",
          "three pre-specified age groups; no model selection was done on the basis of",
          "this figure."),
    paste("Beyond the boundary knots the spline is extrapolated linearly; the number of",
          "patients outside that range is reported in the supplementary methods.")
  ),
  stringsAsFactors = FALSE
)


# -----------------------------------------------------------------------------
# 11. 点検のまとめと保存
# -----------------------------------------------------------------------------

rule("11. Output")

checks <- do.call(rbind, .checks)
say("  REVIEW 項目: ", sum(checks$result == "REVIEW"), " / ", nrow(checks))
if (any(checks$result == "REVIEW")) {
  print(checks[checks$result == "REVIEW", ], row.names = FALSE)
}

w <- function(x, f) {
  p <- file.path(OUT_DIR, f)
  write.csv(x, p, row.names = FALSE, fileEncoding = "UTF-8")
  say("  written: ", p)
}
w(age_or,         "table_supp_fig2_age_or.csv")
w(age_curve_wald, "table_supp_fig2_curve_wald.csv")
w(knot_tab,       "table_supp_fig2_knots.csv")
w(fit_info,       "table_supp_fig2_fit.csv")
w(meta,           "table_supp_fig2_meta.csv")
w(footnotes,      "table_supp_fig2_footnotes.csv")
w(checks,         "table_checks_10.csv")

save(fit_age, age_or, age_curve_wald, age_profile, knot_tab, fit_info,
     meta, footnotes, AGE_KNOTS, AGE_INNER, AGE_BOUND, AGE_TERM_TXT,
     AGE_TERM_LABEL, AGE_COEFS, AGE_INNER_FIT, AGE_BOUND_FIT,
     AGE_FORMULA_TXT, REPAR_FORMULA_TXT, AGE_REF, K, checks,
     file = file.path(OUT_DIR, "10_age_continuous.rda"))
say("  written: ", file.path(OUT_DIR, "10_age_continuous.rda"))

say("\n--- sessionInfo ---")
print(utils::sessionInfo())

say("\nDone. N = ", N_M2, " ; profile grid = ", nrow(age_or), " ages ; reference = ",
    AGE_REF)

.log_close()
