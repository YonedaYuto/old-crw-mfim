# =============================================================================
# 08_po_thresholds.R
#   解析計画書 §8.6(1)「閾値別オッズ比の散らばり」の実装。
#   補足図1 を描くための数値をすべて出力する（作図は別スクリプト）。
#
#   入力 : output/04_ridit_spline.rda     dat_m2 / spline_terms / M2_SETTINGS
#          output/05_m2_clm.rda           M2_FORMULA_TXT / m2_odds_ratio（共通OR）
#          output/06_standardize_boot.rda bin_tab（閾値65の二値モデル由来の標準化値）
#          01_labels.R（表示用の文字列）
#   出力 : output/table_supp_fig1_threshold_or.csv  閾値別ORと profile 95%区間
#          output/table_supp_fig1_spread.csv        散らばりの数値（(i)(ii)）
#          output/table_supp_fig1_binary65.csv      補足図1 に添える標準化値（06 から転記）
#          output/table_supp_fig1_excluded.csv      失敗規則で外した閾値と理由
#          output/table_supp_fig1_fit.csv           各当てはめの収束診断
#          output/table_supp_fig1_meta.csv          補足図1 の脚注に必要な諸元
#          output/table_supp_fig1_footnotes.csv     補足図1 の脚注（下書き）
#          output/table_checks_08.csv               点検結果
#          output/08_po_thresholds.rda
#          output/log_08_po_thresholds.txt          実行ログ
#
#   計画との対応（§8.6(1)）
#     ・対象は M2 の**主解析のみ**。感度分析Sでは再実行しない（§8.8「再実行しないもの」）
#     ・共変量の構成は M2 と同じ（入院時期区分を含む）。ridit とノットは 04 が
#       基準集団（A4）で固定した値をそのまま使う（§6.4-5）
#     ・閾値は 40・65・80 点。各閾値で二値ロジスティック回帰（最尤）を1本当てはめ、
#       年齢群のオッズ比と profile likelihood 95%区間を求める
#     ・**検定（Brant 検定など）による採否は行わない**
#     ・当てはまりは (i) 3閾値のORの最大値・最小値とその比、(ii) 各閾値のORを
#       共通ORで割った値、の2つを数値で示す。閾値を設けた「成り立つ／成り立たない」
#       の判定はしない
#     ・**失敗規則**：(i) 収束しない、(ii) 係数の絶対値または標準誤差が発散する
#       （分離の兆候）のいずれかに当たる閾値は点検から外し、外した事実と理由を
#       補足資料に明記する。罰則付き推定などへの切り替えは行わない
#     ・閾値65の二値モデル由来の標準化 P(Y ≥ 65)・RD は §8.6(2) の値であり、
#       06_standardize_boot.R がブートストラップ区間つきで算出済みである。
#       ここでは補足図1 に添える数値として転記するだけで、再計算しない
#
#   文字コード : UTF-8（Windows の R では source(..., encoding = "UTF-8")）
# =============================================================================


# -----------------------------------------------------------------------------
# 0. 設定
# -----------------------------------------------------------------------------

OUT_DIR   <- "output"
IN_RDA_04 <- file.path(OUT_DIR, "04_ridit_spline.rda")
IN_RDA_05 <- file.path(OUT_DIR, "05_m2_clm.rda")
IN_RDA_06 <- file.path(OUT_DIR, "06_standardize_boot.rda")

CONF_LEVEL <- 0.95            # 未調整95%（§8.5）

## --- 点検する閾値（§8.6(1)）-------------------------------------------------
## 65 は §6.3 の閾値。40・80 はその前後。y >= t を 1 とする。
THRESHOLDS      <- c(40L, 65L, 80L)
MAIN_THRESHOLD  <- 65L        # §8.6(2) の標準化値を添える閾値

## --- 対比（§8.3、§8.4）------------------------------------------------------
CONTRASTS <- list(
  list(key = "G2_vs_G1", ref = "G1", cmp = "G2", role = "secondary"),
  list(key = "G3_vs_G1", ref = "G1", cmp = "G3", role = "primary"),
  list(key = "G3_vs_G2", ref = "G2", cmp = "G3", role = "secondary")
)
PRIMARY_KEY <- "G3_vs_G1"

## --- 失敗規則（§8.6(1)）-----------------------------------------------------
## 05・06 と同じ閾値を用いる。スプライン基底の係数は値域が狭いため除外して判定する。
DIVERGE_COEF   <- 10
DIVERGE_SE     <- 10
IS_SPLINE_TERM <- function(nm) grepl("^ns\\(", nm)
## 分離の兆候をさらに直接に見るための目安（判定には使わず、記録のみ）
MIN_CELL_WARN  <- 5           # 年齢群 × アウトカムのセル度数がこれ未満なら記録


# -----------------------------------------------------------------------------
# 1. 表示用ラベル
# -----------------------------------------------------------------------------

if (!file.exists("01_labels.R")) {
  stop("01_labels.R が見つからない。08_po_thresholds.R は 01_labels.R と同じ",
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

.log_con <- file(file.path(OUT_DIR, "log_08_po_thresholds.txt"), open = "wt",
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

rule(paste0("08_po_thresholds.R  ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
say("R : ", R.version.string)
say("thresholds : ", paste(THRESHOLDS, collapse = ", "))

library(splines)
has_MASS <- requireNamespace("MASS", quietly = TRUE)
check("MASS が利用できる（glm の profile likelihood 区間に必要）", has_MASS,
      if (has_MASS) paste0("version ", as.character(utils::packageVersion("MASS")))
      else "install.packages(\"MASS\") が必要。Wald への差し替えはしない（§8.3）")


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

e04 <- load_env(IN_RDA_04); take(e04, c("dat_m2", "spline_terms", "M2_SETTINGS"), IN_RDA_04)
e05 <- load_env(IN_RDA_05); take(e05, c("M2_FORMULA_TXT", "m2_odds_ratio"), IN_RDA_05)
e06 <- load_env(IN_RDA_06, required = FALSE)

OUTCOME    <- M2_SETTINGS$outcome
AGE_LEVELS <- M2_SETTINGS$age_levels
COVARS_CAT <- M2_SETTINGS$covars_cat

DAT <- as.data.frame(dat_m2, stringsAsFactors = FALSE)
DAT$age_group <- factor(as.character(DAT$age_group), levels = AGE_LEVELS)
N_M2   <- nrow(DAT)
n_by_g <- table(DAT$age_group)

say("dat_m2 (M2 complete-case population) : ", N_M2, " rows")
for (g in AGE_LEVELS) {
  say(sprintf("  %-3s %-14s n = %5d", g, relabel_levels("age_group", g),
              as.integer(n_by_g[[g]])))
}

BIN_VAR <- "y_bin"
BIN_FORMULA_TXT <- sub("^y_ord ~", paste0(BIN_VAR, " ~"), M2_FORMULA_TXT)
say("\nmodel formula (binary) : ", BIN_FORMULA_TXT)

check("共変量の構成が M2 と同じである（§8.6(1)）",
      identical(sub("^y_bin ~", "", BIN_FORMULA_TXT),
                sub("^y_ord ~", "", M2_FORMULA_TXT)),
      "応答だけを二値に差し替えている")
check("ridit とノットを再算出していない（04 の項の文字列をそのまま使う。§6.4-5）",
      all(grepl("Boundary.knots = c\\(", unname(spline_terms))) &&
        all(vapply(unname(spline_terms),
                   function(s) grepl(s, M2_FORMULA_TXT, fixed = TRUE), logical(1))), "")
check("共通オッズ比（05 の出力）が3対比そろっている",
      nrow(m2_odds_ratio) == length(CONTRASTS) &&
        setequal(m2_odds_ratio$contrast, vapply(CONTRASTS, function(x) x$key, "")),
      paste(m2_odds_ratio$contrast, collapse = ", "))


# -----------------------------------------------------------------------------
# 4. 閾値ごとの当てはめ
#    年齢群の参照を G1 と G2 に置いた2本ずつを当てはめる。G3対G2 の区間は
#    参照 G2 の当てはめから取る（§8.3 と同じ扱い。線形対比の Wald は使わない）。
# -----------------------------------------------------------------------------

rule("4. Fitting the threshold-specific binary logistic models")

#' glm の profile likelihood 区間。版や MASS の有無で失敗しても Wald へは
#' 差し替えない（§8.3）。失敗したことを記録し、区間は NA のままにする。
get_profile_ci_glm <- function(fit, parm, level = CONF_LEVEL) {
  if (!has_MASS) return(list(ci = NULL, source = "MASS unavailable"))
  ci <- tryCatch({
    x <- suppressMessages(stats::confint(fit, parm = parm, level = level))
    if (is.null(dim(x))) x <- matrix(x, nrow = 1L, dimnames = list(parm, names(x)))
    as.matrix(x)
  }, error = function(e) {
    say("    confint.glm でエラー: ", conditionMessage(e)); NULL
  })
  if (is.null(ci) || !all(parm %in% rownames(ci))) {
    return(list(ci = NULL, source = "failed"))
  }
  list(ci = ci[parm, , drop = FALSE], source = "MASS::confint.glm (profile likelihood)")
}

get_wald_glm <- function(fit, parm, level = CONF_LEVEL) {
  b  <- stats::coef(fit)[parm]
  se <- tryCatch(sqrt(diag(stats::vcov(fit)))[parm],
                 error = function(e) rep(NA_real_, length(parm)))
  z  <- stats::qnorm(1 - (1 - level) / 2)
  data.frame(term = parm, beta = as.numeric(b), se = as.numeric(se),
             wald_lower = as.numeric(b - z * se), wald_upper = as.numeric(b + z * se),
             stringsAsFactors = FALSE)
}

fit_bin_at <- function(thr, ref) {
  d <- DAT
  d[[BIN_VAR]] <- as.integer(d[[OUTCOME]] >= thr)
  d$age_group  <- stats::relevel(factor(as.character(d$age_group),
                                        levels = AGE_LEVELS), ref = ref)
  t0  <- proc.time()[["elapsed"]]
  fit <- stats::glm(stats::as.formula(BIN_FORMULA_TXT), data = d,
                    family = stats::binomial())
  attr(fit, "seconds") <- proc.time()[["elapsed"]] - t0
  attr(fit, "ref")     <- ref
  attr(fit, "threshold") <- thr
  fit
}

fits      <- list()   # fits[[as.character(thr)]][[ref]]
fit_rows  <- list()
excl_rows <- list()
thr_ok    <- stats::setNames(rep(FALSE, length(THRESHOLDS)), as.character(THRESHOLDS))

for (thr in THRESHOLDS) {
  key_t <- as.character(thr)
  say("\n", strrep("=", 74))
  say(sprintf("Threshold : P(mFIM at discharge >= %d)", thr))

  y1 <- as.integer(DAT[[OUTCOME]] >= thr)
  tab <- table(age_group = DAT$age_group, reached = y1)
  say("  age group x outcome (n):")
  print(tab)
  n_min_cell <- min(tab)
  say(sprintf("  n reaching the threshold : %d / %d (%.1f%%) ; smallest cell = %d",
              sum(y1), length(y1), 100 * mean(y1), n_min_cell))

  ## --- (i) 閾値の両側に観測があるか ------------------------------------------
  if (any(tab == 0L)) {
    reason <- sprintf("empty cell in the age group x outcome table (smallest cell = %d)",
                      n_min_cell)
    say("  => 失敗規則により、この閾値を点検から外す : ", reason)
    excl_rows[[length(excl_rows) + 1L]] <- data.frame(
      threshold = thr, excluded = TRUE, reason = reason,
      smallest_cell = as.integer(n_min_cell), stringsAsFactors = FALSE)
    check(sprintf("[threshold %d] 点検に使える", thr), FALSE, reason)
    next
  }

  ## --- 当てはめ --------------------------------------------------------------
  f1 <- tryCatch(fit_bin_at(thr, "G1"), error = function(e) e)
  f2 <- tryCatch(fit_bin_at(thr, "G2"), error = function(e) e)
  if (inherits(f1, "error") || inherits(f2, "error")) {
    reason <- paste("glm error:",
                    conditionMessage(if (inherits(f1, "error")) f1 else f2))
    say("  => 失敗規則により外す : ", reason)
    excl_rows[[length(excl_rows) + 1L]] <- data.frame(
      threshold = thr, excluded = TRUE, reason = reason,
      smallest_cell = as.integer(n_min_cell), stringsAsFactors = FALSE)
    check(sprintf("[threshold %d] 当てはめができた", thr), FALSE, reason)
    next
  }

  ## --- 収束と発散の判定（失敗規則）------------------------------------------
  diag_one <- function(fit, tag) {
    cf   <- stats::coef(fit)
    cfn  <- cf[names(cf) != "(Intercept)"]
    cfns <- cfn[!IS_SPLINE_TERM(names(cfn))]
    se   <- tryCatch(sqrt(diag(stats::vcov(fit))),
                     error = function(e) stats::setNames(rep(NA_real_, length(cf)),
                                                         names(cf)))
    sen  <- se[names(cf) != "(Intercept)"]
    ## どの項の標準誤差が大きいかを残す。05 と同じく標準誤差の判定は
    ## スプライン基底も含めて当てるため、除外の理由を追えるようにしておく。
    big  <- names(sen)[!is.na(sen) & sen > DIVERGE_SE]
    data.frame(
      threshold   = attr(fit, "threshold"),
      reference   = attr(fit, "ref"),
      model       = tag,
      n           = as.integer(stats::nobs(fit)),
      n_events    = as.integer(sum(fit$y)),
      converged   = isTRUE(fit$converged),
      iterations  = as.integer(fit$iter),
      logLik      = as.numeric(stats::logLik(fit)),
      AIC         = as.numeric(stats::AIC(fit)),
      max_abs_beta = max(abs(cf)),
      max_abs_beta_nonspline = if (length(cfns)) max(abs(cfns)) else NA_real_,
      max_se      = if (all(is.na(sen))) NA_real_ else max(sen, na.rm = TRUE),
      max_se_term = if (all(is.na(sen))) NA_character_
                    else names(sen)[which.max(replace(sen, is.na(sen), -Inf))],
      terms_with_large_se = paste(big, collapse = " | "),
      seconds     = attr(fit, "seconds"),
      stringsAsFactors = FALSE
    )
  }
  d1 <- diag_one(f1, sprintf("binary >= %d (ref G1)", thr))
  d2 <- diag_one(f2, sprintf("binary >= %d (ref G2)", thr))
  fit_rows[[length(fit_rows) + 1L]] <- d1
  fit_rows[[length(fit_rows) + 1L]] <- d2
  print(rbind(d1, d2), row.names = FALSE, digits = 5)

  bad_conv <- !(d1$converged && d2$converged)
  bad_div  <- any(!is.finite(c(d1$max_abs_beta, d2$max_abs_beta))) ||
              any(c(d1$max_abs_beta_nonspline, d2$max_abs_beta_nonspline) > DIVERGE_COEF,
                  na.rm = TRUE) ||
              any(c(d1$max_se, d2$max_se) > DIVERGE_SE, na.rm = TRUE)

  if (bad_conv || bad_div) {
    .big <- unique(unlist(strsplit(c(d1$terms_with_large_se, d2$terms_with_large_se),
                                   " | ", fixed = TRUE)))
    .big <- .big[nzchar(.big)]
    reason <- paste(c(if (bad_conv) "did not converge",
                      if (bad_div) sprintf(
                        "coefficients or standard errors diverged (separation): max|beta(non-spline)| = %.2f, max se = %.2f%s",
                        max(c(d1$max_abs_beta_nonspline, d2$max_abs_beta_nonspline), na.rm = TRUE),
                        max(c(d1$max_se, d2$max_se), na.rm = TRUE),
                        if (length(.big)) paste0(" [terms: ", paste(.big, collapse = ", "), "]")
                        else "")),
                    collapse = "; ")
    say("  => 失敗規則により外す : ", reason)
    say("     §8.6(1)：罰則付き推定などへの切り替えは行わない。")
    excl_rows[[length(excl_rows) + 1L]] <- data.frame(
      threshold = thr, excluded = TRUE, reason = reason,
      smallest_cell = as.integer(n_min_cell), stringsAsFactors = FALSE)
    check(sprintf("[threshold %d] 収束し、発散の兆候がない", thr), FALSE, reason)
    next
  }

  check(sprintf("[threshold %d] 収束し、発散の兆候がない", thr), TRUE,
        sprintf("max|beta(non-spline)| = %.3f, max se = %.3f, smallest cell = %d",
                max(c(d1$max_abs_beta_nonspline, d2$max_abs_beta_nonspline), na.rm = TRUE),
                max(c(d1$max_se, d2$max_se), na.rm = TRUE), n_min_cell))
  if (n_min_cell < MIN_CELL_WARN) {
    check(sprintf("[threshold %d] 最小セル度数が %d 以上", thr, MIN_CELL_WARN),
          FALSE, sprintf("smallest cell = %d。分離に近い可能性を記録する", n_min_cell))
  }

  fits[[key_t]] <- list(G1 = f1, G2 = f2)
  thr_ok[[key_t]] <- TRUE
  excl_rows[[length(excl_rows) + 1L]] <- data.frame(
    threshold = thr, excluded = FALSE, reason = "",
    smallest_cell = as.integer(n_min_cell), stringsAsFactors = FALSE)
}

threshold_excluded <- do.call(rbind, excl_rows)
fit_diag <- if (length(fit_rows)) do.call(rbind, fit_rows) else
  data.frame(threshold = integer(0))

say("\n  thresholds used in the check : ",
    paste(names(thr_ok)[thr_ok], collapse = ", "))
if (any(!thr_ok)) {
  say("  thresholds excluded          : ",
      paste(names(thr_ok)[!thr_ok], collapse = ", "),
      "  （理由は table_supp_fig1_excluded.csv と補足資料に記す）")
}
check("点検に使える閾値が1つ以上ある", any(thr_ok),
      sprintf("used = %s", paste(names(thr_ok)[thr_ok], collapse = ", ")))
check("外した閾値の理由を記録した", TRUE,
      sprintf("excluded = %d", sum(!thr_ok)))


# -----------------------------------------------------------------------------
# 5. 閾値別オッズ比と profile likelihood 95%区間
# -----------------------------------------------------------------------------

rule("5. Threshold-specific odds ratios (profile likelihood 95% CI)")

or_rows <- list()
for (thr in THRESHOLDS) {
  key_t <- as.character(thr)
  if (!isTRUE(thr_ok[[key_t]])) next
  f1 <- fits[[key_t]]$G1
  f2 <- fits[[key_t]]$G2
  pl1 <- get_profile_ci_glm(f1, c("age_groupG2", "age_groupG3"))
  pl2 <- get_profile_ci_glm(f2, c("age_groupG3"))
  say(sprintf("  threshold %d : profile CI source (ref G1) = %s / (ref G2) = %s",
              thr, pl1$source, pl2$source))
  check(sprintf("[threshold %d] profile 区間が得られた（参照 G1）", thr),
        !is.null(pl1$ci), pl1$source)
  check(sprintf("[threshold %d] profile 区間が得られた（参照 G2）", thr),
        !is.null(pl2$ci), pl2$source)

  wd <- rbind(get_wald_glm(f1, c("age_groupG2", "age_groupG3")),
              get_wald_glm(f2, c("age_groupG3")))

  for (ct in CONTRASTS) {
    term <- paste0("age_group", ct$cmp)
    if (identical(ct$ref, "G1")) { f <- f1; pl <- pl1 } else { f <- f2; pl <- pl2 }
    b  <- unname(stats::coef(f)[term])
    ci <- if (!is.null(pl$ci) && term %in% rownames(pl$ci)) pl$ci[term, ] else c(NA, NA)
    wr <- wd[wd$term == term, , drop = FALSE]
    wr <- wr[if (identical(ct$ref, "G1")) 1L else nrow(wr), , drop = FALSE]
    or_rows[[length(or_rows) + 1L]] <- data.frame(
      model            = sprintf("binary logistic (>= %d)", thr),
      threshold        = thr,
      is_main_threshold = identical(as.integer(thr), as.integer(MAIN_THRESHOLD)),
      contrast         = ct$key,
      label            = contrast_label(ct$ref, ct$cmp),
      group_reference  = ct$ref,
      group_comparison = ct$cmp,
      role             = ct$role,
      term             = term,
      n_reference      = as.integer(n_by_g[[ct$ref]]),
      n_comparison     = as.integer(n_by_g[[ct$cmp]]),
      beta             = b,
      odds_ratio       = exp(b),
      ci_lower         = exp(as.numeric(ci[1])),
      ci_upper         = exp(as.numeric(ci[2])),
      conf_level       = CONF_LEVEL,
      ci_method        = "profile likelihood, unadjusted",
      ci_source        = if (identical(ct$ref, "G1")) pl1$source else pl2$source,
      wald_or_lower    = exp(wr$wald_lower[1]),
      wald_or_upper    = exp(wr$wald_upper[1]),
      p_value          = NA_real_,   # §8.5-1：閾値別ORに p 値は付けない
      stringsAsFactors = FALSE
    )
  }
}

## 共通オッズ比（05 の出力）を同じ形で1ブロック足す。補足図1 は両者を並べて描く。
common_rows <- data.frame(
  model            = "cumulative logit (M2), common odds ratio",
  threshold        = NA_integer_,
  is_main_threshold = FALSE,
  contrast         = m2_odds_ratio$contrast,
  label            = m2_odds_ratio$label,
  group_reference  = m2_odds_ratio$group_reference,
  group_comparison = m2_odds_ratio$group_comparison,
  role             = m2_odds_ratio$role,
  term             = m2_odds_ratio$term,
  n_reference      = m2_odds_ratio$n_reference,
  n_comparison     = m2_odds_ratio$n_comparison,
  beta             = m2_odds_ratio$beta,
  odds_ratio       = m2_odds_ratio$odds_ratio,
  ci_lower         = m2_odds_ratio$ci_lower,
  ci_upper         = m2_odds_ratio$ci_upper,
  conf_level       = m2_odds_ratio$conf_level,
  ci_method        = m2_odds_ratio$ci_method,
  ci_source        = "05_m2_clm.R",
  wald_or_lower    = NA_real_,
  wald_or_upper    = NA_real_,
  p_value          = NA_real_,
  stringsAsFactors = FALSE
)

threshold_or <- rbind(do.call(rbind, or_rows), common_rows)
threshold_or <- threshold_or[order(match(threshold_or$contrast,
                                         vapply(CONTRASTS, function(x) x$key, "")),
                                   threshold_or$threshold, na.last = TRUE), ]

say("\n  Supplementary Figure 1 : odds ratios by threshold and the common odds ratio")
print(threshold_or[, c("contrast", "model", "odds_ratio", "ci_lower", "ci_upper")],
      row.names = FALSE, digits = 4)

check("すべての行で区間が点推定値を挟む",
      all(is.na(threshold_or$ci_lower) |
            (threshold_or$ci_lower <= threshold_or$odds_ratio &
               threshold_or$odds_ratio <= threshold_or$ci_upper)), "")
check("閾値別ORに p 値を置いていない（§8.5-1）",
      all(is.na(threshold_or$p_value)), "")
check("Brant 検定などによる採否を行っていない（§8.6(1)）", TRUE,
      "散らばりは数値で示すだけで、成立／不成立の判定はしない")


# -----------------------------------------------------------------------------
# 6. 散らばりの数値（§8.6(1) の (i) と (ii)）
# -----------------------------------------------------------------------------

rule("6. Spread of the threshold-specific odds ratios")

spread_rows <- list()
for (ct in CONTRASTS) {
  sel <- threshold_or$contrast == ct$key & !is.na(threshold_or$threshold)
  ors <- threshold_or$odds_ratio[sel]
  thr <- threshold_or$threshold[sel]
  com <- common_rows$odds_ratio[common_rows$contrast == ct$key][1]
  ok  <- is.finite(ors)
  spread_rows[[length(spread_rows) + 1L]] <- data.frame(
    contrast          = ct$key,
    label             = contrast_label(ct$ref, ct$cmp),
    role              = ct$role,
    n_thresholds_used = sum(ok),
    thresholds_used   = paste(thr[ok], collapse = ", "),
    common_odds_ratio = com,
    or_min            = if (any(ok)) min(ors[ok]) else NA_real_,
    or_max            = if (any(ok)) max(ors[ok]) else NA_real_,
    or_at_min_threshold = if (any(ok)) ors[ok][which.min(thr[ok])] else NA_real_,
    or_at_max_threshold = if (any(ok)) ors[ok][which.max(thr[ok])] else NA_real_,
    ratio_max_over_min = if (any(ok) && min(ors[ok]) > 0) max(ors[ok]) / min(ors[ok])
                         else NA_real_,
    stringsAsFactors = FALSE
  )
}
spread <- do.call(rbind, spread_rows)

## (ii) 各閾値のORを共通ORで割った値
ratio_rows <- list()
for (ct in CONTRASTS) {
  com <- common_rows$odds_ratio[common_rows$contrast == ct$key][1]
  sel <- threshold_or$contrast == ct$key & !is.na(threshold_or$threshold)
  for (i in which(sel)) {
    ratio_rows[[length(ratio_rows) + 1L]] <- data.frame(
      contrast   = ct$key,
      label      = contrast_label(ct$ref, ct$cmp),
      threshold  = threshold_or$threshold[i],
      odds_ratio = threshold_or$odds_ratio[i],
      common_odds_ratio = com,
      ratio_to_common   = threshold_or$odds_ratio[i] / com,
      stringsAsFactors = FALSE
    )
  }
}
spread_ratio <- do.call(rbind, ratio_rows)

say("  (i) max / min of the threshold-specific odds ratios:")
print(spread[, c("contrast", "n_thresholds_used", "thresholds_used",
                 "common_odds_ratio", "or_min", "or_max", "ratio_max_over_min")],
      row.names = FALSE, digits = 4)
say("\n  (ii) each threshold-specific odds ratio divided by the common odds ratio:")
print(spread_ratio[, c("contrast", "threshold", "odds_ratio",
                       "common_odds_ratio", "ratio_to_common")],
      row.names = FALSE, digits = 4)
say("\n  記述の型（§8.6(1)）：「閾値によってオッズ比が x 倍の範囲で動いた」。",
    "閾値を設けた成立／不成立の判定はしない。")

check("散らばりを2種類の数値で示している（(i) 最大/最小の比、(ii) 共通ORとの比）",
      nrow(spread) == length(CONTRASTS) && nrow(spread_ratio) >= 1L, "")
if (sum(thr_ok) < length(THRESHOLDS)) {
  check("外した閾値があるため (i) は残った閾値だけの比である", FALSE,
        sprintf("used = %s / planned = %s",
                paste(names(thr_ok)[thr_ok], collapse = ", "),
                paste(THRESHOLDS, collapse = ", ")))
}


# -----------------------------------------------------------------------------
# 7. 補足図1 に添える標準化値（§8.6(2)。06 の出力を転記するだけ）
# -----------------------------------------------------------------------------

rule("7. Standardised risk and RD from the threshold-65 binary model (from 06)")

binary65 <- NULL
if (!is.null(e06) && exists("bin_tab", envir = e06, inherits = FALSE)) {
  binary65 <- get("bin_tab", envir = e06)
  binary65$source <- "06_standardize_boot.R (bootstrap percentile CI)"
  print(binary65[, c("quantity", "key", "label", "estimate", "ci_lower", "ci_upper")],
        row.names = FALSE, digits = 4)
  check("閾値65の二値モデル由来の標準化値を 06 から転記した", TRUE,
        "ここでは再計算しない（区間はブートストラップで付いている）")
  if (exists("diff_tab", envir = e06, inherits = FALSE)) {
    say("\n  ordinal-model vs binary-model standardised quantities (from 06):")
    print(get("diff_tab", envir = e06), row.names = FALSE, digits = 4)
    say("  本文・Table 2 に載せるのは順序モデル由来の値である（§8.6「事前の規則」）。")
  }
} else {
  binary65 <- data.frame(quantity = character(0), key = character(0),
                         label = character(0), estimate = numeric(0),
                         ci_lower = numeric(0), ci_upper = numeric(0),
                         stringsAsFactors = FALSE)
  check("閾値65の二値モデル由来の標準化値を 06 から転記した", FALSE,
        "06_standardize_boot.rda が無い。先に 06 を実行すること")
}


# -----------------------------------------------------------------------------
# 8. 諸元と脚注
# -----------------------------------------------------------------------------

rule("8. Meta and footnotes")

meta <- data.frame(
  item = c("population", "n", "n_G1", "n_G2", "n_G3", "outcome",
           "thresholds_planned", "thresholds_used", "thresholds_excluded",
           "model_formula", "ci_method_threshold", "ci_method_common",
           "conf_level", "multiplicity_adjustment"),
  value = c("M2 complete-case population", as.character(N_M2),
            as.character(as.integer(n_by_g[["G1"]])),
            as.character(as.integer(n_by_g[["G2"]])),
            as.character(as.integer(n_by_g[["G3"]])),
            relabel_vars(OUTCOME),
            paste(THRESHOLDS, collapse = ", "),
            paste(names(thr_ok)[thr_ok], collapse = ", "),
            if (any(!thr_ok)) paste(names(thr_ok)[!thr_ok], collapse = ", ") else "none",
            BIN_FORMULA_TXT,
            "profile likelihood, unadjusted",
            "profile likelihood, unadjusted (from the cumulative logit model)",
            sprintf("%.2f", CONF_LEVEL), "none (all intervals are unadjusted)"),
  stringsAsFactors = FALSE
)

footnotes <- data.frame(
  n = seq_len(6),
  footnote = c(
    paste("Threshold-specific odds ratios come from binary logistic models with the",
          "same covariates as the cumulative logit model, one model per threshold."),
    paste("No test of the proportional odds assumption was performed. The spread of the",
          "odds ratios is reported as a magnitude, not as a pass/fail judgement."),
    paste("Intervals are unadjusted 95% profile likelihood intervals. The interval for",
          "G3 vs G2 comes from a refit with G2 as the reference level."),
    if (any(!thr_ok))
      sprintf(paste("Threshold(s) %s were excluded because the maximum likelihood fit failed",
                    "(non-convergence or signs of separation). No penalised estimation was used."),
              paste(names(thr_ok)[!thr_ok], collapse = ", "))
    else "All planned thresholds could be fitted.",
    paste("The standardised risk and risk difference shown alongside come from the",
          "threshold-65 binary logistic model with bootstrap percentile intervals;",
          "they do not rely on the proportional odds assumption."),
    paste("The values reported in the main text are those from the cumulative logit",
          "model, whatever the spread seen here (pre-specified rule, plan 8.6).")
  ),
  stringsAsFactors = FALSE
)


# -----------------------------------------------------------------------------
# 9. 点検のまとめと保存
# -----------------------------------------------------------------------------

rule("9. Output")

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
w(threshold_or,       "table_supp_fig1_threshold_or.csv")
w(spread,             "table_supp_fig1_spread.csv")
w(spread_ratio,       "table_supp_fig1_spread_ratio.csv")
w(binary65,           "table_supp_fig1_binary65.csv")
w(threshold_excluded, "table_supp_fig1_excluded.csv")
w(fit_diag,           "table_supp_fig1_fit.csv")
w(meta,               "table_supp_fig1_meta.csv")
w(footnotes,          "table_supp_fig1_footnotes.csv")
w(checks,             "table_checks_08.csv")

save(threshold_or, spread, spread_ratio, binary65, threshold_excluded, fit_diag,
     meta, footnotes, fits, thr_ok, THRESHOLDS, BIN_FORMULA_TXT, checks,
     file = file.path(OUT_DIR, "08_po_thresholds.rda"))
say("  written: ", file.path(OUT_DIR, "08_po_thresholds.rda"))

say("\n--- sessionInfo ---")
print(utils::sessionInfo())

say("\nDone. thresholds used = ", paste(names(thr_ok)[thr_ok], collapse = ", "),
    " ; N = ", N_M2)

.log_close()
