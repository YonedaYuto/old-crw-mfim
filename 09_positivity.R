# =============================================================================
# 09_positivity.R
#   解析計画書 §8.7「周辺標準化の重なり（positivity）」の実装。
#   補足表3 を作るための数値をすべて出力する（作表・作図は別スクリプト）。
#
#   入力 : output/04_ridit_spline.rda  dat_m2 / main / ridit_summary /
#                                      period_by_g / M2_SETTINGS
#          01_labels.R（表示用の文字列）
#   出力 : output/table_supp3_sparse_cells.csv     度数0・10未満のセルの一覧（§8.7-1）
#          output/table_supp3_cell_summary.csv     クロス表の要約（全体は載せない）
#          output/table_supp3_period_by_agegroup.csv 期間区分 × 年齢群の2元表（全体）
#          output/table_supp3_ridit_overlap.csv    群別riditの第5・95百分位と重なり（§8.7-2）
#          output/table_supp3_extrapolation.csv    重なりの範囲の外にある例数（§8.7-3）
#          output/table_supp3.csv                  補足表3（上記を1表にまとめた長形式）
#          output/table_supp3_meta.csv             補足表3 の脚注に必要な諸元
#          output/table_supp3_footnotes.csv        補足表3 の脚注（下書き）
#          output/table_checks_09.csv              点検結果
#          output/09_positivity.rda
#          output/log_09_positivity.txt            実行ログ
#
#   計画との対応（§8.7）
#     1. カテゴリ共変量：年齢群 × 性別 × 疾患区分 × 病前の要介護状態 × 入院時期区分
#        のクロス表を作り、**度数が0または10未満のセルだけ**を一覧にする
#        （クロス表の全体は載せない）。あわせて期間区分 × 年齢群の2元表は全体を示す。
#        疾患区分に「その他」がある場合はその水準も含める（§6.1）
#     2. 連続的な共変量：入院時運動FIM・認知FIM の ridit について、年齢群別の
#        第5・第95百分位を1表に示す。3群すべてが観測されている ridit の範囲を
#        「重なりの範囲」として同じ表に記載する
#     3. 報告：標準化された値のうち、どの程度が重なりの範囲の外への外挿にあたるかを
#        限界の節に書く。ここではその材料となる例数と割合を出力する
#     4. **主解析の設定は変更しない。** 重なりが不十分であることを理由に対象を
#        絞り込む解析は行わない
#
#   対象集団について
#     §8.7 が問題にするのは周辺標準化の対象集団、すなわち M2 の完全ケース集団
#     （dat_m2）である。04_ridit_spline.R が出力した period_by_g と
#     ridit_summary は解析対象（main）にもとづくため、ここでは両方を算出して
#     population 列で区別する。補足表3 に載せるのは完全ケース集団の側である。
#
#   文字コード : UTF-8（Windows の R では source(..., encoding = "UTF-8")）
# =============================================================================


# -----------------------------------------------------------------------------
# 0. 設定
# -----------------------------------------------------------------------------

OUT_DIR   <- "output"
IN_RDA_04 <- file.path(OUT_DIR, "04_ridit_spline.rda")

## --- 疎なセルの基準（§8.7-1）------------------------------------------------
SPARSE_CUT <- 10L          # 度数が0、またはこの値未満のセルを一覧にする

## --- ridit の百分位（§8.7-2）------------------------------------------------
RIDIT_PROBS   <- c(0.05, 0.95)
QUANTILE_TYPE <- 1L        # §8.1 と同じ。04 の設定があればそちらに合わせる


# -----------------------------------------------------------------------------
# 1. 表示用ラベル
# -----------------------------------------------------------------------------

if (!file.exists("01_labels.R")) {
  stop("01_labels.R が見つからない。09_positivity.R は 01_labels.R と同じ",
       "フォルダを作業ディレクトリにして実行すること。")
}
source("01_labels.R", encoding = "UTF-8")


# -----------------------------------------------------------------------------
# 2. 実行の準備（出力先とログ）
# -----------------------------------------------------------------------------

if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

.log_con <- file(file.path(OUT_DIR, "log_09_positivity.txt"), open = "wt",
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

rule(paste0("09_positivity.R  ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
say("R : ", R.version.string)
say("sparse cell cut-off : n < ", SPARSE_CUT, " (zero cells are listed as well)")


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
take(e04, c("dat_m2", "main", "M2_SETTINGS"), IN_RDA_04)
if (exists("ridit_summary", envir = e04, inherits = FALSE)) {
  ridit_summary_04 <- get("ridit_summary", envir = e04)
}
if (exists("period_by_g", envir = e04, inherits = FALSE)) {
  period_by_g_04 <- get("period_by_g", envir = e04)
}

AGE_LEVELS <- M2_SETTINGS$age_levels
COVARS_CAT <- M2_SETTINGS$covars_cat          # sex, class, support_in, period
EXPOSURE   <- M2_SETTINGS$exposure            # age_group
RIDIT_VARS <- M2_SETTINGS$ridit_vars
if (!is.null(M2_SETTINGS$quantile_type)) QUANTILE_TYPE <- M2_SETTINGS$quantile_type

CROSS_VARS <- c(EXPOSURE, COVARS_CAT)         # §8.7-1 のクロス表を作る変数

CC <- as.data.frame(dat_m2, stringsAsFactors = FALSE)
CC$age_group <- factor(as.character(CC$age_group), levels = AGE_LEVELS)
MAIN <- as.data.frame(main, stringsAsFactors = FALSE)

say("dat_m2 (standardisation population, M2 complete cases) : ", nrow(CC), " rows")
say("main   (analysis set)                                  : ", nrow(MAIN), " rows")
say("cross-tabulation variables : ", paste(CROSS_VARS, collapse = " x "))

check("クロス表の変数が M2 のカテゴリ共変量と年齢群である（§8.7-1）",
      setequal(CROSS_VARS, c(EXPOSURE, COVARS_CAT)),
      paste(CROSS_VARS, collapse = ", "))
check("すべてのクロス表変数が因子である",
      all(vapply(CROSS_VARS, function(v) is.factor(CC[[v]]), logical(1))),
      paste(vapply(CROSS_VARS, function(v) sprintf("%s=%s", v, class(CC[[v]])[1]),
                   character(1)), collapse = ", "))
check("疾患区分の水準を確認した（「その他」があれば含める。§6.1）", TRUE,
      paste(levels(CC$class), collapse = " | "))


# -----------------------------------------------------------------------------
# 4. カテゴリ共変量のクロス表（§8.7-1）
#    度数0のセルも並ぶように、全水準の組合せを展開してから度数を数える。
# -----------------------------------------------------------------------------

rule("4. Cross-tabulation of the categorical covariates")

cell_counts <- function(d, vars) {
  tab <- table(d[, vars, drop = FALSE])          # 全水準の組合せを含む
  df  <- as.data.frame(tab, stringsAsFactors = FALSE)
  names(df)[names(df) == "Freq"] <- "n"
  df
}

cells <- cell_counts(CC, CROSS_VARS)
n_cells <- nrow(cells)
n_zero  <- sum(cells$n == 0L)
n_spar  <- sum(cells$n > 0L & cells$n < SPARSE_CUT)

say(sprintf("  total cells                 : %d (= %s)", n_cells,
            paste(vapply(CROSS_VARS, function(v) as.character(nlevels(CC[[v]])),
                         character(1)), collapse = " x ")))
say(sprintf("  cells with n = 0            : %d (%.1f%%)", n_zero,
            100 * n_zero / n_cells))
say(sprintf("  cells with 0 < n < %-2d       : %d (%.1f%%)", SPARSE_CUT, n_spar,
            100 * n_spar / n_cells))
say(sprintf("  patients in sparse cells    : %d (%.2f%% of %d)",
            sum(cells$n[cells$n < SPARSE_CUT]),
            100 * sum(cells$n[cells$n < SPARSE_CUT]) / nrow(CC), nrow(CC)))

## --- 一覧に載せるのは疎なセルだけ（クロス表の全体は載せない。§8.7-1）-------
sparse <- cells[cells$n < SPARSE_CUT, , drop = FALSE]
sparse <- sparse[order(sparse$n,
                       match(sparse[[EXPOSURE]], AGE_LEVELS)), , drop = FALSE]

## 表示用ラベルの列を足す（出力される文字列はすべて英語）
paste_cols <- function(d, cols) {
  if (!nrow(d)) return(character(0))
  apply(d[, cols, drop = FALSE], 1, function(r) paste(as.character(r), collapse = " / "))
}
for (v in CROSS_VARS) {
  sparse[[paste0(v, "_label")]] <- if (nrow(sparse))
    relabel_levels(v, as.character(sparse[[v]]), strict = FALSE) else character(0)
}
sparse$cell       <- paste_cols(sparse, CROSS_VARS)
sparse$cell_label <- paste_cols(sparse, paste0(CROSS_VARS, "_label"))
sparse$cut_off    <- if (nrow(sparse)) SPARSE_CUT else integer(0)
sparse$population <- if (nrow(sparse)) "M2 complete-case population" else character(0)
row.names(sparse) <- NULL

if (nrow(sparse)) {
  say("\n  sparse cells (n = 0 or n < ", SPARSE_CUT, "), the only ones listed in the table:")
  print(head(sparse[, c("cell_label", "n")], 50), row.names = FALSE)
  if (nrow(sparse) > 50) say("    ... (", nrow(sparse), " rows in total; see the csv)")
} else {
  say("\n  no sparse cell: every combination has at least ", SPARSE_CUT, " patients")
}

## --- 年齢群別の疎なセルの要約（外挿がどの群で起きるかを見る）----------------
cell_summary <- do.call(rbind, lapply(AGE_LEVELS, function(g) {
  s <- cells[cells[[EXPOSURE]] == g, , drop = FALSE]
  data.frame(
    age_group       = g,
    age_group_label = relabel_levels("age_group", g),
    n_cells         = nrow(s),
    n_cells_zero    = sum(s$n == 0L),
    n_cells_sparse  = sum(s$n > 0L & s$n < SPARSE_CUT),
    n_patients      = sum(s$n),
    n_patients_in_sparse_cells = sum(s$n[s$n < SPARSE_CUT]),
    percent_in_sparse_cells = if (sum(s$n) > 0)
      100 * sum(s$n[s$n < SPARSE_CUT]) / sum(s$n) else NA_real_,
    cut_off = SPARSE_CUT,
    stringsAsFactors = FALSE
  )
}))
say("\n  summary by age group (the full cross-table is not reported):")
print(cell_summary, row.names = FALSE, digits = 3)

check("疎なセルの一覧を作った（クロス表の全体は出力しない。§8.7-1）", TRUE,
      sprintf("listed %d of %d cells", nrow(sparse), n_cells))
check("度数0のセルがない", n_zero == 0L,
      sprintf("zero cells = %d。標準化はそこで外挿になる（限界の節に書く）", n_zero))
check(sprintf("度数が %d 未満のセルがない", SPARSE_CUT), n_spar == 0L,
      sprintf("sparse cells = %d (cut-off %d)", n_spar, SPARSE_CUT))
check("重なりが不十分でも主解析の設定を変えていない（§8.7-4）", TRUE,
      "対象を絞り込む解析は行わない")


# -----------------------------------------------------------------------------
# 5. 期間区分 × 年齢群の2元表（全体を示す。§8.7-1）
# -----------------------------------------------------------------------------

rule("5. Admission period x age group (the full two-way table)")

two_way <- function(d, label) {
  tb <- table(period = d$period, age_group = d$age_group)
  df <- as.data.frame.matrix(tb)
  out <- data.frame(
    population   = label,
    period       = rownames(df),
    period_label = relabel_levels("period", rownames(df)),
    df, check.names = FALSE, stringsAsFactors = FALSE
  )
  out$total <- rowSums(out[, AGE_LEVELS, drop = FALSE])
  ## 年齢群内の構成比（%）。期間を共変量に加えたことの妥当性が直接読める
  for (g in AGE_LEVELS) {
    out[[paste0(g, "_percent")]] <- 100 * out[[g]] / sum(out[[g]])
  }
  row.names(out) <- NULL
  out
}

period_by_g_cc   <- two_way(CC,   "M2 complete-case population")
period_by_g_main <- two_way(MAIN, "Analysis set")
period_by_agegroup <- rbind(period_by_g_cc, period_by_g_main)

say("  M2 complete-case population:")
print(period_by_g_cc, row.names = FALSE, digits = 4)
say("\n  analysis set (for reference):")
print(period_by_g_main, row.names = FALSE, digits = 4)

check("期間区分 × 年齢群に空セルがない（完全ケース集団。§14-7）",
      all(as.matrix(period_by_g_cc[, AGE_LEVELS, drop = FALSE]) > 0),
      sprintf("zero cells = %d",
              sum(as.matrix(period_by_g_cc[, AGE_LEVELS, drop = FALSE]) == 0)))
if (exists("period_by_g_04", inherits = TRUE)) {
  .a <- as.matrix(period_by_g_main[, AGE_LEVELS, drop = FALSE])
  .b <- as.matrix(period_by_g_04[, AGE_LEVELS, drop = FALSE])
  check("解析対象の2元表が 04 の出力と一致する",
        identical(dim(.a), dim(.b)) && all(.a == .b), "")
}


# -----------------------------------------------------------------------------
# 6. ridit の重なり（§8.7-2）
#    年齢群別の第5・第95百分位と、3群すべてが観測されている範囲。
# -----------------------------------------------------------------------------

rule("6. Overlap of the ridit distributions")

q_of <- function(x, p) {
  x <- x[!is.na(x)]
  if (!length(x)) return(rep(NA_real_, length(p)))
  stats::quantile(x, p, names = FALSE, type = QUANTILE_TYPE)
}

ridit_rows <- list()
overlap_rows <- list()
for (v in names(RIDIT_VARS)) {
  rv <- unname(RIDIT_VARS[v])
  for (g in c(AGE_LEVELS, "total")) {
    x <- if (identical(g, "total")) CC[[rv]] else CC[[rv]][CC$age_group == g]
    x <- x[!is.na(x)]
    qq <- q_of(x, RIDIT_PROBS)
    ridit_rows[[length(ridit_rows) + 1L]] <- data.frame(
      population      = "M2 complete-case population",
      variable        = v,
      variable_label  = relabel_vars(v),
      ridit_variable  = rv,
      age_group       = g,
      age_group_label = if (identical(g, "total")) "Total"
                        else relabel_levels("age_group", g),
      n               = length(x),
      p05             = qq[1],
      p95             = qq[2],
      min             = if (length(x)) min(x) else NA_real_,
      max             = if (length(x)) max(x) else NA_real_,
      stringsAsFactors = FALSE
    )
  }

  ## --- 重なりの範囲 ---------------------------------------------------------
  ## (a) 観測範囲の重なり : [max_g min(x_g), min_g max(x_g)]
  ## (b) 中央90%の重なり  : [max_g p05(x_g), min_g p95(x_g)]
  mins <- vapply(AGE_LEVELS, function(g) {
    x <- CC[[rv]][CC$age_group == g]; x <- x[!is.na(x)]
    if (length(x)) min(x) else NA_real_ }, numeric(1))
  maxs <- vapply(AGE_LEVELS, function(g) {
    x <- CC[[rv]][CC$age_group == g]; x <- x[!is.na(x)]
    if (length(x)) max(x) else NA_real_ }, numeric(1))
  p05s <- vapply(AGE_LEVELS, function(g)
    q_of(CC[[rv]][CC$age_group == g], RIDIT_PROBS[1]), numeric(1))
  p95s <- vapply(AGE_LEVELS, function(g)
    q_of(CC[[rv]][CC$age_group == g], RIDIT_PROBS[2]), numeric(1))

  overlap_rows[[length(overlap_rows) + 1L]] <- data.frame(
    variable        = v,
    variable_label  = relabel_vars(v),
    ridit_variable  = rv,
    overlap_kind    = c("observed range (all three groups)",
                        "central 90% (p05-p95 of all three groups)"),
    lower           = c(max(mins, na.rm = TRUE), max(p05s, na.rm = TRUE)),
    upper           = c(min(maxs, na.rm = TRUE), min(p95s, na.rm = TRUE)),
    stringsAsFactors = FALSE
  )
}
ridit_overlap_q <- do.call(rbind, ridit_rows)
ridit_overlap_r <- do.call(rbind, overlap_rows)
ridit_overlap_r$width <- ridit_overlap_r$upper - ridit_overlap_r$lower
ridit_overlap_r$is_empty <- !(ridit_overlap_r$upper > ridit_overlap_r$lower)

say("  ridit percentiles by age group (standardisation population):")
print(ridit_overlap_q[, c("variable", "age_group", "n", "p05", "p95", "min", "max")],
      row.names = FALSE, digits = 4)
say("\n  overlap ranges:")
print(ridit_overlap_r[, c("variable", "overlap_kind", "lower", "upper", "width")],
      row.names = FALSE, digits = 4)

check("すべての変数・群で百分位を算出できた",
      all(is.finite(ridit_overlap_q$p05) & is.finite(ridit_overlap_q$p95)), "")
check("観測範囲の重なりが空でない",
      all(!ridit_overlap_r$is_empty[
        ridit_overlap_r$overlap_kind == "observed range (all three groups)"]),
      "空であれば標準化はすべて外挿になる")

## --- (3) 重なりの範囲の外にある例数（限界の節の材料）------------------------
extrap_rows <- list()
for (v in names(RIDIT_VARS)) {
  rv <- unname(RIDIT_VARS[v])
  for (kind in unique(ridit_overlap_r$overlap_kind)) {
    rr <- ridit_overlap_r[ridit_overlap_r$ridit_variable == rv &
                            ridit_overlap_r$overlap_kind == kind, , drop = FALSE][1, ]
    for (g in c(AGE_LEVELS, "total")) {
      x <- if (identical(g, "total")) CC[[rv]] else CC[[rv]][CC$age_group == g]
      x <- x[!is.na(x)]
      out <- sum(x < rr$lower | x > rr$upper)
      extrap_rows[[length(extrap_rows) + 1L]] <- data.frame(
        variable        = v,
        variable_label  = relabel_vars(v),
        ridit_variable  = rv,
        overlap_kind    = kind,
        overlap_lower   = rr$lower,
        overlap_upper   = rr$upper,
        age_group       = g,
        age_group_label = if (identical(g, "total")) "Total"
                          else relabel_levels("age_group", g),
        n               = length(x),
        n_outside       = as.integer(out),
        percent_outside = if (length(x)) 100 * out / length(x) else NA_real_,
        stringsAsFactors = FALSE
      )
    }
  }
}
extrapolation <- do.call(rbind, extrap_rows)
say("\n  patients outside the overlap range (material for the limitations section):")
print(extrapolation[extrapolation$age_group != "total",
                    c("variable", "overlap_kind", "age_group", "n",
                      "n_outside", "percent_outside")],
      row.names = FALSE, digits = 3)

check("重なりの範囲の外にある例数を記録した（§8.7-3）", TRUE,
      "標準化のどの程度が外挿かを限界の節に書くための材料")


# -----------------------------------------------------------------------------
# 7. 補足表3（3つの内容を1表にまとめた長形式）
# -----------------------------------------------------------------------------

rule("7. Supplementary Table 3")

t3_sparse <- if (nrow(sparse)) data.frame(
  block = "Sparse cells of the categorical covariates",
  row_key = sparse$cell,
  row_label = sparse$cell_label,
  statistic = "n",
  value = as.numeric(sparse$n),
  detail = sprintf("cut-off n < %d; the full cross-table is not reported", SPARSE_CUT),
  stringsAsFactors = FALSE
) else data.frame(
  block = "Sparse cells of the categorical covariates",
  row_key = "(none)", row_label = "(none)", statistic = "n",
  value = NA_real_,
  detail = sprintf("no cell below the cut-off of %d", SPARSE_CUT),
  stringsAsFactors = FALSE
)

t3_period <- do.call(rbind, lapply(seq_len(nrow(period_by_g_cc)), function(i) {
  do.call(rbind, lapply(AGE_LEVELS, function(g) data.frame(
    block     = "Admission period x age group",
    row_key   = paste(period_by_g_cc$period[i], g, sep = " / "),
    row_label = paste(period_by_g_cc$period_label[i],
                      relabel_levels("age_group", g), sep = " / "),
    statistic = "n",
    value     = as.numeric(period_by_g_cc[[g]][i]),
    detail    = sprintf("%.1f%% of %s", period_by_g_cc[[paste0(g, "_percent")]][i], g),
    stringsAsFactors = FALSE)))
}))

t3_ridit <- do.call(rbind, lapply(seq_len(nrow(ridit_overlap_q)), function(i) {
  r <- ridit_overlap_q[i, ]
  rbind(
    data.frame(block = "Ridit percentiles by age group",
               row_key = paste(r$ridit_variable, r$age_group, sep = " / "),
               row_label = paste(r$variable_label, r$age_group_label, sep = " / "),
               statistic = "5th percentile", value = r$p05,
               detail = sprintf("n = %d", r$n), stringsAsFactors = FALSE),
    data.frame(block = "Ridit percentiles by age group",
               row_key = paste(r$ridit_variable, r$age_group, sep = " / "),
               row_label = paste(r$variable_label, r$age_group_label, sep = " / "),
               statistic = "95th percentile", value = r$p95,
               detail = sprintf("n = %d", r$n), stringsAsFactors = FALSE)
  )
}))

t3_overlap <- do.call(rbind, lapply(seq_len(nrow(ridit_overlap_r)), function(i) {
  r <- ridit_overlap_r[i, ]
  rbind(
    data.frame(block = "Overlap range of the ridit scores",
               row_key = paste(r$ridit_variable, r$overlap_kind, sep = " / "),
               row_label = paste(r$variable_label, r$overlap_kind, sep = " / "),
               statistic = "lower", value = r$lower,
               detail = sprintf("width = %.4f", r$width), stringsAsFactors = FALSE),
    data.frame(block = "Overlap range of the ridit scores",
               row_key = paste(r$ridit_variable, r$overlap_kind, sep = " / "),
               row_label = paste(r$variable_label, r$overlap_kind, sep = " / "),
               statistic = "upper", value = r$upper,
               detail = sprintf("width = %.4f", r$width), stringsAsFactors = FALSE)
  )
}))

table_supp3 <- rbind(t3_sparse, t3_period, t3_ridit, t3_overlap)
table_supp3$population <- "M2 complete-case population"
say("  Supplementary Table 3 assembled: ", nrow(table_supp3), " rows in ",
    length(unique(table_supp3$block)), " blocks")
print(table(table_supp3$block))


# -----------------------------------------------------------------------------
# 8. 諸元と脚注
# -----------------------------------------------------------------------------

meta <- data.frame(
  item = c("population", "n", "cross_tabulation_variables", "n_cells",
           "n_cells_zero", "n_cells_sparse", "sparse_cut_off",
           "n_patients_in_sparse_cells", "ridit_percentiles", "quantile_type",
           "class_levels", "setting_changed_because_of_overlap"),
  value = c("M2 complete-case population (the standardisation population)",
            as.character(nrow(CC)),
            paste(CROSS_VARS, collapse = " x "),
            as.character(n_cells), as.character(n_zero), as.character(n_spar),
            as.character(SPARSE_CUT),
            as.character(sum(cells$n[cells$n < SPARSE_CUT])),
            paste(100 * RIDIT_PROBS, collapse = ", "),
            as.character(QUANTILE_TYPE),
            paste(levels(CC$class), collapse = " | "),
            "no (plan 8.7-4)"),
  stringsAsFactors = FALSE
)

footnotes <- data.frame(
  n = seq_len(5),
  footnote = c(
    sprintf(paste("Only cells with no patients or with fewer than %d patients are listed.",
                  "The full cross-table is not reported."), SPARSE_CUT),
    paste("The two-way table of admission period by age group is shown in full because",
          "it reads directly on the decision to include the period in the model."),
    paste("The overlap range is the interval of ridit scores observed in all three age",
          "groups. Standardised values that rely on covariate patterns outside this",
          "range are extrapolations."),
    paste("Ridit scores were computed in the base population (A4) and were not",
          "recomputed here."),
    paste("The analysis was not restricted because of limited overlap; the main analysis",
          "is reported as pre-specified.")
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
w(sparse,             "table_supp3_sparse_cells.csv")
w(cell_summary,       "table_supp3_cell_summary.csv")
w(period_by_agegroup, "table_supp3_period_by_agegroup.csv")
w(ridit_overlap_q,    "table_supp3_ridit_percentiles.csv")
w(ridit_overlap_r,    "table_supp3_ridit_overlap.csv")
w(extrapolation,      "table_supp3_extrapolation.csv")
w(table_supp3,        "table_supp3.csv")
w(meta,               "table_supp3_meta.csv")
w(footnotes,          "table_supp3_footnotes.csv")
w(checks,             "table_checks_09.csv")

save(cells, sparse, cell_summary, period_by_agegroup,
     ridit_overlap_q, ridit_overlap_r, extrapolation,
     table_supp3, meta, footnotes, SPARSE_CUT, CROSS_VARS, checks,
     file = file.path(OUT_DIR, "09_positivity.rda"))
say("  written: ", file.path(OUT_DIR, "09_positivity.rda"))

say("\n--- sessionInfo ---")
print(utils::sessionInfo())

say("\nDone. cells = ", n_cells, " ; zero = ", n_zero, " ; sparse = ", n_spar,
    " ; N = ", nrow(CC))

.log_close()
