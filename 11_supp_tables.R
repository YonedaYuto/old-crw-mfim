# =============================================================================
# 11_supp_tables.R
#   解析計画書 §8.1 の「Table 1 に載せず、初めから補足表に置くもの」の実装。
#   補足表1（3つの内容を1表にまとめる）と、補足表4（M2の共変量の係数）を作る。
#
#   入力 : output/02_preprocess.rda   dat_base（excl_reason / e3_category を含む）
#                                     e3_tab / flow / flow_g / overlap
#          output/04_ridit_spline.rda base（A4 ＋ 派生変数）/ main / miss_tab /
#                                     M2_SETTINGS
#          output/05_m2_clm.rda       m2_coefficients（補足表4）
#          01_labels.R（表示用の文字列）
#   出力 : output/table_supp1_missing.csv       (A) 変数別の欠測（年齢群別）
#          output/table_supp1_excluded.csv      (B) 除外者の割合と特性（年齢群別）
#          output/table_supp1_e3_reason.csv     (C) E3 の欠測理由の3区分（年齢群別）
#          output/table_supp1.csv               補足表1（A・B・C を1表に）
#          output/table_supp1_meta.csv          補足表1 の脚注に必要な諸元
#          output/table_supp1_footnotes.csv     補足表1 の脚注（下書き）
#          output/table_supp4.csv               補足表4（05 の係数表を整形）
#          output/table_supp4_footnotes.csv     補足表4 の脚注（下書き）
#          output/table_checks_11.csv           点検結果
#          output/11_supp_tables.rda
#          output/log_11_supp_tables.txt        実行ログ
#
#   計画との対応
#     ・補足表1（§8.1、§10.2）= 次の3つを1表にまとめたもの
#         (A) 変数別の欠測数と割合（年齢群別）
#         (B) 除外者（E1〜E3）の割合と特性（入院時の運動FIM・認知FIM、疾患区分。
#             年齢群別）。§14「主解析の対象は退院まで在院した患者に限られる」の
#             材料でもある
#         (C) E3 の欠測理由の3区分（§5.1）の年齢群別集計
#     ・要約の形式は §8.1 に揃える。FIM は中央値［第1四分位, 第3四分位］
#       （`quantile(type = 1)`）、カテゴリは n（%）
#     ・**群間比較の数値（p値・標準化差）は置かない**（§8.1）
#     ・補足表4（§8.3、§10.2）= M2 の共変量の係数。入院時期区分を含む。
#       解釈しない旨を脚注に記す（Table 2 fallacy）。数値は 05 の出力をそのまま使う
#
#   文字コード : UTF-8（Windows の R では source(..., encoding = "UTF-8")）
# =============================================================================


# -----------------------------------------------------------------------------
# 0. 設定
# -----------------------------------------------------------------------------

OUT_DIR   <- "output"
IN_RDA_02 <- file.path(OUT_DIR, "02_preprocess.rda")
IN_RDA_04 <- file.path(OUT_DIR, "04_ridit_spline.rda")
IN_RDA_05 <- file.path(OUT_DIR, "05_m2_clm.rda")

AGE_LEVELS <- c("G1", "G2", "G3")

## --- 要約の形式（§8.1）------------------------------------------------------
QUANTILE_TYPE <- 1L                       # 中央値・四分位に使う type
QUANTILE_PROBS <- c(0.25, 0.50, 0.75)

## --- (A) 欠測を数える変数 ---------------------------------------------------
## 年齢は A2 の適格条件であり欠測がない。アウトカム（退院時運動FIM）は E3 の
## 定義そのものであるため、A4 では欠測が現れ、解析対象では 0 になる。
MISS_VARS <- c("sex", "class", "support_in", "period",
               "mFIM_in", "cFIM_in", "mFIM_out")

## --- (B) 除外者の特性に使う変数（§8.1）--------------------------------------
EXCL_NUM_VARS <- c("mFIM_in", "cFIM_in")   # 中央値［Q1, Q3］
EXCL_CAT_VARS <- c("class")                # n（%）
EXCL_GROUPS   <- c("E1", "E2", "E3", "included")


# -----------------------------------------------------------------------------
# 1. 表示用ラベル
# -----------------------------------------------------------------------------

if (!file.exists("01_labels.R")) {
  stop("01_labels.R が見つからない。11_supp_tables.R は 01_labels.R と同じ",
       "フォルダを作業ディレクトリにして実行すること。")
}
source("01_labels.R", encoding = "UTF-8")


# -----------------------------------------------------------------------------
# 2. 実行の準備（出力先とログ）
# -----------------------------------------------------------------------------

if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

.log_con <- file(file.path(OUT_DIR, "log_11_supp_tables.txt"), open = "wt",
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

rule(paste0("11_supp_tables.R  ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
say("R : ", R.version.string)


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
take_if <- function(e, nm) {
  if (!is.null(e) && exists(nm, envir = e, inherits = FALSE)) get(nm, envir = e)
  else NULL
}

e02 <- load_env(IN_RDA_02); take(e02, c("dat_base"), IN_RDA_02)
e3_tab_02 <- take_if(e02, "e3_tab")
flow_g_02 <- take_if(e02, "flow_g")
overlap_02 <- take_if(e02, "overlap")

e04 <- load_env(IN_RDA_04); take(e04, c("base", "main", "M2_SETTINGS"), IN_RDA_04)
miss_tab_04 <- take_if(e04, "miss_tab")

e05 <- load_env(IN_RDA_05, required = FALSE)
m2_coefficients_05 <- take_if(e05, "m2_coefficients")

if (!is.null(M2_SETTINGS$age_levels)) AGE_LEVELS <- M2_SETTINGS$age_levels
if (!is.null(M2_SETTINGS$quantile_type)) {
  ## 04 の ridit の分位点の type と、§8.1 の要約統計量の type は別の決めごとである。
  ## §8.1 が type = 1 を指定しているため、ここでは上書きしない。
  say("  (note) 04 の quantile type = ", M2_SETTINGS$quantile_type,
      " ; §8.1 の要約統計量は type = ", QUANTILE_TYPE, " を使う")
}

BASE <- as.data.frame(base, stringsAsFactors = FALSE)   # A4 ＋ 派生変数
MAIN <- as.data.frame(main, stringsAsFactors = FALSE)   # 解析対象
BASE$age_group <- factor(as.character(BASE$age_group), levels = AGE_LEVELS)
MAIN$age_group <- factor(as.character(MAIN$age_group), levels = AGE_LEVELS)

say("base (A4)         : ", nrow(BASE), " rows")
say("main (analysis)   : ", nrow(MAIN), " rows")

check("A4 に除外理由の列がある（02 の excl_reason）", "excl_reason" %in% names(BASE),
      paste(setdiff(c("excl_reason", "e3_category"), names(BASE)), collapse = ", "))
check("A4 に E3 の区分の列がある（02 の e3_category）",
      "e3_category" %in% names(BASE), "")
check("A4 と解析対象の例数が整合する（解析対象＝除外理由なし）",
      sum(is.na(BASE$excl_reason)) == nrow(MAIN),
      sprintf("A4 without exclusion = %d / main = %d",
              sum(is.na(BASE$excl_reason)), nrow(MAIN)))
check("年齢群に欠測がない（A4）", !anyNA(BASE$age_group), "")

n_A4_by_g <- table(BASE$age_group)


# -----------------------------------------------------------------------------
# 4. (A) 変数別の欠測数と割合（年齢群別）
# -----------------------------------------------------------------------------

rule("4. (A) Missing values by variable and age group")

miss_of <- function(d, label, vars) {
  vars <- intersect(vars, names(d))
  do.call(rbind, lapply(vars, function(v) {
    do.call(rbind, lapply(c(AGE_LEVELS, "total"), function(g) {
      sel <- if (identical(g, "total")) rep(TRUE, nrow(d)) else d$age_group == g
      n   <- sum(sel)
      nm  <- sum(is.na(d[[v]][sel]))
      data.frame(
        population      = label,
        variable        = v,
        variable_label  = relabel_vars(v),
        age_group       = g,
        age_group_label = if (identical(g, "total")) "Total"
                          else relabel_levels("age_group", g),
        n               = n,
        n_missing       = nm,
        percent_missing = if (n > 0) 100 * nm / n else NA_real_,
        stringsAsFactors = FALSE
      )
    }))
  }))
}

miss_A4   <- miss_of(BASE, "Base population (A4)", MISS_VARS)
miss_main <- miss_of(MAIN, "Analysis set",         MISS_VARS)
supp1_missing <- rbind(miss_A4, miss_main)

say("  base population (A4):")
print(miss_A4[miss_A4$age_group == "total",
              c("variable_label", "n", "n_missing", "percent_missing")],
      row.names = FALSE, digits = 3)
say("\n  analysis set:")
print(miss_main[miss_main$age_group == "total",
                c("variable_label", "n", "n_missing", "percent_missing")],
      row.names = FALSE, digits = 3)

check("すべての対象変数について欠測を数えた",
      all(MISS_VARS %in% supp1_missing$variable),
      paste(setdiff(MISS_VARS, supp1_missing$variable), collapse = ", "))
check("解析対象では退院時運動FIMに欠測がない（E3 で除外済み）",
      isTRUE(sum(miss_main$n_missing[miss_main$variable == "mFIM_out" &
                                       miss_main$age_group == "total"]) == 0L), "")
if (!is.null(miss_tab_04)) {
  .a <- miss_main[miss_main$variable == "mFIM_in" & miss_main$age_group == "total",
                  "n_missing"]
  .b <- miss_tab_04[miss_tab_04$variable == "mFIM_in" &
                      miss_tab_04$age_group == "total", "n_missing"]
  check("解析対象の欠測数が 04 の miss_tab と一致する（共通の変数で照合）",
        length(.a) == 1L && length(.b) == 1L && .a == .b,
        sprintf("11: %s / 04: %s", paste(.a, collapse = ","), paste(.b, collapse = ",")))
}


# -----------------------------------------------------------------------------
# 5. (B) 除外者（E1〜E3）の割合と特性（年齢群別）
#    §8.1 の形式に揃える：FIM は中央値［Q1, Q3］、カテゴリは n（%）。
#    群間比較の数値（p値・標準化差）は置かない。
# -----------------------------------------------------------------------------

rule("5. (B) Excluded patients: proportion and characteristics")

BASE$excl_group <- ifelse(is.na(BASE$excl_reason), "included", BASE$excl_reason)
BASE$excl_group <- factor(BASE$excl_group, levels = EXCL_GROUPS)

excl_group_label <- function(k) {
  switch(k,
         E1 = "Excluded: in-hospital transfer",
         E2 = "Excluded: termination (death etc.)",
         E3 = "Excluded: discharge motor FIM missing",
         included = "Included in the main analysis",
         k)
}

med_iqr <- function(x) {
  x <- x[!is.na(x)]
  if (!length(x)) return(c(n = 0, q1 = NA_real_, med = NA_real_, q3 = NA_real_))
  q <- stats::quantile(x, QUANTILE_PROBS, names = FALSE, type = QUANTILE_TYPE)
  c(n = length(x), q1 = q[1], med = q[2], q3 = q[3])
}

excl_rows <- list()
for (g in c(AGE_LEVELS, "total")) {
  sel_g <- if (identical(g, "total")) rep(TRUE, nrow(BASE)) else BASE$age_group == g
  n_g   <- sum(sel_g)
  for (k in EXCL_GROUPS) {
    sel <- sel_g & BASE$excl_group == k
    n_k <- sum(sel)

    ## 割合（その年齢群の A4 に占める割合）
    excl_rows[[length(excl_rows) + 1L]] <- data.frame(
      age_group = g, excl_group = k, excl_group_label = excl_group_label(k),
      variable = "(count)", variable_label = "Patients",
      level = "", level_label = "",
      statistic = "n (% of A4 in the age group)",
      n = n_k, value = n_k,
      percent = if (n_g > 0) 100 * n_k / n_g else NA_real_,
      q1 = NA_real_, median = NA_real_, q3 = NA_real_,
      denominator = n_g, stringsAsFactors = FALSE)

    ## 連続変数：中央値［Q1, Q3］
    for (v in EXCL_NUM_VARS) {
      s <- med_iqr(BASE[[v]][sel])
      excl_rows[[length(excl_rows) + 1L]] <- data.frame(
        age_group = g, excl_group = k, excl_group_label = excl_group_label(k),
        variable = v, variable_label = relabel_vars(v),
        level = "", level_label = "",
        statistic = "median [Q1, Q3]",
        n = as.integer(s[["n"]]), value = s[["med"]],
        percent = NA_real_,
        q1 = s[["q1"]], median = s[["med"]], q3 = s[["q3"]],
        denominator = n_k, stringsAsFactors = FALSE)
    }

    ## カテゴリ変数：n（%）
    for (v in EXCL_CAT_VARS) {
      lv <- if (is.factor(BASE[[v]])) levels(BASE[[v]])
            else sort(unique(as.character(BASE[[v]][!is.na(BASE[[v]])])))
      xs <- as.character(BASE[[v]][sel])
      for (l in lv) {
        n_l <- sum(!is.na(xs) & xs == l)
        excl_rows[[length(excl_rows) + 1L]] <- data.frame(
          age_group = g, excl_group = k, excl_group_label = excl_group_label(k),
          variable = v, variable_label = relabel_vars(v),
          level = l, level_label = relabel_levels(v, l, strict = FALSE),
          statistic = "n (%)",
          n = n_l, value = n_l,
          percent = if (n_k > 0) 100 * n_l / n_k else NA_real_,
          q1 = NA_real_, median = NA_real_, q3 = NA_real_,
          denominator = n_k, stringsAsFactors = FALSE)
      }
      n_na <- sum(is.na(xs))
      excl_rows[[length(excl_rows) + 1L]] <- data.frame(
        age_group = g, excl_group = k, excl_group_label = excl_group_label(k),
        variable = v, variable_label = relabel_vars(v),
        level = "(missing)", level_label = "Missing",
        statistic = "n (%)",
        n = n_na, value = n_na,
        percent = if (n_k > 0) 100 * n_na / n_k else NA_real_,
        q1 = NA_real_, median = NA_real_, q3 = NA_real_,
        denominator = n_k, stringsAsFactors = FALSE)
    }
  }
}
supp1_excluded <- do.call(rbind, excl_rows)
supp1_excluded$age_group_label <- ifelse(supp1_excluded$age_group == "total", "Total",
                                         relabel_levels("age_group",
                                                        supp1_excluded$age_group,
                                                        strict = FALSE))

say("  number excluded, by reason and age group (% of A4 in the age group):")
.cnt <- supp1_excluded[supp1_excluded$variable == "(count)", ]
print(.cnt[, c("age_group", "excl_group_label", "n", "percent", "denominator")],
      row.names = FALSE, digits = 3)

say("\n  admission FIM of the excluded patients (median [Q1, Q3]):")
.num <- supp1_excluded[supp1_excluded$variable %in% EXCL_NUM_VARS &
                         supp1_excluded$age_group == "total", ]
print(.num[, c("excl_group_label", "variable_label", "n", "q1", "median", "q3")],
      row.names = FALSE, digits = 4)

check("除外群と解析対象の例数の合計が A4 と一致する",
      sum(.cnt$n[.cnt$age_group == "total"]) == nrow(BASE),
      sprintf("sum = %d / A4 = %d", sum(.cnt$n[.cnt$age_group == "total"]),
              nrow(BASE)))
check("解析対象の例数が main と一致する",
      isTRUE(.cnt$n[.cnt$age_group == "total" & .cnt$excl_group == "included"] ==
               nrow(MAIN)), "")
if (!is.null(flow_g_02)) {
  say("\n  (reference) flow by age group from 02_preprocess.R:")
  print(flow_g_02, row.names = FALSE)
}
check("補足表1 に群間比較の数値（p値・標準化差）を置いていない（§8.1）", TRUE,
      "記述統計だけを並べる")


# -----------------------------------------------------------------------------
# 6. (C) E3 の欠測理由の3区分（年齢群別）
# -----------------------------------------------------------------------------

rule("6. (C) Reason for the missing discharge motor FIM (E3)")

E3_CATS <- c("a", "b", "c")
e3_rows <- list()
is_e3 <- !is.na(BASE$excl_reason) & BASE$excl_reason == "E3"
for (k in E3_CATS) {
  sel_k <- is_e3 & !is.na(BASE$e3_category) & BASE$e3_category == k
  for (g in c(AGE_LEVELS, "total")) {
    sel <- if (identical(g, "total")) sel_k else sel_k & BASE$age_group == g
    n_e3 <- if (identical(g, "total")) sum(is_e3) else sum(is_e3 & BASE$age_group == g)
    e3_rows[[length(e3_rows) + 1L]] <- data.frame(
      category        = k,
      age_group       = g,
      age_group_label = if (identical(g, "total")) "Total"
                        else relabel_levels("age_group", g),
      n               = sum(sel),
      n_e3            = n_e3,
      percent_of_e3   = if (n_e3 > 0) 100 * sum(sel) / n_e3 else NA_real_,
      stringsAsFactors = FALSE
    )
  }
}
supp1_e3 <- do.call(rbind, e3_rows)

## 02 が持っているラベル（区分の説明と感度分析Sでの扱い）を引き継ぐ
if (!is.null(e3_tab_02) && all(c("category", "label") %in% names(e3_tab_02))) {
  supp1_e3$label <- e3_tab_02$label[match(supp1_e3$category, e3_tab_02$category)]
  if ("handling_in_S" %in% names(e3_tab_02)) {
    supp1_e3$handling_in_S <-
      e3_tab_02$handling_in_S[match(supp1_e3$category, e3_tab_02$category)]
  }
} else {
  supp1_e3$label <- supp1_e3$category
  check("02 の e3_tab から区分のラベルを引き継げた", FALSE,
        "ラベルが区分の記号のままになっている")
}

print(supp1_e3[supp1_e3$age_group != "total",
               c("category", "label", "age_group", "n", "percent_of_e3")],
      row.names = FALSE, digits = 3)

check("E3 の区分の合計が E3 の例数と一致する",
      sum(supp1_e3$n[supp1_e3$age_group == "total"]) == sum(is_e3),
      sprintf("sum = %d / E3 = %d",
              sum(supp1_e3$n[supp1_e3$age_group == "total"]), sum(is_e3)))
if (!is.null(e3_tab_02)) {
  .a <- supp1_e3$n[supp1_e3$age_group == "total"][match(e3_tab_02$category,
                                                        supp1_e3$category[supp1_e3$age_group == "total"])]
  check("E3 の区分の例数が 02 の e3_tab と一致する",
        all(.a == e3_tab_02$total),
        sprintf("11: %s / 02: %s", paste(.a, collapse = ", "),
                paste(e3_tab_02$total, collapse = ", ")))
}


# -----------------------------------------------------------------------------
# 7. 補足表1（A・B・C を1表に。§8.1「3つを1表にまとめる」）
# -----------------------------------------------------------------------------

rule("7. Supplementary Table 1")

t1_a <- data.frame(
  block     = "(A) Missing values",
  row_key   = paste(supp1_missing$population, supp1_missing$variable, sep = " / "),
  row_label = paste(supp1_missing$population, supp1_missing$variable_label, sep = " / "),
  age_group = supp1_missing$age_group,
  age_group_label = supp1_missing$age_group_label,
  statistic = "n missing (%)",
  n         = supp1_missing$n_missing,
  value     = supp1_missing$n_missing,
  percent   = supp1_missing$percent_missing,
  denominator = supp1_missing$n,
  detail    = "",
  stringsAsFactors = FALSE
)

t1_b <- data.frame(
  block     = "(B) Excluded patients (E1-E3)",
  row_key   = paste(supp1_excluded$excl_group, supp1_excluded$variable,
                    supp1_excluded$level, sep = " / "),
  row_label = trimws(paste(supp1_excluded$excl_group_label,
                           supp1_excluded$variable_label,
                           supp1_excluded$level_label, sep = " / ")),
  age_group = supp1_excluded$age_group,
  age_group_label = supp1_excluded$age_group_label,
  statistic = supp1_excluded$statistic,
  n         = supp1_excluded$n,
  value     = supp1_excluded$value,
  percent   = supp1_excluded$percent,
  denominator = supp1_excluded$denominator,
  detail    = ifelse(is.na(supp1_excluded$median), "",
                     sprintf("%g [%g, %g]", supp1_excluded$median,
                             supp1_excluded$q1, supp1_excluded$q3)),
  stringsAsFactors = FALSE
)

t1_c <- data.frame(
  block     = "(C) Reason for a missing discharge motor FIM",
  row_key   = paste("E3", supp1_e3$category, sep = " / "),
  row_label = paste("E3", supp1_e3$label, sep = " / "),
  age_group = supp1_e3$age_group,
  age_group_label = supp1_e3$age_group_label,
  statistic = "n (% of E3)",
  n         = supp1_e3$n,
  value     = supp1_e3$n,
  percent   = supp1_e3$percent_of_e3,
  denominator = supp1_e3$n_e3,
  detail    = if ("handling_in_S" %in% names(supp1_e3)) supp1_e3$handling_in_S else "",
  stringsAsFactors = FALSE
)

table_supp1 <- rbind(t1_a, t1_b, t1_c)
say("  Supplementary Table 1 assembled: ", nrow(table_supp1), " rows in ",
    length(unique(table_supp1$block)), " blocks")
print(table(table_supp1$block))

check("補足表1 が3つの内容を1表にまとめている（§8.1）",
      length(unique(table_supp1$block)) == 3L,
      paste(unique(table_supp1$block), collapse = " | "))

table_supp1_meta <- data.frame(
  item = c("base_population_A4", "analysis_set", "excluded_E1", "excluded_E2",
           "excluded_E3", "quantile_type", "summary_format_continuous",
           "summary_format_categorical", "group_comparison_statistics"),
  value = c(as.character(nrow(BASE)), as.character(nrow(MAIN)),
            as.character(sum(BASE$excl_group == "E1")),
            as.character(sum(BASE$excl_group == "E2")),
            as.character(sum(BASE$excl_group == "E3")),
            as.character(QUANTILE_TYPE), "median [Q1, Q3]", "n (%)",
            "none (no p-values, no standardised differences)"),
  stringsAsFactors = FALSE
)

table_supp1_footnotes <- data.frame(
  n = seq_len(5),
  footnote = c(
    paste("Continuous variables are summarised as median [first quartile, third",
          "quartile] (quantile type 1); categorical variables as n (%)."),
    paste("No between-group comparison is reported in this table: no p-values and no",
          "standardised differences."),
    paste("E1 is an in-hospital transfer, E2 a termination of the stay (death etc.),",
          "and E3 a missing discharge motor FIM. A patient who meets more than one",
          "criterion is assigned to a single reason in the order E2 > E1 > E3."),
    paste("The characteristics of the excluded patients are shown because the main",
          "analysis is restricted to patients who stayed until discharge; older",
          "patients are expected to be excluded more often."),
    paste("Categories (b) and (c) of a missing discharge motor FIM are removed from",
          "the population of sensitivity analysis S; category (a) is assigned the",
          "worst score there.")
  ),
  stringsAsFactors = FALSE
)


# -----------------------------------------------------------------------------
# 8. 補足表4（M2 の共変量の係数。05 の出力を整形するだけ）
# -----------------------------------------------------------------------------

rule("8. Supplementary Table 4 (covariate coefficients of M2)")

if (!is.null(m2_coefficients_05)) {
  table_supp4 <- m2_coefficients_05
  ## 閾値（切片に当たる母数）は行数が多いため、表では区別できるようにしておく
  table_supp4$block <- ifelse(table_supp4$type == "threshold",
                              "Threshold parameters", "Regression coefficients")
  table_supp4$reported_but_not_interpreted <- !table_supp4$is_age_group
  say("  rows : ", nrow(table_supp4), " (regression coefficients ",
      sum(table_supp4$type == "regression coefficient"), ", thresholds ",
      sum(table_supp4$type == "threshold"), ")")
  print(table_supp4[table_supp4$type == "regression coefficient",
                    c("term_label", "estimate", "std_error", "odds_ratio")],
        row.names = FALSE, digits = 4)
  check("補足表4 に入院時期区分の係数が含まれる（§8.3、§10.2）",
        any(grepl("^period", table_supp4$term)),
        paste(table_supp4$term[grepl("^period", table_supp4$term)], collapse = ", "))
  check("補足表4 の年齢群以外の係数に「解釈しない」旨の印がある",
        all(table_supp4$reported_but_not_interpreted[!table_supp4$is_age_group]), "")
  check("補足表4 の数値を 05 から再計算していない", TRUE,
        "05_m2_clm.R の m2_coefficients をそのまま使う")
} else {
  table_supp4 <- data.frame(term = character(0), stringsAsFactors = FALSE)
  check("補足表4 を 05 から取り込めた", FALSE,
        "05_m2_clm.rda が無い。先に 05 を実行すること")
}

table_supp4_footnotes <- data.frame(
  n = seq_len(3),
  footnote = c(
    paste("Coefficients of the covariates are reported for transparency only and are",
          "not interpreted (Table 2 fallacy). Only the age-group coefficients answer",
          "the question of this study."),
    paste("The admission period is included in the model and its coefficient is shown",
          "here; it is a variable the comparison is conditioned on, not an exposure."),
    paste("Admission motor and cognitive FIM enter as restricted cubic splines of the",
          "ridit scores; the individual basis coefficients have no interpretation on",
          "their own.")
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
w(supp1_missing,          "table_supp1_missing.csv")
w(supp1_excluded,         "table_supp1_excluded.csv")
w(supp1_e3,               "table_supp1_e3_reason.csv")
w(table_supp1,            "table_supp1.csv")
w(table_supp1_meta,       "table_supp1_meta.csv")
w(table_supp1_footnotes,  "table_supp1_footnotes.csv")
w(table_supp4,            "table_supp4.csv")
w(table_supp4_footnotes,  "table_supp4_footnotes.csv")
w(checks,                 "table_checks_11.csv")

save(supp1_missing, supp1_excluded, supp1_e3, table_supp1, table_supp1_meta,
     table_supp1_footnotes, table_supp4, table_supp4_footnotes, checks,
     file = file.path(OUT_DIR, "11_supp_tables.rda"))
say("  written: ", file.path(OUT_DIR, "11_supp_tables.rda"))

say("\n--- sessionInfo ---")
print(utils::sessionInfo())

say("\nDone. A4 = ", nrow(BASE), " ; analysis set = ", nrow(MAIN),
    " ; supplementary table 1 rows = ", nrow(table_supp1))

.log_close()
