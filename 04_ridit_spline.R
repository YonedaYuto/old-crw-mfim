# =============================================================================
# 04_ridit_spline.R
#   解析計画書 §6.4「入院時FIMの表し方」・§6.5「欠測の扱い」の実装
#   （plan.detail §11 のスクリプト 06 と 11 に相当。plan.summary §4-2 の③の前段）
#
#   入力 : output/02_preprocess.rda（dat_base = 基準集団A4、dat_main = 解析対象）
#          01_labels.R（表示用の文字列）
#   出力 : output/table_missing_by_agegroup.csv  変数別の欠測（年齢群別、§6.5-1）
#          output/table_missing_decision.csv     完全ケース／欠測指示変数の分岐判定
#          output/table_ridit_knots.csv          ridit のノット位置（§6.4）
#          output/table_ridit_summary.csv        基準集団の ridit の要約
#          output/table_period_by_agegroup.csv   入院時期区分 × 年齢群（§8.7、§14-7）
#          output/table_covariate_levels.csv     共変量の水準と参照水準の記録
#          output/table_checks_04.csv            点検結果
#          output/04_ridit_spline.rda            dat_m2 ほか（05・06 が読む）
#          output/log_04_ridit_spline.txt        実行ログ
#
#   このスクリプトが決めること
#     (1) 入院年度 fy_in と入院時期区分 period（§6.1）。02 は作っていないためここで作る
#     (2) 年齢群 age_group（§5.3）。02・03 と同じ規則をここに再掲し、02 の出力と照合
#     (3) 変数別の欠測数と合計脱落率、完全ケースとするかの分岐（§6.5、§14-5）
#     (4) 基準集団（A4）を基準とする ridit（§6.4-1）
#     (5) ridit のノット位置（基準集団の第5・35・65・95百分位。退化時は3ノット）
#     (6) M2 の解析データ dat_m2（完全ケース集団）と、モデル式に埋め込む
#         `ns()` の項の文字列
#
#   方針
#     1. ridit とノットは基準集団で一度だけ算出し、主解析・感度分析S・
#        ブートストラップの各回で同じものを用いる（§6.4-5）。このため、ここで
#        算出した数値をそのまま文字列としてモデル式に埋め込み、05・06 は
#        その文字列を使う。再算出の余地を残さない。
#     2. 期待例数・期待結果のハードコードは置かない（§13-4）。分岐判定（§6.5-2）
#        だけは計画が定めた閾値を用いるが、判定結果で解析を自動で切り替えず、
#        判定を出力して記録するにとどめる。
#     3. 出力される文字列はすべて英語（01_labels.R）。
#
#   文字コード : UTF-8（Windows の R では source(..., encoding = "UTF-8")）
# =============================================================================


# -----------------------------------------------------------------------------
# 0. 設定
# -----------------------------------------------------------------------------

OUT_DIR <- "output"
IN_RDA  <- file.path(OUT_DIR, "02_preprocess.rda")

## --- 年齢群（§5.3）。02_preprocess.R・03_m0_ranktest.R と同じ値であること ----
AGE_BREAKS <- c(65, 75, 90)           # G1: 65-74, G2: 75-89, G3: >=90
AGE_LEVELS <- c("G1", "G2", "G3")

## --- 入院時期区分（§6.1）----------------------------------------------------
## 年度は4月始まり。前半 = 2017〜2019年度、後半 = 2020〜2023年度。
## 区切りは COVID-19 の流行開始（2020年4月）と2020年度の診療報酬改定に揃える。
FY_START_MONTH <- 4L
PERIOD_EARLY_FY <- 2017:2019
PERIOD_LATE_FY  <- 2020:2023
PERIOD_LEVELS   <- c("early", "late")  # 参照は early（§6.1）

## --- M2 の変数（§6.1）-------------------------------------------------------
OUTCOME    <- "mFIM_out"                       # 退院時運動FIM（13〜91点）
EXPOSURE   <- "age_group"
COVARS_CAT <- c("sex", "class", "support_in", "period")
COVARS_NUM <- c("mFIM_in", "cFIM_in")          # ridit → スプラインに変換する
M2_VARS    <- c(EXPOSURE, COVARS_CAT, COVARS_NUM)

## --- 参照水準（§6.1）--------------------------------------------------------
## 計画が明記しているのは 年齢群 = G1、疾患区分 = 脳血管、入院時期区分 = early。
## 性別は「男性を参照」（§6.1）。病前の要介護状態の参照水準は計画に明記がないため
## 「なし（要介護認定なし・要支援）」を参照とする。これは表示上の決めごとであり、
## 年齢群のオッズ比にも標準化リスクにも影響しない（点検に記録する）。
REF_LEVELS <- list(
  age_group  = "G1",
  sex        = "男",
  class      = "脳血管",
  support_in = "なし",
  period     = "early"
)

## 疾患区分の水準の並び（§6.1、§14-8）。「その他」は元データにあれば4水準目。
CLASS_ORDER <- c("脳血管", "運動器", "廃用", "その他")

## --- ridit とスプライン（§6.4）----------------------------------------------
RIDIT_VARS     <- c(mFIM_in = "mFIM_in_r", cFIM_in = "cFIM_in_r")
KNOT_PROBS_4   <- c(0.05, 0.35, 0.65, 0.95)   # 4ノット（df = 3）
KNOT_PROBS_3   <- c(0.10, 0.50, 0.90)         # ノットが重なる場合の代替（§6.4-4）
QUANTILE_TYPE  <- 7L                          # R の既定。§6.2 の type = 1 は記述統計用
KNOT_MIN_GAP   <- 1e-8                        # 「重なる」とみなす間隔

## --- 欠測の分岐（§6.5-2、§14-5）---------------------------------------------
MISS_VAR_THRESHOLD   <- 0.05   # 各変数の欠測割合
MISS_TOTAL_THRESHOLD <- 0.10   # 合計脱落率

## --- 値域（§6.2）------------------------------------------------------------
MFIM_RANGE <- c(13, 91)
CFIM_RANGE <- c( 5, 35)


# -----------------------------------------------------------------------------
# 1. 表示用ラベル
# -----------------------------------------------------------------------------

if (!file.exists("01_labels.R")) {
  stop("01_labels.R が見つからない。04_ridit_spline.R は 01_labels.R と同じ",
       "フォルダを作業ディレクトリにして実行すること。")
}
source("01_labels.R", encoding = "UTF-8")


# -----------------------------------------------------------------------------
# 2. 実行の準備（出力先とログ）
# -----------------------------------------------------------------------------

if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

.log_con <- file(file.path(OUT_DIR, "log_04_ridit_spline.txt"), open = "wt",
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

rule(paste0("04_ridit_spline.R  ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
say("R : ", R.version.string)

has_splines <- requireNamespace("splines", quietly = TRUE)
check("splines が利用できる", has_splines, "")
if (has_splines) library(splines)


# -----------------------------------------------------------------------------
# 3. 入力の取り込み
# -----------------------------------------------------------------------------

rule("3. Input")

if (!exists("dat_base", inherits = TRUE) || !exists("dat_main", inherits = TRUE)) {
  if (!file.exists(IN_RDA)) {
    stop(IN_RDA, " が見つからない。先に 02_preprocess.R を実行すること。")
  }
  load(IN_RDA)                  # dat_base / dat_main / dat_sens / flow / flow_g ...
  say("loaded : ", IN_RDA)
} else {
  say("dat_base / dat_main はワークスペース上のものを使う")
}
for (.nm in c("dat_base", "dat_main")) {
  if (!exists(.nm, inherits = TRUE)) {
    stop(.nm, " が見つからない。02_preprocess.R の出力を確認すること。")
  }
}

base <- as.data.frame(dat_base, stringsAsFactors = FALSE)
main <- as.data.frame(dat_main, stringsAsFactors = FALSE)
say("dat_base (A4, ridit の基準集団) : ", nrow(base), " rows")
say("dat_main (analysis set)         : ", nrow(main), " rows")

.need <- c("id", "age", "day_in", "sex", "class", "support_in",
           "mFIM_in", "cFIM_in", OUTCOME)
.miss <- setdiff(.need, names(main))
if (length(.miss)) stop("dat_main に必要な列がない: ", paste(.miss, collapse = ", "))
.miss <- setdiff(c("mFIM_in", "cFIM_in", "age", "day_in"), names(base))
if (length(.miss)) stop("dat_base に必要な列がない: ", paste(.miss, collapse = ", "))


# -----------------------------------------------------------------------------
# 4. 派生変数（年齢群・入院年度・入院時期区分）
#    02_preprocess.R は年齢群を集計にしか使わず、入院時期区分を作っていない。
#    ここで作り、02 の出力（flow_g）と例数を照合する。
# -----------------------------------------------------------------------------

rule("4. Derived variables (age group, fiscal year, admission period)")

make_agegroup <- function(age) {
  factor(ifelse(is.na(age), NA_character_,
         ifelse(age >= AGE_BREAKS[3], "G3",
         ifelse(age >= AGE_BREAKS[2], "G2", "G1"))),
         levels = AGE_LEVELS)
}

#' 入院年度（4月始まり）
make_fy <- function(day) {
  y <- as.integer(format(day, "%Y"))
  m <- as.integer(format(day, "%m"))
  ifelse(is.na(y) | is.na(m), NA_integer_, ifelse(m >= FY_START_MONTH, y, y - 1L))
}

#' 入院時期区分（§6.1）。前半・後半のどちらにも当たらない年度は NA とし、点検に出す。
make_period <- function(fy) {
  factor(ifelse(is.na(fy), NA_character_,
         ifelse(fy %in% PERIOD_EARLY_FY, "early",
         ifelse(fy %in% PERIOD_LATE_FY,  "late", NA_character_))),
         levels = PERIOD_LEVELS)
}

add_derived <- function(d) {
  d$age_group <- make_agegroup(d$age)
  d$fy_in     <- make_fy(d$day_in)
  d$period    <- make_period(d$fy_in)
  d
}

base <- add_derived(base)
main <- add_derived(main)

n_by_g <- table(main$age_group)
for (g in AGE_LEVELS) {
  say(sprintf("  %-3s %-14s n = %5d", g, relabel_levels("age_group", g),
              as.integer(n_by_g[[g]])))
}

## 入院年度の一覧（区分に当たらない年度がないかを見る）
.fy_tab <- table(main$fy_in, useNA = "ifany")
say("\n  fiscal year of admission (analysis set):")
print(.fy_tab)

check("年齢群に欠測がない", !anyNA(main$age_group),
      sprintf("missing = %d", sum(is.na(main$age_group))))
check("入院時期区分に欠測がない（§6.1：入院日から一意に決まる）",
      !anyNA(main$period),
      sprintf("missing = %d ; fiscal years out of FY%d-FY%d = %s",
              sum(is.na(main$period)), min(PERIOD_EARLY_FY), max(PERIOD_LATE_FY),
              paste(sort(unique(main$fy_in[is.na(main$period)])), collapse = ", ")))

## 02 の出力との照合（年齢群の規則が 02 と一致するか）
if (exists("flow_g", inherits = TRUE) && "stage" %in% names(flow_g) &&
    "MAIN" %in% flow_g$stage) {
  .row <- flow_g[flow_g$stage == "MAIN", , drop = FALSE][1, ]
  .exp <- as.integer(c(.row$G1, .row$G2, .row$G3))
  check("年齢群の例数が 02_preprocess.R の flow_g（MAIN）と一致する",
        identical(as.integer(n_by_g), .exp),
        sprintf("04: %s / 02: %s", paste(as.integer(n_by_g), collapse = ", "),
                paste(.exp, collapse = ", ")))
} else {
  check("02 の flow_g と照合できた", FALSE,
        "flow_g が見つからないため照合を行っていない")
}

## 入院時期区分 × 年齢群（§8.7、§14-7。補足表3 の一部でもある）
period_by_g <- as.data.frame.matrix(table(period = main$period,
                                          age_group = main$age_group))
period_by_g <- data.frame(
  period       = rownames(period_by_g),
  period_label = relabel_levels("period", rownames(period_by_g)),
  period_by_g, check.names = FALSE, stringsAsFactors = FALSE
)
period_by_g$total <- rowSums(period_by_g[, AGE_LEVELS, drop = FALSE])
say("\n  Admission period x age group (analysis set):")
print(period_by_g, row.names = FALSE)

check("入院時期区分 × 年齢群に空セルがない（§14-7）",
      all(as.matrix(period_by_g[, AGE_LEVELS, drop = FALSE]) > 0),
      sprintf("zero cells = %d",
              sum(as.matrix(period_by_g[, AGE_LEVELS, drop = FALSE]) == 0)))


# -----------------------------------------------------------------------------
# 5. 共変量の水準と参照水準（§6.1、§14-8）
# -----------------------------------------------------------------------------

rule("5. Covariate levels and reference levels")

#' 水準の並びを決める。参照水準を先頭に置き、既知の並び順（order）に従い、
#' 想定外の水準は末尾に回して点検に出す。
make_factor <- function(x, ref, order = NULL) {
  obs <- sort(unique(as.character(x[!is.na(x)])))
  lv  <- unique(c(ref, intersect(order, obs), setdiff(obs, c(ref, order))))
  factor(as.character(x), levels = lv)
}

main$age_group  <- factor(as.character(main$age_group), levels = AGE_LEVELS)
main$sex        <- make_factor(main$sex,        REF_LEVELS$sex)
main$class      <- make_factor(main$class,      REF_LEVELS$class, CLASS_ORDER)
main$support_in <- make_factor(main$support_in, REF_LEVELS$support_in)
main$period     <- factor(as.character(main$period), levels = PERIOD_LEVELS)

cov_levels <- do.call(rbind, lapply(c(EXPOSURE, COVARS_CAT), function(v) {
  lv <- levels(main[[v]])
  data.frame(
    variable       = v,
    variable_label = relabel_vars(v),
    n_levels       = length(lv),
    reference      = lv[1],
    reference_label = relabel_levels(v, lv[1], strict = FALSE),
    levels         = paste(lv, collapse = " | "),
    levels_label   = paste(relabel_levels(v, lv, strict = FALSE), collapse = " | "),
    stringsAsFactors = FALSE
  )
}))
print(cov_levels[, c("variable", "reference_label", "levels_label")], row.names = FALSE)

check("疾患区分の水準がすべて既知である（§6.1、§14-8：「その他」の有無）",
      all(levels(main$class) %in% CLASS_ORDER),
      sprintf("observed = %s ; 4番目の水準「その他」: %s",
              paste(levels(main$class), collapse = ", "),
              if ("その他" %in% levels(main$class)) "あり（4水準目として保持）" else "なし"))
check("参照水準が計画どおりに置かれている",
      all(vapply(names(REF_LEVELS), function(v)
        identical(levels(main[[v]])[1], REF_LEVELS[[v]]), logical(1))),
      paste(vapply(names(REF_LEVELS), function(v)
        sprintf("%s=%s", v, levels(main[[v]])[1]), character(1)), collapse = ", "))
check("病前の要介護状態の参照水準は計画に明記がなく、このスクリプトで決めている",
      TRUE,
      sprintf("reference = %s (%s)。係数は解釈しない（§8.3、Table 2 fallacy）",
              REF_LEVELS$support_in,
              relabel_levels("support_in", REF_LEVELS$support_in, strict = FALSE)))
check("各水準に少なくとも1例ある",
      all(vapply(c(EXPOSURE, COVARS_CAT), function(v)
        all(table(main[[v]]) > 0), logical(1))), "")


# -----------------------------------------------------------------------------
# 6. 変数別の欠測と完全ケースの分岐（§6.5、§14-5）
#    アウトカム（退院時運動FIM）は E3 で除外済みのため、ここでは共変量のみを見る。
# -----------------------------------------------------------------------------

rule("6. Missing data and the complete-case decision")

miss_rows <- list()
for (v in M2_VARS) {
  for (g in c(AGE_LEVELS, "total")) {
    sel <- if (identical(g, "total")) rep(TRUE, nrow(main)) else main$age_group == g
    n   <- sum(sel)
    nm  <- sum(is.na(main[[v]][sel]))
    miss_rows[[length(miss_rows) + 1L]] <- data.frame(
      variable       = v,
      variable_label = relabel_vars(v),
      age_group      = g,
      age_group_label = if (identical(g, "total")) "Total"
                        else relabel_levels("age_group", g),
      n              = n,
      n_missing      = nm,
      percent_missing = if (n > 0) 100 * nm / n else NA_real_,
      stringsAsFactors = FALSE
    )
  }
}
miss_tab <- do.call(rbind, miss_rows)
say("  Missing by variable and age group:")
print(miss_tab[, c("variable_label", "age_group", "n", "n_missing",
                   "percent_missing")], row.names = FALSE, digits = 3)

## 完全ケース集団（M2 の全変数がそろっている例）
complete <- stats::complete.cases(main[, M2_VARS, drop = FALSE])
n_main   <- nrow(main)
n_cc     <- sum(complete)
dropout  <- if (n_main > 0) (n_main - n_cc) / n_main else NA_real_

## 変数ごとの欠測割合（全体）
p_var <- vapply(M2_VARS, function(v) mean(is.na(main[[v]])), numeric(1))
decision_complete_case <- all(p_var < MISS_VAR_THRESHOLD) &&
                          is.finite(dropout) && dropout < MISS_TOTAL_THRESHOLD

say(sprintf("\n  analysis set N            : %d", n_main))
say(sprintf("  complete cases (M2)       : %d", n_cc))
say(sprintf("  total dropout rate        : %.2f%% (threshold %.0f%%)",
            100 * dropout, 100 * MISS_TOTAL_THRESHOLD))
say(sprintf("  max missing in one var    : %.2f%% (threshold %.0f%%) [%s]",
            100 * max(p_var), 100 * MISS_VAR_THRESHOLD,
            names(p_var)[which.max(p_var)]))
say(sprintf("  decision (§6.5)           : %s",
            if (decision_complete_case) "complete-case analysis"
            else "complete-case for the main analysis + missing-indicator version as a supplement"))

miss_decision <- data.frame(
  item = c("analysis_set_n", "complete_case_n", "total_dropout_rate",
           "max_missing_proportion", "max_missing_variable",
           "threshold_per_variable", "threshold_total",
           "decision", "note"),
  value = c(as.character(n_main), as.character(n_cc),
            sprintf("%.6f", dropout), sprintf("%.6f", max(p_var)),
            names(p_var)[which.max(p_var)],
            sprintf("%.2f", MISS_VAR_THRESHOLD),
            sprintf("%.2f", MISS_TOTAL_THRESHOLD),
            if (decision_complete_case) "complete-case analysis"
            else "complete-case (main) + missing-indicator version (supplement)",
            paste("The main analysis is the complete-case analysis in either case",
                  "(plan 6.5-3). The missing-indicator version, if required,",
                  "is a separate script.")),
  stringsAsFactors = FALSE
)

check("共変量の欠測がすべて 5% 未満（§6.5-2）",
      all(p_var < MISS_VAR_THRESHOLD),
      paste(sprintf("%s=%.2f%%", M2_VARS, 100 * p_var), collapse = ", "))
check("合計脱落率が 10% 未満（§6.5-2）",
      is.finite(dropout) && dropout < MISS_TOTAL_THRESHOLD,
      sprintf("%.2f%%", 100 * dropout))
check("完全ケース集団の各年齢群に例がある",
      all(table(main$age_group[complete]) > 0),
      paste(sprintf("%s=%d", AGE_LEVELS,
                    as.integer(table(main$age_group[complete]))), collapse = ", "))
if (!decision_complete_case) {
  check("§6.5-3 に該当（欠測指示変数の版を補足資料に置く必要がある）", FALSE,
        "主解析は完全ケースのまま。補足の版は別スクリプトで作ること")
}


# -----------------------------------------------------------------------------
# 7. ridit 変換（§6.4-1）
#    基準集団（A4 = dat_base）を基準とし、
#      ridit(x) = (基準集団で x 未満の人数 + 0.5 × x と同点の人数) / 基準集団の人数
#    分母は当該変数が観測されている基準集団の人数とする（欠測は基準に含めない）。
# -----------------------------------------------------------------------------

rule("7. Ridit transformation (reference population = A4)")

#' 基準ベクトルから ridit の変換関数を作る
#' @param ref 基準集団における当該変数の観測値（欠測は取り除く）
make_ridit <- function(ref) {
  ref <- ref[!is.na(ref)]
  n   <- length(ref)
  if (!n) stop("ridit の基準集団に観測値がない")
  u      <- sort(unique(ref))
  cnt    <- as.integer(tabulate(match(ref, u), nbins = length(u)))
  cum_le <- cumsum(cnt)              # その水準以下の人数
  cum_lt <- cum_le - cnt             # その水準未満の人数
  r_u    <- (cum_lt + 0.5 * cnt) / n # 観測水準の ridit
  function(x) {
    out <- rep(NA_real_, length(x))
    ok  <- !is.na(x)
    if (!any(ok)) return(out)
    xi    <- x[ok]
    hit   <- match(xi, u)
    val   <- rep(NA_real_, length(xi))
    inref <- !is.na(hit)
    val[inref] <- r_u[hit[inref]]
    ## 基準集団に現れない値は「未満の人数 ÷ 基準集団の人数」とする（外挿の受け皿）。
    ## 解析対象 ⊆ 基準集団であるため、主解析ではこの枝に入らない（点検で確認する）。
    if (any(!inref)) {
      idx <- findInterval(xi[!inref], u)   # その値より小さい観測水準の数
      val[!inref] <- ifelse(idx > 0L, cum_le[pmax(idx, 1L)], 0) / n
    }
    out[ok] <- val
    out
  }
}

ridit_fun <- list()
ridit_meta_rows <- list()
for (v in names(RIDIT_VARS)) {
  ref_v <- base[[v]]
  ridit_fun[[v]] <- make_ridit(ref_v)
  n_ref <- sum(!is.na(ref_v))
  say(sprintf("  %-8s reference population : n = %d (missing in A4: %d)",
              v, n_ref, sum(is.na(ref_v))))
  ridit_meta_rows[[v]] <- data.frame(
    variable        = v,
    variable_label  = relabel_vars(v),
    ridit_variable  = unname(RIDIT_VARS[v]),
    n_reference     = n_ref,
    n_reference_missing = sum(is.na(ref_v)),
    n_distinct_scores   = length(unique(ref_v[!is.na(ref_v)])),
    stringsAsFactors = FALSE
  )
  ## 基準集団に無い値が解析対象に現れないことの点検（A4 ⊇ 解析対象であるため0のはず）
  .novel <- setdiff(unique(main[[v]][!is.na(main[[v]])]),
                    unique(ref_v[!is.na(ref_v)]))
  check(sprintf("%s：解析対象の値がすべて基準集団に現れる", v),
        length(.novel) == 0,
        if (length(.novel)) paste("not in A4:", paste(sort(.novel), collapse = ", "))
        else "")
}

## 変換の適用（基準集団・解析対象の両方に同じ関数を当てる）
for (v in names(RIDIT_VARS)) {
  rv <- unname(RIDIT_VARS[v])
  base[[rv]] <- ridit_fun[[v]](base[[v]])
  main[[rv]] <- ridit_fun[[v]](main[[v]])
}

## ridit の性質の点検：(0,1) に収まり、元の点数と同順（単調非減少）であること
for (v in names(RIDIT_VARS)) {
  rv <- unname(RIDIT_VARS[v])
  r  <- base[[rv]][!is.na(base[[rv]])]
  check(sprintf("%s の ridit が (0, 1) に収まる", v),
        all(r > 0 & r < 1),
        sprintf("range = [%.6f, %.6f]", min(r), max(r)))
  o  <- order(base[[v]])
  rr <- base[[rv]][o]
  rr <- rr[!is.na(rr)]
  check(sprintf("%s の ridit が点数に対して単調非減少である", v),
        all(diff(rr) >= -1e-12), "")
  check(sprintf("%s の ridit の平均が 0.5 である（基準集団、定義上の帰結）", v),
        abs(mean(r) - 0.5) < 1e-8,
        sprintf("mean = %.10f", mean(r)))
}

ridit_meta <- do.call(rbind, ridit_meta_rows)

## 基準集団と解析対象における ridit の要約（補足表3 の材料にもなる）
ridit_summary <- do.call(rbind, lapply(names(RIDIT_VARS), function(v) {
  rv <- unname(RIDIT_VARS[v])
  do.call(rbind, lapply(c("A4", AGE_LEVELS, "MAIN"), function(g) {
    x <- if (identical(g, "A4")) base[[rv]]
         else if (identical(g, "MAIN")) main[[rv]]
         else main[[rv]][main$age_group == g]
    x <- x[!is.na(x)]
    q <- if (length(x)) stats::quantile(x, c(0.05, 0.50, 0.95), names = FALSE,
                                        type = QUANTILE_TYPE) else rep(NA_real_, 3)
    data.frame(variable = v, ridit_variable = rv, population = g, n = length(x),
               p05 = q[1], p50 = q[2], p95 = q[3],
               min = if (length(x)) min(x) else NA_real_,
               max = if (length(x)) max(x) else NA_real_,
               stringsAsFactors = FALSE)
  }))
}))
say("\n  Ridit distribution (A4 = reference population, G1-G3 = analysis set):")
print(ridit_summary, row.names = FALSE, digits = 4)


# -----------------------------------------------------------------------------
# 8. ノット位置（§6.4-2〜4）
#    基準集団における ridit の第5・35・65・95百分位に4ノット。
#    重なる場合は第10・50・90百分位の3ノットに減らす。
# -----------------------------------------------------------------------------

rule("8. Knot positions for the restricted cubic splines")

#' ノットが狭義単調増加であるか
knots_ok <- function(k) all(is.finite(k)) && all(diff(k) > KNOT_MIN_GAP)

knot_list <- list()
knot_rows <- list()
for (v in names(RIDIT_VARS)) {
  rv <- unname(RIDIT_VARS[v])
  x  <- base[[rv]][!is.na(base[[rv]])]

  k4 <- stats::quantile(x, KNOT_PROBS_4, names = FALSE, type = QUANTILE_TYPE)
  if (knots_ok(k4)) {
    k <- k4; probs <- KNOT_PROBS_4; nknot <- 4L; fb <- FALSE
  } else {
    k3 <- stats::quantile(x, KNOT_PROBS_3, names = FALSE, type = QUANTILE_TYPE)
    if (!knots_ok(k3)) {
      stop("ridit のノットが重なる（", v, "）。4ノットでも3ノットでも狭義単調に",
           "ならない。§6.4-4 の代替を超える事態であり、計画の見直しを要する。")
    }
    k <- k3; probs <- KNOT_PROBS_3; nknot <- 3L; fb <- TRUE
  }

  ## §6.4-3：内側を knots、外側を Boundary.knots に渡す
  inner <- k[-c(1, length(k))]
  bound <- k[c(1, length(k))]
  knot_list[[rv]] <- list(all = k, inner = inner, boundary = bound,
                          probs = probs, n_knots = nknot, fallback = fb)

  say(sprintf("  %-10s %d knots at percentiles %s", rv, nknot,
              paste0(100 * probs, collapse = "/")))
  say(sprintf("             knots           = %s",
              paste(sprintf("%.6f", k), collapse = ", ")))
  say(sprintf("             inner (knots)   = %s",
              paste(sprintf("%.6f", inner), collapse = ", ")))
  say(sprintf("             Boundary.knots  = %s",
              paste(sprintf("%.6f", bound), collapse = ", ")))
  if (fb) say("             （4ノットが重なったため §6.4-4 により3ノットに減らした）")

  for (i in seq_along(k)) {
    knot_rows[[length(knot_rows) + 1L]] <- data.frame(
      variable       = v,
      variable_label = relabel_vars(v),
      ridit_variable = rv,
      n_knots        = nknot,
      knot_index     = i,
      percentile     = 100 * probs[i],
      ridit_value    = k[i],
      passed_as      = if (i == 1L || i == length(k)) "Boundary.knots" else "knots",
      fallback_3knot = fb,
      stringsAsFactors = FALSE
    )
  }
  check(sprintf("%s のノットが狭義単調増加である", rv), knots_ok(k),
        paste(sprintf("%.6f", k), collapse = ", "))
  ## ノットの外側にある観測の数（線形外挿になる範囲。限界の節の材料）
  .out <- sum(main[[rv]] < bound[1] | main[[rv]] > bound[2], na.rm = TRUE)
  check(sprintf("%s：境界ノットの外側にある解析対象の例数を記録した", rv), TRUE,
        sprintf("n = %d (%.1f%%)（ns() は境界の外では線形に外挿する）",
                .out, 100 * .out / nrow(main)))
}
knot_tab <- do.call(rbind, knot_rows)

## --- モデル式に埋め込む項の文字列（§6.4-5：再算出の余地を残さない）-----------
## 数値をそのまま文字列にして式に埋め込む。05・06 はこの文字列を使う。
num_txt <- function(x) paste0("c(", paste(vapply(x, function(z)
  sprintf("%.17g", z), character(1)), collapse = ", "), ")")

spline_terms <- vapply(names(knot_list), function(rv) {
  k <- knot_list[[rv]]
  sprintf("ns(%s, knots = %s, Boundary.knots = %s)",
          rv, num_txt(k$inner), num_txt(k$boundary))
}, character(1))
names(spline_terms) <- names(knot_list)

say("\n  spline terms passed to the model formula:")
for (rv in names(spline_terms)) say("    ", spline_terms[[rv]])

## 基底が作れることをここで確かめる（05・06 で初めて落ちるのを避ける）
if (has_splines) {
  for (rv in names(knot_list)) {
    k  <- knot_list[[rv]]
    xx <- main[[rv]][!is.na(main[[rv]])]     # 欠測は 04 の完全ケース判定で扱う
    B <- tryCatch(splines::ns(xx, knots = k$inner, Boundary.knots = k$boundary),
                  error = function(e) e)
    ok <- !inherits(B, "error") && is.matrix(B) && all(is.finite(B))
    check(sprintf("%s のスプライン基底が作れる（df = %s）", rv,
                  if (ok) ncol(B) else "?"),
          ok, if (ok) sprintf("ncol = %d, rank = %d", ncol(B), qr(B)$rank)
              else paste("error:", conditionMessage(B)))
    if (ok) {
      check(sprintf("%s の基底が列フルランクである", rv), qr(B)$rank == ncol(B),
            sprintf("rank = %d / ncol = %d", qr(B)$rank, ncol(B)))
    }
  }
}


# -----------------------------------------------------------------------------
# 9. M2 の解析データ（完全ケース集団）
# -----------------------------------------------------------------------------

rule("9. M2 analysis data (complete-case population)")

## アウトカムは順序因子にする。水準は観測された点数を昇順に並べたもの。
make_outcome_factor <- function(y, levels_num = NULL) {
  if (is.null(levels_num)) levels_num <- sort(unique(y[!is.na(y)]))
  factor(as.character(y), levels = as.character(levels_num), ordered = TRUE)
}

keep <- c("id", "age", "age_group", "sex", "class", "support_in",
          "fy_in", "period", "mFIM_in", "cFIM_in",
          unname(RIDIT_VARS), OUTCOME)
dat_m2 <- main[complete, intersect(keep, names(main)), drop = FALSE]
dat_m2$y_ord <- make_outcome_factor(dat_m2[[OUTCOME]])
Y_LEVELS <- as.numeric(levels(dat_m2$y_ord))

n_m2_by_g <- table(dat_m2$age_group)
say(sprintf("  M2 complete-case population : %d", nrow(dat_m2)))
for (g in AGE_LEVELS) {
  say(sprintf("    %-3s %-14s n = %5d", g, relabel_levels("age_group", g),
              as.integer(n_m2_by_g[[g]])))
}
say(sprintf("  observed outcome levels     : %d (range %g - %g)",
            length(Y_LEVELS), min(Y_LEVELS), max(Y_LEVELS)))

check("完全ケース集団に欠測がない",
      !anyNA(dat_m2[, c(M2_VARS, unname(RIDIT_VARS), OUTCOME)]), "")
check("完全ケース集団の id が一意（解析単位＝患者、§5.2）",
      !anyDuplicated(dat_m2$id), "")
check("アウトカムが値域に収まる",
      all(Y_LEVELS >= MFIM_RANGE[1] & Y_LEVELS <= MFIM_RANGE[2]),
      sprintf("range = %g - %g", min(Y_LEVELS), max(Y_LEVELS)))
check("アウトカムの順序因子の水準が数値として昇順である",
      !is.unsorted(Y_LEVELS), "")
check("完全ケース集団の例数が n_cc と一致する", nrow(dat_m2) == n_cc,
      sprintf("%d vs %d", nrow(dat_m2), n_cc))

## 解析集団のNの記録（Table 2 と Fig. 2 の脚注に必ず載せる。§6.5「報告の規則」）
analysis_n <- data.frame(
  set   = c("analysis_set_E1_E2", rep("m2_complete_case", 1 + length(AGE_LEVELS))),
  group = c("total", "total", AGE_LEVELS),
  label = c("Analysis set (used by E-1 and E-2)",
            "M2 complete-case population (used by E-3 and E-4)",
            relabel_levels("age_group", AGE_LEVELS)),
  n     = c(n_main, nrow(dat_m2), as.integer(n_m2_by_g[AGE_LEVELS])),
  stringsAsFactors = FALSE
)
say("\n  N by analysis set:")
print(analysis_n, row.names = FALSE)


# -----------------------------------------------------------------------------
# 10. 出来上がりの点検
# -----------------------------------------------------------------------------

rule("10. Verification")

check("ridit とノットは基準集団（A4）だけから決めている（§6.4-5）", TRUE,
      sprintf("A4 N = %d。05・06 は再算出せず、この rda の数値を使う", nrow(base)))
check("ノットの数と `ns()` の自由度の対応",
      all(vapply(names(knot_list), function(rv)
        length(knot_list[[rv]]$inner) == knot_list[[rv]]$n_knots - 2L, logical(1))),
      paste(vapply(names(knot_list), function(rv)
        sprintf("%s: knots=%d, inner=%d", rv, knot_list[[rv]]$n_knots,
                length(knot_list[[rv]]$inner)), character(1)), collapse = ", "))
check("スプラインの項の文字列にノットの数値が埋め込まれている",
      all(grepl("Boundary.knots = c\\(", spline_terms)), "")

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
w(miss_tab,      "table_missing_by_agegroup.csv")
w(miss_decision, "table_missing_decision.csv")
w(knot_tab,      "table_ridit_knots.csv")
w(ridit_summary, "table_ridit_summary.csv")
w(period_by_g,   "table_period_by_agegroup.csv")
w(cov_levels,    "table_covariate_levels.csv")
w(analysis_n,    "table_analysis_n.csv")
w(checks,        "table_checks_04.csv")

## 05・06 が読むもの。設定値も一緒に保存し、後続が再定義しないようにする。
M2_SETTINGS <- list(
  outcome      = OUTCOME,
  exposure     = EXPOSURE,
  covars_cat   = COVARS_CAT,
  covars_num   = COVARS_NUM,
  ridit_vars   = RIDIT_VARS,
  age_levels   = AGE_LEVELS,
  age_breaks   = AGE_BREAKS,
  period_levels = PERIOD_LEVELS,
  ref_levels   = REF_LEVELS,
  y_levels     = Y_LEVELS,
  quantile_type = QUANTILE_TYPE
)

save(dat_m2, base, main, complete,
     knot_list, knot_tab, spline_terms, ridit_fun, ridit_meta, ridit_summary,
     miss_tab, miss_decision, decision_complete_case,
     cov_levels, period_by_g, analysis_n, M2_SETTINGS, checks,
     file = file.path(OUT_DIR, "04_ridit_spline.rda"))
say("  written: ", file.path(OUT_DIR, "04_ridit_spline.rda"))

say("\nDone. M2 complete-case population = ", nrow(dat_m2),
    " / analysis set = ", n_main,
    " ; dropout = ", sprintf("%.2f%%", 100 * dropout))

.log_close()
