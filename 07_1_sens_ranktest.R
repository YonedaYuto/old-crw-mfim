# =============================================================================
# 07_1_sens_ranktest.R
#   解析計画書 §8.8「感度分析S」のうち、**E-2（順位にもとづく比較）**の実装。
#   plan.summary.txt §6「感度分析」。主解析の 03_m0_ranktest.R に対応する。
#
#   入力 : output/02_preprocess.rda    dat_base / dat_sens（mFIM_out_S を含む）
#          output/03_m0_ranktest.rda   m0_relative_effect（主解析との対照。無くても可）
#          01_labels.R（表示用の文字列）
#   出力 : output/table_sens_m0_relative_effect.csv  S の E-2（2群相対効果と95%区間）
#          output/table_sens_perm_meta.csv           順列の諸元（回数・乱数種・所要時間）
#          output/table_sens_m0_compare.csv          主解析との対照（ログ・記録用）
#          output/table_sens_worst_score.csv         最低点（12）と実測の下限（13）の例数（年齢群別）
#          output/table_checks_07_1.csv              点検結果
#          output/07_1_sens_ranktest.rda             07_2 が読む（補足表2 の E-2 行）
#          output/log_07_1_sens_ranktest.txt         実行ログ
#
#   改訂 : 2026-09-22 plan.summary.txt §7-3 により、転院・死亡の患者に与える
#          最低点を 13 点から 12 点に変えた（02 の WORST_SCORE の変更に合わせる）。
#          ・WORST_SCORE を 12 にした。実測の運動FIM の値域（MFIM_RANGE = 13〜91）は
#            そのまま残し、感度分析Sのアウトカムの値域を MFIM_RANGE_S = 12〜91 として
#            別に置いた。値域の点検は MFIM_RANGE_S で行い、12 点を与えた患者と
#            実測値の患者が値として重ならないことも確かめる。
#          ・n_worst（最低点を与えた例数）は、12 点の例数を数える。12 点は実測では
#            ありえない値であるため、実測で 13 点だった患者は数えられない。
#            02 と同じ規則（退院先が転院・終了（死亡等）、または E3 の区分(a)）で
#            数えた例数と一致することを点検する。
#          ・転院・死亡の患者が実測で 13 点の患者より下の順位になり、両者の同順位が
#            なくなる（worst-rank 法、Lachin 1999。査読 W2 の要求に対応）。
#            そのため推定値は旧版と変わる。Brunner–Munzel 検定の手順と乱数種は
#            変えていない。
#
#   計画との対応（§8.8、§8.2）
#     ・対象 : 基準集団（A4）から E3 の区分(b)(c) を除いた集団（02 の dat_sens）。
#       E1（転院）・E2（終了（死亡等））・E3(a) には運動FIM 12点が入っている
#       （§7-3。旧版は 13 点）
#     ・**E-2 を主解析と同じ仕様で繰り返す。** `rankFD::rank.two.samples()`
#       （method = "logit"、permu = TRUE、nperm = 100000、conf.level = 0.95、
#       alternative = "two.sided"）。高齢側を第2水準に置く
#     ・自前の順列分布も生成し、等裾区間の再現と |T*| 区間との差をログに残す（§8.2）
#     ・**検証的な判定には用いず、主対比のp値も示さない**（§8.5、§8.8）。
#       順列p値は算出してログと rda に残すが、出力表の p_value 列は NA にする
#     ・大域検定（ATS）は本文にも補足資料にも報告しない（§8.2、§8.5）。
#       感度分析Sでは算出もしない
#     ・E-1（分布の記述。Fig. 2(a) の階段関数）は §8.8 が繰り返すと定めていない
#       ため、ここでは再掲しない
#
#   乱数種（著者の決定、2026-09-20）
#     主解析（SEED = 20260920）と別の種を割り当てる。SEED_S = SEED + 1000。
#     対比ごとの割り当て方は 03_m0_ranktest.R と同じ。
#
#   実行時間の目安
#     自前の順列は N ≈ 1,500 の対比で 100,000 回あたり約30秒。3対比で約1.5分。
#     rankFD の順列を加えても数分の範囲に収まる見込み。実測値はログに出る。
#
#   文字コード : UTF-8（Windows の R では source(..., encoding = "UTF-8")）
# =============================================================================


# -----------------------------------------------------------------------------
# 0. 設定
# -----------------------------------------------------------------------------

OUT_DIR   <- "output"
IN_RDA_02 <- file.path(OUT_DIR, "02_preprocess.rda")
IN_RDA_03 <- file.path(OUT_DIR, "03_m0_ranktest.rda")

CONF_LEVEL <- 0.95          # 未調整95%（§8.5）
MFIM_RANGE <- c(13, 91)     # 実測の運動FIM の値域（13 項目 × 1〜7 点）

## --- 感度分析Sのアウトカム（§8.8、§7-3）------------------------------------
## 02_preprocess.R が作った複合アウトカム。E1・E2・E3(a) には 12 点が入っている。
## §7-3（2026-09-22）：退院時に実測で 13 点だった患者と区別するため、最低点を
## 実測の下限（13 点）より1点低い 12 点とする（旧設定は 13）。02 の WORST_SCORE と
## 同じ値であること。
OUTCOME_S    <- "mFIM_out_S"
WORST_SCORE  <- 12
MFIM_RANGE_S <- c(WORST_SCORE, MFIM_RANGE[2])   # 感度分析Sのアウトカムの値域（12〜91）

## --- 最低点を与える患者（02_preprocess.R の DISPO_TRANSFER・DISPO_DEATH と同じ値）--
## n_worst の点検に使う。02 の規則：退院先が E1（転院）・E2（終了（死亡等））、
## または E3 の区分(a)（急変・状態悪化による未評価）。
DISPO_WORST <- c("病院・診療所へ転院", "医療機関", "終了（死亡等）")
E3_CATEGORY_WORST <- "a"

## --- 年齢群（§5.3）。02〜06 と同じ値であること ------------------------------
AGE_BREAKS <- c(65, 75, 90)           # G1: 65-74, G2: 75-89, G3: >=90
AGE_LEVELS <- c("G1", "G2", "G3")

## --- 順列 Brunner–Munzel（§8.2、§8.10）--------------------------------------
NPERM       <- 100000L                # rankFD に渡す順列回数
NPERM_OWN   <- NPERM                  # 自前順列の回数（比較のため同数）
BM_METHOD   <- "logit"                # 範囲保存のため（§8.2）
ALTERNATIVE <- "two.sided"
SEED_S      <- 20260920L + 1000L      # 感度分析Sの順列（主解析の SEED と別にする）
## 対比ごとの乱数種。rankFD は大域の乱数状態を使うため、呼び出しの直前に set.seed する。
seed_own_S    <- function(i) SEED_S + 10L * i
seed_rankfd_S <- function(i) SEED_S + 10L * i + 1L

## --- 対比（§8.2）。高齢側を第2水準に置く -----------------------------------
CONTRASTS <- list(
  list(key = "G2_vs_G1", ref = "G1", cmp = "G2", role = "secondary"),
  list(key = "G3_vs_G1", ref = "G1", cmp = "G3", role = "primary"),
  list(key = "G3_vs_G2", ref = "G2", cmp = "G3", role = "secondary")
)
PRIMARY_KEY <- "G3_vs_G1"

## §8.8「報告」：感度分析Sでは主対比のp値も示さない。
## 順列p値は算出してログと rda には残すが、出力表の p_value 列は NA のままにする。
REPORT_P_VALUE_IN_S <- FALSE


# -----------------------------------------------------------------------------
# 1. 表示用ラベル
# -----------------------------------------------------------------------------

if (!file.exists("01_labels.R")) {
  stop("01_labels.R が見つからない。07_1_sens_ranktest.R は 01_labels.R と同じ",
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

.log_con <- file(file.path(OUT_DIR, "log_07_1_sens_ranktest.txt"), open = "wt",
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

rule(paste0("07_1_sens_ranktest.R  ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
say("R : ", R.version.string)
say("permutation seed (base) : ", SEED_S,
    "  [own = SEED_S + 10*i, rankFD = SEED_S + 10*i + 1 ; i は対比の番号]")
say("nperm (rankFD) : ", NPERM, " ; nperm (own) : ", NPERM_OWN)
say("outcome : ", OUTCOME_S, "  （離脱事象と状態悪化による未評価に ", WORST_SCORE,
    " 点を与えた複合アウトカム）")

has_rankFD <- requireNamespace("rankFD", quietly = TRUE)
check("rankFD が利用できる", has_rankFD,
      if (has_rankFD) paste0("version ",
                             as.character(utils::packageVersion("rankFD")))
      else "install.packages(\"rankFD\") が必要。自前の順列で代替する")


# -----------------------------------------------------------------------------
# 3. 入力の取り込み
# -----------------------------------------------------------------------------

rule("3. Input")

## rda は必ず専用の環境に読み込み、必要なオブジェクトだけを取り出す。
## 直接 globalenv() に読み込むと、他のスクリプトの rda に入っている CONTRASTS・
## CONF_LEVEL などがこのスクリプトの設定（§0）を黙って上書きしてしまう。
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
  if (length(miss)) {
    stop(path, " に必要なオブジェクトがない: ", paste(miss, collapse = ", "))
  }
  for (n in nm) assign(n, get(n, envir = e), envir = globalenv())
  invisible(NULL)
}

e02 <- load_env(IN_RDA_02)
take(e02, c("dat_base", "dat_sens"), IN_RDA_02)

## 主解析の E-2（対照としてログと記録に残す。無くても S 単独の表は作れる）
main_m0 <- load_env(IN_RDA_03, required = FALSE)
if (is.null(main_m0)) {
  say("  主解析との対照表は作れない（03_m0_ranktest.R を先に実行すること）")
}

sens <- as.data.frame(dat_sens, stringsAsFactors = FALSE)
say("\ndat_sens (sensitivity S population, A4 minus E3(b)(c)) : ", nrow(sens), " rows")

.need <- c("id", "age", OUTCOME_S)
.miss <- setdiff(.need, names(sens))
if (length(.miss)) stop("dat_sens に必要な列がない: ", paste(.miss, collapse = ", "))

check("Sのアウトカムに欠測がない（§8.8：複合アウトカム）",
      !anyNA(sens[[OUTCOME_S]]),
      sprintf("missing = %d", sum(is.na(sens[[OUTCOME_S]]))))
check(sprintf("Sのアウトカムが値域 %d〜%d に収まる（§7-3：最低点 %d を含む）",
              MFIM_RANGE_S[1], MFIM_RANGE_S[2], WORST_SCORE),
      all(sens[[OUTCOME_S]] >= MFIM_RANGE_S[1] & sens[[OUTCOME_S]] <= MFIM_RANGE_S[2]),
      sprintf("range = %g - %g", min(sens[[OUTCOME_S]]), max(sens[[OUTCOME_S]])))
check(sprintf("最低点（%d）が実測の値域（%d〜%d）の外にある（§7-3：実測の %d 点と区別する）",
              WORST_SCORE, MFIM_RANGE[1], MFIM_RANGE[2], MFIM_RANGE[1]),
      WORST_SCORE < MFIM_RANGE[1],
      sprintf("WORST_SCORE = %d", WORST_SCORE))
## 最低点以外の値は実測の値域に収まる（最低点と実測値が重ならない）
.obs_S <- sens[[OUTCOME_S]][sens[[OUTCOME_S]] != WORST_SCORE]
check("最低点を除く値が実測の値域に収まる",
      all(.obs_S >= MFIM_RANGE[1] & .obs_S <= MFIM_RANGE[2]),
      if (length(.obs_S)) sprintf("range = %g - %g", min(.obs_S), max(.obs_S)) else "no data")
if (exists("DISPO_KNOWN", inherits = TRUE)) {
  check("最低点を与える退院先が 01_labels.R の DISPO_KNOWN と一致する",
        setequal(DISPO_WORST, DISPO_KNOWN),
        sprintf("07_1: %s / 01: %s", paste(DISPO_WORST, collapse = ", "),
                paste(DISPO_KNOWN, collapse = ", ")))
}
check("解析単位が患者である（id が一意。§5.2）", !anyDuplicated(sens$id), "")


# -----------------------------------------------------------------------------
# 4. 年齢群（§5.3）
#    02〜06 と同じ規則をここに再掲し、02 の出力（flow_g）と突き合わせる。
# -----------------------------------------------------------------------------

rule("4. Age group")

make_agegroup <- function(age) {
  factor(ifelse(is.na(age), NA_character_,
         ifelse(age >= AGE_BREAKS[3], "G3",
         ifelse(age >= AGE_BREAKS[2], "G2", "G1"))),
         levels = AGE_LEVELS)
}
sens$age_group <- make_agegroup(sens$age)

n_sens_by_g <- table(sens$age_group)
for (g in AGE_LEVELS) {
  say(sprintf("  %-3s %-14s n = %5d", g, relabel_levels("age_group", g),
              as.integer(n_sens_by_g[[g]])))
}
check("年齢群に欠測がない", !anyNA(sens$age_group),
      sprintf("missing = %d", sum(is.na(sens$age_group))))

if (exists("flow_g", envir = e02, inherits = FALSE)) {
  .fg <- get("flow_g", envir = e02)
  if ("stage" %in% names(.fg) && "SENS" %in% .fg$stage) {
    .row <- .fg[.fg$stage == "SENS", , drop = FALSE][1, ]
    .exp <- as.integer(c(.row$G1, .row$G2, .row$G3))
    check("年齢群の例数が 02_preprocess.R の flow_g（SENS）と一致する",
          identical(as.integer(n_sens_by_g), .exp),
          sprintf("07_1: %s / 02: %s",
                  paste(as.integer(n_sens_by_g), collapse = ", "),
                  paste(.exp, collapse = ", ")))
  }
}

## 最低点を与えられた例の数（補足表2 の脚注と §8.8 の記述の材料）
## §7-3：最低点は 12 点で、実測の値域（13〜91）の外にある。したがって 12 点の例を
## 数えれば、実測で 13 点だった患者を含めずに、最低点を与えた例だけを数えられる。
is_worst <- sens[[OUTCOME_S]] == WORST_SCORE
n_worst <- sum(is_worst)
n_worst_by_g <- table(factor(as.character(sens$age_group[is_worst]), levels = AGE_LEVELS))

## 02 と同じ規則で最低点を与える患者を数え直し、12 点の例と一致することを確かめる
if (all(c("disposition", "e3_category") %in% names(sens))) {
  .rule_worst <- (sens$disposition %in% DISPO_WORST) |
                 (!is.na(sens$e3_category) & sens$e3_category %in% E3_CATEGORY_WORST)
  check(sprintf("%d 点の例が、02 の規則で最低点を与える患者（転院・終了（死亡等）・E3(a)）と一致する（§7-3）",
                WORST_SCORE),
        identical(as.logical(is_worst), as.logical(.rule_worst)),
        sprintf("%d 点 = %d, 規則 = %d, 食い違い = %d", WORST_SCORE, n_worst,
                sum(.rule_worst), sum(is_worst != .rule_worst)))
} else {
  check(sprintf("%d 点の例が、02 の規則で最低点を与える患者と一致する（§7-3）", WORST_SCORE),
        FALSE, "dat_sens に disposition・e3_category の列がなく、照合できなかった")
}
## 実測で 13 点だった患者（最低点とは別に数える。旧版ではここが最低点と同じ値だった）
n_obs_floor <- sum(sens[[OUTCOME_S]] == MFIM_RANGE[1])
say(sprintf("\n  observations at the worst score (%d points) : %d",
            WORST_SCORE, n_worst))
for (g in AGE_LEVELS) {
  say(sprintf("    %-3s n = %4d (%.1f%% of the group)", g,
              as.integer(n_worst_by_g[[g]]),
              100 * as.integer(n_worst_by_g[[g]]) / as.integer(n_sens_by_g[[g]])))
}
say(sprintf("  observed at the measured floor (%d points, not assigned) : %d",
            MFIM_RANGE[1], n_obs_floor))
say("  高齢群ほど最低点が多い場合、主解析との差は「挟み幅」ではなく同方向の",
    "ずれの大きさになる（§8.8）。")

worst_tab <- data.frame(
  age_group       = c(AGE_LEVELS, "total"),
  age_group_label = c(relabel_levels("age_group", AGE_LEVELS), "Total"),
  n               = c(as.integer(n_sens_by_g[AGE_LEVELS]), nrow(sens)),
  n_worst_score   = c(as.integer(n_worst_by_g[AGE_LEVELS]), n_worst),
  stringsAsFactors = FALSE
)
worst_tab$percent_worst_score <- 100 * worst_tab$n_worst_score / worst_tab$n
## §7-3：最低点（12）とは別に、実測で下限（13 点）だった例の数も並べる
.floor_by_g <- table(factor(as.character(sens$age_group[sens[[OUTCOME_S]] == MFIM_RANGE[1]]),
                            levels = AGE_LEVELS))
worst_tab$worst_score           <- WORST_SCORE
worst_tab$n_observed_floor      <- c(as.integer(.floor_by_g[AGE_LEVELS]), n_obs_floor)
worst_tab$observed_floor_score  <- MFIM_RANGE[1]
check("最低点の例数と実測の下限の例数の合計が群別と総数で一致する",
      sum(worst_tab$n_worst_score[worst_tab$age_group != "total"]) == n_worst &&
        sum(worst_tab$n_observed_floor[worst_tab$age_group != "total"]) == n_obs_floor,
      sprintf("worst = %d, observed floor = %d", n_worst, n_obs_floor))


# -----------------------------------------------------------------------------
# 5. Brunner–Munzel の自前実装
#    03_m0_ranktest.R と同一の実装をここに再掲する（03 を読み込むと主解析が
#    再実行されるため）。数式は §8.2 に記した `rank.two.samples()` の仕様。
#      p̂   = P(x < y) + 0.5 × P(x = y)、x は第1水準、y は第2水準
#      se_p = sqrt(n1 σ1² + n2 σ2²) / (n1 n2)
#      f    = logit(p̂)、se_f = se_p / (p̂ (1 − p̂))、T = f / se_f
# -----------------------------------------------------------------------------

logit_f <- function(p) log(p / (1 - p))
expit_f <- function(z) 1 / (1 + exp(-z))

#' 順列p値の表示。順列で観測値以上に極端な並べ替えが1つも出なければ推定値は0に
#' なるが、「P = 0」ではなく分解能 1/nperm を下回ることしか言えない。
format_perm_p <- function(p, nperm) {
  if (!is.finite(p)) return(NA_character_)
  lim <- 1 / nperm
  if (p <= 0) sprintf("P < %s", format(lim, scientific = TRUE, digits = 1))
  else if (p < 0.001) "P < 0.001"
  else sprintf("P = %.3f", p)
}

#' 結合標本の中間順位から Brunner–Munzel の統計量を計算する
bm_logit <- function(ra1, ra2) {
  n1 <- length(ra1); n2 <- length(ra2)
  r1 <- rank(ra1);   r2 <- rank(ra2)
  m1 <- mean(ra1);   m2 <- mean(ra2)
  p  <- (m2 - (n2 + 1) / 2) / n1
  v1 <- sum((ra1 - r1 - m1 + (n1 + 1) / 2)^2) / (n1 - 1)
  v2 <- sum((ra2 - r2 - m2 + (n2 + 1) / 2)^2) / (n2 - 1)
  se_p <- sqrt(n1 * v1 + n2 * v2) / (n1 * n2)
  if (!is.finite(p) || p <= 0 || p >= 1 || !is.finite(se_p) || se_p <= 0) {
    return(c(p = p, se_p = se_p, f = NA_real_, se_f = NA_real_, T = NA_real_))
  }
  f    <- logit_f(p)
  se_f <- se_p / (p * (1 - p))
  c(p = p, se_p = se_p, f = f, se_f = se_f, T = f / se_f)
}

#' 順列分布 T* を作る（応答を並べ替え、群の大きさは固定）
bm_permute <- function(ra, n1, nperm, seed) {
  N  <- length(ra)
  i1 <- seq_len(n1)
  i2 <- (n1 + 1L):N
  set.seed(seed)
  Tp <- numeric(nperm)
  for (b in seq_len(nperm)) {
    s <- sample.int(N)
    Tp[b] <- bm_logit(ra[s[i1]], ra[s[i2]])[["T"]]
  }
  Tp
}

#' 観測統計量と順列分布から、等裾区間と |T*| 対称区間を作る
bm_intervals <- function(obs, Tp, conf.level = CONF_LEVEL) {
  alpha <- 1 - conf.level
  Tp_ok <- Tp[is.finite(Tp)]
  n_bad <- length(Tp) - length(Tp_ok)
  crit1 <- stats::quantile(Tp_ok, 1 - alpha / 2, names = FALSE)
  crit2 <- stats::quantile(Tp_ok,     alpha / 2, names = FALSE)
  lo_et <- expit_f(obs[["f"]] - crit1 * obs[["se_f"]])
  hi_et <- expit_f(obs[["f"]] - crit2 * obs[["se_f"]])
  p_et  <- min(1, 2 * min(mean(Tp_ok <= obs[["T"]]), mean(Tp_ok >= obs[["T"]])))
  cabs   <- stats::quantile(abs(Tp_ok), conf.level, names = FALSE)
  lo_abs <- expit_f(obs[["f"]] - cabs * obs[["se_f"]])
  hi_abs <- expit_f(obs[["f"]] + cabs * obs[["se_f"]])
  p_abs  <- mean(abs(Tp_ok) >= abs(obs[["T"]]))
  list(n_perm_used = length(Tp_ok), n_perm_dropped = n_bad,
       crit_lower = crit2, crit_upper = crit1, crit_abs = cabs,
       et = c(lower = lo_et, upper = hi_et, p_value = p_et),
       abs = c(lower = lo_abs, upper = hi_abs, p_value = p_abs))
}


# -----------------------------------------------------------------------------
# 6. rankFD の返り値から相対効果と区間を取り出す
#    返り値の構造はパッケージの版で変わりうるため、data.frame を総当たりで探し、
#    順列法の行を拾う。拾えなかった場合は str() をログに出し、自前の再現値に
#    切り替えたうえで REVIEW を立てる（黙って別の値を載せないため）。
# -----------------------------------------------------------------------------

.collect_dfs <- function(x) {
  out <- list()
  walk <- function(obj, path) {
    if (is.data.frame(obj)) { out[[path]] <<- obj; return(invisible(NULL)) }
    if (is.matrix(obj) && length(dim(obj)) == 2L && !is.null(colnames(obj))) {
      out[[path]] <<- as.data.frame(obj, stringsAsFactors = FALSE)
      return(invisible(NULL))
    }
    if (is.list(obj)) {
      nm <- names(obj)
      if (is.null(nm) || !length(nm)) nm <- rep("", length(obj))
      nm[!nzchar(nm)] <- paste0("[[", which(!nzchar(nm)), "]]")
      for (i in seq_along(obj)) walk(obj[[i]], paste(path, nm[i], sep = "$"))
    }
    invisible(NULL)
  }
  walk(x, "res")
  out
}
.find_col <- function(df, patterns) {
  nm <- names(df)
  for (p in patterns) {
    hit <- grep(p, nm, ignore.case = TRUE)
    if (length(hit)) return(nm[hit[1]])
  }
  NA_character_
}
.as_num <- function(x) suppressWarnings(as.numeric(as.character(x)))

#' 相対効果・下限・上限・p値をもつ候補行をすべて並べる
extract_two_samples <- function(res) {
  dfs <- .collect_dfs(res)
  rows <- list()
  for (path in names(dfs)) {
    df <- dfs[[path]]
    if (!nrow(df)) next
    c_est <- .find_col(df, c("^Estimator$", "Estimator", "Rel.?Effect",
                             "Relative", "^p$", "Effect"))
    c_lo  <- .find_col(df, c("^Lower$", "Lower", "L\\.?Bound", "CI.?low", "^LCL$"))
    c_hi  <- .find_col(df, c("^Upper$", "Upper", "U\\.?Bound", "CI.?up", "^UCL$"))
    if (is.na(c_est) || is.na(c_lo) || is.na(c_hi)) next
    c_p   <- .find_col(df, c("p\\.?value", "pValue", "^p\\.?val"))
    c_m   <- .find_col(df, c("^Method$", "Method", "^Type$", "Statistic.?name"))
    rn    <- rownames(df); if (is.null(rn)) rn <- rep("", nrow(df))
    for (i in seq_len(nrow(df))) {
      tag <- paste(path, rn[i],
                   if (!is.na(c_m)) as.character(df[[c_m]][i]) else "")
      rows[[length(rows) + 1L]] <- data.frame(
        path      = path,
        row_label = trimws(gsub("\\s+", " ", tag)),
        estimator = .as_num(df[[c_est]][i]),
        lower     = .as_num(df[[c_lo]][i]),
        upper     = .as_num(df[[c_hi]][i]),
        p_value   = if (is.na(c_p)) NA_real_ else .as_num(df[[c_p]][i]),
        stringsAsFactors = FALSE
      )
    }
  }
  if (!length(rows)) {
    return(data.frame(path = character(0), row_label = character(0),
                      estimator = numeric(0), lower = numeric(0),
                      upper = numeric(0), p_value = numeric(0),
                      stringsAsFactors = FALSE))
  }
  do.call(rbind, rows)
}

#' 候補行から「順列法・logit」の行を選ぶ
#' @param p_own 自前に計算した相対効果。これと一致する行だけを候補とすることで、
#'              シフト効果の区間など別の量の行を取り違えない。
pick_permutation_row <- function(cand, p_own) {
  if (!nrow(cand)) return(NULL)
  ok <- is.finite(cand$estimator) & is.finite(cand$lower) & is.finite(cand$upper) &
        cand$lower >= 0 & cand$upper <= 1 & cand$lower <= cand$upper &
        cand$lower <= cand$estimator & cand$estimator <= cand$upper &
        is.finite(p_own) & abs(cand$estimator - p_own) < 1e-6
  cand <- cand[ok, , drop = FALSE]
  if (!nrow(cand)) return(NULL)
  is_perm  <- grepl("perm", cand$row_label, ignore.case = TRUE)
  is_logit <- grepl("logit", cand$row_label, ignore.case = TRUE)
  for (sel in list(is_perm & is_logit, is_perm)) {
    if (any(sel)) return(cand[which(sel)[1], , drop = FALSE])
  }
  NULL
}


# -----------------------------------------------------------------------------
# 7. E-2：3対比の順列 Brunner–Munzel（感度分析S）
#    主解析（03）と同じ仕様。対象集団だけが dat_sens・複合アウトカムに変わる。
# -----------------------------------------------------------------------------

rule("7. E-2 : permuted Brunner-Munzel on the sensitivity S population")
say("  estimand : p = P(older > younger) + 0.5 x P(tie) for the composite outcome")
say("  §8.5・§8.8 : p 値は補足表2 に載せない。ログと rda にのみ残す。")

sens_res_rows  <- list()
sens_perm_meta <- list()

for (ict in seq_along(CONTRASTS)) {
  ct  <- CONTRASTS[[ict]]
  key <- ct$key; ref <- ct$ref; cmp <- ct$cmp
  lab <- contrast_label(ref, cmp)
  say("\n", strrep("=", 74))
  say(sprintf("Contrast %s : %s  [sensitivity S]", key, lab))

  ## --- 因子水準の順序をログに出す（§8.2：高齢側を第2水準に置く）-------------
  lv <- c(ref, cmp)
  say(sprintf("  factor levels passed to rank.two.samples() : c(\"%s\", \"%s\")",
              lv[1], lv[2]))
  say(sprintf("    level 1 (x, younger) = %s = %s", lv[1],
              relabel_levels("age_group", lv[1])))
  say(sprintf("    level 2 (y, older)   = %s = %s", lv[2],
              relabel_levels("age_group", lv[2])))
  say("    estimand : p = P(x < y) + 0.5 x P(x = y) = P(older > younger) + 0.5 P(tie)")

  sub <- sens[sens$age_group %in% lv, , drop = FALSE]
  sub$age_group <- factor(as.character(sub$age_group), levels = lv)
  sub <- sub[order(sub$age_group), , drop = FALSE]      # 第1水準を前に並べる
  n1  <- sum(sub$age_group == lv[1]); n2 <- sum(sub$age_group == lv[2])
  say(sprintf("  n(%s) = %d ; n(%s) = %d ; N = %d", lv[1], n1, lv[2], n2, n1 + n2))

  ## --- 自前の観測統計量と順列分布 -------------------------------------------
  ra  <- rank(sub[[OUTCOME_S]])                 # 結合標本の中間順位
  obs <- bm_logit(ra[seq_len(n1)], ra[(n1 + 1L):(n1 + n2)])
  say(sprintf("  own estimate : p.hat = %.6f ; se(p) = %.6f ; logit f = %.6f ; T = %.4f",
              obs[["p"]], obs[["se_p"]], obs[["f"]], obs[["T"]]))

  t0  <- proc.time()[["elapsed"]]
  Tp  <- bm_permute(ra, n1, NPERM_OWN, seed_own_S(ict))
  ivs <- bm_intervals(obs, Tp, CONF_LEVEL)
  t_own <- proc.time()[["elapsed"]] - t0
  say(sprintf("  own permutation : %d draws in %.1f s, seed = %d (dropped non-finite: %d)",
              NPERM_OWN, t_own, seed_own_S(ict), ivs$n_perm_dropped))
  say(sprintf("    equal-tailed crits : q(%.3f) = %.4f ; q(%.3f) = %.4f",
              (1 - CONF_LEVEL) / 2, ivs$crit_lower,
              1 - (1 - CONF_LEVEL) / 2, ivs$crit_upper))
  say(sprintf("    |T*| crit          : q(%.2f) = %.4f", CONF_LEVEL, ivs$crit_abs))

  ## --- rankFD（補足表2 に載せる値）------------------------------------------
  rf_est <- rf_lo <- rf_hi <- rf_p <- NA_real_
  rf_ok  <- FALSE
  t_rf   <- NA_real_
  if (has_rankFD) {
    .d2 <- data.frame(y = sub[[OUTCOME_S]], age_group = sub$age_group)
    say(sprintf("  rankFD seed : %d", seed_rankfd_S(ict)))
    set.seed(seed_rankfd_S(ict))
    t0  <- proc.time()[["elapsed"]]
    res <- tryCatch({
      call_known_args(rankFD::rank.two.samples,
                      list(formula = y ~ age_group, data = .d2,
                           conf.level = CONF_LEVEL, alternative = ALTERNATIVE,
                           method = BM_METHOD, permu = TRUE, nperm = NPERM,
                           info = FALSE, plot.simci = FALSE))
    }, error = function(e) {
      say("  rankFD::rank.two.samples() でエラー: ", conditionMessage(e)); NULL
    })
    t_rf <- proc.time()[["elapsed"]] - t0
    if (!is.null(res)) {
      cand <- extract_two_samples(res)
      pick <- pick_permutation_row(cand, obs[["p"]])
      if (!is.null(pick)) {
        rf_est <- pick$estimator; rf_lo <- pick$lower
        rf_hi  <- pick$upper;     rf_p  <- pick$p_value
        rf_ok  <- TRUE
        say(sprintf("  picked permutation row : %s", pick$row_label))
      } else {
        say("\n  --- rankFD::rank.two.samples() structure (log only) ---")
        utils::str(res, max.level = 3, give.attr = FALSE)
        if (nrow(cand)) {
          say("\n  candidate rows found in the return value:")
          print(cand, row.names = FALSE, digits = 6)
        }
        say("  順列法の行を特定できなかった。上の構造を確認すること。")
      }
    }
  }

  ## --- 補足表2 に載せる値を決める --------------------------------------------
  if (rf_ok) {
    src <- "rankFD::rank.two.samples (permutation, logit)"
    est <- rf_est; lo <- rf_lo; hi <- rf_hi
    pv  <- if (is.finite(rf_p)) rf_p else ivs$et[["p_value"]]
    if (!is.finite(rf_p)) src <- paste0(src, " + own permutation p-value")
  } else {
    src <- "own reproduction (equal-tailed permutation, logit)"
    est <- obs[["p"]]; lo <- ivs$et[["lower"]]; hi <- ivs$et[["upper"]]
    pv  <- ivs$et[["p_value"]]
  }
  check(sprintf("[S %s] 補足表2 の区間を rankFD から取得できた", key), rf_ok,
        if (rf_ok) "rankFD の順列区間を採用"
        else "自前の再現値に切り替えた。rankFD の返り値の構造をログで確認すること")
  if (rf_ok) {
    d_est <- abs(rf_est - obs[["p"]])
    check(sprintf("[S %s] 採用した行の相対効果が自前の計算と一致する", key),
          is.finite(d_est) && d_est < 1e-6,
          sprintf("rankFD = %.8f, own = %.8f, diff = %.2e", rf_est, obs[["p"]], d_est))
    d_ci <- max(abs(rf_lo - ivs$et[["lower"]]), abs(rf_hi - ivs$et[["upper"]]))
    check(sprintf("[S %s] rankFD の区間と自前の等裾再現の差が 0.01 未満", key),
          is.finite(d_ci) && d_ci < 0.01,
          sprintf("rankFD [%.6f, %.6f] / own [%.6f, %.6f], max diff = %.2e（順列のモンテカルロ誤差の範囲か確認する）",
                  rf_lo, rf_hi, ivs$et[["lower"]], ivs$et[["upper"]], d_ci))
  }
  check(sprintf("[S %s] 区間が (0, 1) に収まる（method = \"logit\" の範囲保存）", key),
        is.finite(lo) && is.finite(hi) && lo > 0 && hi < 1 && lo <= est && est <= hi,
        sprintf("[%.6f, %.6f], p.hat = %.6f", lo, hi, est))

  say(sprintf("  (log only) permutation p-value = %s", format_perm_p(pv, NPERM)))

  sens_res_rows[[length(sens_res_rows) + 1L]] <- data.frame(
    analysis          = "sensitivity S",
    contrast          = key,
    label             = lab,
    group_reference   = ref,
    group_comparison  = cmp,
    reference_label   = relabel_levels("age_group", ref),
    comparison_label  = relabel_levels("age_group", cmp),
    n_reference       = n1,
    n_comparison      = n2,
    role              = ct$role,
    relative_effect   = est,
    ci_lower          = lo,
    ci_upper          = hi,
    conf_level        = CONF_LEVEL,
    ci_method         = "equal-tailed permutation, logit scale, unadjusted",
    ## §8.8「報告」：検証的な判定には用いず、主対比の p 値も示さない
    p_value           = if (REPORT_P_VALUE_IN_S) pv else NA_real_,
    p_value_display   = "",
    p_value_method    = "",
    p_value_log_only  = pv,
    source            = src,
    stringsAsFactors  = FALSE
  )

  sens_perm_meta[[length(sens_perm_meta) + 1L]] <- data.frame(
    contrast = key, label = lab,
    nperm_rankfd = NPERM, nperm_own = NPERM_OWN,
    seed_own = seed_own_S(ict), seed_rankfd = seed_rankfd_S(ict),
    seconds_own = t_own, seconds_rankfd = t_rf,
    n_perm_used = ivs$n_perm_used, n_perm_dropped = ivs$n_perm_dropped,
    p_hat_own = obs[["p"]], p_hat_rankFD = rf_est,
    et_lower_rankFD = rf_lo, et_upper_rankFD = rf_hi,
    et_lower_own = ivs$et[["lower"]], et_upper_own = ivs$et[["upper"]],
    abs_lower_own = ivs$abs[["lower"]], abs_upper_own = ivs$abs[["upper"]],
    p_value_et = ivs$et[["p_value"]], p_value_abs = ivs$abs[["p_value"]],
    crit_lower = ivs$crit_lower, crit_upper = ivs$crit_upper,
    crit_abs = ivs$crit_abs,
    stringsAsFactors = FALSE
  )
}

sens_m0        <- do.call(rbind, sens_res_rows)
sens_perm_meta <- do.call(rbind, sens_perm_meta)

say("\n  E-2 (sensitivity S) : relative effect with unadjusted 95% permutation CI")
print(sens_m0[, c("contrast", "label", "n_reference", "n_comparison", "role",
                  "relative_effect", "ci_lower", "ci_upper")],
      row.names = FALSE, digits = 4)

check("対比が3つある", nrow(sens_m0) == 3L, sprintf("rows = %d", nrow(sens_m0)))
check("感度分析Sの出力表に p 値を入れていない（§8.5、§8.8）",
      all(is.na(sens_m0$p_value)), "順列p値は p_value_log_only 列とログにのみ残す")
check("等裾区間と |T*| 区間の差をログに残した（§8.2）", TRUE,
      "table_sens_perm_meta.csv の et_* と abs_* 列")


# -----------------------------------------------------------------------------
# 8. 主解析との対照（記録用。補足表2 の組み立ては 07_2 が行う）
# -----------------------------------------------------------------------------

rule("8. Comparison with the main analysis (for the record)")

sens_m0_compare <- NULL
if (!is.null(main_m0) && exists("m0_relative_effect", envir = main_m0, inherits = FALSE)) {
  .m <- get("m0_relative_effect", envir = main_m0)
  sens_m0_compare <- data.frame(
    contrast          = sens_m0$contrast,
    label             = sens_m0$label,
    role              = sens_m0$role,
    main_estimate     = .m$relative_effect[match(sens_m0$contrast, .m$contrast)],
    main_ci_lower     = .m$ci_lower[match(sens_m0$contrast, .m$contrast)],
    main_ci_upper     = .m$ci_upper[match(sens_m0$contrast, .m$contrast)],
    sens_estimate     = sens_m0$relative_effect,
    sens_ci_lower     = sens_m0$ci_lower,
    sens_ci_upper     = sens_m0$ci_upper,
    stringsAsFactors  = FALSE
  )
  sens_m0_compare$difference <- sens_m0_compare$sens_estimate -
                                sens_m0_compare$main_estimate
  print(sens_m0_compare[, c("contrast", "main_estimate", "sens_estimate",
                            "difference")], row.names = FALSE, digits = 4)
  say("\n  §8.8：主解析とSは推定対象が違う。両者の間に「真の値」を置かない。",
      "差の向きと大きさを記述するにとどめる。")
  check("主解析との対照表を作れた", TRUE, "")
} else {
  sens_m0_compare <- data.frame(contrast = character(0), stringsAsFactors = FALSE)
  check("主解析との対照表を作れた", FALSE,
        "03_m0_ranktest.rda が無い。補足表2 の主解析列は 07_2 で空になる")
}


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
w(sens_m0,         "table_sens_m0_relative_effect.csv")
w(sens_perm_meta,  "table_sens_perm_meta.csv")
w(worst_tab,       "table_sens_worst_score.csv")
w(sens_m0_compare, "table_sens_m0_compare.csv")
w(checks,          "table_checks_07_1.csv")

## 07_2_sens_clm_boot.R が読む（補足表2 の E-2 行）
save(sens_m0, sens_perm_meta, sens_m0_compare, worst_tab,
     n_worst, n_obs_floor, n_sens_by_g, SEED_S, NPERM, CONF_LEVEL, OUTCOME_S,
     WORST_SCORE, MFIM_RANGE, MFIM_RANGE_S,
     checks,
     file = file.path(OUT_DIR, "07_1_sens_ranktest.rda"))
say("  written: ", file.path(OUT_DIR, "07_1_sens_ranktest.rda"))

say("\n--- sessionInfo ---")
print(utils::sessionInfo())

say("\nDone. sensitivity S population = ", nrow(sens),
    " ; contrasts = ", nrow(sens_m0),
    " ; worst score (", WORST_SCORE, ") assigned = ", n_worst,
    " ; observed at ", MFIM_RANGE[1], " points = ", n_obs_floor)

.log_close()
