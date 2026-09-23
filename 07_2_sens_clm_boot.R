# =============================================================================
# 07_2_sens_clm_boot.R
#   解析計画書 §8.8「感度分析S」のうち、**順序ロジスティック回帰にかかわる部分**
#   （E-3 の当てはめと E-4 の周辺標準化）の実装。
#   plan.summary.txt §6「感度分析」。主解析の 04・05・06 に対応する。
#   あわせて、07_1 の E-2 と並べて**補足表2**を組み立てる。
#
#   入力 : output/02_preprocess.rda      dat_base / dat_sens（mFIM_out_S を含む）
#          output/04_ridit_spline.rda    ridit_fun / spline_terms / M2_SETTINGS /
#                                        cov_levels
#          output/05_m2_clm.rda          M2_FORMULA_TXT / M2_PPO_FORMULA_TXT /
#                                        M2_PPO_NOMINAL_TXT / PPO_THRESHOLDS / PPO_LEVELS
#          output/06_standardize_boot.rda 主解析の E-4（補足表2 の主解析列）
#          output/07_1_sens_ranktest.rda  S の E-2（補足表2 の E-2 行）
#          01_labels.R（表示用の文字列）
#   出力 : output/table_sens_analysis_n.csv        S の集団と完全ケース集団の例数
#          output/table_sens_missing.csv           S の集団の共変量の欠測（年齢群別）
#          output/table_sens_odds_ratio.csv        S の共通オッズ比（比例オッズ。算出するが表に出さない）
#          output/table_sens_ppo_odds_ratio.csv    S の閾値別オッズ比（部分比例オッズ。同上）
#          output/table_sens_fit.csv               当てはめの諸元と収束診断（4本）
#          output/table_sens_standardized_risk.csv S の E-4（年齢群別、70 点・77 点。部分比例オッズ）
#          output/table_sens_risk_difference.csv   S の E-4（3対比のRD、70 点・77 点。同上）
#          output/table_sens_po_standardized.csv   S の比例オッズモデル由来の標準化確率とRD（比較用）
#          output/table_sens_ppo_vs_po.csv         2つのモデルの標準化量の食い違い
#          output/table_sens_or_ci_profile_vs_boot.csv OR の profile 区間とブートストラップ区間
#          output/table_sens_boot_summary.csv      ブートストラップの諸元と失敗回
#          output/table_supp2.csv                  補足表2（主解析とSを並べた表）
#          output/table_supp2_meta.csv             補足表2 の脚注に必要な諸元
#          output/table_supp2_footnotes.csv        補足表2 の脚注（下書き）
#          output/table_checks_07_2.csv            点検結果
#          output/07_2_sens_clm_boot.rda
#          output/log_07_2_sens_clm_boot.txt       実行ログ
#
#   改訂 : 2026-09-22 plan.summary.txt §7-1・§7-2・§7-3 により、次のとおり変更した。
#          §6「解析は主解析と同一」により、主解析（05・06）の変更をそのまま S に適用する。
#          ・§7-3：最低点（WORST_SCORE）を 13 → 12 にした（02 の変更に合わせる）。
#            実測の値域（MFIM_RANGE = 13〜91）は残し、S のアウトカムの値域を
#            MFIM_RANGE_S = 12〜91 として別に置いた。n_worst は 12 点の例数で数える。
#          ・§7-1：閾値を 65 点から 70 点・77 点の2つにした（THRESHOLDS）。
#          ・§7-2：部分比例オッズモデル（M2-PPO：アウトカムは <70・70–76・≥77 の3区分、
#            年齢群を nominal 項）と比例オッズモデル（M2：最低点 12 を含む全水準）の
#            両方を当てはめ、同じ再標本で標準化とブートストラップを行う。補足表2 の
#            標準化確率・RD は M2-PPO 由来とし、M2 由来の値は比較用の表に置く（06 と同じ）。
#          ・閾値 65 の二値ロジスティックモデル（§8.6(2)）を削除した。
#          ・"fast" 経路と plogis(θ_j − x'β) による手計算の点検を、nominal 項を含む形
#            plogis(θ_j + τ_j,g − x'β) に一般化した（06 と同じ関数）。
#          ・E-3（オッズ比）は §8.8 により当てはめるが表に出さない。M2-PPO の閾値別
#            オッズ比は 05 と同じ自前の profile likelihood で区間を求め、ブートストラップ
#            区間とも突き合わせる（06 と同じ手順）。いずれもログと csv にのみ残す。
#          ・補足表2 は、主解析・S とも E-2 と、M2-PPO 由来の 70 点・77 点の標準化確率・
#            RD を並べる形に直した（06 の rda の構造の変更に合わせた。threshold 列を追加）。
#          ・output に残る旧版の table_sens_binary65.csv・table_sens_ordinal_vs_binary.csv
#            は、このスクリプトでは作らない（12 が読まないこと）。
#
#   計画との対応（§8.8）
#     ・対象 : 基準集団（A4）から E3 の区分(b)(c) を除いた集団（02 の dat_sens）。
#       E-4 の標準化の対象は、そのうち M2 の共変量がそろった完全ケース集団とする
#       （主解析で標準化の対象を「解析対象の完全ケース集団」としたのと同じ扱い）
#     ・**E-4（§8.4）を主解析と同じ仕様で繰り返す。** 周辺標準化、年齢群で層化した
#       患者単位のブートストラップ2,000回、パーセンタイル法
#     ・E-3（閾値別オッズ比・共通オッズ比）は当てはめるが表には出さない（reported = FALSE）
#     ・**再実行しないもの** : 閾値別オッズ比の散らばり（§8.6(1)。08_po_thresholds.R）
#     ・報告 : 主解析と並べた表を補足表2 に載せる。p 値は載せない（§8.5）
#
#   ridit とノット（§6.4-5、plan.summary §6 の最終行）
#     ridit の基準集団は主解析と同じく A4 である。04_ridit_spline.R が作った
#     ridit_fun と、ノットの数値が埋め込まれたスプラインの項の文字列
#     （spline_terms / M2_FORMULA_TXT / M2_PPO_FORMULA_TXT）をそのまま使い、再算出しない。
#
#   乱数種（著者の決定、2026-09-20）
#     主解析（BOOT_SEED = 20260920）と別の種を割り当てる。
#     BOOT_SEED_S = BOOT_SEED + 1000。07_1 の SEED_S と同じ値だが、用途も
#     対象データも別であり、同じ乱数列が再利用されることはない。
#     再標本のインデックスは旧版と同じ手順で生成する。
#
#   `ordinal::clm` の符号について（05 と同じ）
#     clm は logit P(Y ≤ j) = θ_j + τ_j,g − x'β と書く（τ は nominal 項。M2 では 0）。
#     位置の係数 β の exp はそのまま Y ≥ y のオッズ比、nominal 項の年齢群の母数 τ は
#     符号を反転した exp(−τ) が Y ≥ t のオッズ比である。
#
#   文字コード : UTF-8（Windows の R では source(..., encoding = "UTF-8")）
# =============================================================================


# -----------------------------------------------------------------------------
# 0. 設定
# -----------------------------------------------------------------------------

OUT_DIR     <- "output"
IN_RDA_02   <- file.path(OUT_DIR, "02_preprocess.rda")
IN_RDA_04   <- file.path(OUT_DIR, "04_ridit_spline.rda")
IN_RDA_05   <- file.path(OUT_DIR, "05_m2_clm.rda")
IN_RDA_06   <- file.path(OUT_DIR, "06_standardize_boot.rda")
IN_RDA_07_1 <- file.path(OUT_DIR, "07_1_sens_ranktest.rda")

CONF_LEVEL <- 0.95          # 未調整95%（§8.5）
## 標準化する到達確率 P(運動FIM ≥ t) の閾値（plan.summary §7-1）。昇順に置く。
## 70 点：入浴と階段昇降以外が見守りレベル（旧設定 65 点を修正）。
## 77 点：入浴と階段昇降以外が修正自立レベル。05 の PPO_THRESHOLDS・06 の THRESHOLDS と
## 一致させる（下の 3 節で点検する）。
THRESHOLDS <- c(70, 77)
MFIM_RANGE <- c(13, 91)     # 実測の運動FIM の値域（13 項目 × 1〜7 点）

## --- 感度分析Sのアウトカム（§8.8、§7-3）------------------------------------
## 02_preprocess.R が作った複合アウトカム。E1・E2・E3(a) には 12 点が入っている。
## §7-3（2026-09-22）：退院時に実測で 13 点だった患者と区別するため、最低点を
## 実測の下限より1点低い 12 点とする（旧設定は 13）。02・07_1 の WORST_SCORE と同じ値。
OUTCOME_S    <- "mFIM_out_S"
WORST_SCORE  <- 12
MFIM_RANGE_S <- c(WORST_SCORE, MFIM_RANGE[2])   # 感度分析Sのアウトカムの値域（12〜91）

## --- 最低点を与える患者（02_preprocess.R の DISPO_TRANSFER・DISPO_DEATH と同じ値）--
DISPO_WORST <- c("病院・診療所へ転院", "医療機関", "終了（死亡等）")
E3_CATEGORY_WORST <- "a"

## --- 年齢群（§5.3）。02〜06 と同じ値であること ------------------------------
AGE_BREAKS <- c(65, 75, 90)
AGE_LEVELS <- c("G1", "G2", "G3")

## --- 入院時期区分（§6.1）。04_ridit_spline.R と同じ規則 ---------------------
FY_START_MONTH  <- 4L
PERIOD_EARLY_FY <- 2017:2019
PERIOD_LATE_FY  <- 2020:2023
PERIOD_LEVELS   <- c("early", "late")

## --- ブートストラップ（§8.4）------------------------------------------------
NBOOT        <- 2000L                 # 主解析（06）と同じ
BOOT_SEED_S  <- 20260920L + 1000L     # 主解析の BOOT_SEED と別にする
CI_TYPE      <- 7L                    # stats::quantile の type。パーセンタイル法
STD_POP_IN_BOOT <- "resample"         # 06 と同じ既定
PAR_CORES      <- 3L                  # 再標本は先に一括生成するため、値によらず結果は同一
PROGRESS_EVERY <- 50L

## --- 標準化リスクの計算方法（06 と同じ）-------------------------------------
##   "predict.clm" : `predict.clm(type = "cum.prob")` を使う（既定）
##   "fast"        : 同じ量を plogis(θ_j(g) − x'β) で直接計算する（06 の説明を参照）。
##                   本体の当てはめで両モデルとも predict.clm と一致しなければ戻す。
PREDICT_METHOD    <- "predict.clm"
PREDICT_AGREE_TOL <- 1e-10

## --- 収束と発散の判定（05・06 と同じ）---------------------------------------
CLM_MAX_ITER     <- 200L
CLM_MAX_LINE     <- 50L
CLM_MAX_MOD_ITER <- 10L
DIVERGE_COEF   <- 10
DIVERGE_SE     <- 10
IS_SPLINE_TERM <- function(nm) grepl("^ns\\(", nm)
FAIL_WARN_RATE <- 0.05

## --- 自前の profile likelihood（05 と同じ設定）-------------------------------
PROFILE_UNIROOT_TOL <- 1e-9
PROFILE_MAX_SE      <- 50
LOGLIK_AGREE_TOL    <- 1e-6
## profile 区間とブートストラップ区間の突き合わせ（06 と同じ許容）
BOOT_PROFILE_TOL <- 0.10

## --- 対比（§8.4）。高齢側を被減数・分子に置く -------------------------------
CONTRASTS <- list(
  list(key = "G2_vs_G1", ref = "G1", cmp = "G2", role = "secondary"),
  list(key = "G3_vs_G1", ref = "G1", cmp = "G3", role = "primary"),
  list(key = "G3_vs_G2", ref = "G2", cmp = "G3", role = "secondary")
)
PRIMARY_KEY <- "G3_vs_G1"

## --- モデルの呼び名 ----------------------------------------------------------
MODEL_PPO_S <- "partial proportional odds, sensitivity"
MODEL_PO_S  <- "proportional odds, sensitivity"


# -----------------------------------------------------------------------------
# 1. 表示用ラベル
# -----------------------------------------------------------------------------

if (!file.exists("01_labels.R")) {
  stop("01_labels.R が見つからない。07_2_sens_clm_boot.R は 01_labels.R と同じ",
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

.log_con <- file(file.path(OUT_DIR, "log_07_2_sens_clm_boot.txt"), open = "wt",
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

rule(paste0("07_2_sens_clm_boot.R  ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
say("R : ", R.version.string)
say("nboot : ", NBOOT, " ; bootstrap seed : ", BOOT_SEED_S,
    " ; conf.level : ", CONF_LEVEL)
say("thresholds : P(composite mFIM at discharge >= t), t = ",
    paste(THRESHOLDS, collapse = ", "))
say("worst score assigned to dropout : ", WORST_SCORE,
    " (measured range ", MFIM_RANGE[1], "-", MFIM_RANGE[2], ")")
say("standardisation population inside the bootstrap : ", STD_POP_IN_BOOT)
say("parallel cores : ", PAR_CORES)

has_ordinal <- requireNamespace("ordinal", quietly = TRUE)
check("ordinal が利用できる", has_ordinal,
      if (has_ordinal) paste0("version ",
                              as.character(utils::packageVersion("ordinal")))
      else "install.packages(\"ordinal\") が必要")
if (!has_ordinal) stop("ordinal が無ければ感度分析Sの E-4 は算出できない。")
library(ordinal)
library(splines)


# -----------------------------------------------------------------------------
# 3. 入力の取り込み
# -----------------------------------------------------------------------------

rule("3. Input")

## rda は必ず専用の環境に読み込み、必要なオブジェクトだけを取り出す。
## 直接 globalenv() に読み込むと、05・06 の rda に入っている CONTRASTS・CONF_LEVEL・
## ctrl・THRESHOLDS などがこのスクリプトの設定（§0）を黙って上書きしてしまう。
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

e04 <- load_env(IN_RDA_04)
take(e04, c("ridit_fun", "spline_terms", "M2_SETTINGS"), IN_RDA_04)
if (exists("cov_levels", envir = e04, inherits = FALSE)) {
  cov_levels <- get("cov_levels", envir = e04)
}

e05 <- load_env(IN_RDA_05)
.need05 <- c("M2_FORMULA_TXT", "M2_PPO_FORMULA_TXT", "M2_PPO_NOMINAL_TXT",
             "PPO_THRESHOLDS", "PPO_LEVELS")
if (!all(vapply(.need05, exists, logical(1), envir = e05, inherits = FALSE))) {
  stop(IN_RDA_05, " に部分比例オッズモデルの式（M2_PPO_FORMULA_TXT ほか）が無い。",
       "§7 版の 05_m2_clm.R を先に実行すること。")
}
take(e05, .need05, IN_RDA_05)

## 補足表2 の材料（無くても S 単独の表は作れる）
main_e4 <- load_env(IN_RDA_06,   required = FALSE)
sens_e2 <- load_env(IN_RDA_07_1, required = FALSE)
if (is.null(main_e4)) say("  補足表2 の主解析の列は空になる（06 を先に実行すること）")
if (is.null(sens_e2)) say("  補足表2 の S の E-2 行は空になる（07_1 を先に実行すること）")

## 06 の rda が §7 版（閾値 70・77、二値モデルなし）であることを確かめる。
## 旧版（閾値 65、bin_tab あり）の数値を補足表2 に混ぜないため、旧版なら主解析の
## E-4 列は空にする。
main_e4_is_s7 <- FALSE
if (!is.null(main_e4)) {
  .thr06 <- if (exists("THRESHOLDS", envir = main_e4, inherits = FALSE))
              as.numeric(get("THRESHOLDS", envir = main_e4)) else NA_real_
  main_e4_is_s7 <- !anyNA(.thr06) &&
    identical(.thr06, as.numeric(THRESHOLDS)) &&
    !exists("bin_tab", envir = main_e4, inherits = FALSE) &&
    exists("std_risk_tab", envir = main_e4, inherits = FALSE) &&
    "threshold" %in% names(get("std_risk_tab", envir = main_e4))
  check("06 の rda が §7 版である（閾値 70・77、二値モデルなし）", main_e4_is_s7,
        sprintf("06 の THRESHOLDS = %s ; bin_tab %s",
                paste(.thr06, collapse = ", "),
                if (exists("bin_tab", envir = main_e4, inherits = FALSE)) "あり（旧版）"
                else "なし"))
  if (!main_e4_is_s7) {
    say("  06 の rda が旧版のため、補足表2 の主解析の E-4 列は空にする。",
        "§7 版の 06 を先に実行すること。")
  }
}

## 07_1 の rda が §7-3 版（最低点 12）であることを確かめる
if (!is.null(sens_e2)) {
  .ws71 <- if (exists("WORST_SCORE", envir = sens_e2, inherits = FALSE))
             get("WORST_SCORE", envir = sens_e2) else NA
  check("07_1 の rda の最低点が 07_2 と一致する（§7-3）",
        isTRUE(.ws71 == WORST_SCORE),
        sprintf("07_1: %s / 07_2: %d", as.character(.ws71), WORST_SCORE))
}

OUTCOME    <- M2_SETTINGS$outcome            # "mFIM_out"（主解析のアウトカム）
EXPOSURE   <- M2_SETTINGS$exposure
COVARS_CAT <- M2_SETTINGS$covars_cat
COVARS_NUM <- M2_SETTINGS$covars_num
RIDIT_VARS <- M2_SETTINGS$ridit_vars
REF_LEVELS <- M2_SETTINGS$ref_levels
M2_VARS    <- c(EXPOSURE, COVARS_CAT, COVARS_NUM)

check("年齢群の規則が 04 の設定と一致する",
      identical(as.character(M2_SETTINGS$age_levels), AGE_LEVELS) &&
        identical(as.numeric(M2_SETTINGS$age_breaks), as.numeric(AGE_BREAKS)),
      sprintf("04: %s / 07_2: %s",
              paste(M2_SETTINGS$age_breaks, collapse = ", "),
              paste(AGE_BREAKS, collapse = ", ")))
check("入院時期区分の水準が 04 の設定と一致する",
      identical(as.character(M2_SETTINGS$period_levels), PERIOD_LEVELS), "")

## --- 部分比例オッズモデルのアウトカム（3区分）。05・06 と同じ規則で作る ------
## 水準名は閾値だけから決まる。最低点 12 は「<70」に入る。
make_ppo_levels <- function(thresholds) {
  th <- sort(thresholds)
  mid <- if (length(th) > 1L) sprintf("%g-%g", head(th, -1L), tail(th, -1L) - 1) else character(0)
  c(sprintf("<%g", th[1]), mid, sprintf(">=%g", th[length(th)]))
}
make_ppo_outcome <- function(y, thresholds) {
  lev <- make_ppo_levels(thresholds)
  factor(lev[findInterval(y, sort(thresholds)) + 1L], levels = lev, ordered = TRUE)
}
PPO_LEVELS_S <- make_ppo_levels(THRESHOLDS)
check("07_2 の閾値が 05 の部分比例オッズモデルの閾値と一致する",
      identical(as.numeric(THRESHOLDS), as.numeric(PPO_THRESHOLDS)) &&
        identical(PPO_LEVELS_S, PPO_LEVELS),
      sprintf("07_2: %s / 05: %s", paste(THRESHOLDS, collapse = ", "),
              paste(PPO_THRESHOLDS, collapse = ", ")))
M2_PPO_NOMINAL_F <- stats::as.formula(M2_PPO_NOMINAL_TXT)

sens <- as.data.frame(dat_sens, stringsAsFactors = FALSE)
say("\ndat_sens (sensitivity S population, A4 minus E3(b)(c)) : ", nrow(sens), " rows")

.need <- c("id", "age", "day_in", "sex", "class", "support_in",
           "mFIM_in", "cFIM_in", OUTCOME_S)
.miss <- setdiff(.need, names(sens))
if (length(.miss)) stop("dat_sens に必要な列がない: ", paste(.miss, collapse = ", "))

check("Sのアウトカムに欠測がない（§8.8：複合アウトカム）",
      !anyNA(sens[[OUTCOME_S]]),
      sprintf("missing = %d", sum(is.na(sens[[OUTCOME_S]]))))
check(sprintf("Sのアウトカムが値域 %d〜%d に収まる（§7-3：最低点 %d を含む）",
              MFIM_RANGE_S[1], MFIM_RANGE_S[2], WORST_SCORE),
      all(sens[[OUTCOME_S]] >= MFIM_RANGE_S[1] & sens[[OUTCOME_S]] <= MFIM_RANGE_S[2]),
      sprintf("range = %g - %g", min(sens[[OUTCOME_S]]), max(sens[[OUTCOME_S]])))
check(sprintf("最低点（%d）が実測の値域（%d〜%d）の外にある（§7-3）",
              WORST_SCORE, MFIM_RANGE[1], MFIM_RANGE[2]),
      WORST_SCORE < MFIM_RANGE[1], sprintf("WORST_SCORE = %d", WORST_SCORE))
.obs_S <- sens[[OUTCOME_S]][sens[[OUTCOME_S]] != WORST_SCORE]
check("最低点を除く値が実測の値域に収まる",
      all(.obs_S >= MFIM_RANGE[1] & .obs_S <= MFIM_RANGE[2]),
      if (length(.obs_S)) sprintf("range = %g - %g", min(.obs_S), max(.obs_S)) else "no data")
if (all(c("disposition", "e3_category") %in% names(sens))) {
  .rule_worst <- (sens$disposition %in% DISPO_WORST) |
                 (!is.na(sens$e3_category) & sens$e3_category %in% E3_CATEGORY_WORST)
  .is_worst <- sens[[OUTCOME_S]] == WORST_SCORE
  check(sprintf("%d 点の例が、02 の規則で最低点を与える患者（転院・終了（死亡等）・E3(a)）と一致する（§7-3）",
                WORST_SCORE),
        identical(as.logical(.is_worst), as.logical(.rule_worst)),
        sprintf("%d 点 = %d, 規則 = %d, 食い違い = %d", WORST_SCORE, sum(.is_worst),
                sum(.rule_worst), sum(.is_worst != .rule_worst)))
}
check("解析単位が患者である（id が一意。§5.2）", !anyDuplicated(sens$id), "")


# -----------------------------------------------------------------------------
# 4. 派生変数と ridit（04_ridit_spline.R に相当）
#    年齢群・入院年度・入院時期区分を作り、共変量の水準の並びを 04 に合わせる。
#    ridit は 04 が A4 で作った変換関数をそのまま当てる（§6.4-5）。
# -----------------------------------------------------------------------------

rule("4. Derived variables and ridit (as in 04_ridit_spline.R)")

make_agegroup <- function(age) {
  factor(ifelse(is.na(age), NA_character_,
         ifelse(age >= AGE_BREAKS[3], "G3",
         ifelse(age >= AGE_BREAKS[2], "G2", "G1"))),
         levels = AGE_LEVELS)
}
make_fy <- function(day) {
  y <- as.integer(format(day, "%Y"))
  m <- as.integer(format(day, "%m"))
  ifelse(is.na(y) | is.na(m), NA_integer_, ifelse(m >= FY_START_MONTH, y, y - 1L))
}
make_period <- function(fy) {
  factor(ifelse(is.na(fy), NA_character_,
         ifelse(fy %in% PERIOD_EARLY_FY, "early",
         ifelse(fy %in% PERIOD_LATE_FY,  "late", NA_character_))),
         levels = PERIOD_LEVELS)
}
#' 04_ridit_spline.R と同じ規則で水準を並べる（参照水準を先頭に置く）
make_factor <- function(x, ref, order = NULL) {
  obs <- sort(unique(as.character(x[!is.na(x)])))
  lv  <- unique(c(ref, intersect(order, obs), setdiff(obs, c(ref, order))))
  factor(as.character(x), levels = lv)
}

sens$age_group <- make_agegroup(sens$age)
sens$fy_in     <- make_fy(sens$day_in)
sens$period    <- make_period(sens$fy_in)

## 水準の並びは 04 が決めたものに合わせる（cov_levels があればそれを正本とする）
if (exists("cov_levels", inherits = TRUE)) {
  for (v in COVARS_CAT) {
    .r <- cov_levels[cov_levels$variable == v, , drop = FALSE]
    if (nrow(.r)) {
      lv  <- strsplit(.r$levels[1], " | ", fixed = TRUE)[[1]]
      obs <- sort(unique(as.character(sens[[v]][!is.na(sens[[v]])])))
      sens[[v]] <- factor(as.character(sens[[v]]),
                          levels = unique(c(lv, setdiff(obs, lv))))
    }
  }
  check("共変量の水準の並びを 04 の cov_levels から引き継いだ", TRUE,
        "主解析と参照水準・水準順が一致する")
} else {
  sens$sex        <- make_factor(sens$sex,        REF_LEVELS$sex)
  sens$class      <- make_factor(sens$class,      REF_LEVELS$class)
  sens$support_in <- make_factor(sens$support_in, REF_LEVELS$support_in)
  check("04 の cov_levels が無いため水準を再構成した", FALSE,
        "参照水準が主解析と一致するかログで確認すること")
}
sens$period    <- factor(as.character(sens$period),    levels = PERIOD_LEVELS)
sens$age_group <- factor(as.character(sens$age_group), levels = AGE_LEVELS)

for (v in c(EXPOSURE, COVARS_CAT)) {
  say(sprintf("  %-12s reference = %-10s levels = %s", v, levels(sens[[v]])[1],
              paste(levels(sens[[v]]), collapse = " | ")))
}
check("参照水準が主解析と同じである",
      all(vapply(names(REF_LEVELS), function(v)
        identical(levels(sens[[v]])[1], REF_LEVELS[[v]]), logical(1))),
      paste(vapply(names(REF_LEVELS), function(v)
        sprintf("%s=%s", v, levels(sens[[v]])[1]), character(1)), collapse = ", "))
check("年齢群に欠測がない", !anyNA(sens$age_group),
      sprintf("missing = %d", sum(is.na(sens$age_group))))
check("入院時期区分に欠測がない", !anyNA(sens$period),
      sprintf("missing = %d", sum(is.na(sens$period))))

## --- ridit（A4 基準。04 の変換関数をそのまま使う。§6.4-5）-------------------
for (v in names(RIDIT_VARS)) {
  rv <- unname(RIDIT_VARS[v])
  sens[[rv]] <- ridit_fun[[v]](sens[[v]])
}
check("ridit を再算出していない（04 の ridit_fun をそのまま当てた。§6.4-5）", TRUE,
      "基準集団は A4。感度分析Sの集団 ⊆ A4 であるため外挿の枝には入らない")
for (v in names(RIDIT_VARS)) {
  rv <- unname(RIDIT_VARS[v])
  r  <- sens[[rv]][!is.na(sens[[rv]])]
  check(sprintf("%s の ridit が (0, 1) に収まる", v),
        length(r) > 0 && all(r > 0 & r < 1),
        if (length(r)) sprintf("range = [%.6f, %.6f]", min(r), max(r)) else "no data")
}

## --- 共変量の欠測と完全ケース集団（§6.5 と同じ数え方）-----------------------
miss_rows <- list()
for (v in M2_VARS) {
  for (g in c(AGE_LEVELS, "total")) {
    sel <- if (identical(g, "total")) rep(TRUE, nrow(sens)) else sens$age_group == g
    n   <- sum(sel)
    nm  <- sum(is.na(sens[[v]][sel]))
    miss_rows[[length(miss_rows) + 1L]] <- data.frame(
      population      = "Sensitivity S population",
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
  }
}
sens_missing <- do.call(rbind, miss_rows)
say("\n  missing covariates in the sensitivity S population:")
print(sens_missing[sens_missing$age_group == "total",
                   c("variable_label", "n", "n_missing", "percent_missing")],
      row.names = FALSE, digits = 3)


# -----------------------------------------------------------------------------
# 5. M2-S の解析データと当てはめ（05_m2_clm.R に相当。E-3 は表に出さない）
#    §7-2：部分比例オッズモデル（M2-PPO）と比例オッズモデル（M2）の2つを、
#    それぞれ年齢群の参照を G1・G2 として当てはめる（計4本）。
# -----------------------------------------------------------------------------

rule("5. M2-PPO and M2 on the sensitivity S population (E-3; fitted but not reported)")

complete_S <- stats::complete.cases(sens[, M2_VARS, drop = FALSE])
keep_cols  <- c("id", "age", "age_group", "sex", "class", "support_in",
                "fy_in", "period", "mFIM_in", "cFIM_in",
                unname(RIDIT_VARS), OUTCOME_S)
DAT <- sens[complete_S, intersect(keep_cols, names(sens)), drop = FALSE]
DAT$age_group <- factor(as.character(DAT$age_group), levels = AGE_LEVELS)

Y_S <- DAT[[OUTCOME_S]]
Y_LEVELS_S <- sort(unique(Y_S))
DAT$y_ord <- factor(as.character(Y_S), levels = as.character(Y_LEVELS_S),
                    ordered = TRUE)
DAT$y_ppo <- make_ppo_outcome(Y_S, THRESHOLDS)     # 3区分（§7-2）
N_S    <- nrow(DAT)
n_by_g <- table(DAT$age_group)

say(sprintf("  sensitivity S population        : %d", nrow(sens)))
say(sprintf("  complete cases for M2 (S)       : %d (dropped %d, %.2f%%)",
            N_S, nrow(sens) - N_S, 100 * (nrow(sens) - N_S) / nrow(sens)))
for (g in AGE_LEVELS) {
  say(sprintf("    %-3s %-14s n = %5d", g, relabel_levels("age_group", g),
              as.integer(n_by_g[[g]])))
}
say(sprintf("  observed outcome levels (M2)    : %d (range %g - %g)",
            length(Y_LEVELS_S), min(Y_LEVELS_S), max(Y_LEVELS_S)))

## §7-3：最低点は 12 点で実測の値域の外にあるため、12 点の例数がそのまま
## 最低点を与えた例数になる（実測で 13 点の例は含まれない）。
n_worst    <- sum(sens[[OUTCOME_S]] == WORST_SCORE)
n_worst_cc <- sum(DAT[[OUTCOME_S]] == WORST_SCORE)
say(sprintf("  observations at the worst score (%d) : %d in S, %d in the M2-S set",
            WORST_SCORE, n_worst, n_worst_cc))
say(sprintf("  observations at the measured floor (%d, not assigned) : %d in S, %d in the M2-S set",
            MFIM_RANGE[1], sum(sens[[OUTCOME_S]] == MFIM_RANGE[1]),
            sum(DAT[[OUTCOME_S]] == MFIM_RANGE[1])))

sens_analysis_n <- data.frame(
  set   = c("sensitivity_S", rep("m2_S_complete_case", 1 + length(AGE_LEVELS))),
  group = c("total", "total", AGE_LEVELS),
  label = c("Sensitivity S population (used by E-2)",
            "M2-S complete-case population (used by E-4)",
            relabel_levels("age_group", AGE_LEVELS)),
  n     = c(nrow(sens), N_S, as.integer(n_by_g[AGE_LEVELS])),
  stringsAsFactors = FALSE
)
say("\n  N by analysis set:")
print(sens_analysis_n, row.names = FALSE)

check("M2-S の完全ケース集団に欠測がない",
      !anyNA(DAT[, c(M2_VARS, unname(RIDIT_VARS), OUTCOME_S)]), "")
check("M2-S の id が一意（解析単位＝患者、§5.2）", !anyDuplicated(DAT$id), "")
check("M2-S の各年齢群に例がある", all(table(DAT$age_group) > 0),
      paste(sprintf("%s=%d", AGE_LEVELS, as.integer(n_by_g)), collapse = ", "))
check("各カテゴリ共変量の各水準に例がある",
      all(vapply(COVARS_CAT, function(v) all(table(DAT[[v]]) > 0), logical(1))), "")
check("アウトカムの順序因子の水準が数値として昇順である",
      !is.unsorted(Y_LEVELS_S), "")
check(sprintf("M2 のアウトカムの最下位の水準が最低点（%d）である（最低点を与えた例がある場合）",
              WORST_SCORE),
      n_worst_cc == 0L || identical(as.numeric(Y_LEVELS_S[1]), as.numeric(WORST_SCORE)),
      sprintf("lowest level = %g", Y_LEVELS_S[1]))

## --- 部分比例オッズモデルのアウトカム（3区分）-------------------------------
ppo_cells <- table(factor(DAT$age_group, levels = AGE_LEVELS), DAT$y_ppo)
say(sprintf("\n  PPO outcome : %s (thresholds %s ; the worst score %d falls in \"%s\")",
            paste(PPO_LEVELS_S, collapse = " / "), paste(THRESHOLDS, collapse = ", "),
            WORST_SCORE, as.character(make_ppo_outcome(WORST_SCORE, THRESHOLDS))))
print(ppo_cells)
check("3区分のアウトカムに欠測がない", !anyNA(DAT$y_ppo), "")
check("年齢群 × 3区分 のすべてのセルに観測がある（nominal 項の推定に必要）",
      all(ppo_cells > 0L), sprintf("最小セル = %d", min(ppo_cells)))

## --- モデル式（05 の文字列をそのまま使う。ノットを再算出しない。§6.4-5）----
say("\n  proportional odds (M2)            : ", M2_FORMULA_TXT)
say("  partial proportional odds (M2-PPO): ", M2_PPO_FORMULA_TXT,
    " ; nominal = ", M2_PPO_NOMINAL_TXT)
check("ridit とノットを再算出していない（05 の式の文字列をそのまま使う）",
      all(grepl("Boundary.knots = c\\(", unname(spline_terms))) &&
        all(vapply(unname(spline_terms), function(s)
          grepl(s, M2_FORMULA_TXT, fixed = TRUE) && grepl(s, M2_PPO_FORMULA_TXT, fixed = TRUE),
          logical(1))),
      "2つのモデル式にノットの数値が literal で埋め込まれている")

ctrl <- call_known_args(ordinal::clm.control,
                        list(maxIter = CLM_MAX_ITER, maxLineIter = CLM_MAX_LINE,
                             maxModIter = CLM_MAX_MOD_ITER))

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
is_nominal_fit    <- function(fit) !is.null(fit[["nom.terms"]]) && length(age_nominal_names(fit)) > 0L

#' 年齢群の参照水準を入れ替えたデータを作る
dat_ref <- function(data, ref) {
  d <- data
  d$age_group <- stats::relevel(factor(as.character(d$age_group),
                                       levels = AGE_LEVELS), ref = ref)
  d
}
#' 比例オッズモデル（M2）
fit_m2S <- function(data, ref) {
  d <- dat_ref(data, ref)
  t0  <- proc.time()[["elapsed"]]
  fit <- ordinal::clm(stats::as.formula(M2_FORMULA_TXT), data = d,
                      link = "logit", control = ctrl)
  attr(fit, "seconds") <- proc.time()[["elapsed"]] - t0
  attr(fit, "ref")     <- ref
  fit
}
#' 部分比例オッズモデル（M2-PPO）。年齢群を nominal 項に置く（§7-2）
fit_m2S_ppo <- function(data, ref) {
  d <- dat_ref(data, ref)
  t0  <- proc.time()[["elapsed"]]
  fit <- ordinal::clm(stats::as.formula(M2_PPO_FORMULA_TXT), nominal = M2_PPO_NOMINAL_F,
                      data = d, link = "logit", control = ctrl)
  attr(fit, "seconds") <- proc.time()[["elapsed"]] - t0
  attr(fit, "ref")     <- ref
  fit
}

say("  fitting M2 with age_group reference = G1 ...")
fitS_G1 <- fit_m2S(DAT, "G1")
say(sprintf("    done in %.1f s", attr(fitS_G1, "seconds")))
say("  fitting M2 with age_group reference = G2 （G3対G2 の profile 区間のため）...")
fitS_G2 <- fit_m2S(DAT, "G2")
say(sprintf("    done in %.1f s", attr(fitS_G2, "seconds")))
say("  fitting M2-PPO with age_group reference = G1 ...")
fitS_ppo_G1 <- fit_m2S_ppo(DAT, "G1")
say(sprintf("    done in %.1f s", attr(fitS_ppo_G1, "seconds")))
say("  fitting M2-PPO with age_group reference = G2 （G3対G2 の profile 区間のため）...")
fitS_ppo_G2 <- fit_m2S_ppo(DAT, "G2")
say(sprintf("    done in %.1f s", attr(fitS_ppo_G2, "seconds")))

#' 発散の判定に使う係数：位置の係数（β）＋ nominal 項の年齢群の母数（M2-PPO のみ）
fit_info <- function(fit, tag) {
  an <- age_nominal_names(fit)
  b  <- c(fit[["beta"]], fit[["alpha"]][an])
  se <- tryCatch(sqrt(diag(vcov(fit)))[names(b)],
                 error = function(e) rep(NA_real_, length(b)))
  mg <- fit[["maxGradient"]]
  yl <- fit[["y.levels"]]
  data.frame(
    model        = tag,
    reference    = attr(fit, "ref"),
    n            = tryCatch(as.integer(stats::nobs(fit))[1],
                            error = function(e) NA_integer_),
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
sens_fit_info <- rbind(fit_info(fitS_G1, "M2-S (ref G1)"),
                       fit_info(fitS_G2, "M2-S (ref G2)"),
                       fit_info(fitS_ppo_G1, "M2-PPO-S (ref G1)"),
                       fit_info(fitS_ppo_G2, "M2-PPO-S (ref G2)"))
say("")
print(sens_fit_info, row.names = FALSE, digits = 6)

for (i in seq_len(nrow(sens_fit_info))) {
  r <- sens_fit_info[i, ]
  check(sprintf("[%s] 収束した", r$model),
        is.na(r$convergence) || r$convergence == 0L,
        sprintf("convergence = %s (%s)", as.character(r$convergence),
                r$convergence_message))
  check(sprintf("[%s] 係数・標準誤差に発散の兆候がない", r$model),
        all(is.finite(c(r$max_abs_beta, r$max_abs_beta_nonspline))) &&
          r$max_abs_beta_nonspline < DIVERGE_COEF &&
          (is.na(r$max_se) || r$max_se < DIVERGE_SE),
        sprintf("max|beta| = %.3f（うちスプライン以外 %.3f）, max se = %s",
                r$max_abs_beta, r$max_abs_beta_nonspline,
                formatC(r$max_se, format = "g", digits = 3)))
}
check("M2-S の2つの当てはめの対数尤度が一致する（参照水準の入れ替えのみであるため）",
      abs(sens_fit_info$logLik[1] - sens_fit_info$logLik[2]) < 1e-6,
      sprintf("%.8f vs %.8f", sens_fit_info$logLik[1], sens_fit_info$logLik[2]))
check("M2-PPO-S の2つの当てはめの対数尤度が一致する（参照水準の入れ替えのみであるため）",
      abs(sens_fit_info$logLik[3] - sens_fit_info$logLik[4]) < 1e-6,
      sprintf("%.8f vs %.8f", sens_fit_info$logLik[3], sens_fit_info$logLik[4]))
check("当てはめに使われた例数が完全ケース集団と一致する",
      all(sens_fit_info$n == N_S),
      sprintf("%s vs %d", paste(sens_fit_info$n, collapse = ", "), N_S))
check("M2-PPO-S の閾値が2つで、年齢群の nominal 母数が 2 × 2 個ある",
      sens_fit_info$n_thresholds[3] == length(THRESHOLDS) &&
        length(age_nominal_names(fitS_ppo_G1)) ==
          length(THRESHOLDS) * (length(AGE_LEVELS) - 1L),
      sprintf("thresholds = %d, age nominal = %s", sens_fit_info$n_thresholds[3],
              paste(age_nominal_names(fitS_ppo_G1), collapse = ", ")))
check("fitS_ppo_G1 が nominal 項をもち、fitS_G1 がもたない",
      is_nominal_fit(fitS_ppo_G1) && !is_nominal_fit(fitS_G1), "")

## nominal 項は閾値の順序を保証しない。年齢群を置き換えた予測で、
## 累積確率が全患者・全年齢群で単調であることを確かめる（05 と同じ点検）。
.nd_ppo <- DAT
.nd_ppo$y_ord <- NULL; .nd_ppo$y_ppo <- NULL
.mono_ppo <- vapply(AGE_LEVELS, function(g) {
  d <- .nd_ppo
  d$age_group <- factor(rep(g, nrow(d)), levels = AGE_LEVELS)
  p  <- stats::predict(fitS_ppo_G1, newdata = d, type = "cum.prob")
  cp <- if (is.list(p) && !is.null(p$cprob1)) p$cprob1 else p
  all(apply(cp, 1L, function(r) all(diff(r) >= -1e-12)))
}, logical(1))
check("M2-PPO-S の累積確率が全患者・全年齢群で単調である（閾値の交差がない）",
      all(.mono_ppo),
      paste(sprintf("%s=%s", AGE_LEVELS, ifelse(.mono_ppo, "ok", "NG")), collapse = ", "))

## --- E-3（M2 の共通オッズ比、profile likelihood 区間）。算出するが補足表2 には載せない
get_profile_ci <- function(fit, parm, level = CONF_LEVEL) {
  as_mat <- function(x, nm) {
    if (is.null(x)) return(NULL)
    if (is.null(dim(x))) x <- matrix(x, nrow = 1L, dimnames = list(nm, names(x)))
    as.matrix(x)
  }
  ci <- tryCatch(as_mat(confint(fit, parm = parm, level = level, type = "profile"),
                        parm),
                 error = function(e) {
                   say("    confint(parm=...) でエラー: ", conditionMessage(e)); NULL
                 })
  src <- "confint(type = \"profile\", parm)"
  if (is.null(ci) || !all(parm %in% rownames(ci))) {
    ci2 <- tryCatch(as_mat(confint(fit, level = level, type = "profile"), parm),
                    error = function(e) NULL)
    if (!is.null(ci2) && all(parm %in% rownames(ci2))) {
      ci  <- ci2[parm, , drop = FALSE]
      src <- "confint(type = \"profile\", all beta)"
    } else {
      return(list(ci = NULL, source = "failed"))
    }
  }
  list(ci = ci[parm, , drop = FALSE], source = src)
}

plS_G1 <- get_profile_ci(fitS_G1, c("age_groupG2", "age_groupG3"))
plS_G2 <- get_profile_ci(fitS_G2, c("age_groupG3"))

sens_or_rows <- list()
for (ct in CONTRASTS) {
  if (identical(ct$ref, "G1")) { f <- fitS_G1; pl <- plS_G1 } else { f <- fitS_G2; pl <- plS_G2 }
  term <- paste0("age_group", ct$cmp)
  b    <- unname(f$beta[term])
  ci   <- if (!is.null(pl$ci) && term %in% rownames(pl$ci)) pl$ci[term, ] else c(NA, NA)
  sens_or_rows[[length(sens_or_rows) + 1L]] <- data.frame(
    analysis         = "sensitivity S",
    contrast         = ct$key,
    label            = contrast_label(ct$ref, ct$cmp),
    group_reference  = ct$ref,
    group_comparison = ct$cmp,
    role             = ct$role,
    term             = term,
    beta             = b,
    odds_ratio       = exp(b),
    ci_lower         = exp(as.numeric(ci[1])),
    ci_upper         = exp(as.numeric(ci[2])),
    conf_level       = CONF_LEVEL,
    ci_method        = "profile likelihood, unadjusted",
    model            = MODEL_PO_S,
    p_value          = NA_real_,
    reported         = FALSE,   # §8.8：Sの E-3 は表に出さない
    note             = "Fitted for comparison only. Not reported (plan 8.8).",
    stringsAsFactors = FALSE
  )
}
sens_odds_ratio <- do.call(rbind, sens_or_rows)
say("\n  E-3 (sensitivity S) : common odds ratio, M2  [NOT reported; log and csv only]")
print(sens_odds_ratio[, c("contrast", "label", "odds_ratio", "ci_lower", "ci_upper")],
      row.names = FALSE, digits = 4)


# -----------------------------------------------------------------------------
# 5b. M2-PPO-S の閾値別オッズ比と profile likelihood 95%区間（算出するが表に出さない）
#     05 の 6b 節と同じ自前の実装（ordinal の profile は nominal 項を扱わない）。
#     実装の正しさは 05 で (i)〜(iii) を確かめている。ここでは S のデータで
#     (i) 対数尤度の一致と (ii) 最大であることを確かめる。
# -----------------------------------------------------------------------------

rule("5b. M2-PPO-S: threshold-specific odds ratios, profile likelihood 95% CI [NOT reported]")

## 閾値 → clm の閾値名（"<70|70-76"、"70-76|>=77"）
PPO_CUTS_S <- stats::setNames(paste(head(PPO_LEVELS_S, -1L), tail(PPO_LEVELS_S, -1L), sep = "|"),
                              as.character(THRESHOLDS))

#' M2-PPO の対数尤度の計算に要る行列をまとめる（05 と同じ）。
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

#' 対数尤度（grad = TRUE なら勾配）。logit P(Y ≤ c) = θ_c + τ_c'n − x'β（05 と同じ）。
ppo_loglik <- function(par, P, grad = FALSE) {
  C <- P$C; q <- P$q
  A <- matrix(par[seq_len(C * q)], nrow = C, ncol = q)
  b <- par[C * q + seq_len(P$p)]
  H <- P$N %*% t(A) - as.vector(P$X %*% b)
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
  G <- matrix(0, n, C)
  G[cbind(iu, k[iu])] <- stats::dlogis(H[cbind(iu, k[iu])]) / pr[iu]
  G[cbind(il, k[il] - 1L)] <- G[cbind(il, k[il] - 1L)] -
    stats::dlogis(H[cbind(il, k[il] - 1L)]) / pr[il]
  c(as.vector(crossprod(G, P$N)), -as.vector(crossprod(P$X, rowSums(G))))
}

#' 母数 idx を val に固定し、残りを最大化する（BFGS、解析的な勾配。05 と同じ）
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

#' profile likelihood 区間（母数 idx の尺度。05 と同じ）
ppo_profile_ci <- function(P, mle, idx, se, level = CONF_LEVEL) {
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

ppo_fit_S <- list(G1 = fitS_ppo_G1, G2 = fitS_ppo_G2)
ppo_P_S   <- list(G1 = ppo_parts(fitS_ppo_G1, dat_ref(DAT, "G1")),
                  G2 = ppo_parts(fitS_ppo_G2, dat_ref(DAT, "G2")))
ppo_mle_S <- list(G1 = c(fitS_ppo_G1$alpha, fitS_ppo_G1$beta),
                  G2 = c(fitS_ppo_G2$alpha, fitS_ppo_G2$beta))

ppo_impl_check_S <- do.call(rbind, lapply(c("G1", "G2"), function(r) {
  own   <- ppo_loglik(ppo_mle_S[[r]], ppo_P_S[[r]])
  refit <- ppo_max_fixed(ppo_P_S[[r]], ppo_mle_S[[r]])
  data.frame(reference = r,
             loglik_clm = as.numeric(stats::logLik(ppo_fit_S[[r]])),
             loglik_own = own,
             loglik_own_remaximised = refit$loglik,
             remax_convergence = refit$convergence,
             stringsAsFactors = FALSE)
}))
print(ppo_impl_check_S, row.names = FALSE, digits = 10)
check("(i) 自前の対数尤度が clm の推定値で logLik(fit) と一致する（S）",
      all(abs(ppo_impl_check_S$loglik_own - ppo_impl_check_S$loglik_clm) < LOGLIK_AGREE_TOL),
      sprintf("max |diff| = %.3e (tol %.0e)",
              max(abs(ppo_impl_check_S$loglik_own - ppo_impl_check_S$loglik_clm)),
              LOGLIK_AGREE_TOL))
check("(ii) clm の推定値から最大化し直しても対数尤度が増えない（S）",
      all(ppo_impl_check_S$loglik_own_remaximised - ppo_impl_check_S$loglik_own <
            LOGLIK_AGREE_TOL),
      sprintf("max gain = %.3e",
              max(ppo_impl_check_S$loglik_own_remaximised - ppo_impl_check_S$loglik_own)))

## clm の nominal 母数 τ は logit P(Y ≤ j) に足される。オッズ比は exp(−τ) であり、
## τ の区間 [τ_L, τ_U] はオッズ比の区間 [exp(−τ_U), exp(−τ_L)] に対応する。
sens_ppo_or_rows <- list()
for (thr in THRESHOLDS) {
  cut <- PPO_CUTS_S[[as.character(thr)]]
  for (ct in CONTRASTS) {
    r    <- ct$ref
    fit  <- ppo_fit_S[[r]]
    term <- paste0(cut, ".age_group", ct$cmp)
    idx  <- match(term, ppo_P_S[[r]]$names)
    if (is.na(idx)) stop("M2-PPO-S の母数が見つからない: ", term)
    tau  <- unname(ppo_mle_S[[r]][idx])
    se   <- tryCatch(unname(sqrt(diag(vcov(fit)))[term]), error = function(e) NA_real_)
    t0   <- proc.time()[["elapsed"]]
    pr   <- ppo_profile_ci(ppo_P_S[[r]], ppo_mle_S[[r]], idx, se)
    sec  <- proc.time()[["elapsed"]] - t0
    say(sprintf("    %-9s y >= %g : %d evaluations, %.1f s", ct$key, thr, pr$n_eval, sec))
    sens_ppo_or_rows[[length(sens_ppo_or_rows) + 1L]] <- data.frame(
      analysis         = "sensitivity S",
      contrast         = ct$key,
      label            = contrast_label(ct$ref, ct$cmp),
      group_reference  = ct$ref,
      group_comparison = ct$cmp,
      role             = ct$role,
      threshold        = thr,
      cut              = cut,
      term             = term,
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
      model            = MODEL_PPO_S,
      p_value          = NA_real_,
      reported         = FALSE,              # §8.8：Sの E-3 は表に出さない
      note             = "Fitted on the S population. Not reported (plan 8.8).",
      stringsAsFactors = FALSE
    )
  }
}
sens_ppo_odds_ratio <- do.call(rbind, sens_ppo_or_rows)
say("\n  E-3 (sensitivity S) : threshold-specific odds ratios, M2-PPO  [NOT reported; log and csv only]")
print(sens_ppo_odds_ratio[, c("threshold", "contrast", "odds_ratio", "ci_lower", "ci_upper")],
      row.names = FALSE, digits = 4)

check("M2-PPO-S のオッズ比が 閾値 × 対比 = 6 行ある",
      nrow(sens_ppo_odds_ratio) == length(THRESHOLDS) * length(CONTRASTS),
      sprintf("rows = %d", nrow(sens_ppo_odds_ratio)))
check("M2-PPO-S のすべてのオッズ比で profile 区間が得られている",
      all(is.finite(sens_ppo_odds_ratio$ci_lower) & is.finite(sens_ppo_odds_ratio$ci_upper)),
      paste(with(sens_ppo_odds_ratio,
                 paste0(contrast, "@", threshold)[!is.finite(ci_lower) | !is.finite(ci_upper)]),
            collapse = ", "))
check("M2-PPO-S の profile の最大化がすべて収束した",
      all(sens_ppo_odds_ratio$profile_not_converged == 0L),
      sprintf("未収束 %d / %d 回", sum(sens_ppo_odds_ratio$profile_not_converged),
              sum(sens_ppo_odds_ratio$profile_evaluations)))
check("M2-PPO-S の区間が点推定値を挟む",
      all(is.na(sens_ppo_odds_ratio$ci_lower) |
          (sens_ppo_odds_ratio$ci_lower <= sens_ppo_odds_ratio$odds_ratio &
           sens_ppo_odds_ratio$odds_ratio <= sens_ppo_odds_ratio$ci_upper)), "")
for (thr in THRESHOLDS) {
  .o <- sens_ppo_odds_ratio[sens_ppo_odds_ratio$threshold == thr, , drop = FALSE]
  .l <- stats::setNames(.o$beta, .o$contrast)
  check(sprintf("M2-PPO-S（y >= %g）：β(G3 vs G1) − β(G2 vs G1) = β(G3 vs G2)", thr),
        abs((.l[["G3_vs_G1"]] - .l[["G2_vs_G1"]]) - .l[["G3_vs_G2"]]) < 1e-6,
        sprintf("%.8f − %.8f = %.8f vs %.8f", .l[["G3_vs_G1"]], .l[["G2_vs_G1"]],
                .l[["G3_vs_G1"]] - .l[["G2_vs_G1"]], .l[["G3_vs_G2"]]))
}
check("SのE-3（閾値別OR・共通OR）を補足表2 に載せない印をつけた（§8.8）",
      all(!sens_odds_ratio$reported) && all(!sens_ppo_odds_ratio$reported), "reported = FALSE")


# -----------------------------------------------------------------------------
# 6. 閾値の水準 j と周辺標準化（06_standardize_boot.R に相当）
#    P(Y ≥ t) = 1 − P(Y ≤ j)。
#    M2（全水準）：j は「t 点未満の最大の観測水準」。
#    M2-PPO（3区分）：j は t 点のすぐ下の区分（<70、70-76）。
# -----------------------------------------------------------------------------

rule("6. Marginal standardisation (point estimates)")

pick_j <- function(levels_num, threshold) {
  lt <- levels_num[levels_num < threshold]
  if (!length(lt)) return(NA_real_)
  max(lt)
}

THR_KEYS <- paste0("t", THRESHOLDS)            # 列名などに使う閾値の鍵（"t70"、"t77"）

J_PO_S <- stats::setNames(vapply(THRESHOLDS, function(t) pick_j(Y_LEVELS_S, t), numeric(1)),
                          THR_KEYS)
J_PO_S_CHR <- stats::setNames(as.character(J_PO_S), THR_KEYS)
J_PPO_CHR  <- stats::setNames(PPO_LEVELS_S[seq_along(THRESHOLDS)], THR_KEYS)

for (i in seq_along(THRESHOLDS)) {
  t <- THRESHOLDS[i]
  say(sprintf("  threshold %d :", t))
  say(sprintf("    M2     (%d levels)  : P(Y >= %d) = 1 - P(Y <= %s)  （%d 点未満の最大の観測水準）",
              length(Y_LEVELS_S), t, J_PO_S_CHR[[i]], t))
  say(sprintf("    M2-PPO (3 categories): P(Y >= %d) = 1 - P(Y <= \"%s\")", t, J_PPO_CHR[[i]]))
}
check("M2-S の j が閾値ごとに決まった", all(is.finite(J_PO_S)),
      paste(sprintf("t = %g: j = %s", THRESHOLDS, J_PO_S_CHR), collapse = ", "))
check("M2-S の j が閾値未満の最大の観測水準である",
      all(vapply(seq_along(THRESHOLDS), function(i) {
        j <- J_PO_S[[i]]; t <- THRESHOLDS[i]
        is.finite(j) && j < t && !any(Y_LEVELS_S > j & Y_LEVELS_S < t)
      }, logical(1))), "")
check("M2-S の j が最大水準ではない（P(Y ≤ j) = 1 にならない）",
      all(is.finite(J_PO_S) & J_PO_S < max(Y_LEVELS_S)), "")
check("M2-PPO-S の j が閾値のすぐ下の区分である",
      all(vapply(seq_along(THRESHOLDS), function(i) {
        lv <- as.character(make_ppo_outcome(c(THRESHOLDS[i] - 1, THRESHOLDS[i]), THRESHOLDS))
        identical(lv[1], J_PPO_CHR[[i]]) && !identical(lv[2], J_PPO_CHR[[i]])
      }, logical(1))),
      paste(sprintf("t = %g: \"%s\"", THRESHOLDS, J_PPO_CHR), collapse = ", "))
if (main_e4_is_s7 && exists("J_PO", envir = main_e4, inherits = FALSE)) {
  .jm <- get("J_PO", envir = main_e4)
  check("M2-S の j が主解析（06）の M2 の j と一致する（閾値ごと）",
        identical(as.numeric(J_PO_S[THR_KEYS]), as.numeric(.jm[THR_KEYS])),
        sprintf("S: %s / main: %s", paste(J_PO_S_CHR, collapse = ", "),
                paste(as.character(.jm[THR_KEYS]), collapse = ", ")))
}

## 素の（調整しない）到達割合（記述。標準化リスクとは別物。査読の要望にも対応）
crude <- vapply(THRESHOLDS, function(t) vapply(AGE_LEVELS, function(g)
  mean(DAT[[OUTCOME_S]][DAT$age_group == g] >= t), numeric(1)),
  numeric(length(AGE_LEVELS)))
crude <- matrix(crude, nrow = length(AGE_LEVELS), dimnames = list(AGE_LEVELS, THR_KEYS))

#' 応答の水準から clm の閾値名を作る（"69|70"、"<70|70-76" など）
cut_names <- function(ylev) paste(head(ylev, -1L), tail(ylev, -1L), sep = "|")

#' 閾値 j・年齢群 g の θ_j(g)。M2 では θ_j、M2-PPO では θ_j + τ_j,g（参照群は 0）。06 と同じ。
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

#' 累積ロジットモデルからの標準化リスク P(Y ≥ t)（年齢群 × 閾値 の行列）。06 と同じ。
std_risk_clm <- function(fit, data, j_chr, groups, ylev) {
  k <- match(j_chr, ylev)
  if (anyNA(k)) stop("j が応答の水準に見つからない: ", paste(j_chr[is.na(k)], collapse = ", "))
  nd0 <- data
  nd0$y_ord <- NULL
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

#' 同じ量を plogis(θ_j(g) − x'β) で直接計算する（PREDICT_METHOD = "fast"）。06 と同じ。
std_risk_clm_fast <- function(fit, data, j_chr, groups, ylev) {
  b  <- fit[["beta"]]
  mm <- stats::model.matrix(stats::delete.response(stats::terms(fit)), data = data)
  if (!all(names(b) %in% colnames(mm))) stop("model.matrix の列が係数名と対応しない")
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

std_risk <- function(fit, data, j_chr, groups, ylev) {
  if (identical(PREDICT_METHOD, "fast")) std_risk_clm_fast(fit, data, j_chr, groups, ylev)
  else std_risk_clm(fit, data, j_chr, groups, ylev)
}

RD_KEYS <- vapply(CONTRASTS, function(ct) ct$key, character(1))

#' 群別リスク（年齢群 × 閾値）から3対比のRDを作る（RD = リスク(高齢側) − リスク(若年側)）
risk_to_rd <- function(risk) {
  out <- t(vapply(CONTRASTS, function(ct) risk[ct$cmp, ] - risk[ct$ref, ],
                  numeric(ncol(risk))))
  if (ncol(risk) == 1L) out <- t(out)
  dimnames(out) <- list(RD_KEYS, colnames(risk))
  out
}

#' 行列（行 = 年齢群または対比、列 = 閾値）を名前つきベクトルにする（"G1_t70" など）
flat <- function(m) stats::setNames(as.vector(m),
                                    as.vector(outer(rownames(m), colnames(m), paste, sep = "_")))
RISK_KEYS <- as.vector(outer(AGE_LEVELS, THR_KEYS, paste, sep = "_"))
RDT_KEYS  <- as.vector(outer(RD_KEYS,    THR_KEYS, paste, sep = "_"))

#' M2-PPO の閾値別の対数オッズ比（Y ≥ t）。log OR = −(τ_cmp − τ_ref)。06 と同じ。
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
#' M2 の共通対数オッズ比（3対比）
logor_po <- function(fit) {
  b <- fit[["beta"]]
  bg <- function(g) { tm <- paste0("age_group", g); if (tm %in% names(b)) unname(b[[tm]]) else 0 }
  stats::setNames(vapply(CONTRASTS, function(ct) bg(ct$cmp) - bg(ct$ref), numeric(1)), RD_KEYS)
}

## --- 本体の推定値 -----------------------------------------------------------
Y_LEVELS_S_CHR <- levels(DAT$y_ord)
PPO_LEVELS_CHR <- levels(DAT$y_ppo)

FITS  <- list(ppo = fitS_ppo_G1, po = fitS_G1)
JCHR  <- list(ppo = J_PPO_CHR,   po = J_PO_S_CHR)
YLEVS <- list(ppo = PPO_LEVELS_CHR, po = Y_LEVELS_S_CHR)
MODEL_NAME <- c(ppo = MODEL_PPO_S, po = MODEL_PO_S)

## --- 2つの計算方法の一致を確かめる（06 と同じ手順。両モデル）----------------
risk_pkg <- list(); agree <- logical(0)
for (m in names(FITS)) {
  risk_pkg[[m]] <- std_risk_clm(FITS[[m]], DAT, JCHR[[m]], AGE_LEVELS, YLEVS[[m]])
  rf <- tryCatch(std_risk_clm_fast(FITS[[m]], DAT, JCHR[[m]], AGE_LEVELS, YLEVS[[m]]),
                 error = function(e) { say("  fast path error: ", conditionMessage(e)); NULL })
  agree[m] <- !is.null(rf) && max(abs(risk_pkg[[m]] - rf)) < PREDICT_AGREE_TOL
  check(sprintf("[%s] 2つの計算方法（predict.clm と fast）が一致する", MODEL_NAME[[m]]),
        agree[[m]],
        if (is.null(rf)) "fast の計算に失敗した"
        else sprintf("max |diff| = %.3e", max(abs(risk_pkg[[m]] - rf))))
}
if (identical(PREDICT_METHOD, "fast") && !all(agree)) {
  PREDICT_METHOD <- "predict.clm"
  check("PREDICT_METHOD を predict.clm に戻した", FALSE,
        "fast が一致しなかったため、§8.4 の記載どおりの実装に戻す")
}
say("  PREDICT_METHOD in the bootstrap : ", PREDICT_METHOD)

risk_point    <- risk_pkg$ppo                  # 補足表2（M2-PPO）
rd_point      <- risk_to_rd(risk_point)
risk_point_po <- risk_pkg$po                   # 比較用（M2）
rd_point_po   <- risk_to_rd(risk_point_po)
logor_point_ppo <- logor_ppo(fitS_ppo_G1, PPO_LEVELS_CHR)
logor_point_po  <- logor_po(fitS_G1)

for (m in names(FITS)) {
  rk <- risk_pkg[[m]]; rd <- risk_to_rd(rk)
  say(sprintf("\n  standardised risk, %s (point estimates):", MODEL_NAME[[m]]))
  for (g in AGE_LEVELS) {
    say(sprintf("    %-3s %-14s %s", g, relabel_levels("age_group", g),
                paste(sprintf("P(Y >= %d) = %.6f (crude %.4f)", THRESHOLDS, rk[g, ], crude[g, ]),
                      collapse = " ; ")))
  }
  say("  risk differences (older minus younger):")
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

## --- predict.clm の返り値を手計算で裏づける（06 と同じ。両モデル・両閾値）----
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
        if (!is.finite(manual_max_diff[[m]])) "手計算の照合を行えなかった"
        else sprintf("max |diff| = %.3e", manual_max_diff[[m]]))
}

## --- 5 節・5b 節のオッズ比と、θ_j(g) から求めた値が一致すること --------------
.d5 <- vapply(seq_len(nrow(sens_ppo_odds_ratio)), function(i)
  abs(sens_ppo_odds_ratio$beta[i] -
        logor_point_ppo[sens_ppo_odds_ratio$contrast[i],
                        paste0("t", sens_ppo_odds_ratio$threshold[i])]), numeric(1))
check("M2-PPO-S の閾値別対数オッズ比が 5b 節の表と一致する（参照 G2 の当てはめを含む）",
      nrow(sens_ppo_odds_ratio) == length(RDT_KEYS) && all(.d5 < 1e-6),
      sprintf("max |diff| = %.3e", max(.d5)))
.d5po <- abs(sens_odds_ratio$beta - logor_point_po[sens_odds_ratio$contrast])
check("M2-S の共通対数オッズ比が 5 節の表と一致する（参照 G2 の当てはめを含む）",
      all(.d5po < 1e-6), sprintf("max |diff| = %.3e", max(.d5po)))


# -----------------------------------------------------------------------------
# 7. ブートストラップ（§8.4 と同じ仕様。年齢群で層化した患者単位の復元抽出）
#    各回で M2-PPO と M2 の2本を当てはめ、両閾値の標準化リスクと、
#    オッズ比（M2-PPO の閾値別、M2 の共通）を記録する（06 と同じ）。
# -----------------------------------------------------------------------------

rule("7. Stratified nonparametric bootstrap (sensitivity S)")

idx_by_g <- lapply(AGE_LEVELS, function(g) which(DAT$age_group == g))
names(idx_by_g) <- AGE_LEVELS

set.seed(BOOT_SEED_S)
BOOT_IDX <- vector("list", NBOOT)
for (b in seq_len(NBOOT)) {
  BOOT_IDX[[b]] <- unlist(lapply(AGE_LEVELS, function(g)
    sample(idx_by_g[[g]], length(idx_by_g[[g]]), replace = TRUE)), use.names = FALSE)
}
say(sprintf("  %d resamples generated (stratified by age group; sizes fixed at %s)",
            NBOOT, paste(sprintf("%s=%d", AGE_LEVELS,
                                 vapply(idx_by_g, length, integer(1))), collapse = ", ")))
check("再標本の大きさが各回で観測値に固定されている",
      all(vapply(BOOT_IDX, length, integer(1)) == N_S), sprintf("N = %d", N_S))

#' 当てはめの成否（収束と発散）を判定する（06 と同じ）
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

## --- 1回分の処理（06 の boot_one と同じ。アウトカムだけが S の複合アウトカム）--
## 失敗の条件（§8.4）：収束しない／係数・標準誤差の発散／共変量の水準の欠落／
## いずれかの閾値の上または下に観測が存在しない。M2-PPO では、3区分のいずれかに
## 観測がない回も失敗とする。
boot_one <- function(idx) {
  out <- list(ok_ppo = FALSE, ok_po = FALSE,
              risk_ppo = stats::setNames(rep(NA_real_, length(RISK_KEYS)), RISK_KEYS),
              risk_po  = stats::setNames(rep(NA_real_, length(RISK_KEYS)), RISK_KEYS),
              logor_ppo = stats::setNames(rep(NA_real_, length(RDT_KEYS)), RDT_KEYS),
              logor_po  = stats::setNames(rep(NA_real_, length(RD_KEYS)), RD_KEYS),
              j_used = stats::setNames(rep(NA_real_, length(THRESHOLDS)), THR_KEYS),
              msg_ppo = "", msg_po = "")

  d <- DAT[idx, , drop = FALSE]
  y <- d[[OUTCOME_S]]
  if (!all(vapply(THRESHOLDS, function(t) any(y >= t) && any(y < t), logical(1)))) {
    out$msg_ppo <- out$msg_po <- "no observation on one side of a threshold"
    return(out)
  }
  for (v in c("age_group", "sex", "class", "support_in", "period")) {
    if (is.factor(DAT[[v]])) {
      d[[v]] <- factor(as.character(d[[v]]), levels = levels(DAT[[v]]))
      if (any(table(d[[v]]) == 0L)) {
        out$msg_ppo <- out$msg_po <- paste0("empty level in ", v)
        return(out)
      }
    }
  }
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

  ## --- M2（全水準、比例オッズ）-----------------------------------------------
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
      r <- tryCatch(std_risk(f, pop, j_chr, AGE_LEVELS, ylev_chr), error = function(e) e)
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

t_start <- proc.time()[["elapsed"]]
if (PAR_CORES > 1L && requireNamespace("parallel", quietly = TRUE)) {
  say(sprintf("  running in parallel on %d cores ...", PAR_CORES))
  cl <- parallel::makePSOCKcluster(PAR_CORES)
  boot_res <- tryCatch({
    parallel::clusterEvalQ(cl, { library(ordinal); library(splines); NULL })
    parallel::clusterExport(cl, varlist = c(
      "DAT", "AGE_LEVELS", "OUTCOME_S", "THRESHOLDS", "THR_KEYS", "M2_FORMULA_TXT",
      "M2_PPO_FORMULA_TXT", "M2_PPO_NOMINAL_F", "PPO_LEVELS_CHR", "J_PPO_CHR",
      "ctrl", "CONTRASTS", "RD_KEYS", "RISK_KEYS", "RDT_KEYS", "STD_POP_IN_BOOT",
      "DIVERGE_COEF", "DIVERGE_SE", "IS_SPLINE_TERM", "pick_j", "conv_code",
      "age_nominal_names", "is_nominal_fit", "cut_names", "theta_jg", "fit_status",
      "make_ppo_levels", "make_ppo_outcome", "flat", "logor_ppo", "logor_po",
      "PREDICT_METHOD", "std_risk", "std_risk_clm", "std_risk_clm_fast",
      "risk_to_rd", "boot_one"), envir = environment())
    parallel::parLapply(cl, BOOT_IDX, boot_one)
  }, error = function(e) { say("  parallel error: ", conditionMessage(e)); NULL })
  try(parallel::stopCluster(cl), silent = TRUE)
  if (is.null(boot_res)) stop("並列実行に失敗した。PAR_CORES <- 1L で再実行すること。")
} else {
  if (PAR_CORES > 1L) say("  parallel が使えないため逐次で実行する。")
  say("  running sequentially ...")
  boot_res <- vector("list", NBOOT)
  for (b in seq_len(NBOOT)) {
    boot_res[[b]] <- boot_one(BOOT_IDX[[b]])
    if (b %% PROGRESS_EVERY == 0L || b == NBOOT) {
      el <- proc.time()[["elapsed"]] - t_start
      say(sprintf("    %5d / %d  (%5.1f s elapsed, ETA %5.1f s)",
                  b, NBOOT, el, el / b * (NBOOT - b)))
    }
  }
}
t_boot <- proc.time()[["elapsed"]] - t_start
say(sprintf("  bootstrap finished in %.1f s (%.3f s per resample)",
            t_boot, t_boot / max(NBOOT, 1L)))


# -----------------------------------------------------------------------------
# 8. 区間の算出（パーセンタイル法、§8.4）
# -----------------------------------------------------------------------------

rule("8. Percentile confidence intervals")

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
say(sprintf("  successful resamples : M2-PPO-S %d / %d (failed %d, %.2f%%)",
            n_ok_ppo, NBOOT, NBOOT - n_ok_ppo, 100 * (NBOOT - n_ok_ppo) / NBOOT))
say(sprintf("                         M2-S     %d / %d (failed %d, %.2f%%)",
            n_ok_po, NBOOT, NBOOT - n_ok_po, 100 * (NBOOT - n_ok_po) / NBOOT))

fail_reason <- function(msgs, ok) {
  m <- msgs[!ok]; m <- m[nzchar(m)]
  if (!length(m)) return(data.frame(reason = character(0), n = integer(0),
                                    stringsAsFactors = FALSE))
  m  <- sub(":.*$", "", m)
  tb <- sort(table(m), decreasing = TRUE)
  data.frame(reason = names(tb), n = as.integer(tb), stringsAsFactors = FALSE)
}
fr_ppo <- fail_reason(vapply(boot_res, function(x) x$msg_ppo, character(1)), ok_ppo)
fr_po  <- fail_reason(vapply(boot_res, function(x) x$msg_po,  character(1)), ok_po)
if (nrow(fr_ppo)) { say("\n  failure reasons (M2-PPO-S):"); print(fr_ppo, row.names = FALSE) }
if (nrow(fr_po))  { say("\n  failure reasons (M2-S):");     print(fr_po,  row.names = FALSE) }

say("\n  j used in the resamples (M2-S):")
for (tk in THR_KEYS) {
  say("  ", tk, ":")
  print(table(j_used[, tk], useNA = "ifany"))
}
check("再標本の j が S の本体の j と一致している（M2-S、閾値ごと）",
      all(vapply(THR_KEYS, function(tk)
        all(is.na(j_used[, tk]) | j_used[, tk] == J_PO_S[[tk]]), logical(1))),
      paste(vapply(THR_KEYS, function(tk)
        sprintf("%s: j = %s ; different in %d resamples", tk, J_PO_S_CHR[[tk]],
                sum(!is.na(j_used[, tk]) & j_used[, tk] != J_PO_S[[tk]])), character(1)),
        collapse = " / "))
check("M2-PPO-S の失敗が 5% 以下（§8.4）",
      (NBOOT - n_ok_ppo) / NBOOT <= FAIL_WARN_RATE,
      sprintf("failed %d / %d = %.2f%%。超えた場合は区間の解釈を保留する旨を記す",
              NBOOT - n_ok_ppo, NBOOT, 100 * (NBOOT - n_ok_ppo) / NBOOT))
check("M2-S の失敗が 5% 以下（§8.4）",
      (NBOOT - n_ok_po) / NBOOT <= FAIL_WARN_RATE,
      sprintf("failed %d / %d = %.2f%%", NBOOT - n_ok_po, NBOOT,
              100 * (NBOOT - n_ok_po) / NBOOT))
check("失敗した回を補充していない（§8.4-2）", TRUE, "追加の再標本化は行っていない")


# -----------------------------------------------------------------------------
# 9. 出力表の組み立て
# -----------------------------------------------------------------------------

rule("9. Assembling the output tables")

## --- S の E-4：年齢群別の標準化リスク（M2-PPO、閾値 × 年齢群）-----------------
sens_std_risk <- do.call(rbind, lapply(seq_along(THRESHOLDS), function(i) {
  t <- THRESHOLDS[i]; tk <- THR_KEYS[i]
  do.call(rbind, lapply(AGE_LEVELS, function(g) {
    ci <- pct_ci(risk_boot_ppo[ok_ppo, paste(g, tk, sep = "_")])
    data.frame(
      analysis        = "sensitivity S",
      threshold       = t,
      age_group       = g,
      age_group_label = relabel_levels("age_group", g),
      n               = as.integer(n_by_g[[g]]),
      estimand        = sprintf("Standardised P(composite mFIM at discharge >= %d)", t),
      estimate        = risk_point[g, tk],
      ci_lower        = unname(ci[["lower"]]),
      ci_upper        = unname(ci[["upper"]]),
      percent         = 100 * risk_point[g, tk],
      percent_lower   = 100 * unname(ci[["lower"]]),
      percent_upper   = 100 * unname(ci[["upper"]]),
      conf_level      = CONF_LEVEL,
      ci_method       = "bootstrap percentile, unadjusted",
      n_boot_used     = n_ok_ppo,
      model           = MODEL_PPO_S,
      crude_proportion = unname(crude[g, tk]),
      stringsAsFactors = FALSE
    )
  }))
}))

## --- S の E-4：3対比のリスク差（M2-PPO、閾値 × 対比）--------------------------
sens_rd <- do.call(rbind, lapply(seq_along(THRESHOLDS), function(i) {
  t <- THRESHOLDS[i]; tk <- THR_KEYS[i]
  do.call(rbind, lapply(CONTRASTS, function(ct) {
    k  <- ct$key
    ci <- pct_ci(rd_boot_ppo[ok_ppo, paste(k, tk, sep = "_")])
    data.frame(
      analysis         = "sensitivity S",
      threshold        = t,
      contrast         = k,
      label            = contrast_label(ct$ref, ct$cmp),
      group_reference  = ct$ref,
      group_comparison = ct$cmp,
      role             = ct$role,
      estimand         = sprintf("Risk difference in standardised P(composite mFIM >= %d), older minus younger",
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
      model            = MODEL_PPO_S,
      p_value          = NA_real_,
      stringsAsFactors = FALSE
    )
  }))
}))

## --- S の E-3：閾値別オッズ比にブートストラップ区間を添える（表に出さない）-----
sens_ppo_odds_ratio$boot_ci_lower <- NA_real_
sens_ppo_odds_ratio$boot_ci_upper <- NA_real_
for (i in seq_len(nrow(sens_ppo_odds_ratio))) {
  key <- paste(sens_ppo_odds_ratio$contrast[i], paste0("t", sens_ppo_odds_ratio$threshold[i]),
               sep = "_")
  ci  <- pct_ci(logor_boot_ppo[ok_ppo, key])
  sens_ppo_odds_ratio$boot_ci_lower[i] <- exp(unname(ci[["lower"]]))
  sens_ppo_odds_ratio$boot_ci_upper[i] <- exp(unname(ci[["upper"]]))
}
sens_ppo_odds_ratio$boot_ci_method <- "bootstrap percentile, unadjusted"
sens_ppo_odds_ratio$n_boot_used    <- n_ok_ppo
sens_ppo_odds_ratio <- sens_ppo_odds_ratio[order(sens_ppo_odds_ratio$threshold,
                                                 match(sens_ppo_odds_ratio$contrast, RD_KEYS)), ,
                                           drop = FALSE]
rownames(sens_ppo_odds_ratio) <- NULL

## --- profile 区間とブートストラップ区間の突き合わせ（06 と同じ。記録用）--------
rel_shift <- function(pl, pu, bl, bu) {
  pmax(abs(log(bl) - log(pl)), abs(log(bu) - log(pu))) / (log(pu) - log(pl))
}
.po_boot <- t(vapply(sens_odds_ratio$contrast, function(k) {
  ci <- pct_ci(logor_boot_po[ok_po, k]); exp(c(ci[["lower"]], ci[["upper"]]))
}, numeric(2)))
sens_or_ci_compare <- rbind(
  data.frame(model = MODEL_PPO_S, threshold = sens_ppo_odds_ratio$threshold,
             contrast = sens_ppo_odds_ratio$contrast, role = sens_ppo_odds_ratio$role,
             odds_ratio = sens_ppo_odds_ratio$odds_ratio,
             profile_lower = sens_ppo_odds_ratio$ci_lower,
             profile_upper = sens_ppo_odds_ratio$ci_upper,
             profile_source = "own implementation (as in 05)",
             boot_lower = sens_ppo_odds_ratio$boot_ci_lower,
             boot_upper = sens_ppo_odds_ratio$boot_ci_upper,
             n_boot_used = n_ok_ppo, used_for_check = TRUE,
             stringsAsFactors = FALSE),
  data.frame(model = MODEL_PO_S, threshold = NA_real_,
             contrast = sens_odds_ratio$contrast, role = sens_odds_ratio$role,
             odds_ratio = sens_odds_ratio$odds_ratio,
             profile_lower = sens_odds_ratio$ci_lower, profile_upper = sens_odds_ratio$ci_upper,
             profile_source = "ordinal::confint",
             boot_lower = unname(.po_boot[, 1]), boot_upper = unname(.po_boot[, 2]),
             n_boot_used = n_ok_po, used_for_check = FALSE,
             stringsAsFactors = FALSE)
)
sens_or_ci_compare$relative_shift <- rel_shift(sens_or_ci_compare$profile_lower,
                                               sens_or_ci_compare$profile_upper,
                                               sens_or_ci_compare$boot_lower,
                                               sens_or_ci_compare$boot_upper)
sens_or_ci_compare$within_tolerance <- sens_or_ci_compare$relative_shift <= BOOT_PROFILE_TOL
sens_or_ci_compare$reported <- FALSE
say("  odds ratios (sensitivity S, NOT reported): profile likelihood vs bootstrap percentile")
print(sens_or_ci_compare[, c("model", "threshold", "contrast", "odds_ratio", "profile_lower",
                             "boot_lower", "profile_upper", "boot_upper", "relative_shift")],
      row.names = FALSE, digits = 4)
.cmp_ppo <- sens_or_ci_compare[sens_or_ci_compare$used_for_check, , drop = FALSE]
check(sprintf(paste0("M2-PPO-S の閾値別オッズ比で、profile 区間とブートストラップ区間に",
                     "大きな差がない（端のずれ ≤ profile 区間の幅の %.0f%%、対数尺度。記録用）"),
              100 * BOOT_PROFILE_TOL),
      all(is.finite(.cmp_ppo$relative_shift)) && all(.cmp_ppo$within_tolerance),
      if (any(is.finite(.cmp_ppo$relative_shift)))
        sprintf("max = %.1f%%（%s）", 100 * max(.cmp_ppo$relative_shift, na.rm = TRUE),
                with(.cmp_ppo[which.max(.cmp_ppo$relative_shift), ],
                     sprintf("%s, y >= %g", contrast, threshold)))
      else "比較できなかった")

## --- M2-S（比例オッズ）由来の標準化量（比較用。補足表2 には載せない）----------
sens_po_tab <- rbind(
  do.call(rbind, lapply(seq_along(THRESHOLDS), function(i) {
    t <- THRESHOLDS[i]; tk <- THR_KEYS[i]
    do.call(rbind, lapply(AGE_LEVELS, function(g) {
      ci <- pct_ci(risk_boot_po[ok_po, paste(g, tk, sep = "_")])
      data.frame(quantity = "standardised risk", threshold = t, key = g,
                 label = relabel_levels("age_group", g),
                 estimate = risk_point_po[g, tk],
                 ci_lower = unname(ci[["lower"]]), ci_upper = unname(ci[["upper"]]),
                 j_used = J_PO_S_CHR[[tk]],
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
                 j_used = J_PO_S_CHR[[tk]],
                 stringsAsFactors = FALSE)
    }))
  }))
)
sens_po_tab$analysis    <- "sensitivity S"
sens_po_tab$conf_level  <- CONF_LEVEL
sens_po_tab$ci_method   <- "bootstrap percentile, unadjusted"
sens_po_tab$n_boot_used <- n_ok_po
sens_po_tab$model       <- MODEL_PO_S

## --- M2-PPO-S と M2-S の食い違い（§8.6「事前の規則」に準じて数値で示す）--------
sens_diff <- do.call(rbind, lapply(seq_along(THRESHOLDS), function(i) {
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
sens_diff$difference <- sens_diff$po_model - sens_diff$ppo_model
say("\n  M2-PPO-S vs M2-S (proportional odds) standardised quantities:")
print(sens_diff, row.names = FALSE, digits = 4)
say("  補足表2 に載せるのは M2-PPO 由来の値である（plan.summary §7-2、主解析と同じ）。")

boot_summary <- data.frame(
  model = c(MODEL_PPO_S, MODEL_PO_S),
  n_boot = c(NBOOT, NBOOT),
  n_used = c(n_ok_ppo, n_ok_po),
  n_failed = c(NBOOT - n_ok_ppo, NBOOT - n_ok_po),
  percent_failed = c(100 * (NBOOT - n_ok_ppo) / NBOOT, 100 * (NBOOT - n_ok_po) / NBOOT),
  exceeds_5_percent = c((NBOOT - n_ok_ppo) / NBOOT > FAIL_WARN_RATE,
                        (NBOOT - n_ok_po) / NBOOT > FAIL_WARN_RATE),
  seed = c(BOOT_SEED_S, BOOT_SEED_S),
  stratified_by = c("age group", "age group"),
  standardisation_population = c(STD_POP_IN_BOOT, STD_POP_IN_BOOT),
  ci_type = c(CI_TYPE, CI_TYPE),
  replacement_of_failures = c("none (failed resamples are dropped)",
                              "none (failed resamples are dropped)"),
  seconds = c(t_boot, NA_real_),   # 2つのモデルを同じ回で当てはめた合計
  stringsAsFactors = FALSE
)


# -----------------------------------------------------------------------------
# 10. 補足表2（主解析とSを並べる。§8.8「報告」）
#     E-2 は 07_1 の出力、E-4 は 06（主解析）とこのスクリプト（S）の出力を使う。
#     §7-1・§7-2：E-4 は M2-PPO 由来の 70 点・77 点の標準化確率とRD。
#     二値モデルの行は削除した。threshold 列で閾値を見分ける（E-2 の行は NA）。
# -----------------------------------------------------------------------------

rule("10. Supplementary Table 2 (main analysis vs sensitivity S)")

blk <- function(analysis, block, threshold, key, label, role, n, quantity, est, lo, hi,
                scale, ci_method, model) {
  data.frame(analysis = analysis, block = block, threshold = threshold,
             row_key = key, row_label = label,
             role = role, n = n, quantity = quantity, estimate = est,
             ci_lower = lo, ci_upper = hi, scale = scale,
             ci_method = ci_method, model = model, conf_level = CONF_LEVEL,
             p_value = NA_real_, stringsAsFactors = FALSE)
}

supp2_rows <- list()

## --- 主解析（03・06 の出力から転記）----------------------------------------
main_m0_tab <- NULL
if (!is.null(sens_e2) && exists("sens_m0_compare", envir = sens_e2, inherits = FALSE)) {
  .c <- get("sens_m0_compare", envir = sens_e2)
  if (nrow(.c)) main_m0_tab <- .c
}
if (!is.null(main_m0_tab)) {
  supp2_rows[[length(supp2_rows) + 1L]] <- blk(
    "main", "Relative effect (E-2)", NA_real_, main_m0_tab$contrast, main_m0_tab$label,
    main_m0_tab$role, NA_integer_,
    "P(older > younger) + 0.5 x P(tie)", main_m0_tab$main_estimate,
    main_m0_tab$main_ci_lower, main_m0_tab$main_ci_upper, "probability",
    "equal-tailed permutation, logit scale, unadjusted", "Brunner-Munzel")
  check("主解析の E-2 を補足表2 に取り込めた", TRUE, "07_1 の対照表から転記した")
} else {
  check("主解析の E-2 を補足表2 に取り込めた", FALSE,
        if (is.null(sens_e2)) "07_1_sens_ranktest.rda が無い。03 → 07_1 の順に実行すること"
        else "07_1 が 03_m0_ranktest.rda を読めていない。03 → 07_1 の順に実行すること")
}

if (main_e4_is_s7 && exists("rd_tab", envir = main_e4, inherits = FALSE)) {
  .r <- get("std_risk_tab", envir = main_e4)
  .d <- get("rd_tab",       envir = main_e4)
  .r <- .r[order(.r$threshold, match(.r$age_group, AGE_LEVELS)), , drop = FALSE]
  .d <- .d[order(.d$threshold, match(.d$contrast, RD_KEYS)), , drop = FALSE]
  supp2_rows[[length(supp2_rows) + 1L]] <- blk(
    "main", "Standardised risk (E-4)", .r$threshold, .r$age_group, .r$age_group_label, "",
    .r$n, .r$estimand, .r$estimate, .r$ci_lower, .r$ci_upper, "probability",
    .r$ci_method, .r$model)
  supp2_rows[[length(supp2_rows) + 1L]] <- blk(
    "main", "Risk difference (E-4)", .d$threshold, .d$contrast, .d$label, .d$role,
    NA_integer_, .d$estimand, .d$estimate, .d$ci_lower, .d$ci_upper,
    "probability difference", .d$ci_method, .d$model)
  check("主解析の E-4（M2-PPO、70 点・77 点）を補足表2 に取り込めた",
        setequal(unique(.r$threshold), THRESHOLDS) && setequal(unique(.d$threshold), THRESHOLDS),
        sprintf("標準化リスク %d 行・RD %d 行", nrow(.r), nrow(.d)))
} else {
  check("主解析の E-4（M2-PPO、70 点・77 点）を補足表2 に取り込めた", FALSE,
        "§7 版の 06_standardize_boot.rda が無い。先に §7 版の 06 を実行すること")
}

## --- 感度分析S -------------------------------------------------------------
if (!is.null(sens_e2) && exists("sens_m0", envir = sens_e2, inherits = FALSE)) {
  .s <- get("sens_m0", envir = sens_e2)
  supp2_rows[[length(supp2_rows) + 1L]] <- blk(
    "sensitivity S", "Relative effect (E-2)", NA_real_, .s$contrast, .s$label,
    .s$role, .s$n_reference + .s$n_comparison,
    "P(older > younger) + 0.5 x P(tie)", .s$relative_effect,
    .s$ci_lower, .s$ci_upper, "probability", .s$ci_method, "Brunner-Munzel")
  check("S の E-2 を補足表2 に取り込めた", TRUE,
        sprintf("07_1_sens_ranktest.rda から %d 行", nrow(.s)))
} else {
  check("S の E-2 を補足表2 に取り込めた", FALSE,
        "07_1_sens_ranktest.rda が無い。先に 07_1 を実行すること")
}

supp2_rows[[length(supp2_rows) + 1L]] <- blk(
  "sensitivity S", "Standardised risk (E-4)", sens_std_risk$threshold, sens_std_risk$age_group,
  sens_std_risk$age_group_label, "", sens_std_risk$n, sens_std_risk$estimand,
  sens_std_risk$estimate, sens_std_risk$ci_lower, sens_std_risk$ci_upper,
  "probability", sens_std_risk$ci_method, sens_std_risk$model)
supp2_rows[[length(supp2_rows) + 1L]] <- blk(
  "sensitivity S", "Risk difference (E-4)", sens_rd$threshold, sens_rd$contrast, sens_rd$label,
  sens_rd$role, NA_integer_, sens_rd$estimand, sens_rd$estimate,
  sens_rd$ci_lower, sens_rd$ci_upper, "probability difference",
  sens_rd$ci_method, sens_rd$model)

table_supp2 <- do.call(rbind, supp2_rows)
rownames(table_supp2) <- NULL
say("\n  Supplementary Table 2 (main vs sensitivity S):")
print(table_supp2[, c("analysis", "block", "threshold", "row_key", "estimate",
                      "ci_lower", "ci_upper")], row.names = FALSE, digits = 4)

check("補足表2 に p 値を置いていない（§8.5、§8.8）",
      all(is.na(table_supp2$p_value)), "")
check("補足表2 に S のオッズ比（E-3）を入れていない（§8.8）",
      !any(grepl("odds ratio", table_supp2$quantity, ignore.case = TRUE)), "")
check("補足表2 に二値モデル・閾値65 の行がない（§7-1、§7-2）",
      !any(grepl("binary|65", paste(table_supp2$block, table_supp2$quantity),
                 ignore.case = TRUE)) &&
        all(is.na(table_supp2$threshold) | table_supp2$threshold %in% THRESHOLDS), "")
.e4 <- table_supp2[table_supp2$block != "Relative effect (E-2)", , drop = FALSE]
check("補足表2 の E-4 の行が S で 閾値 × (3群 + 3対比) = 12 行ある",
      sum(.e4$analysis == "sensitivity S") == length(THRESHOLDS) * (length(AGE_LEVELS) + length(RD_KEYS)),
      sprintf("S の E-4 = %d 行", sum(.e4$analysis == "sensitivity S")))
if (any(.e4$analysis == "main")) {
  .km <- sort(paste(.e4$block, .e4$threshold, .e4$row_key)[.e4$analysis == "main"])
  .ks <- sort(paste(.e4$block, .e4$threshold, .e4$row_key)[.e4$analysis == "sensitivity S"])
  check("補足表2 で主解析と S の E-4 の行（ブロック・閾値・行）が1対1に対応する",
        identical(.km, .ks), sprintf("main %d 行 / S %d 行", length(.km), length(.ks)))
}

table_supp2_meta <- data.frame(
  item = c("base_population_A4", "sensitivity_S_population",
           "removed_E3_b_c", "n_worst_score", "worst_score", "measured_range",
           "m2_S_complete_case_n", "n_worst_score_m2_S", "outcome", "threshold",
           "ppo_categories", "j_used_ppo", "j_used_po",
           "nboot", "nboot_used_ppo", "nboot_used_po", "bootstrap_seed",
           "permutation_settings",
           "ridit_reference_population", "model", "model_nominal", "model_common_or",
           "odds_ratios_S"),
  value = c(as.character(nrow(dat_base)), as.character(nrow(sens)),
            as.character(nrow(dat_base) - nrow(sens)),
            as.character(n_worst), as.character(WORST_SCORE),
            paste(MFIM_RANGE, collapse = "-"),
            as.character(N_S), as.character(n_worst_cc), OUTCOME_S,
            paste(THRESHOLDS, collapse = ";"), paste(PPO_LEVELS_CHR, collapse = ";"),
            paste(J_PPO_CHR, collapse = ";"), paste(J_PO_S_CHR, collapse = ";"),
            as.character(NBOOT), as.character(n_ok_ppo), as.character(n_ok_po),
            as.character(BOOT_SEED_S),
            if (!is.null(sens_e2) && exists("NPERM", envir = sens_e2, inherits = FALSE))
              sprintf("nperm = %s, seed base = %s",
                      as.character(get("NPERM", envir = sens_e2)),
                      as.character(get("SEED_S", envir = sens_e2)))
            else "see 07_1_sens_ranktest.R",
            "A4 (base population), identical to the main analysis",
            M2_PPO_FORMULA_TXT, M2_PPO_NOMINAL_TXT, M2_FORMULA_TXT,
            "fitted but not reported (plan 8.8)"),
  stringsAsFactors = FALSE
)

table_supp2_footnotes <- data.frame(
  n = seq_len(8),
  footnote = c(
    sprintf(paste("Sensitivity analysis S assigns a motor FIM score of %d, one point below",
                  "the lowest possible observed score (%d), to in-hospital transfer,",
                  "termination (death etc.) and to discharge scores that were not recorded",
                  "because of acute deterioration, so that these patients rank below every",
                  "patient with an observed score (worst-rank)."),
            WORST_SCORE, MFIM_RANGE[1]),
    paste("The population is the base population (A4) after removing patients whose",
          "discharge motor FIM was missing for reasons unrelated to the patient's",
          "condition (categories (b) and (c) of the pre-specified rule)."),
    paste("The estimand of S differs from that of the main analysis. The two sets of",
          "numbers are two scenarios for the handling of dropout events, not a lower",
          "and an upper bound: the true value does not lie between them."),
    paste0("Standardised probabilities of reaching ", paste(THRESHOLDS, collapse = " and "),
           " points and risk differences are from a partial proportional odds model in ",
           "which the effect of age group was allowed to differ between the thresholds ",
           "(outcome grouped as ", paste(PPO_LEVELS_CHR, collapse = ", "), " points; ",
           "the assigned score of ", WORST_SCORE, " falls in the lowest category) and ",
           "covariate effects were common to both thresholds, as in the main analysis."),
    paste("Intervals are unadjusted 95% intervals of the same construction as the main",
          "analysis: equal-tailed permutation intervals on the logit scale for the",
          "relative effect, and bootstrap percentile intervals for the standardised",
          "probabilities and risk differences."),
    paste("No p-value is reported for the sensitivity analysis. The single confirmatory",
          "hypothesis is the main analysis of G3 vs G1."),
    paste("The threshold-specific and common odds ratios were fitted on the S",
          "population but are not reported (plan 8.8)."),
    paste("Ridit scores and spline knots were fixed in the base population (A4)",
          "and were not recomputed here.")
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
w(sens_analysis_n,       "table_sens_analysis_n.csv")
w(sens_missing,          "table_sens_missing.csv")
w(sens_odds_ratio,       "table_sens_odds_ratio.csv")
w(sens_ppo_odds_ratio,   "table_sens_ppo_odds_ratio.csv")
w(sens_fit_info,         "table_sens_fit.csv")
w(sens_std_risk,         "table_sens_standardized_risk.csv")
w(sens_rd,               "table_sens_risk_difference.csv")
w(sens_po_tab,           "table_sens_po_standardized.csv")
w(sens_diff,             "table_sens_ppo_vs_po.csv")
w(sens_or_ci_compare,    "table_sens_or_ci_profile_vs_boot.csv")
w(boot_summary,          "table_sens_boot_summary.csv")
w(table_supp2,           "table_supp2.csv")
w(table_supp2_meta,      "table_supp2_meta.csv")
w(table_supp2_footnotes, "table_supp2_footnotes.csv")
w(checks,                "table_checks_07_2.csv")

## risk_point / rd_point は M2-PPO-S（補足表2）、*_po は M2-S（比較用）。
## いずれも 行 = 年齢群または対比、列 = 閾値（"t70"、"t77"）の行列（06 と同じ形）。
dat_m2_S <- DAT
save(sens_analysis_n, sens_missing, sens_odds_ratio, sens_ppo_odds_ratio, sens_fit_info,
     sens_std_risk, sens_rd, sens_po_tab, sens_diff, sens_or_ci_compare, boot_summary,
     table_supp2, table_supp2_meta, table_supp2_footnotes,
     dat_m2_S, fitS_G1, fitS_G2, fitS_ppo_G1, fitS_ppo_G2, ppo_impl_check_S, ppo_cells,
     risk_point, rd_point, risk_point_po, rd_point_po, logor_point_ppo, logor_point_po,
     risk_boot_ppo, rd_boot_ppo, logor_boot_ppo, risk_boot_po, rd_boot_po, logor_boot_po,
     ok_ppo, ok_po, j_used, J_PO_S, J_PO_S_CHR, J_PPO_CHR, THRESHOLDS, THR_KEYS,
     NBOOT, BOOT_SEED_S, fr_ppo, fr_po, n_worst, n_worst_cc, WORST_SCORE,
     MFIM_RANGE, MFIM_RANGE_S, checks,
     file = file.path(OUT_DIR, "07_2_sens_clm_boot.rda"))
say("  written: ", file.path(OUT_DIR, "07_2_sens_clm_boot.rda"))

say("\n--- sessionInfo ---")
print(utils::sessionInfo())

say("\nDone. sensitivity S population = ", nrow(sens),
    " ; M2-S complete cases = ", N_S,
    " ; worst score (", WORST_SCORE, ") assigned = ", n_worst,
    " ; bootstrap M2-PPO-S ", n_ok_ppo, " / ", NBOOT, " usable, M2-S ", n_ok_po,
    " / ", NBOOT, " usable ; thresholds = ", paste(THRESHOLDS, collapse = ", "),
    " ; M2-S j = ", paste(J_PO_S_CHR, collapse = ", "))

.log_close()
