# =============================================================================
# 02_preprocess.R
#   解析計画書 §5.1「対象の選抜」（plan.summary.txt §2）の実装
#   入力 : Alldata（ワークスペース上に存在する前提。2,919 行 74 列）
#   出力 : output/table_cleaning.csv            クリーニングの記録（§4）
#          output/table_flow_selection.csv      フロー図の数値（A1〜A4, E1〜E3）
#          output/table_flow_by_agegroup.csv    同・年齢群別（Fig.1 用、§5.1）
#          output/table_disposition_base.csv    基準集団の退院先の内訳
#          output/table_exclusion_overlap.csv   E1〜E3 の重なり（§5.1）
#          output/table_e3_reason.csv           E3 の欠測理由の区分（補足表1、§5.1）
#          output/02_preprocess.rda             dat_base / dat_main / dat_sens / flow
#          output/log_02_preprocess.txt         実行ログ（点検結果を含む）
#   依存 : 01_labels.R（同じフォルダに置く。表示用の文字列はすべてそちらに置く）
#   文字コード : UTF-8（Windows の R では source(..., encoding = "UTF-8")）
#   注意 : §13-4 により、期待例数のハードコードは置かない。
#          例数は検証せず、集計して報告するだけにとどめる。
#   改訂 : 2026-09-21 plan.summary.txt §7-3 により、感度分析の最低点
#          （WORST_SCORE）を 13 → 12 に変更した。変わるのは dat_sens の
#          mFIM_out_S だけで、dat_base・dat_main は変わらない。
# =============================================================================


# -----------------------------------------------------------------------------
# 0. 設定
# -----------------------------------------------------------------------------

OUT_DIR <- "output"

## --- 対象期間（§4、§5.1 A1）------------------------------------------------
PERIOD_START <- as.Date("2017-04-01")
PERIOD_END   <- as.Date("2024-03-31")

## --- 適格年齢（§5.1 A2）。age は入院時年齢（著者確認済み）------------------
AGE_MIN <- 65

## --- A4：抽出時点で在院中の患者（§5.1、§14-4）------------------------------
## 著者確認：Alldata は全例が退院確定済みであるため A4 = A3 となる。
## 将来「在院中」を表す値が現れた場合に備え、判定規則は設定として外に出す。
IN_HOSPITAL_LEVELS <- character(0)  # 在院中を表す disposition の値（無ければ空）
IN_HOSPITAL_NA     <- FALSE         # disposition 欠測を在院中とみなすか

## --- 除外の定義（§5.1 E1・E2）----------------------------------------------
## §14-1 の未決事項への対応：著者の確認により、「病院・診療所へ転院」と
## 「医療機関」の両方を E1（在院中の転院）として扱う。
DISPO_TRANSFER <- c("病院・診療所へ転院", "医療機関")  # E1
DISPO_DEATH    <- c("終了（死亡等）")                   # E2

## --- 除外の優先順位（§5.1：E2 ＞ E1 ＞ E3）----------------------------------
EXCL_PRIORITY <- c("E2", "E1", "E3")

## --- 年齢群（§5.3）。フロー図の年齢群構成の集計にのみ用いる ----------------
## 英語名（pre-old / old / oldest-old 等）は原典未確認（§14-11）のため、
## ラベルには年齢の範囲だけを書く。
AGE_BREAKS <- c(65, 75, 90)  # G1: 65-74, G2: 75-89, G3: >=90

## --- 値域（§6.2）------------------------------------------------------------
MFIM_RANGE <- c(13, 91)  # 運動FIM 13 項目 × 1〜7 点
CFIM_RANGE <- c( 5, 35)  # 認知FIM  5 項目 × 1〜7 点

## --- 範囲外の値のクリーニング（§4、2026-09-20 に確定）------------------------
## 著者の判断：運動FIM の 12 点は 13 点の入力ミスとみなし、13 に置換する。
## 運動FIM は 13 項目それぞれの下限が 1 点であり、合計 12 点は定義上ありえない。
## 置換の件数は table_cleaning.csv とログに残し、論文の§4（クリーニング）に記す。
## §7-3：この置換は感度分析Sの集団（本スクリプトの 6-2）を作る前に行う。
## そのため WORST_SCORE = 12 としても、dat_sens で 12 点をもつのは最低点を
## 与えた転院・死亡の患者だけになる（本スクリプトの 8 で点検する）。
FIX_MFIM_12   <- TRUE
FIX_MFIM_FROM <- 12
FIX_MFIM_TO   <- 13
MFIM_VARS     <- c("mFIM_in", "mFIM_out")

## --- E3 の欠測理由の分類規則（§5.1、§14-3。2026-09-20 に確定）---------------
## 著者の決定：「転院ならば全て急変」。退院時運動FIMが欠測した例のうち、
##   (a) 急変・状態悪化による未評価   : 転院した例    → 感度分析Sで最低点を与える
##   (b) 記録漏れ・通常退院時の未評価 : 転院でない例  → 感度分析Sの集団から除く
##   (c) その他                       : 置かない
## 除外の優先順位（E2 ＞ E1 ＞ E3）により、転院かつ欠測の例は E1 に割り当てられ、
## §8.8 により E1 にも最低点が与えられる。よって区分(a)と扱いが一致する。
E3_CATEGORY_EXCLUDE <- c("b", "c")  # 感度分析Sの集団から除く区分（§8.8）
## §7-3（2026-09-21）：転院・死亡の患者に 12 点を代入し、退院時に実測で
## 13 点だった患者と区別する（旧設定は 13）。
WORST_SCORE         <- 12           # 複合アウトカムの最低点（§8.8、§7-3）

## 欠測例の区分を決める規則。転院なら (a)、そうでなければ (b)。
classify_e3 <- function(is_transfer) ifelse(is_transfer, "a", "b")


# -----------------------------------------------------------------------------
# 1. 表示用ラベル
#    表示用の文字列は 01_labels.R に一元化する（このスクリプトでは定義しない）。
# -----------------------------------------------------------------------------

if (!file.exists("01_labels.R")) {
  stop("01_labels.R が見つからない。02_preprocess.R は 01_labels.R と同じ",
       "フォルダを作業ディレクトリにして実行すること。")
}
source("01_labels.R", encoding = "UTF-8")

## 退院先の表示。欠測は "(missing)" とする。
## 未定義の水準は 01_labels.R 側で警告が出る（§14-1 の追記漏れの検知）。
show_dispo <- function(x, strict = TRUE) {
  out <- relabel_levels("disposition", x, strict = strict)
  out[is.na(x)] <- "(missing)"
  out
}

## 年齢群（§5.3）。フロー図と補足表1の年齢群構成の集計にのみ用いる。
make_agegroup <- function(age) {
  factor(ifelse(is.na(age), NA_character_,
         ifelse(age >= AGE_BREAKS[3], "G3",
         ifelse(age >= AGE_BREAKS[2], "G2", "G1"))),
         levels = c("G1", "G2", "G3"))
}
count_by_g <- function(x) {
  a <- if (is.data.frame(x)) x$age else x
  if (!length(a)) return(c(G1 = 0L, G2 = 0L, G3 = 0L))
  stats::setNames(as.integer(table(make_agegroup(a))), c("G1", "G2", "G3"))
}


# -----------------------------------------------------------------------------
# 2. 実行の準備（出力先とログ）
# -----------------------------------------------------------------------------

if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

.log_con <- file(file.path(OUT_DIR, "log_02_preprocess.txt"), open = "wt", encoding = "UTF-8")
sink(.log_con, split = TRUE)

## 途中で停止した場合にもコンソールへの出力を必ず戻す
.old_error <- getOption("error")
.log_close <- function() {
  while (sink.number() > 0) sink()
  try(close(.log_con), silent = TRUE)
  options(error = .old_error)
}
options(error = function() .log_close())

say <- function(...) cat(..., "\n", sep = "")
rule <- function(title) say("\n", strrep("-", 74), "\n", title, "\n", strrep("-", 74))

## 点検結果の収集
.checks <- list()
check <- function(name, ok, detail = "") {
  .checks[[length(.checks) + 1L]] <<-
    data.frame(check = name, result = if (isTRUE(ok)) "OK" else "REVIEW",
               detail = detail, stringsAsFactors = FALSE)
  say(sprintf("  [%-6s] %s%s", if (isTRUE(ok)) "OK" else "REVIEW", name,
              if (nzchar(detail)) paste0(" : ", detail) else ""))
  invisible(ok)
}

rule(paste0("02_preprocess.R  ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
say("R : ", R.version.string)


# -----------------------------------------------------------------------------
# 3. 入力の取り込みと型の整備
# -----------------------------------------------------------------------------

rule("3. Input")

if (!exists("Alldata", inherits = TRUE)) {
  stop("Alldata がワークスペースに見つからない。先にデータを読み込むこと。")
}

d0 <- as.data.frame(Alldata, stringsAsFactors = FALSE)
say("Alldata : ", nrow(d0), " rows x ", ncol(d0), " cols")

REQUIRED_VARS <- c("id", "age", "sex", "day_in", "disposition", "class",
                   "support_in", "mFIM_in", "cFIM_in", "mFIM_out")
.miss <- setdiff(REQUIRED_VARS, names(d0))
if (length(.miss)) {
  stop("Alldata に必要な列がない: ", paste(.miss, collapse = ", "))
}

d <- d0[, REQUIRED_VARS, drop = FALSE]

## 型の整備 -------------------------------------------------------------------
as_date_safe <- function(x) {
  if (inherits(x, "Date"))   return(x)
  if (inherits(x, "POSIXt")) return(as.Date(x))
  as.Date(as.character(x), format = "%Y-%m-%d")
}
as_chr <- function(x) {
  x <- as.character(x)
  x[!is.na(x) & !nzchar(trimws(x))] <- NA_character_  # 空白のみは欠測とみなす
  trimws(x)
}

d$id          <- as_chr(d$id)
d$day_in      <- as_date_safe(d$day_in)
d$age         <- suppressWarnings(as.numeric(d$age))
d$sex         <- as_chr(d$sex)
d$disposition <- as_chr(d$disposition)
d$class       <- as_chr(d$class)
d$support_in  <- as_chr(d$support_in)
d$mFIM_in     <- suppressWarnings(as.numeric(d$mFIM_in))
d$cFIM_in     <- suppressWarnings(as.numeric(d$cFIM_in))
d$mFIM_out    <- suppressWarnings(as.numeric(d$mFIM_out))
d$.row        <- seq_len(nrow(d))  # 並べ替えを再現可能にするための元の行番号


# -----------------------------------------------------------------------------
# 3-2. クリーニング（§4）
#      (1) 運動FIM の 12 点 → 13 点（入力ミスとみなす）
#      (2) 退院先の集約（計画で名指しされた3つ以外は「自宅相当に退院」）
#      いずれも選抜と点検の前に行い、件数をログと表に残す。
# -----------------------------------------------------------------------------

rule("3-2. Cleaning")

clean_log <- list()
.rec <- function(target, action, n) {
  clean_log[[length(clean_log) + 1L]] <<-
    data.frame(target = target, action = action, n = n, stringsAsFactors = FALSE)
}

## --- (1) 運動FIM の 12 点を 13 点に置換 -------------------------------------
if (isTRUE(FIX_MFIM_12)) {
  for (v in MFIM_VARS) {
    hit <- !is.na(d[[v]]) & d[[v]] == FIX_MFIM_FROM
    n   <- sum(hit)
    d[[v]][hit] <- FIX_MFIM_TO
    .rec(v, sprintf("%d -> %d (data entry error)", FIX_MFIM_FROM, FIX_MFIM_TO), n)
    say(sprintf("  %-10s %d -> %d : %4d replaced", v, FIX_MFIM_FROM, FIX_MFIM_TO, n))
  }
  check("運動FIM に 13 点未満の値が残っていない",
        all(vapply(MFIM_VARS, function(v)
          all(is.na(d[[v]]) | d[[v]] >= MFIM_RANGE[1]), logical(1))),
        paste(vapply(MFIM_VARS, function(v)
          sprintf("%s: %d", v, sum(!is.na(d[[v]]) & d[[v]] < MFIM_RANGE[1])),
          character(1)), collapse = ", "))
}

## --- (2) 退院先の集約（01_labels.R の collapse_disposition）-----------------
## 在院中を表す値があれば集約の対象から外す（A4 の判定を壊さないため）
.keep_dispo <- unique(c(DISPO_KNOWN, IN_HOSPITAL_LEVELS))
.dispo_raw  <- d$disposition
d$disposition <- collapse_disposition(.dispo_raw, keep = .keep_dispo)

.collapsed <- !is.na(.dispo_raw) & .dispo_raw != d$disposition
.rec("disposition",
     sprintf("collapsed to \"%s\" (%s)", DISPO_COMMUNITY_CODE,
             unname(LEVELS_BY_VAR$disposition[DISPO_COMMUNITY_CODE])),
     sum(.collapsed))
say(sprintf("\n  disposition collapsed to \"%s\" (%s) : %d",
            DISPO_COMMUNITY_CODE,
            unname(LEVELS_BY_VAR$disposition[DISPO_COMMUNITY_CODE]),
            sum(.collapsed)))
if (any(.collapsed)) {
  .from <- sort(table(.dispo_raw[.collapsed]), decreasing = TRUE)
  for (i in seq_along(.from)) {
    say(sprintf("    %-24s %6d", names(.from)[i], as.integer(.from[i])))
  }
}

check("E1・E2 の設定値が 01_labels.R の DISPO_KNOWN と一致する",
      setequal(c(DISPO_TRANSFER, DISPO_DEATH), DISPO_KNOWN),
      paste0("02: ", paste(sort(c(DISPO_TRANSFER, DISPO_DEATH)), collapse = ", "),
             " / 01: ", paste(sort(DISPO_KNOWN), collapse = ", ")))

clean_tab <- do.call(rbind, clean_log)


# -----------------------------------------------------------------------------
# 4. 入力の点検（アウトカムの分布には触れない）
#    ※ check_alldata_labels() はクリーニング前の Alldata を見るため、
#      disposition の水準は集約前のもの（自宅、介護施設など）が並ぶ。
# -----------------------------------------------------------------------------

rule("4. Input checks")

## ラベルの被覆（01_labels.R）。未定義の列・水準があればここで分かる。
.cov <- check_alldata_labels(d0)
check("Alldata の全列がラベル済みか使用しない列として登録されている",
      nrow(.cov$variables) ==
        sum(.cov$variables$status %in% c("labelled", "not used")),
      paste0("UNDEFINED: ",
             paste(.cov$variables$variable[.cov$variables$status == "UNDEFINED"],
                   collapse = ", ")))

check("id has no missing", !anyNA(d$id),
      sprintf("missing = %d", sum(is.na(d$id))))
check("age has no missing", !anyNA(d$age),
      sprintf("missing = %d", sum(is.na(d$age))))
check("day_in parsed as Date", !anyNA(d$day_in),
      sprintf("unparsed/missing = %d", sum(is.na(d$day_in))))

.lv <- function(x) paste(sort(unique(x[!is.na(x)])), collapse = " | ")

check("sex levels are {男, 女}",
      setequal(unique(d$sex[!is.na(d$sex)]), c("男", "女")),
      .lv(d$sex))
check("support_in levels are {あり, なし}",
      setequal(unique(d$support_in[!is.na(d$support_in)]), c("あり", "なし")),
      .lv(d$support_in))

## §6.1・§14-8：疾患区分に「その他」があるかを確認する（あれば 4 水準目として保持）
.class_lv <- sort(unique(d$class[!is.na(d$class)]))
check("class levels are {脳血管, 運動器, 廃用} only",
      setequal(.class_lv, c("脳血管", "運動器", "廃用")),
      paste0(.lv(d$class), "  <- 4 水準目があれば §6.1 に従い保持する"))

## §14-1：退院先の水準の一覧（集約後）と、英語ラベルの有無の確認
.dispo_lv  <- sort(unique(d$disposition[!is.na(d$disposition)]))
.dispo_lab <- show_dispo(.dispo_lv)   # 未定義の水準があればここで警告が出る
say("\n  disposition levels observed (n = ", length(.dispo_lv), "):")
for (i in seq_along(.dispo_lv)) {
  say(sprintf("    %-24s %6d   -> %s", .dispo_lv[i],
              sum(d$disposition == .dispo_lv[i], na.rm = TRUE), .dispo_lab[i]))
}
check("退院先の全水準に英語ラベルがある（§14-1）",
      all(.dispo_lab != .dispo_lv),
      paste0("unlabelled: ",
             paste(.dispo_lv[.dispo_lab == .dispo_lv], collapse = ", ")))
if (anyNA(d$disposition)) say(sprintf("    %-24s %6d", "(missing)", sum(is.na(d$disposition))))

check("E1/E2 の設定値がすべて元データに存在する",
      all(c(DISPO_TRANSFER, DISPO_DEATH) %in% .dispo_lv),
      paste0("not found: ",
             paste(setdiff(c(DISPO_TRANSFER, DISPO_DEATH), .dispo_lv), collapse = ", ")))

check("FIM in range (mFIM_in 13-91)",
      all(is.na(d$mFIM_in)  | (d$mFIM_in  >= MFIM_RANGE[1] & d$mFIM_in  <= MFIM_RANGE[2])),
      sprintf("out of range = %d",
              sum(!is.na(d$mFIM_in) & (d$mFIM_in < MFIM_RANGE[1] | d$mFIM_in > MFIM_RANGE[2]))))
check("FIM in range (cFIM_in 5-35)",
      all(is.na(d$cFIM_in)  | (d$cFIM_in  >= CFIM_RANGE[1] & d$cFIM_in  <= CFIM_RANGE[2])),
      sprintf("out of range = %d",
              sum(!is.na(d$cFIM_in) & (d$cFIM_in < CFIM_RANGE[1] | d$cFIM_in > CFIM_RANGE[2]))))
check("FIM in range (mFIM_out 13-91)",
      all(is.na(d$mFIM_out) | (d$mFIM_out >= MFIM_RANGE[1] & d$mFIM_out <= MFIM_RANGE[2])),
      sprintf("out of range = %d",
              sum(!is.na(d$mFIM_out) & (d$mFIM_out < MFIM_RANGE[1] | d$mFIM_out > MFIM_RANGE[2]))))


# -----------------------------------------------------------------------------
# 5. 対象の選抜（§5.1）
#    年齢 → 初回入院 → 退院確定 の順。この順序により ridit の基準集団と
#    感度分析 S の母集団が、ともに A4 に一本化される。
# -----------------------------------------------------------------------------

rule("5. Selection (A1 - A4)")

n_input <- nrow(d)

## --- A1：2017年4月〜2024年3月の連続入院 -------------------------------------
keep_A1 <- !is.na(d$day_in) & d$day_in >= PERIOD_START & d$day_in <= PERIOD_END
dA1 <- d[keep_A1, , drop = FALSE]
say(sprintf("A1  in-period admissions        : %5d   (dropped from input: %d)",
            nrow(dA1), n_input - nrow(dA1)))
check("入力がすでに対象期間に限定されている", nrow(dA1) == n_input,
      sprintf("out-of-period or missing day_in = %d", n_input - nrow(dA1)))

## --- A2：入院時の年齢が65歳以上 ---------------------------------------------
keep_A2 <- !is.na(dA1$age) & dA1$age >= AGE_MIN
dA2 <- dA1[keep_A2, , drop = FALSE]
n_age_na <- sum(is.na(dA1$age))
say(sprintf("A2  age >= %d at admission       : %5d   (removed: %d, of which age missing: %d)",
            AGE_MIN, nrow(dA2), nrow(dA1) - nrow(dA2), n_age_na))
check("年齢の欠測による除外がない（欠測があれば A2 の定義を確認する）",
      n_age_na == 0, sprintf("age missing = %d", n_age_na))

## --- A3：期間内の初回入院（A2 のうち、同一患者の2回目以降を除く）-----------
## 解析単位は患者（§5.2）。同一 id・同一 day_in が複数ある場合は、
## 元の行順で先頭を採る（規則を固定するための決めごと。件数を点検に出す）。
ord   <- order(dA2$id, dA2$day_in, dA2$.row)
dA2o  <- dA2[ord, , drop = FALSE]
first <- !duplicated(dA2o$id)
dA3   <- dA2o[first, , drop = FALSE]
n_repeat <- nrow(dA2) - nrow(dA3)
say(sprintf("A3  first admission per patient  : %5d   (removed repeat admissions: %d)",
            nrow(dA3), n_repeat))

.tie <- duplicated(dA2o[, c("id", "day_in")]) & !is.na(dA2o$day_in)
check("同一患者・同一入院日の重複レコードがない", sum(.tie) == 0,
      sprintf("ties = %d（元の行順で先頭を採用）", sum(.tie)))

## --- A4：抽出時点で退院が確定している ---------------------------------------
in_hosp <- (dA3$disposition %in% IN_HOSPITAL_LEVELS) |
           (IN_HOSPITAL_NA & is.na(dA3$disposition))
dA4 <- dA3[!in_hosp, , drop = FALSE]
n_inhosp <- sum(in_hosp)
say(sprintf("A4  discharge confirmed          : %5d   (still hospitalised: %d)",
            nrow(dA4), n_inhosp))

## §5.1：未退院例が A3 の 3% を超える場合は対象期間の終端を前倒しする
.p_inhosp <- if (nrow(dA3) > 0) n_inhosp / nrow(dA3) else 0
check("未退院例が A3 の 3% 以下（§5.1）", .p_inhosp <= 0.03,
      sprintf("%d / %d = %.2f%%", n_inhosp, nrow(dA3), 100 * .p_inhosp))

## day_out があれば、A4 の判定を実測で裏づける（§5.1、§14-4）。
## 退院日が入っていない例は在院中とみなせるため、両者の数が一致するはずである。
if ("day_out" %in% names(d0)) {
  .dayout_A3 <- as_date_safe(d0[["day_out"]])[dA3$.row]
  .n_nodate  <- sum(is.na(.dayout_A3))
  check("day_out の欠測数が未退院例の数と一致する（A4 の判定の裏づけ）",
        .n_nodate == n_inhosp,
        sprintf("day_out missing = %d, judged in-hospital = %d",
                .n_nodate, n_inhosp))
  .neg_los <- sum(!is.na(.dayout_A3) & .dayout_A3 < dA3$day_in)
  check("day_out が day_in より前になっている例がない",
        .neg_los == 0, sprintf("day_out < day_in : %d", .neg_los))
}

## --- 基準集団 ---------------------------------------------------------------
dat_base <- dA4
say(sprintf("\nBase population (A4)             : %5d", nrow(dat_base)))


# -----------------------------------------------------------------------------
# 6. 除外（§5.1 E1〜E3）。重なりは上位から一つだけ割り当てる（E2 > E1 > E3）
# -----------------------------------------------------------------------------

rule("6. Exclusions (E1 - E3), priority E2 > E1 > E3")

flag <- data.frame(
  E1 = dat_base$disposition %in% DISPO_TRANSFER,
  E2 = dat_base$disposition %in% DISPO_DEATH,
  E3 = is.na(dat_base$mFIM_out)
)

excl <- rep(NA_character_, nrow(dat_base))
for (k in EXCL_PRIORITY) excl[is.na(excl) & flag[[k]]] <- k
dat_base$excl_reason <- excl

n_flag <- vapply(flag, sum, integer(1))
n_excl <- c(E1 = sum(excl == "E1", na.rm = TRUE),
            E2 = sum(excl == "E2", na.rm = TRUE),
            E3 = sum(excl == "E3", na.rm = TRUE))

for (k in c("E1", "E2", "E3")) {
  say(sprintf("%s  %-46s flagged: %4d   assigned: %4d",
              k, STAGE_LABELS[[k]], n_flag[[k]], n_excl[[k]]))
}

## 重なりの一覧（§5.1：重なりを確認する）
ov_key <- apply(flag, 1, function(r) paste(c("E1", "E2", "E3")[r], collapse = "+"))
ov_key[ov_key == ""] <- "(none)"
overlap <- as.data.frame(table(pattern = ov_key), stringsAsFactors = FALSE)
names(overlap) <- c("flag_pattern", "n")
overlap <- overlap[order(-overlap$n), , drop = FALSE]
say("\n  Overlap of exclusion flags:")
print(overlap, row.names = FALSE)

## --- 解析対象（主解析）------------------------------------------------------
dat_main <- dat_base[is.na(dat_base$excl_reason), , drop = FALSE]
say(sprintf("\nAnalysis set (main analysis)     : %5d", nrow(dat_main)))


# -----------------------------------------------------------------------------
# 6-2. E3 の欠測理由の区分（§5.1、§14-3）と感度分析Sの集団（§8.8）
#      規則：「転院ならば全て急変」。区分は優先順位の割り当て後の E3
#      （＝転院でも死亡等でもない欠測例）に対して行う。
# -----------------------------------------------------------------------------

rule("6-2. E3 missing reason and sensitivity analysis S")

is_e3  <- !is.na(dat_base$excl_reason) & dat_base$excl_reason == "E3"
e3_cat <- rep(NA_character_, nrow(dat_base))
e3_cat[is_e3] <- classify_e3(flag$E1[is_e3])
dat_base$e3_category <- e3_cat

for (k in c("a", "b", "c")) {
  say(sprintf("  %s  %-54s n = %4d", k, E3_REASON_LABELS[[k]],
              sum(!is.na(e3_cat) & e3_cat == k)))
}

## 参考：退院時運動FIMが欠測した例のうち、転院・死亡等として先に除外された数。
## 規則により転院例は急変とみなすが、これらは E1・E2 として §8.8 で最低点を
## 与えられるため、区分(a)と扱いは一致する。
n_miss_e1 <- sum(flag$E3 & flag$E1)
n_miss_e2 <- sum(flag$E3 & flag$E2 & !flag$E1)
say(sprintf("\n  (reference) missing mFIM_out assigned to E1 (transfer) : %4d", n_miss_e1))
say(sprintf("  (reference) missing mFIM_out assigned to E2 (death etc.): %4d", n_miss_e2))
say("  いずれも §8.8 により最低点が与えられる（区分(a)と同じ扱い）")

## --- 感度分析Sの解析集団（§8.8）--------------------------------------------
## 基準集団（A4）から、E3 の区分(b)(c) に該当する患者だけを除く。
drop_s   <- is_e3 & !is.na(e3_cat) & e3_cat %in% E3_CATEGORY_EXCLUDE
dat_sens <- dat_base[!drop_s, , drop = FALSE]
n_sens   <- nrow(dat_sens)

## 複合アウトカム（§8.8）：E1・E2・E3(a) には最低点を与え、他は実測値を使う。
worst <- (dat_sens$disposition %in% c(DISPO_TRANSFER, DISPO_DEATH)) |
         (!is.na(dat_sens$e3_category) & dat_sens$e3_category == "a")
dat_sens$mFIM_out_S <- ifelse(worst, WORST_SCORE, dat_sens$mFIM_out)

say(sprintf("\nAnalysis set (sensitivity S)     : %5d   (E3(b)(c) removed: %d)",
            n_sens, sum(drop_s)))
say(sprintf("  of which assigned the worst score (%d) : %d", WORST_SCORE, sum(worst)))

check("感度分析Sの集団に欠測アウトカムが残っていない",
      !anyNA(dat_sens$mFIM_out_S),
      sprintf("missing = %d", sum(is.na(dat_sens$mFIM_out_S))))
check("区分(c) は置かない（規則どおり0件）",
      sum(!is.na(e3_cat) & e3_cat == "c") == 0)
check("区分(a) は優先順位により0件（転院例は E1 に割り当てられる）",
      sum(!is.na(e3_cat) & e3_cat == "a") == 0,
      sprintf("E3(a) = %d, 転院かつ欠測 = %d",
              sum(!is.na(e3_cat) & e3_cat == "a"), n_miss_e1))

## --- 補足表1 用：E3 の区分の年齢群別集計（§5.1、§8.8）----------------------
e3_tab <- do.call(rbind, lapply(c("a", "b", "c"), function(k) {
  sel <- !is.na(e3_cat) & e3_cat == k
  cnt <- count_by_g(dat_base$age[sel])
  data.frame(category = k,
             label = unname(E3_REASON_LABELS[k]),
             handling_in_S = unname(E3_HANDLING_LABELS[k]),
             G1 = cnt[["G1"]], G2 = cnt[["G2"]], G3 = cnt[["G3"]],
             total = sum(cnt), stringsAsFactors = FALSE)
}))
say("\n  E3 missing reason by age group:")
print(e3_tab, row.names = FALSE)


# -----------------------------------------------------------------------------
# 7. フロー図の数値（英語ラベル）
# -----------------------------------------------------------------------------

rule("7. Flow table")

flow <- data.frame(
  stage     = c("A1", "A2", "A3", "A4", "BASE", "E1", "E2", "E3", "MAIN", "SENS"),
  n         = c(nrow(dA1), nrow(dA2), nrow(dA3), nrow(dA4), nrow(dat_base),
                NA, NA, NA, nrow(dat_main), n_sens),
  n_removed = c(n_input - nrow(dA1), nrow(dA1) - nrow(dA2), n_repeat, n_inhosp, NA,
                n_excl[["E1"]], n_excl[["E2"]], n_excl[["E3"]], NA, NA),
  n_flagged = c(NA, NA, NA, NA, NA,
                n_flag[["E1"]], n_flag[["E2"]], n_flag[["E3"]], NA, NA),
  stringsAsFactors = FALSE
)
flow$label <- relabel_stage(flow$stage)
flow <- flow[, c("stage", "label", "n", "n_removed", "n_flagged")]
print(flow, row.names = FALSE)

## --- 年齢群別の内訳（Fig.1 の年齢群構成、§5.1）------------------------------
sets <- list(A3 = dA3, A4 = dA4, BASE = dat_base)
sets$E1 <- dat_base[!is.na(dat_base$excl_reason) & dat_base$excl_reason == "E1", ]
sets$E2 <- dat_base[!is.na(dat_base$excl_reason) & dat_base$excl_reason == "E2", ]
sets$E3 <- dat_base[!is.na(dat_base$excl_reason) & dat_base$excl_reason == "E3", ]
sets$MAIN   <- dat_main
sets$SENS   <- dat_sens
sets$INHOSP <- dA3[in_hosp, ]

flow_g <- do.call(rbind, lapply(names(sets), function(k) {
  cnt <- count_by_g(sets[[k]])
  data.frame(stage = k, label = relabel_stage(k),
             G1 = cnt[["G1"]], G2 = cnt[["G2"]], G3 = cnt[["G3"]],
             total = sum(cnt), stringsAsFactors = FALSE)
}))
attr(flow_g, "agegroup_labels") <- LEVELS_BY_VAR$age_group
say("\n  By age group (",
    paste(names(LEVELS_BY_VAR$age_group), unname(LEVELS_BY_VAR$age_group),
          sep = " = ", collapse = ", "), "):")
print(flow_g, row.names = FALSE)

## --- 基準集団の退院先の内訳 -------------------------------------------------
dispo_tab <- as.data.frame(table(disposition = dat_base$disposition, useNA = "ifany"),
                           stringsAsFactors = FALSE)
names(dispo_tab) <- c("disposition_ja", "n")
dispo_tab$disposition <- show_dispo(dispo_tab$disposition_ja, strict = FALSE)
dispo_tab$assigned_to <- ifelse(dispo_tab$disposition_ja %in% DISPO_DEATH, "E2",
                         ifelse(dispo_tab$disposition_ja %in% DISPO_TRANSFER, "E1", "-"))
dispo_tab <- dispo_tab[order(-dispo_tab$n), c("disposition", "assigned_to", "n")]
say("\n  Disposition in the base population:")
print(dispo_tab, row.names = FALSE)


# -----------------------------------------------------------------------------
# 8. 出来上がりの点検
# -----------------------------------------------------------------------------

rule("8. Verification")

check("A1 = A2 + age ineligible", nrow(dA1) == nrow(dA2) + sum(!keep_A2))
check("A2 = A3 + repeats",   nrow(dA2) == nrow(dA3) + n_repeat)
check("A3 = A4 + in-hospital", nrow(dA3) == nrow(dA4) + n_inhosp)
check("BASE = MAIN + E1 + E2 + E3",
      nrow(dat_base) == nrow(dat_main) + sum(n_excl),
      sprintf("%d = %d + %d", nrow(dat_base), nrow(dat_main), sum(n_excl)))
check("解析対象の id が一意（解析単位＝患者、§5.2）",
      !anyDuplicated(dat_main$id))
check("解析対象に退院時運動FIMの欠測がない", !anyNA(dat_main$mFIM_out))
check("解析対象に転院・死亡等が残っていない",
      !any(dat_main$disposition %in% c(DISPO_TRANSFER, DISPO_DEATH)))
check("解析対象の年齢がすべて 65 歳以上", all(dat_main$age >= AGE_MIN))
check("解析対象の入院日がすべて対象期間内",
      all(dat_main$day_in >= PERIOD_START & dat_main$day_in <= PERIOD_END))
check("BASE = SENS + E3(b)(c)（§8.8）",
      nrow(dat_base) == n_sens + sum(drop_s),
      sprintf("%d = %d + %d", nrow(dat_base), n_sens, sum(drop_s)))
check("感度分析Sの集団が解析対象を含む（S ⊇ MAIN）",
      all(dat_main$id %in% dat_sens$id))
check("感度分析Sの最低点付与が E1・E2・E3(a) と一致する（§8.8）",
      sum(worst) == n_excl[["E1"]] + n_excl[["E2"]] +
                    sum(!is.na(e3_cat) & e3_cat == "a"),
      sprintf("worst = %d, E1 + E2 + E3(a) = %d", sum(worst),
              n_excl[["E1"]] + n_excl[["E2"]] +
                sum(!is.na(e3_cat) & e3_cat == "a")))
## §7-3：最低点（12）が実測の値と重ならないこと。12→13 の置換が先に
## 行われていれば、最低点を与えていない患者に 12 点は残らない。
.n_worst_val <- sum(dat_sens$mFIM_out_S == WORST_SCORE)
check(sprintf("感度分析Sで %d 点をもつのは最低点を与えた患者だけ（§7-3）", WORST_SCORE),
      .n_worst_val == sum(worst) &&
        !any(dat_sens$mFIM_out_S[!worst] == WORST_SCORE),
      sprintf("mFIM_out_S = %d : %d, worst = %d, 実測で %d 点 = %d",
              WORST_SCORE, .n_worst_val, sum(worst), WORST_SCORE,
              sum(dat_sens$mFIM_out_S[!worst] == WORST_SCORE)))

checks <- do.call(rbind, .checks)
say("\n  REVIEW 項目: ", sum(checks$result == "REVIEW"), " / ", nrow(checks))


# -----------------------------------------------------------------------------
# 9. 保存
# -----------------------------------------------------------------------------

rule("9. Output")

w <- function(x, f) {
  p <- file.path(OUT_DIR, f)
  write.csv(x, p, row.names = FALSE, fileEncoding = "UTF-8")
  say("  written: ", p)
}
w(clean_tab, "table_cleaning.csv")
w(flow,      "table_flow_selection.csv")
w(flow_g,    "table_flow_by_agegroup.csv")
w(dispo_tab, "table_disposition_base.csv")
w(overlap,   "table_exclusion_overlap.csv")
w(e3_tab,    "table_e3_reason.csv")
w(checks,    "table_checks_02.csv")

dat_base$.row <- NULL
dat_main$.row <- NULL
dat_sens$.row <- NULL
save(dat_base, dat_main, dat_sens, flow, flow_g, dispo_tab, overlap, e3_tab,
     clean_tab, checks,
     file = file.path(OUT_DIR, "02_preprocess.rda"))
say("  written: ", file.path(OUT_DIR, "02_preprocess.rda"))

say("\nDone. dat_base (base population, A4) = ", nrow(dat_base),
    " ; dat_main (analysis set) = ", nrow(dat_main),
    " ; dat_sens (sensitivity S) = ", nrow(dat_sens))

## ログを閉じ、コンソールへの出力を戻す
.log_close()
