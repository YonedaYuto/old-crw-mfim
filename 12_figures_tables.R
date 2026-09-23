# =============================================================================
# 12_figures_tables.R
#   解析計画書 §10「図表」の作図・作表（plan.detail §11 のスクリプト16 に相当）。
#   Fig. 1（フロー図）は描かない（著者の指示、2026-09-21）。
#
#   入力 : output/03_m0_ranktest.rda      Fig. 2（E-1・E-2）、Table 1 の到達割合の照合
#          output/04_ridit_spline.rda     Table 1 の集計元（main）、参照水準
#          output/05_m2_clm.rda           補足表4（11 が無い場合の代替）、閾値の照合
#          output/06_standardize_boot.rda Table 2（E-3・E-4）
#          output/07_2_sens_clm_boot.rda  補足表2
#          output/09_positivity.rda       補足表3
#          output/10_age_continuous.rda   補足図1
#          output/11_supp_tables.rda      補足表1・補足表4
#          Alldata（ワークスペース。Table 1 の在棟日数にだけ使う。無ければ空欄にする）
#          01_labels.R（表示用の文字列）
#   出力 : output/figures/Fig2.tiff        Fig. 2（(a) 階段関数、(b) 2群相対効果）
#          output/figures/SuppFig1.tiff    補足図1（連続年齢のオッズ比曲線。旧版の補足図2）
#          output/tables/Table1.docx       Table 1（年齢群別の記述統計）
#          output/tables/Table2.docx       Table 2（70 点・77 点の標準化確率・RD・閾値別オッズ比、
#                                          共通オッズ比）
#          output/tables/SuppTable1.docx〜SuppTable4.docx
#          output/tables/AllTables.docx    上の表を1ファイルにまとめたもの（任意）
#          output/table1.csv               Table 1 の表示値（他の表の数値は 03〜11 の CSV にある）
#          output/table_outputs_12.csv     出力の一覧
#          output/table_checks_12.csv      点検結果
#          output/12_figures_tables.rda    図のオブジェクトと表の表示値
#          output/log_12_figures_tables.txt 実行ログ
#
#   方針
#     1. 数値は再計算しない。03〜11 の rda をそのまま使う。例外は Table 1 で、これを
#        作るスクリプトが無いため、04 の解析対象（main）からここで集計する（§8.1）。
#     2. すべてグレースケール（黒・白・灰）で描く。色は §0 の設定にまとめ、灰色だけで
#        あることを点検する。magick があれば TIFF 自体も 8 bit グレースケールに変換する。
#     3. 図は TIFF（600 dpi、LZW 圧縮）。文字の大きさは「数値の列と作図領域が図の幅に
#        収まる範囲で最大（上限 FONT_PT_MAX）」を図ごとに自動で決める。
#     4. 表は Word（.docx）。flextable の表として書き出すため、Word 上で編集できる。
#     5. タイトル・サブタイトル・脚注は SHOW（全体）と SHOW_BY_ITEM（図表ごと）で
#        表示と非表示を切り替える。
#     6. 出力される文字列はすべて英語（01_labels.R）。計画の決めごと（p値は主対比の
#        E-2 にだけ付ける、Table 2 に p値と RR を置かない等。§8.5）を表示でも守る。
#
#   改訂 : 2026-09-22 plan.summary.txt §7-1〜§7-3 により、次のとおり変更した。
#          ・§7-1：閾値を 65 点から 70 点・77 点の2つにした。Table 1 の「≥ 65」の行を
#            70 点・77 点の2行に、Fig. 2(a) の縦線を2本にした。閾値の値はここでは
#            決めず、03・06・07_2 の rda の諸元（"70;77" の形）から読む。
#            THRESHOLDS_EXPECTED と 05 の PPO_THRESHOLDS に突き合わせて点検する。
#          ・§7-2：Table 2 を、部分比例オッズモデル（M2-PPO）の閾値別オッズ比と
#            2つの閾値での標準化確率・RD に、比例オッズモデル（M2）の共通オッズ比を
#            併記する形にした。閾値を列に置く（著者の決定、2026-09-22）。
#          ・§7-2：閾値別オッズ比と共通オッズ比の図（旧版の補足図1。08 の出力）を
#            外した。08_po_thresholds.rda は読まない。連続年齢の図（旧版の補足図2）を
#            補足図1 に繰り上げ（著者の決定、2026-09-22）、識別子・設定名・ファイル名も
#            SuppFig1 に揃えた。
#          ・§7-2：補足表2 から閾値65の二値モデルの行を外し、E-4 を 70 点・77 点に
#            分けて示す（07_2 の table_supp2 の threshold 列を使う）。
#          ・§7-2：補足表4 で、年齢群の閾値別の係数（M2-PPO の nominal 母数。05 で
#            符号を反転済み）を閾値ごとに示す。
#          ・§7-3：感度分析Sの最低点（12 点）の説明は 07_2 の脚注にあり、ここでは変えない。
#          ・Fig. 2(a) の縦線の上端に点数（70、77）を書く（FIG2A_VLINE_LABELS）。横軸の
#            目盛（13 点刻み）に 70・77 が無いため。凡例は2つの閾値の到達割合を並べる。
#          ・凡例が長くなって右端が切れたため、文字幅の計測（text_width_mm）を 720 dpi
#            相当に改め（72 dpi では 4% ほど短く出る）、凡例の幅の見積もりで差し引く
#            余白を 30 mm から 34 mm に広げた。他の図の文字の大きさがわずかに変わりうる。
#          ・§7 以前の版の rda（閾値65・二値モデルを含むもの）を読んだ場合は、その図表を
#            作らずに止め、実行し直すスクリプトを点検結果に出す。output に残る旧版の
#            csv（table_binary65_standardized.csv、table_sens_binary65.csv など）と
#            08_po_thresholds.rda は読まない。
#
#   実行 : 01〜11 と同じフォルダを作業ディレクトリにし、Alldata を読み込んだ R で
#          source("12_figures_tables.R", encoding = "UTF-8")
#          （Alldata は Table 1 の在棟日数にだけ使う。無くても他の図表はすべて作られる）
#
#   必要なパッケージ : ggplot2（3.4 以上）, patchwork, flextable, officer
#   あれば使うもの   : ragg（TIFF の描画）, systemfonts（文字幅の計測）,
#                      magick（TIFF のグレースケール化）
#   文字コード : UTF-8（Windows の R では source(..., encoding = "UTF-8")）
# =============================================================================


# -----------------------------------------------------------------------------
# 0. 設定
# -----------------------------------------------------------------------------

OUT_DIR <- "output"
FIG_DIR <- file.path(OUT_DIR, "figures")
TAB_DIR <- file.path(OUT_DIR, "tables")

RDA_FILES <- c(
  "03" = "03_m0_ranktest.rda",
  "04" = "04_ridit_spline.rda",
  "05" = "05_m2_clm.rda",
  "06" = "06_standardize_boot.rda",
  "07_2" = "07_2_sens_clm_boot.rda",
  "09" = "09_positivity.rda",
  "10" = "10_age_continuous.rda",
  "11" = "11_supp_tables.rda"
)

## --- 作る図表（FALSE にすると作らない）--------------------------------------
##   SuppFig1 は連続年齢のオッズ比曲線（旧版の補足図2。§7-2 で旧版の補足図1 を外した）
MAKE <- c(Fig2 = TRUE, SuppFig1 = TRUE,
          Table1 = TRUE, Table2 = TRUE,
          SuppTable1 = TRUE, SuppTable2 = TRUE, SuppTable3 = TRUE, SuppTable4 = TRUE)
MAKE_ALL_TABLES_DOCX <- TRUE          # 表を1つの docx にまとめたものも出す

## --- タイトル・サブタイトル・脚注の表示（全体の既定）------------------------
##   投稿用の図（本文に図の説明を別に書く場合）は title/subtitle/footnote を
##   FALSE にする。表は Word 上の段落として出力されるので、残して編集してもよい。
SHOW <- list(title = FALSE, subtitle = FALSE, footnote = FALSE)

## --- 図表ごとの上書き（書かなければ SHOW に従う）----------------------------
##   例：SHOW_BY_ITEM$Fig2 <- list(title = FALSE, footnote = FALSE)
SHOW_BY_ITEM <- list(
  Fig2       = list(),
  SuppFig1   = list(),
  Table1     = list(),
  Table2     = list(),
  SuppTable1 = list(),
  SuppTable2 = list(),
  SuppTable3 = list(),
  SuppTable4 = list()
)

## --- タイトル（英語）---------------------------------------------------------
TITLES <- list(
  Fig2       = "Figure 2. Motor FIM at discharge by age group and two-sample relative effects",
  SuppFig1   = "Supplementary Figure 1. Odds ratio by age as a continuous variable",
  Table1     = "Table 1. Characteristics of the patients by age group",
  ## 閾値（70・77）が rda と食い違えば、点検（15 節）が知らせる
  Table2     = paste("Table 2. Standardised probabilities of a motor FIM at discharge",
                     "\u2265 70 and \u2265 77, risk differences and odds ratios by age group"),
  SuppTable1 = paste("Supplementary Table 1. Missing values, excluded patients and reasons",
                     "for a missing motor FIM at discharge, by age group"),
  SuppTable2 = "Supplementary Table 2. Main analysis and sensitivity analysis S",
  SuppTable3 = paste("Supplementary Table 3. Overlap of the covariate distributions",
                     "across age groups (positivity)"),
  SuppTable4 = paste("Supplementary Table 4. Coefficients of the partial proportional odds",
                     "model (M2-PPO)")
)

## --- 脚注を足す場合（図表ごと。文字列のベクトル）----------------------------
EXTRA_FOOTNOTES <- list()             # 例：EXTRA_FOOTNOTES$Table1 <- "..."

## --- 図の共通設定 -----------------------------------------------------------
FIG_DPI          <- 600
FIG_COMPRESSION  <- "lzw"
FIG_WIDTH_MM     <- 174               # 2段抜きの幅（1段なら 84 前後）
FIG_MAX_HEIGHT_MM <- 234              # 1ページに収まる高さの目安（図の本体が超えたら点検に出す）
FONT_FAMILY      <- "sans"            # Windows では Arial に対応する
FONT_PT_MAX      <- 16                # 文字の大きさの上限（pt、印刷寸法での値）
FONT_PT_MIN      <- 8                 # これより小さくはしない（収まらなければ警告）
MIN_PLOT_MM      <- 60                # 数値の列を並べた図で作図領域に残す最小の幅
COL_PAD_MM       <- 3                 # 文字の列どうしの間隔
TITLE_REL        <- 1.10              # タイトルの大きさ（本文の文字に対する比）
CAPTION_REL      <- 0.85              # 脚注の大きさ（同上）
LINEHEIGHT       <- 1.0
ROW_PAD          <- 0.35              # 数値の表の行間（行の高さに対する比）
LW               <- 0.8               # 主な線の太さ（mm）
TIFF_GRAYSCALE   <- TRUE              # magick があれば 8 bit グレースケールに変換する

## --- グレースケールの描き分け（年齢群）--------------------------------------
##   線種・灰色の濃さ・記号の3つで区別する。値はすべて灰色（R = G = B）であること。
GROUP_STYLE <- list(
  colour   = c(G1 = "black", G2 = "grey55", G3 = "black"),
  linetype = c(G1 = "solid", G2 = "solid",  G3 = "42"),
  shape    = c(G1 = 21,      G2 = 24,       G3 = 22),
  fill     = c(G1 = "black", G2 = "grey55", G3 = "white")
)
REF_COLOUR   <- "grey40"              # 基準線（0.5、1）
GUIDE_COLOUR <- "grey50"              # 補助線（閾値の 70・77 点、年齢群の境界）
BAND_FILL    <- "grey80"              # 区間の帯
POINT_MAIN   <- list(shape = 23, fill = "black")   # 主対比・共通オッズ比
POINT_SUB    <- list(shape = 21, fill = "white")   # 副次対比・閾値別オッズ比

## --- 閾値（plan.summary §7-1）-----------------------------------------------
##   図表に使う閾値は 03・06・07_2 の rda の諸元から読む（ここでは決めない）。
##   次の値はその点検にだけ使う（旧版の rda を読んでいないかを確かめる）。
##   70 点：入浴と階段昇降以外が見守りレベル（旧設定 65 点を修正）。
##   77 点：入浴と階段昇降以外が修正自立レベル。
THRESHOLDS_EXPECTED <- c(70, 77)

## --- Fig. 2 -----------------------------------------------------------------
FIG2A_PLOT_MM    <- 80                # (a) の作図領域の高さ
FIG2A_VLINE_LABELS <- TRUE            # (a) の閾値の縦線の上端に点数（70、77）を書く
FIG2A_LEGEND     <- "bottom"          # "bottom"、または作図領域内の位置 c(x, y)（0〜1）
FIG2B_XLIM       <- NULL              # (b) の横軸。NULL なら区間と 0.5 が入る範囲。全域なら c(0, 1)
FIG2_AXIS_X      <- "Motor FIM at discharge, y (points)"
FIG2_AXIS_Y      <- "P(motor FIM at discharge \u2265 y)"
FIG2B_AXIS_X     <- "Relative effect"

## --- 補足図1（連続年齢。旧版の補足図2）-------------------------------------
SUPPFIG1_PLOT_MM     <- 100
SUPPFIG1_GROUP_LINES <- TRUE          # 年齢群の境界（75・90歳）に補助線と群名を置く
SUPPFIG1_LABEL_AGES  <- NULL          # 数値を添える年齢（例：c(75, 90)）。NULL なら添えない
SUPPFIG1_AXIS_X      <- "Age at admission (years)"

## --- 表の共通設定 -----------------------------------------------------------
TABLE_FONT      <- "Arial"
TABLE_PT        <- 10
TABLE_TITLE_PT  <- 11
TABLE_NOTE_PT   <- 9
TABLE_INDENT_PT <- 12                 # 字下げ1段の幅（pt）
TABLE_MAX_WIDTH_IN <- 6.3             # 表の幅の上限（インチ。A4 縦・余白 25 mm で約 6.3）
TABLE1_TOTAL    <- TRUE               # Table 1 に全体の列を置く
DIGITS <- list(percent = 1, re = 2, or = 2, risk = 1, rd = 1,
               coef = 3, ridit = 3, ratio = 2)
NA_TEXT       <- "\u2014"             # 値が無い（算出できない）セル
UNICODE_MINUS <- TRUE                 # 負の数に U+2212 の負号を使う

## --- 在棟日数（Table 1、§6.1・§8.1）----------------------------------------
##   著者の決定（2026-09-21）：在棟日数 = 退院日 − 入院日（日）。
##   02_preprocess.R は退院日を dat_main に残していないため、ワークスペースの
##   Alldata から id と入院日で照合して取る（02 と同じ規則：同一 id・同一入院日が
##   複数あれば元の行順で先頭）。Alldata が無ければ Table 1 の該当行は空欄にする。
LOS_OFFSET <- 0L

## --- 補足表3・補足表4 --------------------------------------------------------
SUPP3_SHOW_EXTRAPOLATION <- FALSE     # 重なりの範囲の外にある例数（§8.7-3）も載せる
SUPP4_INCLUDE_THRESHOLDS <- FALSE     # 閾値（切片）の母数も載せる
SUPP4_OR_FOR_SPLINE      <- FALSE     # スプライン基底の係数にも exp(β) を載せる

## --- 対比の並び（§8.2）------------------------------------------------------
CONTRAST_KEYS <- c("G3_vs_G1", "G2_vs_G1", "G3_vs_G2")
PRIMARY_KEY   <- "G3_vs_G1"


# -----------------------------------------------------------------------------
# 1. 表示用ラベル
# -----------------------------------------------------------------------------

if (!file.exists("01_labels.R")) {
  stop("01_labels.R が見つからない。12_figures_tables.R は 01_labels.R と同じ",
       "フォルダを作業ディレクトリにして実行すること。")
}
source("01_labels.R", encoding = "UTF-8")


# -----------------------------------------------------------------------------
# 2. 実行の準備（出力先とログ）
# -----------------------------------------------------------------------------

for (.d in c(OUT_DIR, FIG_DIR, TAB_DIR)) {
  if (!dir.exists(.d)) dir.create(.d, recursive = TRUE)
}

.log_con <- file(file.path(OUT_DIR, "log_12_figures_tables.txt"), open = "wt",
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

.outputs <- list()
record_output <- function(id, file, status, detail = "") {
  .outputs[[length(.outputs) + 1L]] <<-
    data.frame(item = id, file = file, status = status, detail = detail,
               stringsAsFactors = FALSE)
}

rule(paste0("12_figures_tables.R  ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
say("R : ", R.version.string)


# -----------------------------------------------------------------------------
# 3. パッケージ
# -----------------------------------------------------------------------------

rule("3. Packages")

.need_pkg <- c("ggplot2", "patchwork", "flextable", "officer")
.miss_pkg <- .need_pkg[!vapply(.need_pkg, requireNamespace, logical(1), quietly = TRUE)]
if (length(.miss_pkg)) {
  .log_close()
  stop("必要なパッケージがない: ", paste(.miss_pkg, collapse = ", "),
       "\n  install.packages(c(\"", paste(.miss_pkg, collapse = "\", \""),
       "\")) を実行すること。")
}
suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})
HAS_RAGG        <- requireNamespace("ragg", quietly = TRUE)
HAS_SYSTEMFONTS <- requireNamespace("systemfonts", quietly = TRUE)
HAS_MAGICK      <- requireNamespace("magick", quietly = TRUE)
GG_35           <- utils::packageVersion("ggplot2") >= "3.5.0"

for (.p in c(.need_pkg, "ragg", "systemfonts", "magick")) {
  say(sprintf("  %-12s %s", .p,
              if (requireNamespace(.p, quietly = TRUE))
                as.character(utils::packageVersion(.p)) else "(not installed)"))
}
check("ggplot2 が 3.4 以上（linewidth を使う）",
      utils::packageVersion("ggplot2") >= "3.4.0",
      as.character(utils::packageVersion("ggplot2")))
check("TIFF を ragg で描く（無ければ grDevices::tiff）", HAS_RAGG,
      if (HAS_RAGG) "ragg::agg_tiff" else "grDevices::tiff に切り替えた")
check("文字幅を systemfonts で計る（無ければ文字数から近似）", HAS_SYSTEMFONTS, "")
if (TIFF_GRAYSCALE) {
  check("TIFF を 8 bit グレースケールに変換できる（magick）", HAS_MAGICK,
        if (HAS_MAGICK) "" else "RGB のまま保存する（描画は灰色のみ）")
}


# -----------------------------------------------------------------------------
# 4. 入力の取り込み
#    rda は環境ごとに読み、オブジェクト名の衝突（meta、footnotes など）を避ける。
# -----------------------------------------------------------------------------

rule("4. Input")

E <- list()
for (.k in names(RDA_FILES)) {
  .path <- file.path(OUT_DIR, RDA_FILES[[.k]])
  if (file.exists(.path)) {
    .e <- new.env(parent = emptyenv())
    load(.path, envir = .e)
    E[[.k]] <- .e
    say("  loaded    : ", .path)
  } else {
    say("  not found : ", .path)
  }
}

#' rda の中のオブジェクトを取る（無ければ NULL）
getx <- function(key, nm) {
  e <- E[[key]]
  if (!is.null(e) && exists(nm, envir = e, inherits = FALSE)) get(nm, envir = e) else NULL
}
#' 必要なオブジェクトが揃っていなければ止める（図表ごとの tryCatch で受ける）
need <- function(key, nm) {
  x <- getx(key, nm)
  if (is.null(x)) stop(file.path(OUT_DIR, RDA_FILES[[key]]), " に ", nm,
                       " が無い。先に該当のスクリプトを実行すること。")
  x
}
or_else <- function(x, y) if (is.null(x) || !length(x)) y else x

.set04     <- getx("04", "M2_SETTINGS")
AGE_LEVELS <- or_else(.set04$age_levels, c("G1", "G2", "G3"))
AGE_BREAKS <- or_else(.set04$age_breaks, c(65, 75, 90))
say("  age groups : ", paste(AGE_LEVELS, relabel_levels("age_group", AGE_LEVELS),
                             sep = " = ", collapse = ", "))


# -----------------------------------------------------------------------------
# 5. 共通の道具
# -----------------------------------------------------------------------------

## --- 5-1. 数値の書式 ---------------------------------------------------------

.minus <- function(s) {
  if (isTRUE(UNICODE_MINUS)) sub("^-", "\u2212", s) else s
}

fmt_num <- function(x, digits) {
  x   <- suppressWarnings(as.numeric(x))
  out <- rep(NA_character_, length(x))
  ok  <- is.finite(x)
  r   <- round(x[ok], digits)
  r[r == 0] <- 0                                  # "-0.0" を避ける
  out[ok] <- .minus(formatC(r, format = "f", digits = digits))
  out
}

#' 点推定値（区間）。どちらかの端が負なら区切りを " to " にする
fmt_ci <- function(est, lo, hi, digits, scale = 1) {
  e <- fmt_num(scale * est, digits)
  l <- fmt_num(scale * lo,  digits)
  h <- fmt_num(scale * hi,  digits)
  neg <- (is.finite(lo) & round(scale * lo, digits) < 0) |
         (is.finite(hi) & round(scale * hi, digits) < 0)
  sep <- ifelse(neg, " to ", "\u2013")
  ifelse(is.na(e), NA_TEXT,
         ifelse(is.na(l) | is.na(h), e, sprintf("%s (%s%s%s)", e, l, sep, h)))
}

fmt_n_pct <- function(n, d, digits = DIGITS$percent) {
  n <- as.numeric(n); d <- as.numeric(d)
  ifelse(is.finite(d) & d > 0,
         sprintf("%d (%s)", as.integer(n), fmt_num(100 * n / d, digits)),
         sprintf("%d", as.integer(n)))
}

#' 中央値［第1四分位, 第3四分位］（§8.1：quantile(type = 1)）
fmt_med_iqr <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[!is.na(x)]
  if (!length(x)) return(NA_TEXT)
  q  <- stats::quantile(x, c(0.50, 0.25, 0.75), type = 1, names = FALSE)
  dg <- if (all(abs(q - round(q)) < 1e-9)) 0L else 1L
  sprintf("%s [%s, %s]", fmt_num(q[1], dg), fmt_num(q[2], dg), fmt_num(q[3], dg))
}

#' 3つの数（中央値・Q1・Q3）から同じ書式を作る（11 の出力を表示するとき）
fmt_med_iqr3 <- function(med, q1, q3) {
  ok <- is.finite(med) & is.finite(q1) & is.finite(q3)
  dg <- ifelse(ok & abs(med - round(med)) < 1e-9 & abs(q1 - round(q1)) < 1e-9 &
                 abs(q3 - round(q3)) < 1e-9, 0L, 1L)
  out <- rep(NA_TEXT, length(med))
  for (i in which(ok)) {
    out[i] <- sprintf("%s [%s, %s]", fmt_num(med[i], dg[i]), fmt_num(q1[i], dg[i]),
                      fmt_num(q3[i], dg[i]))
  }
  out
}

#' 順列p値の表示。順列で観測値以上に極端な並べ替えが1つも無ければ 1/nperm 未満とだけ書く
fmt_p_perm <- function(p, nperm) {
  p <- suppressWarnings(as.numeric(p))
  if (!length(p) || !is.finite(p)) return("")
  if (p <= 0) {
    return(sprintf("P < %s", format(1 / nperm, scientific = FALSE, drop0trailing = TRUE)))
  }
  if (p < 0.001) return("P < 0.001")
  sprintf("P = %s", formatC(p, format = "f", digits = 3))
}

meta_get <- function(meta, item, default = NA_character_) {
  if (is.null(meta) || !all(c("item", "value") %in% names(meta))) return(default)
  v <- meta$value[meta$item == item]
  if (length(v)) as.character(v[1]) else default
}

#' 脚注の表（03・06 は tag/text、07_2・08〜11 は n/footnote）を文字列に揃える
footnote_text <- function(fn, roman = FALSE) {
  if (is.null(fn) || !nrow(fn)) return(character(0))
  txt <- if ("text" %in% names(fn)) fn$text else if ("footnote" %in% names(fn)) fn$footnote
         else fn[[ncol(fn)]]
  txt <- gsub(">=", "\u2265", gsub("<=", "\u2264", as.character(txt), fixed = TRUE),
              fixed = TRUE)
  tag <- if ("tag" %in% names(fn)) fn$tag else NULL
  if (roman && !is.null(tag)) txt <- sprintf("(%s) %s", tag, txt)
  as.character(txt)
}

## --- 5-2. 名前 ---------------------------------------------------------------

age_label <- function(g) relabel_levels("age_group", g, strict = FALSE)

#' 対比の短い表示名。"≥90 vs 65-74 years" の形（高齢側を先に書く。§8.4）
contrast_short <- function(ref, cmp) {
  sprintf("%s vs %s", sub(" years$", "", age_label(cmp)), age_label(ref))
}
contrast_parts <- function(key) {
  p <- strsplit(key, "_vs_", fixed = TRUE)[[1]]
  list(cmp = p[1], ref = p[2])
}
contrast_name <- function(key, primary_mark = TRUE, sep = " ") {
  p <- contrast_parts(key)
  paste0(contrast_short(p$ref, p$cmp),
         if (primary_mark && identical(key, PRIMARY_KEY)) paste0(sep, "(primary)") else "")
}

## --- 5-3. 表示の切り替え ------------------------------------------------------

show_part <- function(id, part) {
  v <- SHOW_BY_ITEM[[id]][[part]]
  if (is.null(v)) isTRUE(SHOW[[part]]) else isTRUE(v)
}

## --- 5-4. 文字の寸法 ----------------------------------------------------------

pt_mm   <- function(pt) pt * 25.4 / 72
line_mm <- function(pt) pt_mm(pt) * 1.2 * LINEHEIGHT
n_lines <- function(x) {
  if (is.null(x) || !length(x)) return(0L)
  x <- as.character(x); x[is.na(x)] <- ""
  max(lengths(regmatches(x, gregexpr("\n", x, fixed = TRUE))) + 1L)
}

#' 文字列の幅（mm）。複数行は最も長い行。systemfonts が無ければ文字数から近似する
#' §7 改訂（2026-09-22）：10 pt・72 dpi（10 px）で測ると字送りが整数の画素に丸められ、
#' 600 dpi で描いた幅より 4% ほど短く出る（Fig. 2(a) の凡例が長くなって右端が切れた）。
#' 720 dpi で測って pt に直す。
text_width_mm <- function(x, pt, bold = FALSE) {
  x <- as.character(x); x <- x[!is.na(x)]
  ln <- unlist(strsplit(x, "\n", fixed = TRUE))
  ln <- ln[nzchar(ln)]
  if (!length(ln)) return(0)
  w <- NULL
  if (HAS_SYSTEMFONTS) {
    w <- tryCatch(systemfonts::string_width(ln, family = FONT_FAMILY, bold = bold,
                                            size = 10, res = 720) / 10 * pt / 10,
                  error = function(e) NULL)
  }
  if (is.null(w) || !all(is.finite(w))) {
    w <- nchar(ln, type = "width") * pt * (if (bold) 0.60 else 0.55)
  }
  max(w) * 25.4 / 72
}

#' 文字の列（通常と太字が混在）の幅
col_width_mm <- function(text, bold, pt) {
  bold <- rep_len(as.logical(bold), length(text))
  max(text_width_mm(text[!bold], pt, FALSE), text_width_mm(text[bold], pt, TRUE))
}

## --- 5-5. テーマと保存 -------------------------------------------------------

theme_fig <- function(pt) {
  ggplot2::theme_classic(base_size = pt, base_family = FONT_FAMILY) +
    ggplot2::theme(
      text              = ggplot2::element_text(colour = "black"),
      axis.text         = ggplot2::element_text(size = pt, colour = "black"),
      axis.title        = ggplot2::element_text(size = pt, colour = "black"),
      axis.line         = ggplot2::element_line(colour = "black", linewidth = LW * 0.6),
      axis.ticks        = ggplot2::element_line(colour = "black", linewidth = LW * 0.6),
      axis.ticks.length = grid::unit(pt_mm(pt) * 0.35, "mm"),
      legend.text       = ggplot2::element_text(size = pt, colour = "black"),
      legend.title      = ggplot2::element_text(size = pt, colour = "black"),
      legend.background = ggplot2::element_blank(),
      legend.key        = ggplot2::element_blank(),
      panel.background  = ggplot2::element_blank(),
      plot.background   = ggplot2::element_blank(),
      plot.margin       = ggplot2::margin(2, 3, 1, 1, "mm")
    )
}

theme_tag <- function(pt) {
  ggplot2::theme(plot.tag = ggplot2::element_text(size = pt * TITLE_REL, face = "bold",
                                                  family = FONT_FAMILY, colour = "black"),
                 plot.tag.position = "topleft")
}

#' 図の幅に合わせて折り返す（段落ごと）。実際の文字幅を測って単語単位で折る
wrap_text <- function(s, pt, width_mm, bold = FALSE) {
  if (is.null(s) || !length(s)) return(NULL)
  s <- as.character(s[!is.na(s) & nzchar(s)])
  if (!length(s)) return(NULL)
  avail <- width_mm - 10
  one <- function(z) {
    words <- strsplit(z, " ", fixed = TRUE)[[1]]
    lines <- character(0); cur <- ""
    for (w in words) {
      cand <- if (nzchar(cur)) paste(cur, w) else w
      if (nzchar(cur) && text_width_mm(cand, pt, bold) > avail) {
        lines <- c(lines, cur); cur <- w
      } else cur <- cand
    }
    paste(c(lines, cur), collapse = "\n")
  }
  paste(vapply(s, one, character(1)), collapse = "\n")
}

#' タイトル・サブタイトル・脚注を付ける（表示の切り替えに従う）。
#' 付けた分だけ必要になる高さ（mm）も返す。
annotate_fig <- function(p, id, title, subtitle, notes, pt, width_mm) {
  cap_pt <- pt * CAPTION_REL
  ttl <- if (show_part(id, "title"))    wrap_text(title,    pt * TITLE_REL, width_mm, TRUE) else NULL
  sub <- if (show_part(id, "subtitle")) wrap_text(subtitle, pt,             width_mm) else NULL
  cap <- if (show_part(id, "footnote")) wrap_text(notes,    cap_pt,         width_mm) else NULL
  extra <- n_lines(ttl) * line_mm(pt * TITLE_REL) + n_lines(sub) * line_mm(pt) +
           n_lines(cap) * line_mm(cap_pt) * 1.05 +
           3 * sum(!vapply(list(ttl, sub, cap), is.null, logical(1)))
  q <- p + patchwork::plot_annotation(
    title = ttl, subtitle = sub, caption = cap,
    theme = ggplot2::theme(
      plot.title    = ggplot2::element_text(size = pt * TITLE_REL, face = "bold",
                                            family = FONT_FAMILY, hjust = 0,
                                            margin = ggplot2::margin(0, 0, 2, 0, "mm")),
      plot.subtitle = ggplot2::element_text(size = pt, family = FONT_FAMILY, hjust = 0,
                                            margin = ggplot2::margin(0, 0, 2, 0, "mm")),
      plot.caption  = ggplot2::element_text(size = cap_pt, family = FONT_FAMILY, hjust = 0,
                                            lineheight = 1.05,
                                            margin = ggplot2::margin(3, 0, 0, 0, "mm")),
      plot.title.position   = "plot",
      plot.caption.position = "plot",
      plot.background = ggplot2::element_rect(fill = "white", colour = NA),
      plot.margin     = ggplot2::margin(3, 3, 3, 3, "mm")))
  list(plot = q, extra_mm = extra + 6,
       annotated = !all(vapply(list(ttl, sub, cap), is.null, logical(1))))
}

#' TIFF（600 dpi、LZW）で保存し、可能なら 8 bit グレースケールに変換する
save_tiff <- function(p, file, width_mm, height_mm, annotated = FALSE) {
  if (HAS_RAGG) {
    ragg::agg_tiff(filename = file, width = width_mm, height = height_mm, units = "mm",
                   res = FIG_DPI, compression = FIG_COMPRESSION, background = "white")
  } else {
    .args <- list(filename = file, width = width_mm, height = height_mm, units = "mm",
                  res = FIG_DPI, compression = FIG_COMPRESSION, bg = "white")
    if (isTRUE(capabilities("cairo"))) .args$type <- "cairo"
    do.call(grDevices::tiff, .args)
  }
  .dev <- grDevices::dev.cur()
  tryCatch(print(p), finally = grDevices::dev.off(.dev))

  gray <- FALSE
  if (isTRUE(TIFF_GRAYSCALE) && HAS_MAGICK) {
    gray <- tryCatch({
      img <- magick::image_read(file)
      img <- magick::image_convert(img, type = "Grayscale", colorspace = "gray")
      tmp <- paste0(file, ".tmp.tiff")
      magick::image_write(img, path = tmp, format = "tiff",
                          density = sprintf("%dx%d", FIG_DPI, FIG_DPI),
                          compression = "LZW")
      file.copy(tmp, file, overwrite = TRUE)
      unlink(tmp)
      TRUE
    }, error = function(e) {
      say("    グレースケールへの変換に失敗（RGB のまま）: ", conditionMessage(e))
      FALSE
    })
  }
  info <- sprintf("%.0f x %.0f mm, %d dpi, %s", width_mm, height_mm, FIG_DPI,
                  if (gray) "8-bit grayscale" else "RGB (grey colours only)")
  if (HAS_MAGICK) {
    ii <- tryCatch(magick::image_info(magick::image_read(file)), error = function(e) NULL)
    if (!is.null(ii)) {
      info <- paste0(info, sprintf(" ; %d x %d px, colorspace %s", ii$width[1],
                                   ii$height[1], ii$colorspace[1]))
    }
  }
  say("    written: ", file, "  (", info, ")")
  ## 図の本体の高さを点検する。タイトル・脚注を表示している場合は、その分で
  ## 1ページを超えることがあるため、点検にはせず注記だけを出す。
  if (annotated) {
    if (height_mm > FIG_MAX_HEIGHT_MM) {
      say(sprintf(paste0("    （注）タイトル・脚注を含めた高さ %.0f mm は目安の %g mm を超える。",
                        "投稿用には SHOW / SHOW_BY_ITEM で非表示にする"),
                  height_mm, FIG_MAX_HEIGHT_MM))
    }
  } else {
    check(sprintf("%s：高さが %g mm 以下", basename(file), FIG_MAX_HEIGHT_MM),
          height_mm <= FIG_MAX_HEIGHT_MM,
          sprintf("%.0f mm。FIG_WIDTH_MM・行間（ROW_PAD）・作図領域の高さを見直す", height_mm))
  }
  info
}

## --- 5-6. 軸の目盛 -----------------------------------------------------------

#' 目盛のラベルが重ならない範囲で最も細かい目盛を選ぶ
choose_breaks <- function(lim, log_x, panel_mm, pt) {
  lab_fun <- function(b) {
    out <- format(b, trim = TRUE, drop0trailing = TRUE, scientific = FALSE)
    .minus(out)
  }
  if (log_x) {
    dec  <- 10^(-3:3)
    sets <- list(sort(as.vector(outer(c(1, 1.5, 2, 3, 5, 7), dec))),
                 sort(as.vector(outer(c(1, 2, 5), dec))),
                 sort(as.vector(outer(c(1, 3), dec))),
                 dec)
    tr <- log10
  } else {
    span <- diff(lim)
    steps <- c(0.01, 0.02, 0.025, 0.05, 0.1, 0.2, 0.25, 0.5, 1, 2, 2.5, 5, 10, 20, 25, 50)
    steps <- steps[steps >= span / 12]
    sets <- lapply(steps, function(s) seq(floor(lim[1] / s) * s, ceiling(lim[2] / s) * s, by = s))
    tr <- identity
  }
  best <- NULL
  for (b in sets) {
    b <- b[b >= lim[1] - 1e-9 & b <= lim[2] + 1e-9]
    if (length(b) < 2) next
    gap_mm <- panel_mm * min(diff(tr(b))) / diff(tr(lim))
    if (gap_mm >= 1.35 * text_width_mm(lab_fun(b), pt)) { best <- b; break }
    best <- b
  }
  if (is.null(best)) best <- lim
  list(breaks = best, labels = lab_fun(best))
}


# -----------------------------------------------------------------------------
# 5-7. 数値の列を並べた図（フォレストプロットと文字の表）
#      行ごとに「左の見出し」「作図（任意）」「右の数値の列」を横に並べる。
#      行の高さは行内の最大の行数で決め、全パネルで同じ縦座標を使う。
#      span = TRUE の行（見出し行）は左の見出しが右隣へはみ出してもよい。
# -----------------------------------------------------------------------------

#' 文字の大きさを決める（左の見出しと数値の列が収まり、作図領域が MIN_PLOT_MM 以上）
#' span_text   : 見出し行の文字。span_all = TRUE なら右の列すべてにはみ出してよい
#'               （その行に値が無い場合）。FALSE なら左の列と作図領域の幅に収める。
fit_font_pt <- function(cols, span_text = character(0), has_plot = TRUE,
                        width_mm = FIG_WIDTH_MM, span_all = FALSE) {
  w10 <- sum(vapply(cols, function(cc) col_width_mm(cc$text, cc$bold, 10), numeric(1)))
  fixed <- length(cols) * COL_PAD_MM + 6
  lim1 <- if (w10 > 0) (width_mm - fixed - if (has_plot) MIN_PLOT_MM else 0) / w10 * 10 else Inf
  lim2 <- Inf
  if (length(span_text)) {
    ws  <- text_width_mm(span_text, 10, TRUE)
    wr  <- if (length(cols) > 1 && !span_all) sum(vapply(cols[-1], function(cc)
      col_width_mm(cc$text, cc$bold, 10), numeric(1))) else 0
    lim2 <- (width_mm - fixed) / (ws + wr) * 10
  }
  pt <- floor(min(FONT_PT_MAX, lim1, lim2) * 2) / 2
  if (pt < FONT_PT_MIN) {
    check("文字の大きさが下限以上で収まる", FALSE,
          sprintf("必要な大きさ %.1f pt < 下限 %g pt。FIG_WIDTH_MM を広げるか MIN_PLOT_MM を小さくする",
                  pt, FONT_PT_MIN))
    pt <- FONT_PT_MIN
  }
  pt
}

#' rows : data.frame。列 cols[1] が左の見出し、cols[-1] が右の数値の列。
#'        .face（左の見出しの書体）、.span（見出し行）を持つ。
#'        forest を描く場合は est, lo, hi, shape, fill, psize を持つ。
#' forest: NULL（文字だけ）または list(log, ref, xlim, xlab, vsegs)
build_row_figure <- function(rows, cols, headers, pt, forest = NULL,
                             width_mm = FIG_WIDTH_MM) {
  n  <- nrow(rows)
  rows$.face <- or_else(rows$.face, rep("plain", n))
  rows$.span <- or_else(rows$.span, rep(FALSE, n))
  for (cc in cols) { v <- as.character(rows[[cc]]); v[is.na(v)] <- ""; rows[[cc]] <- v }

  ## 行の高さ（行数の単位）と縦座標
  h  <- vapply(seq_len(n), function(i) max(vapply(cols, function(cc)
          n_lines(rows[[cc]][i]), integer(1))), integer(1)) + ROW_PAD
  hh <- max(vapply(headers, n_lines, integer(1)), 1L) + ROW_PAD
  body_total <- sum(h)
  total <- body_total + hh
  rows$y <- body_total - (cumsum(h) - h / 2)
  hy     <- body_total + hh / 2
  ylim   <- c(0, total)

  ## 列の幅（見出し行ははみ出してよいので幅の計算から外す）
  bold1 <- rows$.face == "bold"
  wl <- col_width_mm(c(rows[[cols[1]]][!rows$.span], headers[1]),
                     c(bold1[!rows$.span], TRUE), pt) + COL_PAD_MM
  wr <- vapply(seq_along(cols)[-1], function(k)
          col_width_mm(c(rows[[cols[k]]], headers[k]),
                       c(rep(FALSE, n), TRUE), pt) + COL_PAD_MM, numeric(1))
  wp <- if (is.null(forest)) 0 else width_mm - wl - sum(wr) - 6

  text_panel <- function(label, face, header) {
    d <- data.frame(y = rows$y, label = label, face = face, stringsAsFactors = FALSE)
    d <- d[nzchar(d$label), , drop = FALSE]
    p <- ggplot2::ggplot() +
      ggplot2::geom_blank(data = data.frame(x = c(0, 1), y = ylim),
                          ggplot2::aes(x = x, y = y))
    if (nrow(d)) {
      p <- p + ggplot2::geom_text(data = d,
                                  ggplot2::aes(x = 0, y = y, label = label, fontface = face),
                                  hjust = 0, vjust = 0.5, size = pt / ggplot2::.pt,
                                  family = FONT_FAMILY, lineheight = LINEHEIGHT,
                                  colour = "black")
    }
    if (nzchar(header)) {
      p <- p + ggplot2::annotate("text", x = 0, y = hy, label = header, hjust = 0,
                                 vjust = 0.5, fontface = "bold", size = pt / ggplot2::.pt,
                                 family = FONT_FAMILY, lineheight = LINEHEIGHT,
                                 colour = "black")
    }
    p + ggplot2::coord_cartesian(xlim = c(0, 1), ylim = ylim, expand = FALSE, clip = "off") +
      ggplot2::theme_void(base_family = FONT_FAMILY) +
      ggplot2::theme(plot.margin = ggplot2::margin(0, 0, 0, 1, "mm"),
                     plot.background = ggplot2::element_blank(),
                     panel.background = ggplot2::element_blank())
  }

  plist  <- list(text_panel(rows[[cols[1]]], rows$.face, headers[1]))
  widths <- wl

  axis_mm <- 0
  if (!is.null(forest)) {
    it  <- rows[is.finite(rows$est), , drop = FALSE]
    rng <- range(c(it$lo, it$hi, it$est, forest$ref), na.rm = TRUE, finite = TRUE)
    lim <- forest$xlim
    if (is.null(lim)) {
      if (isTRUE(forest$log)) {
        f <- (rng[2] / rng[1])^0.08
        lim <- c(rng[1] / f, rng[2] * f)
      } else {
        lim <- rng + c(-1, 1) * 0.08 * diff(rng)
      }
    }
    br <- choose_breaks(lim, isTRUE(forest$log), wp, pt)
    ref_span <- or_else(forest$ref_span, c(0, body_total))

    p <- ggplot2::ggplot() +
      ggplot2::geom_blank(data = data.frame(x = lim, y = ylim), ggplot2::aes(x = x, y = y)) +
      ggplot2::annotate("segment", x = forest$ref, xend = forest$ref,
                        y = ref_span[1], yend = ref_span[2],
                        linetype = "dashed", colour = REF_COLOUR, linewidth = LW * 0.6)
    if (!is.null(forest$vsegs) && nrow(forest$vsegs)) {
      p <- p + ggplot2::geom_segment(data = forest$vsegs,
                                     ggplot2::aes(x = x, xend = x, y = y0, yend = y1),
                                     linetype = "dotted", colour = GUIDE_COLOUR,
                                     linewidth = LW * 0.6)
    }
    if (nrow(it)) {
      p <- p +
        ggplot2::geom_segment(data = it, ggplot2::aes(x = lo, xend = hi, y = y, yend = y),
                              linewidth = LW, colour = "black", lineend = "butt") +
        ggplot2::geom_point(data = it,
                            ggplot2::aes(x = est, y = y, shape = shape, fill = fill,
                                         size = psize),
                            colour = "black", stroke = LW * 0.9)
    }
    p <- p + ggplot2::scale_shape_identity() + ggplot2::scale_fill_identity() +
      ggplot2::scale_size_identity() +
      (if (isTRUE(forest$log))
         ggplot2::scale_x_log10(breaks = br$breaks, labels = br$labels)
       else ggplot2::scale_x_continuous(breaks = br$breaks, labels = br$labels)) +
      ggplot2::coord_cartesian(xlim = lim, ylim = ylim, expand = FALSE, clip = "on") +
      ggplot2::labs(x = forest$xlab, y = NULL) +
      theme_fig(pt) +
      ggplot2::theme(axis.line.y = ggplot2::element_blank(),
                     axis.text.y = ggplot2::element_blank(),
                     axis.ticks.y = ggplot2::element_blank(),
                     axis.title.y = ggplot2::element_blank(),
                     plot.margin = ggplot2::margin(0, 3, 0, 3, "mm"))
    plist  <- c(plist, list(p))
    widths <- c(widths, wp)
    axis_mm <- line_mm(pt) * (1 + n_lines(forest$xlab)) + pt_mm(pt) * 0.35 + 3
  }
  for (k in seq_along(cols)[-1]) {
    plist  <- c(plist, list(text_panel(rows[[cols[k]]], rep("plain", n), headers[k])))
    widths <- c(widths, wr[k - 1])
  }
  ## 文字の列は mm で固定し、作図領域（文字だけの表では最後の列）に残りの幅を回す
  flex <- if (is.null(forest)) length(widths) else 2L
  wunit <- rep("mm", length(widths)); wunit[flex] <- "null"
  wval <- widths; wval[flex] <- 1
  fig <- patchwork::wrap_plots(plist, nrow = 1, widths = grid::unit(wval, wunit))
  list(plot = fig, height_mm = total * line_mm(pt) + axis_mm + 2, pt = pt,
       widths = widths, plot_mm = wp)
}


# -----------------------------------------------------------------------------
# 5-8. Word の表
#      rows は data.frame（1列目が見出し、2列目以降が値）に、書式の列
#      .indent（字下げの段数）、.bold、.italic、.merge（横に結合する見出し行）を持つ。
# -----------------------------------------------------------------------------

new_rows <- function(values) {
  data.frame(values, .indent = integer(0), .bold = logical(0), .italic = logical(0),
             .merge = logical(0), check.names = FALSE, stringsAsFactors = FALSE)
}

#' 表の1行を作る
trow <- function(label, values, cols, indent = 0L, bold = FALSE, italic = FALSE,
                 merge = FALSE) {
  v <- as.character(values)
  if (length(v) <= 1L) v <- rep_len(if (length(v)) v else "", length(cols))
  v[is.na(v)] <- NA_TEXT
  d <- as.data.frame(as.list(c(label, v)), stringsAsFactors = FALSE)
  names(d) <- c("label", cols)
  d$.indent <- as.integer(indent); d$.bold <- bold; d$.italic <- italic; d$.merge <- merge
  d
}

make_ft <- function(tab, header) {
  meta <- tab[, c(".indent", ".bold", ".italic", ".merge")]
  df   <- tab[, setdiff(names(tab), c(".indent", ".bold", ".italic", ".merge")), drop = FALSE]
  df[] <- lapply(df, function(z) { z <- as.character(z); z[is.na(z)] <- ""; z })
  keys <- paste0("c", seq_along(df))
  names(df) <- keys
  nc <- ncol(df)

  ft <- flextable::flextable(df)
  ft <- flextable::set_header_labels(ft, values = stats::setNames(as.list(header), keys))
  ft <- flextable::font(ft, fontname = TABLE_FONT, part = "all")
  ft <- flextable::fontsize(ft, size = TABLE_PT, part = "all")
  ft <- flextable::color(ft, color = "black", part = "all")
  ft <- flextable::bold(ft, part = "header")
  ft <- flextable::align(ft, j = 1, align = "left", part = "all")
  if (nc > 1) ft <- flextable::align(ft, j = 2:nc, align = "center", part = "all")
  ft <- flextable::valign(ft, valign = "bottom", part = "header")

  for (i in which(meta$.bold))   ft <- flextable::bold(ft, i = i, j = 1, part = "body")
  for (i in which(meta$.italic)) ft <- flextable::italic(ft, i = i, j = 1, part = "body")
  for (i in which(meta$.indent > 0)) {
    ft <- flextable::padding(ft, i = i, j = 1, part = "body",
                             padding.left = 3 + TABLE_INDENT_PT * meta$.indent[i])
  }
  if (nc > 1) {
    for (i in which(meta$.merge)) {
      ft <- flextable::merge_at(ft, i = i, j = seq_len(nc), part = "body")
      ft <- flextable::align(ft, i = i, j = 1, align = "left", part = "body")
    }
  }
  ft <- flextable::padding(ft, padding.top = 1, padding.bottom = 1, part = "all")

  ## 罫線：表の上下と見出しの下だけ（黒）
  b <- officer::fp_border(color = "black", width = 1)
  ft <- flextable::border_remove(ft)
  ft <- flextable::hline_top(ft, border = b, part = "header")
  ft <- flextable::hline_bottom(ft, border = b, part = "header")
  ft <- flextable::hline_bottom(ft, border = b, part = "body")
  ## 列幅は内容に合わせて決め（autofit）、固定レイアウトで書き出す。
  ## 本文の幅を超える場合は文字を小さくせず、1列目（見出し）の幅を詰めて折り返す。
  ## それでも収まらなければ全列の幅を同じ比で詰める。Word 上で列幅は自由に変えられる。
  ft <- flextable::autofit(ft, add_w = 0.05, add_h = 0)
  wd <- dim(ft)$widths                 # dim.flextable（列幅、インチ）
  if (sum(wd) > TABLE_MAX_WIDTH_IN) {
    rest <- sum(wd[-1])
    if (nc > 1 && rest <= TABLE_MAX_WIDTH_IN - 1.5) {
      ft <- flextable::width(ft, j = 1, width = TABLE_MAX_WIDTH_IN - rest)
    } else {
      ft <- flextable::width(ft, j = seq_len(nc), width = wd * TABLE_MAX_WIDTH_IN / sum(wd))
    }
  }
  ft <- flextable::set_table_properties(ft, layout = "fixed")
  ft
}

fp_title <- function() officer::fp_text(font.size = TABLE_TITLE_PT, bold = TRUE,
                                        font.family = TABLE_FONT, color = "black")
fp_sub   <- function() officer::fp_text(font.size = TABLE_PT, italic = TRUE,
                                        font.family = TABLE_FONT, color = "black")
fp_head  <- function() officer::fp_text(font.size = TABLE_PT, bold = TRUE,
                                        font.family = TABLE_FONT, color = "black")
fp_note  <- function() officer::fp_text(font.size = TABLE_NOTE_PT,
                                        font.family = TABLE_FONT, color = "black")

add_par <- function(doc, text, prop) {
  officer::body_add_fpar(doc, officer::fpar(officer::ftext(text, prop = prop)))
}

#' 表の仕様（spec）を docx に書き足す
#' spec: list(id, title, subtitle, blocks = list(list(heading, ft)), notes)
docx_add_spec <- function(doc, spec) {
  if (show_part(spec$id, "title") && length(spec$title) && nzchar(spec$title)) {
    doc <- add_par(doc, spec$title, fp_title())
  }
  if (show_part(spec$id, "subtitle")) {
    for (s in spec$subtitle) doc <- add_par(doc, s, fp_sub())
  }
  for (k in seq_along(spec$blocks)) {
    b <- spec$blocks[[k]]
    if (!is.null(b$heading) && nzchar(b$heading)) doc <- add_par(doc, b$heading, fp_head())
    doc <- flextable::body_add_flextable(doc, b$ft, align = "left")
    if (k < length(spec$blocks)) doc <- officer::body_add_par(doc, "")
  }
  notes <- c(spec$notes, EXTRA_FOOTNOTES[[spec$id]])
  if (show_part(spec$id, "footnote") && length(notes)) {
    doc <- officer::body_add_par(doc, "")
    for (s in notes) doc <- add_par(doc, s, fp_note())
  }
  doc
}

TABLE_SPECS <- list()

write_table_docx <- function(spec, file) {
  doc <- officer::read_docx()
  doc <- docx_add_spec(doc, spec)
  print(doc, target = file)
  say("    written: ", file)
  TABLE_SPECS[[spec$id]] <<- spec
  invisible(file)
}

#' 図表1つ分を実行する。失敗しても他の図表は続ける。
run_part <- function(id, fun) {
  if (!isTRUE(MAKE[[id]])) {
    say("  skipped (MAKE$", id, " = FALSE)")
    record_output(id, "", "skipped")
    return(invisible(NULL))
  }
  res <- tryCatch(fun(), error = function(e) e)
  if (inherits(res, "error")) {
    msg <- conditionMessage(res)
    say("  エラー : ", msg)
    check(sprintf("[%s] を作成できた", id), FALSE, msg)
    record_output(id, "", "failed", msg)
  } else {
    check(sprintf("[%s] を作成できた", id), TRUE, basename(res))
    record_output(id, res, "written")
  }
  invisible(res)
}

FIG_DATA <- list()


# -----------------------------------------------------------------------------
# 5-9. 閾値（plan.summary §7-1）
#      図表に使う閾値は各 rda の諸元から読む。複数の閾値は ";" で区切った1つの値
#      （例 "70;77"）で保存されている（03・06・07_2 の改訂版）。旧版の rda（"65"）を
#      読んでいないことを THRESHOLDS_EXPECTED と 05 の PPO_THRESHOLDS で確かめる。
# -----------------------------------------------------------------------------

#' 諸元の閾値（"70;77"）を数値のベクトル（昇順）にする。無ければ NULL
parse_thresholds <- function(meta, item = "threshold") {
  v <- meta_get(meta, item, NA_character_)
  if (is.na(v) || !nzchar(v)) return(NULL)
  out <- suppressWarnings(as.numeric(strsplit(v, ";", fixed = TRUE)[[1]]))
  if (!length(out) || anyNA(out)) return(NULL)
  sort(out)
}
#' 閾値の表示（例 "70 and 77"）
thr_text <- function(t, sep = " and ") paste(format(t, trim = TRUE), collapse = sep)

rule("5-9. Thresholds (plan.summary 7-1)")

.ppo_thr05 <- getx("05", "PPO_THRESHOLDS")
THR_SOURCES <- list(
  "03 fig2_meta"          = parse_thresholds(getx("03", "fig2_meta")),
  "05 PPO_THRESHOLDS"     = if (is.null(.ppo_thr05)) NULL else sort(as.numeric(.ppo_thr05)),
  "06 table2_meta"        = parse_thresholds(getx("06", "table2_meta")),
  "07_2 table_supp2_meta" = parse_thresholds(getx("07_2", "table_supp2_meta"))
)
for (.k in names(THR_SOURCES)) {
  say(sprintf("  %-22s : %s", .k,
              if (is.null(THR_SOURCES[[.k]])) "(not available)"
              else paste(THR_SOURCES[[.k]], collapse = ", ")))
}
.thr_found <- !vapply(THR_SOURCES, is.null, logical(1))
.thr_ok    <- vapply(THR_SOURCES, function(x)
  is.null(x) || isTRUE(all.equal(as.numeric(x), as.numeric(sort(THRESHOLDS_EXPECTED)))),
  logical(1))
check(sprintf("各 rda の閾値が §7-1 の %s 点と一致する（旧版の rda を読んでいない）",
              thr_text(sort(THRESHOLDS_EXPECTED), "・")),
      all(.thr_ok) && all(.thr_found),
      paste0(if (any(!.thr_ok)) paste0("一致しない: ", paste(names(THR_SOURCES)[!.thr_ok],
                                                             collapse = ", "),
                                        "。script_section7 の版で実行し直すこと。") else "",
             if (any(!.thr_found)) paste0("読めない（rda が無いか旧版）: ",
                                         paste(names(THR_SOURCES)[!.thr_found],
                                               collapse = ", ")) else ""))


# -----------------------------------------------------------------------------
# 6. Fig. 2（E-1・E-2、§10.1）
#    (a) 年齢群別の P(退院時運動FIM ≥ y)、y = 13〜91 の階段関数。70点と77点に縦線
#        （§7-1。閾値は 03 の fig2_meta から読む）。
#    (b) 3対比の2群相対効果と未調整95%順列区間（基準線0.5）。p値は主対比のみ。
# -----------------------------------------------------------------------------

rule("6. Figure 2")

make_fig2 <- function() {
  id <- "Fig2"
  fig2a <- need("03", "fig2a")
  re    <- need("03", "m0_relative_effect")
  meta  <- getx("03", "fig2_meta")
  fn    <- getx("03", "fig2_footnotes")

  THR   <- parse_thresholds(meta)                  # §7-1：70, 77（昇順）
  if (is.null(THR)) stop("03 の fig2_meta に閾値（threshold）が無い。§7 版の 03 を実行すること")
  NPERM <- as.numeric(meta_get(meta, "nperm_rankFD", "100000"))
  N_ALL <- as.integer(meta_get(meta, "analysis_set_n", as.character(sum(unique(fig2a[, c("age_group", "n")])$n))))

  ## --- (b) の行 --------------------------------------------------------------
  re <- re[match(intersect(CONTRAST_KEYS, re$contrast), re$contrast), , drop = FALSE]
  is_pri <- re$contrast == PRIMARY_KEY
  pv <- re$p_value[is_pri]
  pstr <- if (length(pv) && is.finite(pv[1])) fmt_p_perm(pv[1], NPERM) else ""
  rows_b <- data.frame(
    left  = vapply(re$contrast, contrast_name, character(1), sep = "\n"),
    value = paste0(fmt_ci(re$relative_effect, re$ci_lower, re$ci_upper, DIGITS$re),
                   ifelse(is_pri & nzchar(pstr), paste0("\n", pstr), "")),
    est = re$relative_effect, lo = re$ci_lower, hi = re$ci_upper,
    shape = ifelse(is_pri, POINT_MAIN$shape, POINT_SUB$shape),
    fill  = ifelse(is_pri, POINT_MAIN$fill,  POINT_SUB$fill),
    .face = ifelse(is_pri, "bold", "plain"),
    stringsAsFactors = FALSE
  )
  headers_b <- c("Contrast", "Relative effect\n(95% CI)")

  .has_p <- grepl("P [<=]", rows_b$value)
  check("Fig. 2(b)：p値を付けたのは主対比だけ（§8.5-2）",
        all(.has_p == (is_pri & nzchar(pstr))) && all(is.na(re$p_value[!is_pri])),
        sprintf("primary : %s", if (nzchar(pstr)) pstr else "(no p-value)"))
  check("Fig. 2(b)：3対比とも同じ構成の区間（未調整95%）",
        length(unique(re$ci_method)) == 1L && length(unique(re$conf_level)) == 1L,
        paste(unique(re$ci_method), collapse = " / "))

  ## --- (a) の凡例（各閾値での到達割合を添える）------------------------------
  d <- fig2a
  d$age_group <- factor(as.character(d$age_group), levels = AGE_LEVELS)
  thr_pts <- d[d$y %in% THR, , drop = FALSE]
  check("Fig. 2(a)：各閾値の行が年齢群ごとに1本ずつある",
        nrow(thr_pts) == length(THR) * length(AGE_LEVELS) &&
          !anyDuplicated(paste(thr_pts$y, thr_pts$age_group)),
        sprintf("rows = %d (thresholds %s)", nrow(thr_pts), thr_text(THR, ", ")))
  .t1   <- thr_pts[thr_pts$y == THR[1], , drop = FALSE]
  n_g   <- .t1$n[match(AGE_LEVELS, .t1$age_group)]
  pct_txt <- vapply(AGE_LEVELS, function(g) {
    paste(vapply(THR, function(t) {
      z <- thr_pts[thr_pts$y == t & as.character(thr_pts$age_group) == g, , drop = FALSE]
      sprintf("%s%% \u2265 %g", if (nrow(z)) fmt_num(100 * z$prob[1], DIGITS$percent) else NA_TEXT, t)
    }, character(1)), collapse = ", ")
  }, character(1))
  leg_lab <- stats::setNames(
    sprintf("%s (n = %d): %s", age_label(AGE_LEVELS), n_g, pct_txt), AGE_LEVELS)

  ## --- 文字の大きさ（(b) の数値の列と (a) の凡例が収まる最大）---------------
  pt_b <- fit_font_pt(list(list(text = c(rows_b$left, headers_b[1]),
                                bold = c(rows_b$.face == "bold", TRUE)),
                           list(text = c(rows_b$value, headers_b[2]),
                                bold = c(rep(FALSE, nrow(rows_b)), TRUE))))
  ## 凡例は作図領域の左端（縦軸の目盛と見出しの右、約 27 mm）から始まり、右に余白が
  ## 約 6 mm ある。§7 改訂で凡例が2つの閾値の値をもち長くなったため、差し引く幅を
  ## 30 mm から 34 mm に広げた。
  w_leg10 <- text_width_mm(leg_lab, 10) + 18 * 10 / FONT_PT_MAX
  pt_a <- floor(min(FONT_PT_MAX, (FIG_WIDTH_MM - 34) / w_leg10 * 10) * 2) / 2
  pt <- max(FONT_PT_MIN, min(pt_a, pt_b))
  say(sprintf("  font size : %.1f pt (limit from panel (b) %.1f pt, legend of (a) %.1f pt)",
              pt, pt_b, pt_a))

  ## --- (a) ------------------------------------------------------------------
  ps <- pt / ggplot2::.pt * 0.75
  ## 閾値の縦線。点数を上端に書く場合は、文字と重ならないよう線を文字の下で止める
  ## （作図領域の高さ FIG2A_PLOT_MM から文字の高さを縦軸の単位に換算する）。
  vl_top <- if (isTRUE(FIG2A_VLINE_LABELS)) 1 - 1.15 * line_mm(pt) / FIG2A_PLOT_MM else Inf
  p_a <- ggplot2::ggplot(d, ggplot2::aes(x = y, y = prob, colour = age_group,
                                         linetype = age_group)) +
    ggplot2::annotate("segment", x = THR, xend = THR, y = -Inf, yend = vl_top,
                      linetype = "dotted", colour = GUIDE_COLOUR, linewidth = LW * 0.7)
  if (isTRUE(FIG2A_VLINE_LABELS)) {
    p_a <- p_a + ggplot2::annotate("text", x = THR, y = 1, label = format(THR, trim = TRUE),
                                   vjust = 1, hjust = 0.5, size = pt * CAPTION_REL / ggplot2::.pt,
                                   family = FONT_FAMILY, colour = "black")
  }
  p_a <- p_a +
    ggplot2::geom_step(direction = "vh", linewidth = LW) +
    ggplot2::geom_point(data = thr_pts,
                        ggplot2::aes(x = y, y = prob, shape = age_group, fill = age_group),
                        size = ps, stroke = LW * 0.9) +
    ggplot2::scale_colour_manual(values = GROUP_STYLE$colour,   labels = leg_lab, name = NULL) +
    ggplot2::scale_linetype_manual(values = GROUP_STYLE$linetype, labels = leg_lab, name = NULL) +
    ggplot2::scale_shape_manual(values = GROUP_STYLE$shape,     labels = leg_lab, name = NULL) +
    ggplot2::scale_fill_manual(values = GROUP_STYLE$fill,       labels = leg_lab, name = NULL) +
    ggplot2::scale_x_continuous(breaks = seq(13, 91, by = 13), limits = c(13, 91),
                                expand = ggplot2::expansion(mult = c(0.01, 0.02))) +
    ggplot2::scale_y_continuous(breaks = seq(0, 1, by = 0.25),
                                labels = c("0", "0.25", "0.50", "0.75", "1"),
                                limits = c(0, 1),
                                expand = ggplot2::expansion(mult = c(0.01, 0.03))) +
    ggplot2::labs(x = FIG2_AXIS_X, y = FIG2_AXIS_Y) +
    theme_fig(pt) +
    ggplot2::theme(plot.margin = ggplot2::margin(line_mm(pt * TITLE_REL) + 1, 3, 1, 1, "mm"),
                   legend.key.width = grid::unit(pt_mm(pt) * 3.2, "mm"),
                   legend.key.height = grid::unit(line_mm(pt), "mm"),
                   legend.spacing.y = grid::unit(0, "mm"),
                   legend.margin = ggplot2::margin(0, 0, 0, 0))
  leg_mm <- 0
  if (is.numeric(FIG2A_LEGEND)) {
    p_a <- p_a + if (GG_35) ggplot2::theme(legend.position = "inside",
                                           legend.position.inside = FIG2A_LEGEND,
                                           legend.justification = FIG2A_LEGEND)
                 else ggplot2::theme(legend.position = FIG2A_LEGEND,
                                     legend.justification = FIG2A_LEGEND)
  } else {
    p_a <- p_a + ggplot2::theme(legend.position = "bottom", legend.direction = "vertical",
                                legend.justification = "left")
    leg_mm <- length(AGE_LEVELS) * line_mm(pt) * 1.1 + 3
  }
  p_a <- p_a + ggplot2::guides(colour = ggplot2::guide_legend(ncol = 1),
                               linetype = ggplot2::guide_legend(ncol = 1),
                               shape = ggplot2::guide_legend(ncol = 1),
                               fill = ggplot2::guide_legend(ncol = 1))

  ## --- (b) ------------------------------------------------------------------
  rows_b$psize <- ifelse(rows_b$.face == "bold", ps * 1.25, ps)
  fb <- build_row_figure(rows_b, cols = c("left", "value"), headers = headers_b, pt = pt,
                         forest = list(log = FALSE, ref = 0.5, xlim = FIG2B_XLIM,
                                       xlab = FIG2B_AXIS_X))

  ## --- 組み立て --------------------------------------------------------------
  tag_mm <- line_mm(pt * TITLE_REL) + 2
  h_a <- FIG2A_PLOT_MM + line_mm(pt) * (1 + n_lines(FIG2_AXIS_X)) + 5 + leg_mm + tag_mm
  h_b <- fb$height_mm + tag_mm
  wa <- patchwork::wrap_elements(full = p_a) + ggplot2::labs(tag = "(a)") + theme_tag(pt)
  wb <- patchwork::wrap_elements(full = fb$plot) + ggplot2::labs(tag = "(b)") + theme_tag(pt)
  fig <- (wa / wb) + patchwork::plot_layout(heights = c(h_a, h_b))

  n_txt <- paste(sprintf("%s, n\u00a0=\u00a0%d", age_label(AGE_LEVELS), n_g), collapse = "; ")
  subtitle <- sprintf("Analysis set, N = %d (%s)", N_ALL, n_txt)
  notes <- c(
    sprintf(paste("(a) Proportion of patients whose motor FIM at discharge is y points or more,",
                  "by age group (y = 13 to 91); the dotted vertical lines mark %s points and",
                  "the legend gives the proportions reaching %s points or more."),
            thr_text(THR), thr_text(THR)),
    paste("(b) Two-sample relative effect p = P(older > younger) + 0.5 \u00d7 P(tie);",
          "values below 0.5 indicate that the older group tends to have a lower motor FIM",
          "at discharge. The dashed line marks 0.5."),
    sprintf(paste("P is from the two-sided permuted Brunner\u2013Munzel test (%s permutations)",
                  "and is shown for the primary contrast only."),
            format(NPERM, big.mark = ",", scientific = FALSE)),
    footnote_text(fn, roman = TRUE),
    "FIM, Functional Independence Measure; CI, confidence interval."
  )
  ann <- annotate_fig(fig, id, TITLES[[id]], subtitle, notes, pt, FIG_WIDTH_MM)

  file <- file.path(FIG_DIR, "Fig2.tiff")
  height <- h_a + h_b + ann$extra_mm
  save_tiff(ann$plot, file, FIG_WIDTH_MM, height, ann$annotated)
  FIG_DATA$Fig2 <<- list(plot = ann$plot, panel_a = p_a, rows_b = rows_b, pt = pt,
                         thresholds = THR, width_mm = FIG_WIDTH_MM, height_mm = height)
  file
}
run_part("Fig2", make_fig2)


# -----------------------------------------------------------------------------
# 7. 補足図1（§8.9、§10.2。旧版の補足図2。§7-2 で旧版の補足図1 を外して繰り上げた）
#    70歳を参照とした年齢の共通オッズ比の曲線。線は細かい格子の点推定値、
#    帯は profile likelihood 95%区間（10 の格子）。
#    10 の出力（table_supp_fig2_*.csv など）は旧番号のファイル名のまま。数値は同じ。
# -----------------------------------------------------------------------------

rule("7. Supplementary Figure 1 (age as a continuous variable)")

make_suppfig1 <- function() {
  id <- "SuppFig1"
  aor  <- need("10", "age_or")
  acw  <- getx("10", "age_curve_wald")
  meta <- getx("10", "meta")
  fn   <- getx("10", "footnotes")

  ref_age <- as.numeric(aor$reference_age[1])
  band <- aor[is.finite(aor$ci_lower) & is.finite(aor$ci_upper), , drop = FALSE]
  band <- band[order(band$age), , drop = FALSE]
  line <- if (!is.null(acw)) acw[acw$age >= min(band$age) & acw$age <= max(band$age), ,
                                 drop = FALSE] else aor
  line <- line[order(line$age), , drop = FALSE]
  check("補足図1：profile 区間が得られていない年齢がない",
        all(aor$ci_ok | aor$age == ref_age),
        sprintf("failed at: %s", paste(aor$age[!aor$ci_ok & aor$age != ref_age],
                                       collapse = ", ")))
  check("補足図1：帯は profile likelihood 区間である（Wald は描かない）",
        all(grepl("profile", unique(aor$ci_method))), paste(unique(aor$ci_method)))

  pt <- FONT_PT_MAX
  yr <- range(c(band$ci_lower, band$ci_upper, line$odds_ratio, 1), finite = TRUE)
  ylim <- c(yr[1] / 1.15, yr[2] * (if (SUPPFIG1_GROUP_LINES) 1.9 else 1.15))
  xr <- range(band$age)
  br <- choose_breaks(ylim, TRUE, SUPPFIG1_PLOT_MM, pt)

  p <- ggplot2::ggplot() +
    ggplot2::geom_ribbon(data = band, ggplot2::aes(x = age, ymin = ci_lower, ymax = ci_upper),
                         fill = BAND_FILL, colour = NA) +
    ggplot2::geom_hline(yintercept = 1, linetype = "dashed", colour = REF_COLOUR,
                        linewidth = LW * 0.6)
  if (isTRUE(SUPPFIG1_GROUP_LINES)) {
    cut <- AGE_BREAKS[-1] - 0.5
    lo  <- c(xr[1], cut); hi <- c(cut, xr[2])
    lab <- sub(" years$", "", age_label(AGE_LEVELS))
    p <- p + ggplot2::geom_vline(xintercept = cut, linetype = "dotted", colour = GUIDE_COLOUR,
                                 linewidth = LW * 0.7) +
      ggplot2::annotate("text", x = (lo + hi) / 2, y = ylim[2] / 1.04, label = lab,
                        vjust = 1, size = pt / ggplot2::.pt, family = FONT_FAMILY,
                        colour = "black")
  }
  p <- p +
    ggplot2::geom_line(data = line, ggplot2::aes(x = age, y = odds_ratio),
                       linewidth = LW, colour = "black") +
    ggplot2::annotate("point", x = ref_age, y = 1, shape = 21, fill = "black",
                      colour = "black", size = pt / ggplot2::.pt * 0.75)
  if (length(SUPPFIG1_LABEL_AGES)) {
    la <- aor[aor$age %in% SUPPFIG1_LABEL_AGES, , drop = FALSE]
    if (nrow(la)) {
      la$lab <- fmt_ci(la$odds_ratio, la$ci_lower, la$ci_upper, DIGITS$or)
      p <- p + ggplot2::geom_point(data = la, ggplot2::aes(x = age, y = odds_ratio),
                                   shape = 21, fill = "white", size = pt / ggplot2::.pt * 0.7) +
        ggplot2::geom_text(data = la, ggplot2::aes(x = age, y = ci_lower, label = lab),
                           vjust = 1.3, size = pt * 0.85 / ggplot2::.pt, family = FONT_FAMILY)
    }
  }
  p <- p +
    ggplot2::scale_y_log10(breaks = br$breaks, labels = br$labels) +
    ggplot2::scale_x_continuous(breaks = seq(ceiling(xr[1] / 5) * 5, floor(xr[2] / 5) * 5, 5)) +
    ggplot2::coord_cartesian(xlim = xr, ylim = ylim,
                             expand = TRUE) +
    ggplot2::labs(x = SUPPFIG1_AXIS_X,
                  y = sprintf("Odds ratio (reference: %g years)", ref_age)) +
    theme_fig(pt)

  h_plot <- SUPPFIG1_PLOT_MM + line_mm(pt) * 2 + 6
  n_m2 <- meta_get(meta, "n", NA_character_)
  subtitle <- sprintf("M2 complete-case population, N = %s", n_m2)
  notes <- c(
    sprintf(paste("Solid line: odds ratio relative to %g years; shaded band: 95%% profile",
                  "likelihood interval (computed at %s-year steps). The dashed line marks 1%s."),
            ref_age, meta_get(meta, "profile_grid_step", "1"),
            if (SUPPFIG1_GROUP_LINES) "; dotted vertical lines mark the age-group boundaries" else ""),
    footnote_text(fn),
    "FIM, Functional Independence Measure."
  )
  ann <- annotate_fig(p, id, TITLES[[id]], subtitle, notes, pt, FIG_WIDTH_MM)
  file <- file.path(FIG_DIR, "SuppFig1.tiff")
  height <- h_plot + ann$extra_mm
  save_tiff(ann$plot, file, FIG_WIDTH_MM, height, ann$annotated)
  check("補足図1：p値を置いていない（§8.5）", all(is.na(aor$p_value)), "")
  FIG_DATA$SuppFig1 <<- list(plot = ann$plot, band = band, line = line, pt = pt,
                             width_mm = FIG_WIDTH_MM, height_mm = height)
  file
}
run_part("SuppFig1", make_suppfig1)


# -----------------------------------------------------------------------------
# 8. Table 1（E-1、§8.1）
#    解析対象（04 の main）から集計する。群間比較の数値（p値・標準化差）は置かない。
#    中央値・四分位は quantile(type = 1)。カテゴリは n（%）、分母は欠測を除いた人数。
#    §7-1：退院時運動FIM の到達割合は 70 点・77 点の2行（閾値は 03 の fig2_meta から読む）。
# -----------------------------------------------------------------------------

rule("8. Table 1")

.as_chr <- function(x) {
  x <- as.character(x)
  x[!is.na(x) & !nzchar(trimws(x))] <- NA_character_
  trimws(x)
}
.as_date <- function(x) {
  if (inherits(x, "Date"))   return(x)
  if (inherits(x, "POSIXt")) return(as.Date(x))
  as.Date(as.character(x), format = "%Y-%m-%d")
}

#' 在棟日数（退院日 − 入院日 + LOS_OFFSET）。退院日は main にあればそれを、
#' 無ければワークスペースの Alldata から id と入院日で照合して取る。
compute_los <- function(M) {
  if ("los" %in% names(M)) return(list(los = as.numeric(M$los), source = "main$los"))
  if ("day_out" %in% names(M)) {
    return(list(los = as.numeric(.as_date(M$day_out) - .as_date(M$day_in)) + LOS_OFFSET,
                source = "main$day_out"))
  }
  if (!exists("Alldata", inherits = TRUE)) {
    return(list(los = rep(NA_real_, nrow(M)), source = "unavailable"))
  }
  A <- as.data.frame(get("Alldata", inherits = TRUE), stringsAsFactors = FALSE)
  if (!all(c("id", "day_in", "day_out") %in% names(A))) {
    return(list(los = rep(NA_real_, nrow(M)), source = "Alldata lacks id/day_in/day_out"))
  }
  key_a <- paste(.as_chr(A$id), format(.as_date(A$day_in)))
  key_m <- paste(.as_chr(M$id), format(.as_date(M$day_in)))
  idx <- match(key_m, key_a)                  # 同一キーが複数あれば元の行順で先頭（02 と同じ）
  dout <- .as_date(A$day_out)[idx]
  list(los = as.numeric(dout - .as_date(M$day_in)) + LOS_OFFSET,
       source = "Alldata (matched on id and date of admission)",
       n_unmatched = sum(is.na(idx)))
}

make_table1 <- function() {
  id <- "Table1"
  M  <- as.data.frame(need("04", "main"), stringsAsFactors = FALSE)
  M$age_group <- factor(as.character(M$age_group), levels = AGE_LEVELS)
  THR <- parse_thresholds(getx("03", "fig2_meta"))  # §7-1：70, 77
  if (is.null(THR)) {
    THR <- sort(THRESHOLDS_EXPECTED)
    check("Table 1：到達割合の閾値を 03 の fig2_meta から読めた", FALSE,
          sprintf("03_m0_ranktest.rda が無いか旧版。THRESHOLDS_EXPECTED（%s）を使った",
                  thr_text(THR, ", ")))
  }

  ## 在棟日数
  L <- compute_los(M)
  M$los <- L$los
  say("  length of stay source : ", L$source)
  if (identical(L$source, "unavailable") || all(is.na(M$los))) {
    check("Table 1：在棟日数を算出できた", FALSE,
          "Alldata がワークスペースに無い。Alldata を読み込んでから実行すると埋まる")
  } else {
    check("Table 1：在棟日数を算出できた（全例で照合できた）", sum(is.na(M$los)) == 0L,
          sprintf("missing = %d / %d ; definition = day_out - day_in%s", sum(is.na(M$los)),
                  nrow(M), if (LOS_OFFSET != 0) sprintf(" + %d", LOS_OFFSET) else ""))
    check("Table 1：在棟日数に負の値がない", !any(M$los < 0, na.rm = TRUE),
          sprintf("negative = %d", sum(M$los < 0, na.rm = TRUE)))
  }

  groups <- c(AGE_LEVELS, if (TABLE1_TOTAL) "Total")
  sel <- lapply(groups, function(g) if (g == "Total") rep(TRUE, nrow(M))
                else !is.na(M$age_group) & M$age_group == g)
  n_g <- vapply(sel, sum, integer(1))

  rows <- list()
  add <- function(...) rows[[length(rows) + 1L]] <<- trow(..., cols = groups)
  cont <- function(v, label) {
    add(label, vapply(sel, function(s) fmt_med_iqr(M[[v]][s]), character(1)))
  }
  catv <- function(v, label) {
    x  <- M[[v]]
    lv <- if (is.factor(x)) levels(x) else sort(unique(as.character(x[!is.na(x)])))
    add(paste0(label, ", n (%)"), rep("", length(groups)), merge = TRUE)
    for (l in lv) {
      add(relabel_levels(v, l, strict = FALSE),
          vapply(sel, function(s) {
            xs <- as.character(x[s])
            fmt_n_pct(sum(xs == l, na.rm = TRUE), sum(!is.na(xs)))
          }, character(1)), indent = 1L)
    }
  }

  cont("age", paste0(relabel_vars("age"), ", years, median [Q1, Q3]"))
  catv("sex",        relabel_vars("sex"))
  catv("class",      relabel_vars("class"))
  catv("support_in", relabel_vars("support_in"))
  catv("period",     relabel_vars("period"))
  cont("mFIM_in",  paste0(relabel_vars("mFIM_in"),  ", median [Q1, Q3]"))
  cont("cFIM_in",  paste0(relabel_vars("cFIM_in"),  ", median [Q1, Q3]"))
  cont("mFIM_out", paste0(relabel_vars("mFIM_out"), ", median [Q1, Q3]"))
  for (t in THR) {
    add(sprintf("%s \u2265 %g, n (%%)", relabel_vars("mFIM_out"), t),
        vapply(sel, function(s) {
          y <- M$mFIM_out[s]
          fmt_n_pct(sum(y >= t, na.rm = TRUE), sum(!is.na(y)))
        }, character(1)))
  }
  cont("los", paste0(relabel_vars("los"), ", median [Q1, Q3]"))

  tab <- do.call(rbind, rows)
  header <- c("", sprintf("%s\n(n = %d)", ifelse(groups == "Total", "Total",
                                                  age_label(groups)), n_g))

  ## --- 点検 -------------------------------------------------------------------
  an <- getx("04", "analysis_n")
  if (!is.null(an)) {
    n_e1 <- an$n[an$set == "analysis_set_E1_E2" & an$group == "total"]
    check("Table 1：N が解析対象（04 の analysis_n）と一致する",
          length(n_e1) == 1L && n_e1 == nrow(M), sprintf("%d vs %s", nrow(M), n_e1))
  }
  thr3 <- getx("03", "fig2a_thr")
  if (!is.null(thr3)) {
    ## fig2a_thr は 閾値 × 年齢群 の行をもつ（§7 版の 03）。閾値ごとに突き合わせる
    for (t in THR) {
      mine <- vapply(AGE_LEVELS, function(g) sum(M$mFIM_out[M$age_group == g] >= t,
                                                 na.rm = TRUE), numeric(1))
      z    <- thr3[thr3$y == t, , drop = FALSE]
      ref3 <- z$n_ge_threshold[match(AGE_LEVELS, z$age_group)]
      check(sprintf("Table 1：運動FIM \u2265 %g の例数が Fig. 2(a)（03）と一致する", t),
            isTRUE(all(mine == ref3)),
            sprintf("Table 1: %s / 03: %s", paste(mine, collapse = ", "),
                    paste(ref3, collapse = ", ")))
    }
  }
  miss_vars <- c("sex", "class", "support_in", "period", "mFIM_in", "cFIM_in")
  n_miss <- vapply(miss_vars, function(v) sum(is.na(M[[v]])), integer(1))
  say("  missing values in the analysis set : ",
      paste(sprintf("%s=%d", miss_vars, n_miss), collapse = ", "))
  check("Table 1：群間比較の数値（p値・標準化差）を置いていない（§8.1）", TRUE, "")

  ## --- 脚注 -------------------------------------------------------------------
  notes <- c(
    paste("Values are median [first quartile, third quartile] (quantile type 1, which",
          "returns observed values) or n (%). No between-group comparison is shown."),
    paste("Percentages are calculated among patients with non-missing values;",
          "the numbers of missing values are given in Supplementary Table 1."),
    paste("Pre-admission care need: Independent, no certification or certified as",
          "requiring support; Needed, certified as requiring long-term care (care levels",
          "1 to 5)."),
    "Admission period: fiscal years (April to March) of admission.",
    sprintf("Length of stay = date of discharge \u2212 date of admission%s (days).",
            if (LOS_OFFSET != 0) sprintf(" + %d", LOS_OFFSET) else ""),
    paste("CVD, cerebrovascular disease; MSD, musculoskeletal disorder; DS, disuse",
          "syndrome; FIM, Functional Independence Measure; mFIM, motor FIM;",
          "cFIM, cognitive FIM; Q1, first quartile; Q3, third quartile.")
  )
  spec <- list(id = id, title = TITLES[[id]],
               subtitle = sprintf("Analysis set (main analysis), N = %d", nrow(M)),
               blocks = list(list(heading = NULL, ft = make_ft(tab, header))),
               notes = notes)
  file <- file.path(TAB_DIR, "Table1.docx")
  write_table_docx(spec, file)

  out <- tab[, c("label", groups)]
  names(out) <- c("row", ifelse(groups == "Total", "Total", age_label(groups)))
  utils::write.csv(out, file.path(OUT_DIR, "table1.csv"), row.names = FALSE,
                   fileEncoding = "UTF-8")
  say("    written: ", file.path(OUT_DIR, "table1.csv"))
  FIG_DATA$Table1 <<- list(table = out, los_source = L$source, thresholds = THR)
  file
}
run_part("Table1", make_table1)


# -----------------------------------------------------------------------------
# 9. Table 2（E-3・E-4、§10.1、plan.summary §7-1・§7-2）
#    部分比例オッズモデル（M2-PPO）による年齢群別の標準化 P(Y ≥ t)、3対比の RD、
#    3対比の閾値別オッズ比（t = 70、77。閾値は 06 の table2 から読む）。
#    閾値を列に置く（著者の決定、2026-09-22）。比較として比例オッズモデル（M2）の
#    共通オッズ比を併記する。共通オッズ比は閾値によらない1つの値なので、
#    閾値の列を結合したセルに置く。
#    p値と RR は載せない（§8.5-1、§8.4-3）。
# -----------------------------------------------------------------------------

rule("9. Table 2")

make_table2 <- function() {
  id <- "Table2"
  t2   <- need("06", "table2")
  meta <- getx("06", "table2_meta")
  fn   <- getx("06", "table2_footnotes")
  if (!"threshold" %in% names(t2) || !any(t2$block == "Odds ratio") ||
      !is.null(getx("06", "bin_tab"))) {
    stop("06_standardize_boot.rda が §7 以前の版（閾値65・共通オッズ比のみ）である。",
         "script_section7 の 05 → 06 を実行してから作り直すこと")
  }
  THR      <- sort(unique(stats::na.omit(t2$threshold)))
  THR_META <- parse_thresholds(meta)
  N        <- meta_get(meta, "analysis_n", as.character(t2$analysis_n[1]))

  check("Table 2：p値を置いていない（§8.5-1）", all(is.na(t2$p_value)), "")
  check("Table 2：RR を置いていない（§8.4-3）",
        !any(grepl("risk ratio", c(t2$quantity, t2$scale), ignore.case = TRUE)), "")
  check("Table 2：表の閾値が 06 の諸元（table2_meta）と一致する",
        !is.null(THR_META) && isTRUE(all.equal(as.numeric(THR), as.numeric(THR_META))),
        sprintf("table2: %s / table2_meta: %s", thr_text(THR, ", "),
                if (is.null(THR_META)) "(none)" else thr_text(THR_META, ", ")))
  check("Table 2：閾値65の値を置いていない（§7-1）",
        !(65 %in% THR) && !any(grepl("65", t2$quantity, fixed = TRUE)), "")
  .n_or  <- sum(t2$block == "Odds ratio")
  .n_com <- sum(t2$block == "Common odds ratio")
  check("Table 2：閾値別オッズ比が 閾値 × 3対比、共通オッズ比が 3対比 ある",
        .n_or == length(THR) * length(CONTRAST_KEYS) && .n_com == length(CONTRAST_KEYS),
        sprintf("threshold-specific OR = %d, common OR = %d", .n_or, .n_com))

  cols <- paste0("t", THR)
  rows <- list()
  add <- function(...) rows[[length(rows) + 1L]] <<- trow(..., cols = cols)
  #' 1つの閾値の値（M2-PPO の行）
  cell <- function(block, key, t, digits, scale = 1) {
    z <- t2[t2$block == block & t2$row_key == key & !is.na(t2$threshold) &
              t2$threshold == t, , drop = FALSE]
    if (nrow(z) != 1L) return(NA_TEXT)
    fmt_ci(z$estimate, z$ci_lower, z$ci_upper, digits, scale)
  }
  by_thr <- function(block, key, digits, scale = 1) {
    vapply(THR, function(t) cell(block, key, t, digits, scale), character(1))
  }
  n_of <- function(g) {
    z <- t2$n[t2$block == "Standardised risk" & t2$row_key == g]
    if (length(z)) z[1] else NA
  }

  add("Standardised probability, %", "", merge = TRUE)
  for (g in AGE_LEVELS) {
    add(sprintf("%s (n = %s)", age_label(g), n_of(g)),
        by_thr("Standardised risk", g, DIGITS$risk, 100), indent = 1L)
  }
  add("Risk difference, percentage points", "", merge = TRUE)
  for (k in CONTRAST_KEYS) {
    add(contrast_name(k), by_thr("Risk difference", k, DIGITS$rd, 100), indent = 1L)
  }
  add("Odds ratio, partial proportional odds model", "", merge = TRUE)
  for (k in CONTRAST_KEYS) {
    add(contrast_name(k), by_thr("Odds ratio", k, DIGITS$or), indent = 1L)
  }
  add("Common odds ratio, proportional odds model", "", merge = TRUE)
  common_rows <- integer(0)
  for (k in CONTRAST_KEYS) {
    z <- t2[t2$block == "Common odds ratio" & t2$row_key == k, , drop = FALSE]
    v <- if (nrow(z) == 1L) fmt_ci(z$estimate, z$ci_lower, z$ci_upper, DIGITS$or) else NA_TEXT
    add(contrast_name(k), c(v, rep("", length(THR) - 1L)), indent = 1L)
    common_rows <- c(common_rows, length(rows))
  }
  tab <- do.call(rbind, rows)

  ## 値のセルがすべて埋まっていること（共通オッズ比の行は結合する前の2列目以降を除く）
  .vals <- as.matrix(tab[!tab$.merge, cols, drop = FALSE])
  .vals[match(common_rows, which(!tab$.merge)), -1L] <- "x"
  check("Table 2：すべての値のセルが埋まっている", !any(.vals == NA_TEXT | .vals == ""),
        sprintf("empty = %d", sum(.vals == NA_TEXT | .vals == "")))

  header <- c("", sprintf("%s ≥ %g\nEstimate (95%% CI)", relabel_vars("mFIM_out"), THR))
  ft <- make_ft(tab, header)
  if (length(THR) > 1L) {
    .j <- seq_along(THR) + 1L
    for (i in common_rows) ft <- flextable::merge_at(ft, i = i, j = .j, part = "body")
    ft <- flextable::align(ft, i = common_rows, j = .j, align = "center", part = "body")
  }

  notes <- c(footnote_text(fn, roman = TRUE),
             paste("Contrasts are written older vs younger; a negative risk difference and",
                   "an odds ratio below 1 indicate a lower motor FIM at discharge in the",
                   "older group. No p-values are reported for these estimates."),
             paste("The common odds ratio does not depend on the threshold and is shown once",
                   "across the columns."),
             "CI, confidence interval; FIM, Functional Independence Measure; mFIM, motor FIM.")
  spec <- list(id = id, title = TITLES[[id]],
               subtitle = sprintf("M2 complete-case population, N = %s", N),
               blocks = list(list(heading = NULL, ft = ft)),
               notes = notes)
  file <- file.path(TAB_DIR, "Table2.docx")
  write_table_docx(spec, file)
  FIG_DATA$Table2 <<- list(table = tab, thresholds = THR, common_rows = common_rows)
  file
}
run_part("Table2", make_table2)


# -----------------------------------------------------------------------------
# 10. 補足表1（§8.1、§10.2）
#     (A) 変数別の欠測、(B) 除外者（E1〜E3）の割合と特性、(C) E3 の欠測理由を1表に。
# -----------------------------------------------------------------------------

rule("10. Supplementary Table 1")

make_supptable1 <- function() {
  id <- "SuppTable1"
  sm <- need("11", "supp1_missing")
  se <- need("11", "supp1_excluded")
  s3 <- need("11", "supp1_e3")
  fn <- getx("11", "table_supp1_footnotes")
  M04 <- getx("04", "main")
  gs <- c(AGE_LEVELS, "total")

  rows <- list()
  add <- function(...) rows[[length(rows) + 1L]] <<- trow(..., cols = gs)

  ## (A)
  add("(A) Missing values, n (%)", "", bold = TRUE, merge = TRUE)
  for (pop in unique(sm$population)) {
    add(pop, "", italic = TRUE, merge = TRUE, indent = 0L)
    s <- sm[sm$population == pop, , drop = FALSE]
    for (v in unique(s$variable)) {
      z <- s[s$variable == v, , drop = FALSE]
      z <- z[match(gs, z$age_group), , drop = FALSE]
      add(z$variable_label[1],
          ifelse(z$n > 0, sprintf("%d (%s)", as.integer(z$n_missing),
                                  fmt_num(z$percent_missing, DIGITS$percent)), "0"),
          indent = 1L)
    }
  }

  ## (B)
  add("(B) Excluded and included patients", "", bold = TRUE, merge = TRUE)
  for (k in unique(se$excl_group)) {
    s  <- se[se$excl_group == k, , drop = FALSE]
    ct <- s[s$variable == "(count)", , drop = FALSE]
    ct <- ct[match(gs, ct$age_group), , drop = FALSE]
    add(paste0(ct$excl_group_label[1], ", n (% of A4)"),
        fmt_n_pct(ct$n, ct$denominator), italic = TRUE)
    for (v in setdiff(unique(s$variable[s$statistic == "median [Q1, Q3]"]), "(count)")) {
      z <- s[s$variable == v, , drop = FALSE]
      z <- z[match(gs, z$age_group), , drop = FALSE]
      add(paste0(z$variable_label[1], ", median [Q1, Q3]"),
          ifelse(z$n > 0, fmt_med_iqr3(z$median, z$q1, z$q3), NA_TEXT), indent = 1L)
    }
    for (v in unique(s$variable[s$statistic == "n (%)"])) {
      z0 <- s[s$variable == v, , drop = FALSE]
      add(paste0(z0$variable_label[1], ", n (%)"), "", indent = 1L, merge = TRUE)
      ## 水準の並びは 04 の因子の水準（参照水準が先頭）に揃える。欠測は最後
      lv  <- unique(z0$level)
      ref <- if (!is.null(M04) && v %in% names(M04) && is.factor(M04[[v]])) levels(M04[[v]])
             else character(0)
      lv  <- lv[order(is.na(match(lv, ref)), match(lv, ref), lv == "(missing)")]
      for (l in lv) {
        z <- z0[z0$level == l, , drop = FALSE]
        z <- z[match(gs, z$age_group), , drop = FALSE]
        if (identical(l, "(missing)") && all(z$n == 0)) next
        add(z$level_label[1], fmt_n_pct(z$n, z$denominator), indent = 2L)
      }
    }
  }

  ## (C)
  add("(C) Reason for a missing motor FIM at discharge (E3), n (% of E3)", "",
      bold = TRUE, merge = TRUE)
  for (k in unique(s3$category)) {
    z <- s3[s3$category == k, , drop = FALSE]
    z <- z[match(gs, z$age_group), , drop = FALSE]
    lab <- if ("label" %in% names(z)) z$label[1] else k
    add(lab, fmt_n_pct(z$n, z$n_e3), indent = 1L)   # 感度分析Sでの扱いは脚注（11）にある
  }
  tab <- do.call(rbind, rows)
  header <- c("", age_label(AGE_LEVELS), "Total")

  spec <- list(id = id, title = TITLES[[id]],
               subtitle = "Base population (A4) and analysis set (main analysis)",
               blocks = list(list(heading = NULL, ft = make_ft(tab, header))),
               notes = c(footnote_text(fn),
                         "FIM, Functional Independence Measure; Q1, first quartile; Q3, third quartile."))
  file <- file.path(TAB_DIR, "SuppTable1.docx")
  write_table_docx(spec, file)
  FIG_DATA$SuppTable1 <<- list(table = tab)
  file
}
run_part("SuppTable1", make_supptable1)


# -----------------------------------------------------------------------------
# 11. 補足表2（§8.8、§10.2、plan.summary §7-1〜§7-3）
#     主解析と感度分析S の E-2 と、M2-PPO 由来の E-4（70 点・77 点の標準化確率と
#     RD）を並べる。閾値65の二値モデルの行は外した（§7-2）。閾値は 07_2 の
#     table_supp2 の threshold 列で見分ける（E-2 の行は NA）。
#     p値は載せない。S のオッズ比（E-3）は載せない（§8.8）。
# -----------------------------------------------------------------------------

rule("11. Supplementary Table 2")

make_supptable2 <- function() {
  id <- "SuppTable2"
  s2    <- need("07_2", "table_supp2")
  meta2 <- getx("07_2", "table_supp2_meta")
  fn2   <- getx("07_2", "table_supp2_footnotes")
  if (!"threshold" %in% names(s2) || any(grepl("binary", s2$block, ignore.case = TRUE))) {
    stop("07_2_sens_clm_boot.rda が §7 以前の版（閾値65の二値モデルを含む）である。",
         "script_section7 の 07_1 → 07_2 を実行してから作り直すこと")
  }
  THR      <- sort(unique(stats::na.omit(s2$threshold)))
  THR_META <- parse_thresholds(meta2)
  n_main_re <- meta_get(getx("03", "fig2_meta"), "analysis_set_n", NA_TEXT)
  n_main_cc <- meta_get(getx("06", "table2_meta"), "analysis_n", NA_TEXT)
  n_s_re    <- meta_get(meta2, "sensitivity_S_population", NA_TEXT)
  n_s_cc    <- meta_get(meta2, "m2_S_complete_case_n", NA_TEXT)

  check("補足表2：p値を置いていない（§8.5、§8.8）", all(is.na(s2$p_value)), "")
  check("補足表2：S のオッズ比を置いていない（§8.8）",
        !any(grepl("odds ratio", s2$quantity, ignore.case = TRUE)), "")
  check("補足表2：二値モデル・閾値65の行がない（§7-1、§7-2）",
        !(65 %in% THR) &&
          !any(grepl("binary|65", paste(s2$block, s2$quantity), ignore.case = TRUE)), "")
  check("補足表2：表の閾値が 07_2 の諸元（table_supp2_meta）と一致する",
        !is.null(THR_META) && isTRUE(all.equal(as.numeric(THR), as.numeric(THR_META))),
        sprintf("table_supp2: %s / meta: %s", thr_text(THR, ", "),
                if (is.null(THR_META)) "(none)" else thr_text(THR_META, ", ")))
  check("補足表2：主解析の E-4 がある（07_2 が §7 版の 06 を読んだ）",
        any(s2$analysis == "main" & s2$block == "Standardised risk (E-4)"),
        "無ければ §7 版の 06 → 07_2 の順に実行し直すこと")

  an <- c("main", "sensitivity S")
  cols <- c("main", "sens")
  #' t = NA は閾値をもたない行（E-2）
  val <- function(block, key, t = NA_real_, digits, scale) {
    vapply(an, function(a) {
      z <- s2[s2$analysis == a & s2$block == block & s2$row_key == key, , drop = FALSE]
      z <- if (is.na(t)) z[is.na(z$threshold), , drop = FALSE]
           else z[!is.na(z$threshold) & z$threshold == t, , drop = FALSE]
      if (!nrow(z)) return(NA_TEXT)
      fmt_ci(z$estimate[1], z$ci_lower[1], z$ci_upper[1], digits, scale)
    }, character(1))
  }
  rows <- list()
  add <- function(...) rows[[length(rows) + 1L]] <<- trow(..., cols = cols)

  add("Two-sample relative effect", "", merge = TRUE, bold = TRUE)
  add("N", c(n_main_re, n_s_re), indent = 1L, italic = TRUE)
  for (k in CONTRAST_KEYS) add(contrast_name(k), val("Relative effect (E-2)", k,
                                                      digits = DIGITS$re, scale = 1), indent = 1L)
  for (i in seq_along(THR)) {
    t <- THR[i]
    add(sprintf("Standardised P(%s ≥ %g), %%", relabel_vars("mFIM_out"), t), "",
        merge = TRUE, bold = TRUE)
    if (i == 1L) add("N (complete cases)", c(n_main_cc, n_s_cc), indent = 1L, italic = TRUE)
    for (g in AGE_LEVELS) add(age_label(g), val("Standardised risk (E-4)", g, t,
                                                digits = DIGITS$risk, scale = 100), indent = 1L)
    add(sprintf("Risk difference in P(%s ≥ %g), percentage points",
                relabel_vars("mFIM_out"), t), "", merge = TRUE, bold = TRUE)
    for (k in CONTRAST_KEYS) add(contrast_name(k), val("Risk difference (E-4)", k, t,
                                                        digits = DIGITS$rd, scale = 100),
                                 indent = 1L)
  }
  tab <- do.call(rbind, rows)
  .vals <- as.matrix(tab[!tab$.merge & !tab$.italic, cols, drop = FALSE])
  check("補足表2：主解析・S のすべての値のセルが埋まっている", !any(.vals == NA_TEXT),
        sprintf("empty = %d（主解析の列が空なら 06 → 07_2 の順に実行し直す）",
                sum(.vals == NA_TEXT)))
  header <- c("", "Main analysis\nestimate (95% CI)", "Sensitivity analysis S\nestimate (95% CI)")

  spec <- list(id = id, title = TITLES[[id]],
               subtitle = sprintf(paste("Main analysis: analysis set; sensitivity analysis S:",
                                        "base population (A4) without E3 categories (b) and (c),",
                                        "N = %s"), n_s_re),
               blocks = list(list(heading = NULL, ft = make_ft(tab, header))),
               notes = c(footnote_text(fn2),
                         "CI, confidence interval; FIM, Functional Independence Measure; mFIM, motor FIM."))
  file <- file.path(TAB_DIR, "SuppTable2.docx")
  write_table_docx(spec, file)
  FIG_DATA$SuppTable2 <<- list(table = tab, thresholds = THR)
  file
}
run_part("SuppTable2", make_supptable2)


# -----------------------------------------------------------------------------
# 12. 補足表3（§8.7、§10.2）
#     (A) 度数10未満のセルの一覧、(B) 期間区分 × 年齢群の2元表、
#     (C) 群別 ridit の第5・95百分位と重なりの範囲。
# -----------------------------------------------------------------------------

rule("12. Supplementary Table 3")

make_supptable3 <- function() {
  id <- "SuppTable3"
  sparse <- need("09", "sparse")
  pbg    <- need("09", "period_by_agegroup")
  rq     <- need("09", "ridit_overlap_q")
  ro     <- need("09", "ridit_overlap_r")
  ex     <- getx("09", "extrapolation")
  meta   <- getx("09", "meta")
  fn     <- getx("09", "footnotes")
  cut    <- or_else(getx("09", "SPARSE_CUT"), 10L)
  xvars  <- or_else(getx("09", "CROSS_VARS"), c("age_group", "sex", "class", "support_in", "period"))
  gs     <- c(AGE_LEVELS, "total")
  blocks <- list()

  ## (A) 疎なセル
  lab_cols <- paste0(xvars, "_label")
  if (nrow(sparse) && all(lab_cols %in% names(sparse))) {
    ta <- sparse[, c(lab_cols, "n"), drop = FALSE]
    ta[] <- lapply(ta, as.character)
    names(ta) <- c("label", paste0("v", seq_along(xvars)[-1]), "n")
    ta$.indent <- 0L; ta$.bold <- FALSE; ta$.italic <- FALSE; ta$.merge <- FALSE
  } else {
    ta <- data.frame(label = "(none)", stringsAsFactors = FALSE)
    for (k in seq_along(xvars)[-1]) ta[[paste0("v", k)]] <- ""
    ta$n <- ""
    ta$.indent <- 0L; ta$.bold <- FALSE; ta$.italic <- FALSE; ta$.merge <- TRUE
  }
  blocks[[1]] <- list(
    heading = sprintf("(A) Combinations of %s with no patients or fewer than %d patients",
                      paste(tolower(relabel_vars(xvars)), collapse = " \u00d7 "), cut),
    ft = make_ft(ta, c(relabel_vars(xvars), "n")))

  ## (B) 期間区分 × 年齢群
  pb <- pbg[pbg$population == "M2 complete-case population", , drop = FALSE]
  rows <- list()
  for (i in seq_len(nrow(pb))) {
    v <- c(vapply(AGE_LEVELS, function(g)
             sprintf("%d (%s)", as.integer(pb[[g]][i]),
                     fmt_num(pb[[paste0(g, "_percent")]][i], DIGITS$percent)), character(1)),
           sprintf("%d (%s)", as.integer(pb$total[i]),
                   fmt_num(100 * pb$total[i] / sum(pb$total), DIGITS$percent)))
    rows[[i]] <- trow(pb$period_label[i], v, cols = gs)
  }
  blocks[[2]] <- list(heading = "(B) Admission period \u00d7 age group, n (% of the age group)",
                      ft = make_ft(do.call(rbind, rows),
                                   c(relabel_vars("period"), age_label(AGE_LEVELS), "Total")))

  ## (C) ridit の百分位と重なりの範囲
  rows <- list()
  merge_rows <- integer(0)
  for (v in unique(rq$variable)) {
    z <- rq[rq$variable == v, , drop = FALSE]
    z <- z[match(gs, z$age_group), , drop = FALSE]
    rows[[length(rows) + 1L]] <- trow(paste0(z$variable_label[1], " (ridit)"), "",
                                      cols = gs, merge = TRUE)
    rows[[length(rows) + 1L]] <- trow("5th percentile", fmt_num(z$p05, DIGITS$ridit),
                                      cols = gs, indent = 1L)
    rows[[length(rows) + 1L]] <- trow("95th percentile", fmt_num(z$p95, DIGITS$ridit),
                                      cols = gs, indent = 1L)
    o <- ro[ro$variable == v, , drop = FALSE]
    for (j in seq_len(nrow(o))) {
      rng <- if (isTRUE(o$is_empty[j])) "empty" else
        sprintf("%s\u2013%s", fmt_num(o$lower[j], DIGITS$ridit), fmt_num(o$upper[j], DIGITS$ridit))
      rows[[length(rows) + 1L]] <- trow(sprintf("Overlap range, %s", o$overlap_kind[j]),
                                        c(rng, rep("", length(gs) - 1L)), cols = gs,
                                        indent = 1L)
      merge_rows <- c(merge_rows, length(rows))
    }
  }
  tc <- do.call(rbind, rows)
  ftc <- make_ft(tc, c("", age_label(AGE_LEVELS), "Total"))
  for (i in merge_rows) ftc <- flextable::merge_at(ftc, i = i, j = 2:(length(gs) + 1L),
                                                   part = "body")
  blocks[[3]] <- list(heading = "(C) Ridit scores of the admission FIM by age group",
                      ft = ftc)

  ## (D) 任意：重なりの範囲の外にある例数
  if (isTRUE(SUPP3_SHOW_EXTRAPOLATION) && !is.null(ex)) {
    rows <- list()
    for (v in unique(ex$variable)) {
      for (kd in unique(ex$overlap_kind[ex$variable == v])) {
        z <- ex[ex$variable == v & ex$overlap_kind == kd, , drop = FALSE]
        z <- z[match(gs, z$age_group), , drop = FALSE]
        rows[[length(rows) + 1L]] <- trow(sprintf("%s, outside the %s", z$variable_label[1], kd),
                                          fmt_n_pct(z$n_outside, z$n), cols = gs)
      }
    }
    blocks[[4]] <- list(heading = "(D) Patients outside the overlap range, n (%)",
                        ft = make_ft(do.call(rbind, rows),
                                     c("", age_label(AGE_LEVELS), "Total")))
  }

  spec <- list(id = id, title = TITLES[[id]],
               subtitle = sprintf("M2 complete-case population (standardisation population), N = %s",
                                  meta_get(meta, "n", NA_TEXT)),
               blocks = blocks,
               notes = c(footnote_text(fn), "FIM, Functional Independence Measure."))
  file <- file.path(TAB_DIR, "SuppTable3.docx")
  write_table_docx(spec, file)
  FIG_DATA$SuppTable3 <<- list(sparse = ta, ridit = tc)
  file
}
run_part("SuppTable3", make_supptable3)


# -----------------------------------------------------------------------------
# 13. 補足表4（§8.3、§10.2、plan.summary §7-2）
#     M2-PPO の係数。年齢群は閾値別の係数（nominal 母数。05 で符号を反転済みで、
#     exp(estimate) が Y ≥ t のオッズ比）を閾値ごとに示す。共変量の係数は2つの
#     閾値で共通。年齢群以外は解釈しない（Table 2 fallacy）。
#     §7 以前の比例オッズモデルの係数表（閾値別の係数が無いもの）も従来どおり作れるが、
#     Table 2 と食い違うため点検に出す。
# -----------------------------------------------------------------------------

rule("13. Supplementary Table 4")

make_supptable4 <- function() {
  id <- "SuppTable4"
  co11 <- getx("11", "table_supp4")
  co05 <- getx("05", "m2_coefficients")
  co <- or_else(co11, co05)
  if (is.null(co) || !nrow(co)) stop("補足表4の係数（11 の table_supp4 または 05 の m2_coefficients）が無い")
  fn <- getx("11", "table_supp4_footnotes")
  cl <- getx("04", "cov_levels")

  is_ts  <- co$type == "threshold-specific coefficient"
  is_ppo <- any(is_ts)
  check("補足表4：部分比例オッズモデル（M2-PPO、§7-2）の係数である", is_ppo,
        if (is_ppo) sprintf("threshold-specific = %d, common = %d, thresholds = %d",
                            sum(is_ts), sum(co$type == "regression coefficient"),
                            sum(co$type == "threshold"))
        else "年齢群の閾値別の係数が無い（§7 以前の 05 の出力）。§7 版の 05 → 11 を実行し直すこと")
  if (!is.null(co11) && nrow(co11) && !is.null(co05)) {
    .same <- nrow(co11) == nrow(co05) && identical(as.character(co11$term),
                                                   as.character(co05$term)) &&
      isTRUE(all.equal(as.numeric(co11$estimate), as.numeric(co05$estimate)))
    check("補足表4：11 の係数表が 05 の m2_coefficients と一致する（11 を 05 の後に実行した）",
          .same, if (.same) "" else "11 が古い 05 の出力を転記している。05 の後に 11 を実行し直すこと")
  }

  ## 閾値別の係数の閾値：05 の PPO_CUTS（閾値 → clm の閾値名）で term から引く。
  ## 無ければ term_label の末尾（"P(... >= 70)"）から読む。
  thr_of <- rep(NA_real_, nrow(co))
  if (is_ppo) {
    cuts05 <- getx("05", "PPO_CUTS")
    if (!is.null(cuts05)) {
      cut_of <- sub("\\.age_group.*$", "", co$term)
      thr_of[is_ts] <- suppressWarnings(as.numeric(names(cuts05)[match(cut_of[is_ts], cuts05)]))
    }
    .miss <- is_ts & is.na(thr_of)
    thr_of[.miss] <- suppressWarnings(as.numeric(
      sub("^.*(>=|≥)\\s*([0-9.]+)\\)?\\s*$", "\\2", co$term_label[.miss])))
    check("補足表4：閾値別の係数の閾値をすべて特定できた", !anyNA(thr_of[is_ts]),
          sprintf("thresholds = %s", thr_text(sort(unique(thr_of[is_ts])), ", ")))
  }
  ts_thr <- sort(unique(stats::na.omit(thr_of[is_ts])))

  cols <- c("b", "se", "or")
  is_sp <- grepl("^ns\\(", co$term)
  rows <- list()
  add <- function(...) rows[[length(rows) + 1L]] <<- trow(..., cols = cols)
  one <- function(i, label = co$term_label[i], indent = 1L) {
    show_or <- co$type[i] %in% c("regression coefficient", "threshold-specific coefficient") &&
      (!is_sp[i] || SUPP4_OR_FOR_SPLINE)
    add(gsub(">=", "≥", label, fixed = TRUE),
        c(fmt_num(co$estimate[i], DIGITS$coef), fmt_num(co$std_error[i], DIGITS$coef),
          if (show_or) fmt_num(co$odds_ratio[i], DIGITS$or) else NA_TEXT),
        indent = indent)
  }
  reg <- which(co$type == "regression coefficient")
  if (is_ppo) {
    add("Age group, threshold-specific", "", merge = TRUE, bold = TRUE)
    for (t in ts_thr) {
      add(sprintf("%s ≥ %g", relabel_vars("mFIM_out"), t), "", merge = TRUE,
          italic = TRUE, indent = 1L)
      for (i in which(is_ts & !is.na(thr_of) & thr_of == t)) {
        one(i, label = sub(",\\s*P\\(.*$", "", co$term_label[i]), indent = 2L)
      }
    }
    add("Covariates, common to both thresholds (reported for transparency, not interpreted)",
        "", merge = TRUE, bold = TRUE)
  } else {
    add("Age group", "", merge = TRUE, bold = TRUE)
    for (i in reg[co$is_age_group[reg]]) one(i)
    add("Covariates (reported for transparency, not interpreted)", "", merge = TRUE, bold = TRUE)
  }
  for (i in reg[!co$is_age_group[reg]]) one(i)
  if (isTRUE(SUPP4_INCLUDE_THRESHOLDS)) {
    add("Threshold parameters", "", merge = TRUE, bold = TRUE)
    for (i in which(co$type == "threshold")) one(i)
  }
  tab <- do.call(rbind, rows)

  ref_note <- NULL
  if (!is.null(cl)) {
    ref_note <- paste0("Reference levels: ",
                       paste(sprintf("%s, %s", tolower(cl$variable_label), cl$reference_label),
                             collapse = "; "), ".")
  }
  notes <- c(
    if (is_ppo)
      paste0("Partial proportional odds model: the coefficients of age group differ between ",
             "the thresholds (", thr_text(ts_thr), " points) and are log odds ratios for ",
             "a motor FIM at discharge at or above each threshold; the coefficients of the ",
             "covariates are common to both thresholds. A positive coefficient indicates a ",
             "higher motor FIM at discharge. Standard errors are shown; OR = exp(coefficient).")
    else paste("Estimates are regression coefficients on the log odds scale of P(Y ≥ y)",
               "with standard errors; OR = exp(coefficient)."),
    if (is_ppo)
      paste("The age-group coefficients are the nominal parameters of ordinal::clm with the",
            "sign reversed (clm models logit P(Y ≤ j) = θj + τj − x′β)."),
    ref_note,
    if (!SUPP4_OR_FOR_SPLINE) "OR is not shown for the spline basis coefficients.",
    if (!SUPP4_INCLUDE_THRESHOLDS)
      sprintf("The %d threshold parameters are not shown.", sum(co$type == "threshold")),
    footnote_text(fn),
    "OR, odds ratio; SE, standard error; FIM, Functional Independence Measure; mFIM, motor FIM; cFIM, cognitive FIM."
  )
  spec <- list(id = id, title = TITLES[[id]],
               subtitle = sprintf("M2 complete-case population, N = %s",
                                  meta_get(getx("06", "table2_meta"), "analysis_n", NA_TEXT)),
               blocks = list(list(heading = NULL,
                                  ft = make_ft(tab, c("", "Coefficient", "SE", "OR")))),
               notes = notes)
  file <- file.path(TAB_DIR, "SuppTable4.docx")
  write_table_docx(spec, file)
  FIG_DATA$SuppTable4 <<- list(table = tab, thresholds = ts_thr)
  file
}
run_part("SuppTable4", make_supptable4)


# -----------------------------------------------------------------------------
# 14. 表を1つの docx にまとめる（任意）
# -----------------------------------------------------------------------------

rule("14. All tables in one document")

if (isTRUE(MAKE_ALL_TABLES_DOCX) && length(TABLE_SPECS)) {
  .ord  <- intersect(c("Table1", "Table2", "SuppTable1", "SuppTable2", "SuppTable3",
                       "SuppTable4"), names(TABLE_SPECS))
  .file <- file.path(TAB_DIR, "AllTables.docx")
  .res <- tryCatch({
    doc <- officer::read_docx()
    for (i in seq_along(.ord)) {
      doc <- docx_add_spec(doc, TABLE_SPECS[[.ord[i]]])
      if (i < length(.ord)) doc <- officer::body_add_break(doc)
    }
    print(doc, target = .file)
    say("    written: ", .file)
    TRUE
  }, error = function(e) { say("  エラー : ", conditionMessage(e)); FALSE })
  record_output("AllTables", if (.res) .file else "", if (.res) "written" else "failed")
}


# -----------------------------------------------------------------------------
# 15. 出来上がりの点検と保存
# -----------------------------------------------------------------------------

rule("15. Verification")

## グレースケールであること（設定した色がすべて R = G = B）
.cols <- unique(c(GROUP_STYLE$colour, GROUP_STYLE$fill, REF_COLOUR, GUIDE_COLOUR, BAND_FILL,
                  POINT_MAIN$fill, POINT_SUB$fill, "black", "white"))
.rgb  <- grDevices::col2rgb(.cols)
check("図の色がすべて灰色（R = G = B）である",
      all(.rgb[1, ] == .rgb[2, ] & .rgb[2, ] == .rgb[3, ]),
      paste(.cols, collapse = ", "))
check("年齢群が線種・濃さ・記号の組で区別できる",
      !anyDuplicated(paste(GROUP_STYLE$colour, GROUP_STYLE$linetype, GROUP_STYLE$shape)), "")

## §7-1：表題の閾値が rda の閾値と合っているか（TITLES は手で書くため）
.thr_used <- or_else(FIG_DATA$Table2$thresholds, THR_SOURCES[["06 table2_meta"]])
if (isTRUE(MAKE[["Table2"]]) && length(.thr_used)) {
  check("Table 2 の表題に rda の閾値がすべて書かれている（§7-1）",
        all(vapply(.thr_used, function(t) grepl(sprintf("\u2265 %g\\b", t), TITLES$Table2),
                   logical(1))),
        sprintf("thresholds = %s ; title = %s", thr_text(.thr_used, ", "), TITLES$Table2))
}

## §7-2：旧版で作った図のうち、この版では作らないファイルが残っていないか。
## 補足図の番号を繰り上げたため、旧版の SuppFig2.tiff（連続年齢の図）は不要になった。
.stale <- file.path(FIG_DIR, "SuppFig2.tiff")
.stale <- .stale[file.exists(.stale)]
check("旧版の図のファイルが output/figures に残っていない（§7-2）", !length(.stale),
      if (length(.stale)) paste0(paste(basename(.stale), collapse = ", "),
                                 " は旧版の出力（この版では作らない）。削除するか別の場所へ移すこと")
      else "")

outputs <- do.call(rbind, .outputs)
for (i in seq_len(nrow(outputs))) {
  if (outputs$status[i] == "written") {
    f <- outputs$file[i]
    check(sprintf("[%s] ファイルがある", outputs$item[i]),
          file.exists(f) && file.info(f)$size > 0,
          sprintf("%s (%s bytes)", f, format(file.info(f)$size, big.mark = ",")))
  }
}
say("\n  display switches : title = ", SHOW$title, ", subtitle = ", SHOW$subtitle,
    ", footnote = ", SHOW$footnote)
print(outputs, row.names = FALSE)

checks <- do.call(rbind, .checks)
say("\n  REVIEW 項目: ", sum(checks$result == "REVIEW"), " / ", nrow(checks))
if (any(checks$result == "REVIEW")) print(checks[checks$result == "REVIEW", ], row.names = FALSE)

rule("16. Output")

utils::write.csv(outputs, file.path(OUT_DIR, "table_outputs_12.csv"), row.names = FALSE,
                 fileEncoding = "UTF-8")
utils::write.csv(checks, file.path(OUT_DIR, "table_checks_12.csv"), row.names = FALSE,
                 fileEncoding = "UTF-8")
say("  written: ", file.path(OUT_DIR, "table_outputs_12.csv"))
say("  written: ", file.path(OUT_DIR, "table_checks_12.csv"))

save(FIG_DATA, TABLE_SPECS, outputs, checks,
     file = file.path(OUT_DIR, "12_figures_tables.rda"))
say("  written: ", file.path(OUT_DIR, "12_figures_tables.rda"))

say("\n--- sessionInfo ---")
print(utils::sessionInfo())

say("\nDone. figures -> ", FIG_DIR, " ; tables -> ", TAB_DIR)

.log_close()
