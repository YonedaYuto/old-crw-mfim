# =============================================================================
# 03_m0_ranktest.R
#   解析計画書 §8.2「M0：順位にもとづく比較」（plan.summary.txt §4-1）の実装。
#   あわせて Fig. 2 を描くための数値をすべて出力する（作図は別スクリプト）。
#
#   入力 : output/02_preprocess.rda（dat_main。02_preprocess.R の出力）
#          01_labels.R（表示用の文字列）
#   出力 : output/table_m0_relative_effect.csv  Fig. 2(b) の数値（E-2）
#          output/fig2a_prob_by_agegroup.csv    Fig. 2(a) の階段関数（E-1）
#          output/fig2_meta.csv                 Fig. 2 の脚注に必要な諸元
#          output/fig2_footnotes.csv            Fig. 2 の脚注（下書き）
#          output/table_checks_03.csv           点検結果
#          output/03_m0_ranktest.rda            上記＋ログのみの記録
#          output/log_03_m0_ranktest.txt        実行ログ
#
#   計画との対応（plan.detail §11 のスクリプト10）
#     (a) rankFD の大域検定（ATS）と群別 pseudo-rank 相対効果 → ログのみ（§8.2、§8.5）
#     (b) 3対比の順列 Brunner–Munzel（method = "logit"、nperm = 100000、
#         conf.level = 0.95）と未調整95%順列信頼区間 → E-2（出力表）。
#         主対比 G3対G1 の順列p値だけを出力表に入れる
#     (c) 順列分布を自前で生成し、|T*| の95%分位点による対称区間と
#         p = mean(|T*| ≥ |T|) を算出して (b) との差をログに記録する（§8.2）
#
#   方針
#     1. 期待例数・期待結果のハードコードは置かない（§13-4）。点検は「集計して
#        報告する」にとどめ、閾値で解析を分岐させない。
#     2. 本文・図表に載せる区間は `rankFD::rank.two.samples()` の等裾95%順列区間
#        （logit スケール）である（§8.2、付録 AA10）。自前の順列分布は (c) の
#        比較と、パッケージの返り値の構造が版によって変わった場合の受け皿に使う。
#     3. 出力される文字列はすべて英語（01_labels.R）。
#     4. 年齢群は 02_preprocess.R と同じ規則でこのスクリプト内で作る（02 は変更
#        しない）。規則が一致することを 02 の出力（flow_g）と突き合わせて点検する。
#
#   実行時間の目安（§8.10、README に記録する）
#     自前の順列は N ≈ 1,600 の対比で 100,000 回あたり約30秒。3対比で約1.5分。
#     rankFD の順列を加えても数分の範囲に収まる見込み。実測値はログに出る。
#
#   文字コード : UTF-8（Windows の R では source(..., encoding = "UTF-8")）
#
#   改訂 : 2026-09-21 plan.summary.txt §7-1 により、Fig. 2(a) の縦線（閾値）を
#          65 点から 70 点・77 点の2本に変更した。あわせて is_threshold、
#          fig2a_thr、閾値の行の点検、fig2_meta の threshold を2つの閾値に
#          合わせた。Brunner–Munzel 検定（5〜8 節）は変更していない。乱数の
#          使い方も変わらないため、同じ乱数種で同じ結果になる。
# =============================================================================


# -----------------------------------------------------------------------------
# 0. 設定
# -----------------------------------------------------------------------------

OUT_DIR <- "output"
IN_RDA  <- file.path(OUT_DIR, "02_preprocess.rda")

## --- 年齢群（§5.3）----------------------------------------------------------
## 02_preprocess.R の AGE_BREAKS と同じ値であること。§3-2 で突き合わせる。
AGE_BREAKS  <- c(65, 75, 90)          # G1: 65-74, G2: 75-89, G3: >=90
AGE_LEVELS  <- c("G1", "G2", "G3")

## --- アウトカム（§6.1、§6.2）------------------------------------------------
OUTCOME     <- "mFIM_out"             # 退院時運動FIM。13〜91点の順序尺度
MFIM_RANGE  <- c(13, 91)

## --- Fig. 2(a)（§10.1）------------------------------------------------------
## 年齢群別の P(退院時運動FIM ≥ y) を y = 13〜91 の階段関数で描き、70点と77点に縦線を引く
## §7-1（2026-09-21）：入浴と階段昇降以外が見守りレベルとなる 70 点（旧設定の 65 点を
## 修正）と、入浴と階段昇降以外が修正自立レベルとなる 77 点の2つ。
FIG2A_GRID  <- seq(MFIM_RANGE[1], MFIM_RANGE[2], by = 1)
FIG2_VLINE  <- c(70, 77)              # plan.summary §7-1 の閾値（昇順に置く）

## --- 順列 Brunner–Munzel（§8.2、§8.10）--------------------------------------
CONF_LEVEL  <- 0.95                   # 未調整95%（§8.5）
NPERM       <- 100000L                # rankFD に渡す順列回数（§8.2）
NPERM_OWN   <- NPERM                  # (c) の自前順列の回数。比較のため同数にする
BM_METHOD   <- "logit"                # 範囲保存のため（§8.2）
ALTERNATIVE <- "two.sided"
SEED        <- 20260920L              # 乱数種。対比ごとに SEED + offset を使う（下記）
## 対比ごとの乱数種の割り当て。rankFD は大域の乱数状態を使うため、呼び出しの直前に
## 明示的に set.seed() する。こうしておかないと、先に走る自前順列の回数を変えた
## だけで本文の数値が変わる。
seed_own    <- function(i) SEED + 10L * i          # (c) 自前順列
seed_rankfd <- function(i) SEED + 10L * i + 1L     # (b) rankFD

## --- 3対比（§8.2）-----------------------------------------------------------
## 高齢側を第2水準に置く。第1水準を x、第2水準を y として
## p = P(x < y) + 0.5 × P(x = y) が推定されるため、この順序で
## p = P(高齢側 > 若年側) + 0.5 × P(同点) となる（§8.2）。
CONTRASTS <- list(
  list(key = "G2_vs_G1", ref = "G1", cmp = "G2", role = "secondary"),
  list(key = "G3_vs_G1", ref = "G1", cmp = "G3", role = "primary"),
  list(key = "G3_vs_G2", ref = "G2", cmp = "G3", role = "secondary")
)

## 検証的主要仮説はこの1つだけ。順列p値を出力表に入れるのもこれだけ（§1.2、§8.5）
PRIMARY_KEY <- "G3_vs_G1"


# -----------------------------------------------------------------------------
# 1. 表示用ラベル
# -----------------------------------------------------------------------------

if (!file.exists("01_labels.R")) {
  stop("01_labels.R が見つからない。03_m0_ranktest.R は 01_labels.R と同じ",
       "フォルダを作業ディレクトリにして実行すること。")
}
source("01_labels.R", encoding = "UTF-8")

## 対比の表示名。"G3 vs G1 (≥90 years vs 65-74 years)" の形にする。
contrast_label <- function(ref, cmp) {
  sprintf("%s vs %s (%s vs %s)", cmp, ref,
          relabel_levels("age_group", cmp), relabel_levels("age_group", ref))
}


# -----------------------------------------------------------------------------
# 2. 実行の準備（出力先とログ）
# -----------------------------------------------------------------------------

if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

.log_con <- file(file.path(OUT_DIR, "log_03_m0_ranktest.txt"), open = "wt",
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

rule(paste0("03_m0_ranktest.R  ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
say("R : ", R.version.string)
say("seed (base) : ", SEED,
    "  [own = SEED + 10*i, rankFD = SEED + 10*i + 1 ; i は対比の番号]")
say("nperm (rankFD) : ", NPERM, " ; nperm (own) : ", NPERM_OWN)
say("method : ", BM_METHOD, " ; conf.level : ", CONF_LEVEL,
    " ; alternative : ", ALTERNATIVE)


# -----------------------------------------------------------------------------
# 3. 入力の取り込み
# -----------------------------------------------------------------------------

rule("3. Input")

if (!exists("dat_main", inherits = TRUE)) {
  if (!file.exists(IN_RDA)) {
    stop(IN_RDA, " が見つからない。先に 02_preprocess.R を実行すること。")
  }
  load(IN_RDA)                      # dat_base / dat_main / dat_sens / flow / flow_g ...
  say("loaded : ", IN_RDA)
} else {
  say("dat_main はワークスペース上のものを使う")
}
if (!exists("dat_main", inherits = TRUE)) {
  stop("dat_main が見つからない。02_preprocess.R の出力を確認すること。")
}

dat <- as.data.frame(dat_main, stringsAsFactors = FALSE)
say("dat_main : ", nrow(dat), " rows x ", ncol(dat), " cols")

.need <- c("id", "age", OUTCOME)
.miss <- setdiff(.need, names(dat))
if (length(.miss)) stop("dat_main に必要な列がない: ", paste(.miss, collapse = ", "))


# -----------------------------------------------------------------------------
# 3-2. 年齢群（§5.3）
#      02_preprocess.R の make_agegroup() と同じ規則をここに再掲する。
#      02 は年齢群を集計にしか使わず dat_main に列として残していないため。
#      規則が一致していることは、02 の出力 flow_g の MAIN 行と突き合わせて点検する。
# -----------------------------------------------------------------------------

rule("3-2. Age group")

make_agegroup <- function(age) {
  factor(ifelse(is.na(age), NA_character_,
         ifelse(age >= AGE_BREAKS[3], "G3",
         ifelse(age >= AGE_BREAKS[2], "G2", "G1"))),
         levels = AGE_LEVELS)
}

dat$age_group <- make_agegroup(dat$age)

n_by_g <- table(dat$age_group)
for (g in AGE_LEVELS) {
  say(sprintf("  %-3s %-14s n = %5d", g,
              relabel_levels("age_group", g), as.integer(n_by_g[[g]])))
}
say(sprintf("  %-3s %-14s n = %5d", "", "total", nrow(dat)))

check("年齢群に欠測がない", !anyNA(dat$age_group),
      sprintf("missing = %d", sum(is.na(dat$age_group))))
check("年齢群の各群に少なくとも2例ある（BM検定の前提）",
      all(as.integer(n_by_g) >= 2),
      paste(sprintf("%s=%d", AGE_LEVELS, as.integer(n_by_g)), collapse = ", "))

## 02 の出力との突き合わせ（年齢群の規則が 02 と一致するか）
if (exists("flow_g", inherits = TRUE) && "stage" %in% names(flow_g) &&
    "MAIN" %in% flow_g$stage) {
  .row <- flow_g[flow_g$stage == "MAIN", , drop = FALSE][1, ]
  .exp <- as.integer(c(.row$G1, .row$G2, .row$G3))
  check("年齢群の例数が 02_preprocess.R の flow_g（MAIN）と一致する",
        identical(as.integer(n_by_g), .exp),
        sprintf("03: %s / 02: %s",
                paste(as.integer(n_by_g), collapse = ", "),
                paste(.exp, collapse = ", ")))
} else {
  check("02 の flow_g と突き合わせできた", FALSE,
        "flow_g が見つからないため突き合わせを行っていない")
}


# -----------------------------------------------------------------------------
# 3-3. アウトカムの点検（分布の記述は §8.1／Fig. 2(a) で行う）
# -----------------------------------------------------------------------------

rule("3-3. Outcome checks")

y_all <- suppressWarnings(as.numeric(dat[[OUTCOME]]))
dat[[OUTCOME]] <- y_all

check("退院時運動FIMに欠測がない（E3 で除外済み）", !anyNA(y_all),
      sprintf("missing = %d", sum(is.na(y_all))))
check(sprintf("退院時運動FIMが %d〜%d の範囲にある", MFIM_RANGE[1], MFIM_RANGE[2]),
      all(is.na(y_all) | (y_all >= MFIM_RANGE[1] & y_all <= MFIM_RANGE[2])),
      sprintf("out of range = %d",
              sum(!is.na(y_all) & (y_all < MFIM_RANGE[1] | y_all > MFIM_RANGE[2]))))
check("退院時運動FIMが整数である",
      all(is.na(y_all) | abs(y_all - round(y_all)) < 1e-8),
      sprintf("non-integer = %d",
              sum(!is.na(y_all) & abs(y_all - round(y_all)) >= 1e-8)))
check("解析対象の id が一意（解析単位＝患者、§5.2）", !anyDuplicated(dat$id))

say(sprintf("\n  observed levels of %s : %d (range %d - %d)",
            OUTCOME, length(unique(y_all[!is.na(y_all)])),
            as.integer(min(y_all, na.rm = TRUE)), as.integer(max(y_all, na.rm = TRUE))))

## 同点の多さ（順位法の性質に関わるため記録する。判定には使わない）
.tie_tab <- table(y_all)
say(sprintf("  ties : largest tied group = %d (%.1f%% of N)",
            as.integer(max(.tie_tab)), 100 * max(.tie_tab) / length(y_all)))

## 解析に使う集団を固定する（欠測があれば落とす。上の点検で0件のはず）
ok_row <- !is.na(dat$age_group) & !is.na(y_all)
if (any(!ok_row)) {
  say(sprintf("\n  dropped rows with missing age group or outcome : %d", sum(!ok_row)))
}
dat_m0 <- dat[ok_row, , drop = FALSE]
N_M0   <- nrow(dat_m0)
say(sprintf("\nM0 analysis set : %d", N_M0))


# -----------------------------------------------------------------------------
# 4. Fig. 2(a) の数値（E-1、§8.1・§10.1）
#    年齢群別の P(退院時運動FIM ≥ y)、y = 13〜91 の階段関数
# -----------------------------------------------------------------------------

rule("4. Fig. 2(a)  P(mFIM at discharge >= y) by age group")

fig2a <- do.call(rbind, lapply(AGE_LEVELS, function(g) {
  v <- dat_m0[[OUTCOME]][dat_m0$age_group == g]
  n <- length(v)
  k <- vapply(FIG2A_GRID, function(y) sum(v >= y), integer(1))
  data.frame(
    age_group       = g,
    age_group_label = relabel_levels("age_group", g),
    n               = n,
    y               = FIG2A_GRID,
    n_at_or_above   = k,
    prob            = if (n > 0) k / n else NA_real_,
    is_threshold    = FIG2A_GRID %in% FIG2_VLINE,
    stringsAsFactors = FALSE
  )
}))

## 点検：y = 13 で 1、単調非増加、各閾値の行が群ごとに1本
.p13 <- fig2a$prob[fig2a$y == MFIM_RANGE[1]]
check(sprintf("P(Y >= %d) = 1（全例が下限以上）", MFIM_RANGE[1]),
      all(abs(.p13 - 1) < 1e-12), paste(sprintf("%.6f", .p13), collapse = ", "))
.mono <- vapply(AGE_LEVELS, function(g) {
  p <- fig2a$prob[fig2a$age_group == g]
  all(diff(p) <= 1e-12)
}, logical(1))
check("階段関数が単調非増加である", all(.mono),
      paste(sprintf("%s=%s", AGE_LEVELS, ifelse(.mono, "ok", "NG")), collapse = ", "))
## 年齢群 × 閾値 の各セルにちょうど1行（閾値が y の格子にない場合もここで分かる）
.thr_cnt <- table(factor(fig2a$age_group[fig2a$is_threshold], levels = AGE_LEVELS),
                  factor(fig2a$y[fig2a$is_threshold], levels = FIG2_VLINE))
check("各閾値の行が年齢群ごとに1本ずつある",
      all(.thr_cnt == 1L) &&
        sum(fig2a$is_threshold) == length(AGE_LEVELS) * length(FIG2_VLINE),
      sprintf("rows = %d (expected %d = %d groups x %d thresholds; y = %s)",
              sum(fig2a$is_threshold), length(AGE_LEVELS) * length(FIG2_VLINE),
              length(AGE_LEVELS), length(FIG2_VLINE),
              paste(FIG2_VLINE, collapse = ", ")))

## 閾値 70 点・77 点での到達割合（Fig. 2(a) の注記・Table 1 と同じ量。記述であって比較ではない）
## 行は 閾値 × 年齢群 の 6 行。閾値の昇順に並べ、その中を G1 → G3 とする。
## どの閾値の行かは y 列で見分ける（列の構成は変えていない）。
fig2a_thr <- fig2a[fig2a$is_threshold,
                   c("age_group", "age_group_label", "n", "y",
                     "n_at_or_above", "prob")]
fig2a_thr <- fig2a_thr[order(fig2a_thr$y, match(fig2a_thr$age_group, AGE_LEVELS)), ,
                       drop = FALSE]
rownames(fig2a_thr) <- NULL
names(fig2a_thr)[names(fig2a_thr) == "n_at_or_above"] <- "n_ge_threshold"
fig2a_thr$percent <- 100 * fig2a_thr$prob
say(sprintf("\n  P(Y >= y) at y = %s by age group (descriptive, not a between-group comparison):",
            paste(FIG2_VLINE, collapse = ", ")))
print(fig2a_thr, row.names = FALSE)


# -----------------------------------------------------------------------------
# 5. Brunner–Munzel の自前実装（(c) の順列分布と、(b) の受け皿）
#    §8.2 に記した `rankFD::rank.two.samples()` の仕様をそのまま実装する。
#      p̂   = P(x < y) + 0.5 × P(x = y)、x は第1水準、y は第2水準
#      se_p = sqrt(n1 σ1² + n2 σ2²) / (n1 n2)      （Brunner–Munzel 2000）
#      f    = logit(p̂)、se_f = se_p / (p̂ (1 − p̂))  （デルタ法）
#      T    = f / se_f
#    順列は「群ラベルを固定し応答を並べ替える」。応答の並べ替えは結合標本の
#    中間順位の並べ替えと同値であるため（順位は値の単調変換で同点構造も保つ）、
#    中間順位を一度だけ計算して使い回す。
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
#' @param ra1 第1水準（若年側）の観測値がもつ結合標本での中間順位
#' @param ra2 第2水準（高齢側）の観測値がもつ結合標本での中間順位
bm_logit <- function(ra1, ra2) {
  n1 <- length(ra1); n2 <- length(ra2)
  r1 <- rank(ra1);   r2 <- rank(ra2)          # 群内の中間順位
  m1 <- mean(ra1);   m2 <- mean(ra2)
  p  <- (m2 - (n2 + 1) / 2) / n1              # = P(x < y) + 0.5 P(x = y)
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

#' 観測統計量と順列分布から、(b) の等裾区間と (c) の |T*| 対称区間を作る
bm_intervals <- function(obs, Tp, conf.level = CONF_LEVEL) {
  alpha <- 1 - conf.level
  Tp_ok <- Tp[is.finite(Tp)]
  n_bad <- length(Tp) - length(Tp_ok)

  ## (b) の再現：符号付き統計量の等裾分位点による検定の反転（§8.2）
  crit1 <- stats::quantile(Tp_ok, 1 - alpha / 2, names = FALSE)   # type 7（既定）
  crit2 <- stats::quantile(Tp_ok,     alpha / 2, names = FALSE)
  lo_et <- expit_f(obs[["f"]] - crit1 * obs[["se_f"]])
  hi_et <- expit_f(obs[["f"]] - crit2 * obs[["se_f"]])
  p_et  <- min(1, 2 * min(mean(Tp_ok <= obs[["T"]]), mean(Tp_ok >= obs[["T"]])))

  ## (c)：|T*| の (1 − α) 分位点による対称区間
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

#' 関数に実在する引数だけを渡して呼ぶ（パッケージの版で引数名が変わっても落ちない）
call_known_args <- function(fun, args) {
  fml  <- setdiff(names(formals(fun)), "...")
  drop <- setdiff(names(args), fml)
  if (length(drop)) {
    say("    （この版にない引数は渡さない: ", paste(drop, collapse = ", "), "）")
  }
  do.call(fun, args[names(args) %in% fml])
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
# 7. (a) rankFD の大域検定（ATS）と群別 pseudo-rank 相対効果
#        §8.2・§8.5 により、本文にも補足資料にも載せない。ログにのみ残す。
# -----------------------------------------------------------------------------

rule("7. (a) Global ATS and group-wise pseudo-rank relative effects  [LOG ONLY]")
say("  §8.2・§8.5：本文にも補足資料にも報告しない。記録としてログにのみ残す。")

has_rankFD <- requireNamespace("rankFD", quietly = TRUE)
check("rankFD が利用できる", has_rankFD,
      if (has_rankFD) paste0("version ",
                             as.character(utils::packageVersion("rankFD")))
      else "install.packages(\"rankFD\") が必要")

m0_global <- NULL
if (has_rankFD) {
  .d <- data.frame(y = dat_m0[[OUTCOME]],
                   age_group = factor(as.character(dat_m0$age_group),
                                      levels = AGE_LEVELS))
  m0_global <- tryCatch({
    call_known_args(rankFD::rankFD,
                    list(formula = y ~ age_group, data = .d,
                         effect = "unweighted", hypothesis = "H0p",
                         CI.method = "logit", alpha = 1 - CONF_LEVEL,
                         info = FALSE))
  }, error = function(e) {
    say("  rankFD::rankFD() でエラー: ", conditionMessage(e))
    NULL
  })
  if (!is.null(m0_global)) {
    say("\n  --- rankFD::rankFD() output (log only) ---")
    print(m0_global)
  }
}
check("大域検定の結果を出力表・図表に入れていない", TRUE,
      "ログと rda にのみ保存する（§8.5）")


# -----------------------------------------------------------------------------
# 8. (b) 3対比の順列 Brunner–Munzel ＝ E-2（Fig. 2(b)）
#        あわせて (c) 自前の順列分布による比較を行う
# -----------------------------------------------------------------------------

rule("8. (b) Permuted Brunner-Munzel for the three contrasts  [E-2, Fig. 2(b)]")

res_rows  <- list()   # 出力表（Fig. 2(b)）
cmp_rows  <- list()   # (c) の比較（ログのみ）
perm_meta <- list()

for (ict in seq_along(CONTRASTS)) {
  ct  <- CONTRASTS[[ict]]
  key <- ct$key; ref <- ct$ref; cmp <- ct$cmp
  lab <- contrast_label(ref, cmp)
  say("\n", strrep("=", 74))
  say(sprintf("Contrast %s : %s", key, lab))

  ## --- 因子水準の順序をログに出す（§8.2：高齢側を第2水準に置く）-------------
  lv <- c(ref, cmp)
  say(sprintf("  factor levels passed to rank.two.samples() : c(\"%s\", \"%s\")",
              lv[1], lv[2]))
  say(sprintf("    level 1 (x, younger) = %s = %s", lv[1],
              relabel_levels("age_group", lv[1])))
  say(sprintf("    level 2 (y, older)   = %s = %s", lv[2],
              relabel_levels("age_group", lv[2])))
  say("    estimand : p = P(x < y) + 0.5 x P(x = y) = P(older > younger) + 0.5 P(tie)")

  sub <- dat_m0[dat_m0$age_group %in% lv, , drop = FALSE]
  sub$age_group <- factor(as.character(sub$age_group), levels = lv)
  sub <- sub[order(sub$age_group), , drop = FALSE]      # 第1水準を前に並べる
  n1  <- sum(sub$age_group == lv[1]); n2 <- sum(sub$age_group == lv[2])
  say(sprintf("  n(%s) = %d ; n(%s) = %d ; N = %d", lv[1], n1, lv[2], n2, n1 + n2))

  ## --- 自前の観測統計量と順列分布（(c)、および (b) の受け皿）----------------
  ra  <- rank(sub[[OUTCOME]])                 # 結合標本の中間順位
  obs <- bm_logit(ra[seq_len(n1)], ra[(n1 + 1L):(n1 + n2)])
  say(sprintf("  own estimate : p.hat = %.6f ; se(p) = %.6f ; logit f = %.6f ; T = %.4f",
              obs[["p"]], obs[["se_p"]], obs[["f"]], obs[["T"]]))

  t0 <- proc.time()[["elapsed"]]
  Tp <- bm_permute(ra, n1, NPERM_OWN, seed_own(ict))
  ivs <- bm_intervals(obs, Tp, CONF_LEVEL)
  t_own <- proc.time()[["elapsed"]] - t0
  say(sprintf("  own permutation : %d draws in %.1f s, seed = %d (dropped non-finite: %d)",
              NPERM_OWN, t_own, seed_own(ict), ivs$n_perm_dropped))
  say(sprintf("    equal-tailed crits : q(%.3f) = %.4f ; q(%.3f) = %.4f",
              (1 - CONF_LEVEL) / 2, ivs$crit_lower,
              1 - (1 - CONF_LEVEL) / 2, ivs$crit_upper))
  say(sprintf("    |T*| crit          : q(%.2f) = %.4f", CONF_LEVEL, ivs$crit_abs))

  ## --- rankFD（本文・図表に載せる値）----------------------------------------
  rf_est <- rf_lo <- rf_hi <- rf_p <- NA_real_
  rf_ok  <- FALSE
  t_rf   <- NA_real_
  if (has_rankFD) {
    .d2 <- data.frame(y = sub[[OUTCOME]], age_group = sub$age_group)
    say(sprintf("  rankFD seed : %d", seed_rankfd(ict)))
    set.seed(seed_rankfd(ict))
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
      say("\n  --- rankFD::rank.two.samples() structure (log only) ---")
      utils::str(res, max.level = 3, give.attr = FALSE)
      cand <- extract_two_samples(res)
      if (nrow(cand)) {
        say("\n  candidate rows found in the return value:")
        print(cand, row.names = FALSE, digits = 6)
      }
      pick <- pick_permutation_row(cand, obs[["p"]])
      if (!is.null(pick)) {
        rf_est <- pick$estimator; rf_lo <- pick$lower
        rf_hi  <- pick$upper;     rf_p  <- pick$p_value
        rf_ok  <- TRUE
        say(sprintf("\n  picked permutation row : %s", pick$row_label))
      } else {
        say("\n  順列法の行を特定できなかった。上の構造を確認すること。")
      }
    }
  }

  ## --- 本文・図表に載せる値を決める ------------------------------------------
  if (rf_ok) {
    src <- "rankFD::rank.two.samples (permutation, logit)"
    est <- rf_est; lo <- rf_lo; hi <- rf_hi
    pv  <- if (is.finite(rf_p)) rf_p else ivs$et[["p_value"]]
    if (!is.finite(rf_p)) {
      src <- paste0(src, " + own permutation p-value")
    }
  } else {
    src <- "own reproduction (equal-tailed permutation, logit)"
    est <- obs[["p"]]; lo <- ivs$et[["lower"]]; hi <- ivs$et[["upper"]]
    pv  <- ivs$et[["p_value"]]
  }
  check(sprintf("[%s] 本文の区間を rankFD から取得できた", key), rf_ok,
        if (rf_ok) "rankFD の順列区間を採用"
        else "自前の再現値に切り替えた。rankFD の返り値の構造をログで確認すること")

  ## --- 点検：採用した行の相対効果（行の選択条件でもある。値をログに残す）-----
  if (rf_ok) {
    d_est <- abs(rf_est - obs[["p"]])
    check(sprintf("[%s] 採用した行の相対効果が自前の計算と一致する", key),
          is.finite(d_est) && d_est < 1e-6,
          sprintf("rankFD = %.8f, own = %.8f, diff = %.2e", rf_est, obs[["p"]], d_est))
  }
  ## --- 点検：rankFD の区間と自前の等裾再現が同じ構成法とみて矛盾しないか -----
  ## 両者は独立に引いた順列であるため、分位点のモンテカルロ誤差の分だけずれる。
  ## 桁違いにずれる場合は構成法が想定と違う（版の変更など）ことを疑う。
  if (rf_ok) {
    d_ci <- max(abs(rf_lo - ivs$et[["lower"]]), abs(rf_hi - ivs$et[["upper"]]))
    check(sprintf("[%s] rankFD の区間と自前の等裾再現の差が 0.01 未満", key),
          is.finite(d_ci) && d_ci < 0.01,
          sprintf("rankFD [%.6f, %.6f] / own [%.6f, %.6f], max diff = %.2e（順列のモンテカルロ誤差の範囲か確認する）",
                  rf_lo, rf_hi, ivs$et[["lower"]], ivs$et[["upper"]], d_ci))
  }

  check(sprintf("[%s] 区間が (0, 1) に収まる（method = \"logit\" の範囲保存）", key),
        is.finite(lo) && is.finite(hi) && lo > 0 && hi < 1 && lo <= est && est <= hi,
        sprintf("[%.6f, %.6f], p.hat = %.6f", lo, hi, est))

  ## --- 出力表の1行（Fig. 2(b)）----------------------------------------------
  ## 順列p値は主対比だけに入れる（§8.5-2、§8.5-3）
  res_rows[[length(res_rows) + 1L]] <- data.frame(
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
    p_value           = if (identical(key, PRIMARY_KEY)) pv else NA_real_,
    ## 順列p値が 0 になった場合に図表へそのまま「P = 0」と書かないための表示用文字列。
    ## 値そのものは p_value 列にある（§8.5-3：p値を記すのは本文1か所と Fig. 2(b)）。
    p_value_display   = if (identical(key, PRIMARY_KEY)) format_perm_p(pv, NPERM) else "",
    p_value_method    = if (identical(key, PRIMARY_KEY))
                          "permuted Brunner-Munzel, two-sided" else "",
    source            = src,
    stringsAsFactors  = FALSE
  )

  ## --- (c) の比較（ログのみ、§8.2）------------------------------------------
  cmp_rows[[length(cmp_rows) + 1L]] <- data.frame(
    contrast            = key,
    label               = lab,
    p_hat_own           = obs[["p"]],
    p_hat_rankFD        = rf_est,
    et_lower_rankFD     = rf_lo,
    et_upper_rankFD     = rf_hi,
    et_lower_own        = ivs$et[["lower"]],
    et_upper_own        = ivs$et[["upper"]],
    abs_lower_own       = ivs$abs[["lower"]],
    abs_upper_own       = ivs$abs[["upper"]],
    p_value_rankFD      = rf_p,
    p_value_et_own      = ivs$et[["p_value"]],
    p_value_abs_own     = ivs$abs[["p_value"]],
    diff_lower_et_vs_abs = ivs$abs[["lower"]] - ivs$et[["lower"]],
    diff_upper_et_vs_abs = ivs$abs[["upper"]] - ivs$et[["upper"]],
    diff_p_et_vs_abs     = ivs$abs[["p_value"]] - ivs$et[["p_value"]],
    crit_lower           = ivs$crit_lower,
    crit_upper           = ivs$crit_upper,
    crit_abs             = ivs$crit_abs,
    perm_asymmetry       = ivs$crit_upper + ivs$crit_lower,  # 0 なら対称
    seed_own             = seed_own(ict),
    seed_rankFD          = seed_rankfd(ict),
    seconds_own          = t_own,
    seconds_rankFD       = t_rf,
    stringsAsFactors     = FALSE
  )

  perm_meta[[key]] <- list(obs = obs, intervals = ivs,
                           n1 = n1, n2 = n2, levels = lv,
                           perm_quantiles = stats::quantile(
                             Tp[is.finite(Tp)],
                             c(0.005, 0.025, 0.25, 0.5, 0.75, 0.975, 0.995),
                             names = TRUE))
}

m0_relative_effect <- do.call(rbind, res_rows)
m0_perm_compare    <- do.call(rbind, cmp_rows)

say("\n", strrep("=", 74))
say("E-2 : relative effect and unadjusted 95% permutation CI (Fig. 2(b))")
print(m0_relative_effect[, c("contrast", "label", "n_reference", "n_comparison",
                             "role", "relative_effect", "ci_lower", "ci_upper",
                             "p_value")],
      row.names = FALSE, digits = 4)


# -----------------------------------------------------------------------------
# 8-2. (c) 等裾区間と |T*| 対称区間の差  [LOG ONLY]
#      差が無視できない大きさであった場合は、順列分布の非対称性を限界の節に
#      1文で記す（§8.2）。記録を見てから本文の構成法を選び直すことはしない。
# -----------------------------------------------------------------------------

rule("8-2. (c) Equal-tailed vs |T*| symmetric interval  [LOG ONLY]")

print(m0_perm_compare[, c("contrast", "et_lower_own", "et_upper_own",
                          "abs_lower_own", "abs_upper_own",
                          "diff_lower_et_vs_abs", "diff_upper_et_vs_abs",
                          "p_value_et_own", "p_value_abs_own",
                          "diff_p_et_vs_abs", "perm_asymmetry")],
      row.names = FALSE, digits = 4)

say("\n  perm_asymmetry = q(0.975) + q(0.025)。順列分布が0のまわりで対称なら0に近い。")
say("  この表は §8.2 によりログにのみ残す。本文・図表には載せない。")

for (i in seq_len(nrow(m0_perm_compare))) {
  r <- m0_perm_compare[i, ]
  say(sprintf("\n  %s : 区間端の差 lower %+.4f / upper %+.4f、p値の差 %+.4f",
              r$contrast, r$diff_lower_et_vs_abs, r$diff_upper_et_vs_abs,
              r$diff_p_et_vs_abs))
}
check("(c) の比較を出力表・図表に入れていない", TRUE,
      "ログと rda にのみ保存する（§8.2）")


# -----------------------------------------------------------------------------
# 9. Fig. 2 の諸元と脚注
# -----------------------------------------------------------------------------

rule("9. Fig. 2 metadata and footnotes")

.pv <- m0_relative_effect$p_value[m0_relative_effect$contrast == PRIMARY_KEY]
.pv <- if (length(.pv)) .pv[1] else NA_real_

## threshold は §7-1 により2つ。";" で区切った1つの値にする（例 "70;77"）。
## 読む側は strsplit(value, ";") で分ける。項目を1行に保つのは、12 の meta_get() が
## 同名の項目の先頭しか返さないため（行を重ねると黙って 70 点だけが使われる）。
fig2_meta <- data.frame(
  item = c("analysis_set_n",
           paste0("n_", AGE_LEVELS),
           "outcome", "outcome_label", "outcome_min", "outcome_max",
           "threshold", "conf_level", "ci_method", "nperm_rankFD",
           "nperm_own", "seed", "primary_contrast", "primary_p_value",
           "primary_p_value_display", "test", "panel_a", "panel_b"),
  value = c(as.character(N_M0),
            as.character(as.integer(n_by_g[AGE_LEVELS])),
            OUTCOME, relabel_vars(OUTCOME),
            as.character(MFIM_RANGE[1]), as.character(MFIM_RANGE[2]),
            paste(FIG2_VLINE, collapse = ";"), as.character(CONF_LEVEL),
            "equal-tailed permutation, logit scale, unadjusted",
            as.character(NPERM), as.character(NPERM_OWN), as.character(SEED),
            PRIMARY_KEY,
            if (is.finite(.pv)) formatC(.pv, format = "g", digits = 6) else NA_character_,
            format_perm_p(.pv, NPERM),
            "permuted Brunner-Munzel test, two-sided 5%",
            sprintf("P(%s >= y), y = %d to %d, vertical lines at %s",
                    relabel_vars(OUTCOME), MFIM_RANGE[1], MFIM_RANGE[2],
                    paste(FIG2_VLINE, collapse = " and ")),
            "Relative effect with unadjusted 95% CI, reference line at 0.5"),
  stringsAsFactors = FALSE
)
print(fig2_meta, row.names = FALSE)

## 脚注（§10.1 の (i)〜(iv)）。作図スクリプトがそのまま使えるように英語で置く。
fig2_footnotes <- data.frame(
  tag = c("i", "ii", "iii", "iv"),
  text = c(
    "The two-sample relative effect need not be transitive across the three groups.",
    paste("The single confirmatory hypothesis is the G3 vs G1 contrast",
          "(>=90 years vs 65-74 years); the other two contrasts are estimates."),
    paste0("Intervals are equal-tailed permutation intervals on the logit scale ",
           "(inversion of the permutation test, unadjusted ",
           format(100 * CONF_LEVEL), "%)."),
    paste0("Analysis set: N = ", N_M0, " (",
           paste(sprintf("%s %s, n = %d", AGE_LEVELS,
                         relabel_levels("age_group", AGE_LEVELS),
                         as.integer(n_by_g[AGE_LEVELS])), collapse = "; "), ").")
  ),
  stringsAsFactors = FALSE
)
say("")
for (i in seq_len(nrow(fig2_footnotes))) {
  say(sprintf("  (%s) %s", fig2_footnotes$tag[i], fig2_footnotes$text[i]))
}


# -----------------------------------------------------------------------------
# 10. 出来上がりの点検
# -----------------------------------------------------------------------------

rule("10. Verification")

check("対比が3つある", nrow(m0_relative_effect) == 3L,
      sprintf("rows = %d", nrow(m0_relative_effect)))
check("順列p値が主対比の1つだけに入っている",
      sum(!is.na(m0_relative_effect$p_value)) == 1L &&
        !is.na(m0_relative_effect$p_value[m0_relative_effect$contrast == PRIMARY_KEY]),
      sprintf("p値のある行 = %s",
              paste(m0_relative_effect$contrast[!is.na(m0_relative_effect$p_value)],
                    collapse = ", ")))
check("3対比とも同じ構成法・同じ水準の区間である",
      length(unique(m0_relative_effect$ci_method)) == 1L &&
        length(unique(m0_relative_effect$conf_level)) == 1L,
      sprintf("method = %s, level = %s",
              paste(unique(m0_relative_effect$ci_method), collapse = " / "),
              paste(unique(m0_relative_effect$conf_level), collapse = " / ")))
check("各対比の例数が年齢群の例数と一致する",
      all(m0_relative_effect$n_reference ==
            as.integer(n_by_g[m0_relative_effect$group_reference])) &&
      all(m0_relative_effect$n_comparison ==
            as.integer(n_by_g[m0_relative_effect$group_comparison])), "")

## 主対比：区間が0.5を含まないことと、順列p値が5%未満であることの対応（§8.2）
.pr <- m0_relative_effect[m0_relative_effect$contrast == PRIMARY_KEY, ]
.excl_half <- (.pr$ci_lower > 0.5) || (.pr$ci_upper < 0.5)
.reject    <- is.finite(.pr$p_value) && .pr$p_value < (1 - CONF_LEVEL)
check("主対比で「区間が0.5を含まない」と「順列p値 < 0.05」が一致する",
      identical(.excl_half, .reject),
      sprintf("CI = [%.4f, %.4f], p = %s（分位点の補間による境界例は §8.2 の想定内）",
              .pr$ci_lower, .pr$ci_upper,
              if (is.finite(.pr$p_value)) formatC(.pr$p_value, format = "g", digits = 4)
              else "NA"))

check("Fig. 2(a) の行数が 年齢群 × y の数と一致する",
      nrow(fig2a) == length(AGE_LEVELS) * length(FIG2A_GRID),
      sprintf("%d rows (expected %d)", nrow(fig2a),
              length(AGE_LEVELS) * length(FIG2A_GRID)))
## fig2a_thr は 閾値 × 年齢群 の行をもつため、閾値ごとに合計して突き合わせる
.n_thr <- vapply(FIG2_VLINE, function(t) as.numeric(sum(fig2a_thr$n[fig2a_thr$y == t])),
                 numeric(1))
check("Fig. 2(a) の n が M0 の解析集団と一致する（閾値ごと）",
      all(.n_thr == N_M0),
      paste0(paste(sprintf("y = %g: %d", FIG2_VLINE, as.integer(.n_thr)),
                   collapse = ", "), " vs ", N_M0))
check("Fig. 2(a) に欠測がない", !anyNA(fig2a$prob), "")

checks <- do.call(rbind, .checks)
say("\n  REVIEW 項目: ", sum(checks$result == "REVIEW"), " / ", nrow(checks))
if (any(checks$result == "REVIEW")) {
  print(checks[checks$result == "REVIEW", ], row.names = FALSE)
}


# -----------------------------------------------------------------------------
# 11. 保存
#     §8.2・§8.5 により、(a) の大域検定と (c) の比較は CSV にしない
#     （ログと rda にのみ残す）。
# -----------------------------------------------------------------------------

rule("11. Output")

w <- function(x, f) {
  p <- file.path(OUT_DIR, f)
  write.csv(x, p, row.names = FALSE, fileEncoding = "UTF-8")
  say("  written: ", p)
}
w(m0_relative_effect, "table_m0_relative_effect.csv")
w(fig2a,              "fig2a_prob_by_agegroup.csv")
w(fig2a_thr,          "fig2a_threshold.csv")
w(fig2_meta,          "fig2_meta.csv")
w(fig2_footnotes,     "fig2_footnotes.csv")
w(checks,             "table_checks_03.csv")

save(m0_relative_effect, m0_perm_compare, m0_global, perm_meta,
     fig2a, fig2a_thr, fig2_meta, fig2_footnotes, checks,
     file = file.path(OUT_DIR, "03_m0_ranktest.rda"))
say("  written: ", file.path(OUT_DIR, "03_m0_ranktest.rda"))

say("\n--- sessionInfo ---")
print(utils::sessionInfo())

say("\nDone. M0 analysis set = ", N_M0,
    " ; contrasts = ", nrow(m0_relative_effect),
    " ; primary contrast = ", PRIMARY_KEY)

.log_close()
