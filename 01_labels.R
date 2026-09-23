# =============================================================================
# 01_labels.R  （改訂版 2026-09-20）
#   表示用ラベル。出力される文字列はすべて英語。
#
#   対象 : Alldata の17列と、02 以降のスクリプトが作る派生変数
#   方針 :
#     1. 本計画（plan.detail §6.1）で使う変数にだけラベルを与える。使わない列は
#        VARS_EXCLUDED に理由つきで列挙し、ラベルは与えない。
#     2. 水準ラベルは変数ごとに持つ（LEVELS_BY_VAR）。同じ文字列が変数によって
#        別の意味をもつ場合に取り違えない。旧版の平坦な LEVEL_LABELS も
#        後方互換のために残す（変数を指定しない呼び出しの受け皿）。
#     3. factor を渡しても正しく動く。旧版は factor の整数コードで添字を引くため、
#        エラーを出さずに別のラベルを返していた。
#     4. 未定義の名前は黙って素通しせず警告を出す（LABEL_STRICT）。
#        英語の表に和文が混じるのを、表を作る時点で検知するため。
#
#   文字コード : UTF-8。source("01_labels.R", encoding = "UTF-8") で読むこと。
# =============================================================================

## 未定義の変数名・水準に警告を出すか。表の作成中は TRUE のままにする。
LABEL_STRICT <- TRUE


# -----------------------------------------------------------------------------
# 1. 変数ラベル
# -----------------------------------------------------------------------------

VAR_LABELS <- c(
  ## --- Alldata の列（§6.1 で使うもの）---------------------------------------
  id          = "Patient ID",
  age         = "Age",
  sex         = "Sex",
  day_in      = "Date of admission",
  day_out     = "Date of discharge",
  disposition = "Discharge destination",
  class       = "Disease category",
  support_in  = "Pre-admission care need",
  mFIM_in     = "mFIM at admission",
  cFIM_in     = "cFIM at admission",
  mFIM_out    = "mFIM at discharge",

  ## --- 02 以降で作る派生変数 -------------------------------------------------
  age_group   = "Age group",
  period      = "Admission period",
  fy_in       = "Fiscal year of admission",
  los         = "Length of stay, days",
  mFIM_in_r   = "mFIM at admission (ridit)",
  cFIM_in_r   = "cFIM at admission (ridit)"
)

## 本計画で使わない Alldata の列。理由を残し、ラベルは与えない。
## check_alldata_labels() がこの一覧を使い、「ラベルの付け忘れ」と
## 「意図して使わない列」を区別する。
VARS_EXCLUDED <- c(
  key         = "【要確認】レコードの識別子とみられる。解析単位は患者であり id を使う（§5.2）",
  birth       = "生年月日。入院時年齢 age を使う（§5.1、著者確認済み）",
  cFIM_out    = "退院時認知FIM。変数一覧（§6.1）に含まれない",
  mFIM_diff   = "FIM利得。本計画では扱わない（§15）",
  cFIM_diff   = "FIM利得。本計画では扱わない（§15）",
  totFIM_diff = "FIM利得。本計画では扱わない（§15）"
)


# -----------------------------------------------------------------------------
# 2. 水準ラベル（変数ごと）
# -----------------------------------------------------------------------------

LEVELS_BY_VAR <- list(

  sex = c(
    "男" = "Male",
    "女" = "Female"
  ),

  ## 「その他」は元データにあれば4水準目として保持する（§6.1、§14-8）
  class = c(
    "脳血管" = "CVD",
    "運動器" = "MSD",
    "廃用"   = "DS",
    "その他" = "Other"
  ),

  support_in = c(
    "あり" = "Needed",
    "なし" = "Independent"
  ),

  ## 計画で名指しされている3水準と、それ以外をまとめた「地域社会への退院」
  ## （著者の決定、2026-09-20。§14-1 を解決）。集約は collapse_disposition() が行う。
  disposition = c(
    "医療機関"           = "Medical institution",
    "病院・診療所へ転院" = "Transfer to hospital/clinic",
    "終了（死亡等）"     = "Termination (death etc.)",
    "地域社会への退院"   = "Discharge to the community"
  ),

  ## 年齢群（§5.3）。pre-old / old / oldest-old などの英語名は原典未確認のため
  ## （§14-11）、確認できるまで年齢の範囲だけを書く。
  age_group = c(
    "G1" = "65-74 years",
    "G2" = "75-89 years",
    "G3" = "≥90 years"
  ),

  ## 入院時期区分（§6.1）。02 以降はこの水準コードを使うこと。
  period = c(
    "early" = "FY2017-2019",
    "late"  = "FY2020-2023"
  )
)

## 後方互換：変数を指定しない呼び出しのための平坦な対応表。
## 同じ文字列が複数の変数で使われている場合は先に定義されたものが優先される。
LEVEL_LABELS <- unlist(unname(LEVELS_BY_VAR))
LEVEL_LABELS <- LEVEL_LABELS[!duplicated(names(LEVEL_LABELS))]


# -----------------------------------------------------------------------------
# 2-2. 退院先の集約（§5.1、§14-1。2026-09-20 に確定）
#      元データの disposition は「自宅」「特養」など多様な値をとる。
#      本計画が区別するのは E1（転院）と E2（終了（死亡等））だけであり、
#      それ以外はいずれも解析対象として残る。そこで、計画で名指しされた3つ
#      （DISPO_KNOWN）を除いた残りすべてを「地域社会への退院」にまとめる。
#      値を列挙するのではなく「3つ以外の全て」として定義するので、元データに
#      未知の退院先があっても取りこぼさない。
#      欠測は集約しない。欠測の意味（抽出時点で在院中かどうか）は 02 が決める。
# -----------------------------------------------------------------------------

## E1・E2 の判定に使うため、元の表記のまま残す退院先
DISPO_KNOWN <- c("医療機関", "病院・診療所へ転院", "終了（死亡等）")

## 上記3つ以外をまとめる先のコード
DISPO_COMMUNITY_CODE <- "地域社会への退院"

#' 退院先を4区分に集約する（DISPO_KNOWN の3つ ＋ 地域社会への退院）
#' @param x    退院先のベクトル（factor 可）
#' @param keep 元の表記のまま残す値。既定は DISPO_KNOWN。02 では在院中を表す
#'             値があればこれに加えて渡す（集約で消さないため）
#' @return keep にも NA にも当たらない値をすべて DISPO_COMMUNITY_CODE に置換した
#'         文字列ベクトル。NA は NA のまま返す。
collapse_disposition <- function(x, keep = DISPO_KNOWN) {
  x <- as.character(x)
  if (!length(x)) return(character(0))
  out <- x
  chg <- !is.na(x) & !(x %in% keep)
  out[chg] <- DISPO_COMMUNITY_CODE
  out
}


# -----------------------------------------------------------------------------
# 3. 選抜の段階（§5.1）と解析のパイプライン
# -----------------------------------------------------------------------------

STAGE_LABELS <- c(
  A1     = "A1. Consecutive admissions, Apr 2017 - Mar 2024",
  A2     = "A2. Age ≥ 65 years at admission",
  A3     = "A3. First admission within the period",
  A4     = "A4. Discharge confirmed at data extraction",
  BASE   = "Base population (A4)",
  INHOSP = "Still hospitalised at data extraction",
  E1     = "E1. Excluded: transfer during the stay",
  E2     = "E2. Excluded: termination of stay (death etc.)",
  E3     = "E3. Excluded: missing motor FIM at discharge",
  MAIN   = "Analysis set (main analysis)",
  SENS   = "Analysis set (sensitivity analysis S)"
)

## 本計画のパイプラインは主解析と感度分析Sの2つだけである。
PIPE_LABELS <- c(
  main = "Main analysis",
  sens = "Sensitivity analysis S"
)

## E3（退院時運動FIMの欠測）の理由区分（§5.1）。
## 分類規則は「転院ならば全て急変」
E3_REASON_LABELS <- c(
  a = "(a) Not assessed owing to deterioration",
  b = "(b) Not assessed at routine discharge / not recorded",
  c = "(c) Other"
)

## 感度分析Sでの扱い（§8.8）
E3_HANDLING_LABELS <- c(
  a = "Assigned the worst score",
  b = "Excluded from the population",
  c = "Excluded from the population"
)

relabel_e3_reason <- function(x) {
  x <- as.character(x)
  if (!length(x)) return(character(0))
  hit <- x %in% names(E3_REASON_LABELS)
  out <- x
  out[hit] <- unname(E3_REASON_LABELS[x[hit]])
  out
}


# -----------------------------------------------------------------------------
# 4. 変換関数
# -----------------------------------------------------------------------------

.label_warn <- function(what, x) {
  x <- unique(x[!is.na(x)])
  if (isTRUE(LABEL_STRICT) && length(x)) {
    warning("01_labels: ラベル未定義の", what, ": ",
            paste(x, collapse = ", "), call. = FALSE)
  }
  invisible(NULL)
}

#' 水準を英語ラベルに変換する
#' @param var 変数名（LEVELS_BY_VAR の要素名）。NULL なら平坦な対応表を使う
#' @param lv  水準のベクトル。factor でも character でもよい
relabel_levels <- function(var = NULL, lv, strict = LABEL_STRICT) {
  lv <- as.character(lv)                       # factor の整数添字を避ける
  if (!length(lv)) return(character(0))
  map <- if (!is.null(var) && length(var) == 1L && !is.na(var) &&
             var %in% names(LEVELS_BY_VAR)) LEVELS_BY_VAR[[var]] else LEVEL_LABELS
  hit <- lv %in% names(map)
  out <- lv
  out[hit] <- unname(map[lv[hit]])
  if (strict) {
    .label_warn(paste0("水準（", if (is.null(var)) "変数指定なし" else var, "）"),
                lv[!hit])
  }
  out
}

#' 変数名・係数名を英語ラベルに変換する
#' 対応する形：
#'   age                       -> "Age"
#'   class:運動器              -> "Disease category: MSD"
#'   age_groupG2               -> "Age group: 75-89 years"      （clm/lm の係数名）
#'   ns(mFIM_in_r, knots = ...)1 -> "mFIM at admission (ridit), spline 1"
#'   mFIM_in_measurable / mFIM_in_c
relabel_vars <- function(x, strict = LABEL_STRICT) {
  x <- as.character(x)
  if (!length(x)) return(character(0))
  unknown <- character(0)
  cand <- names(VAR_LABELS)[order(-nchar(names(VAR_LABELS)))]  # 長い名前から照合

  one <- function(s) {
    if (is.na(s)) return(NA_character_)

    ## (1) 完全一致
    if (s %in% names(VAR_LABELS)) return(unname(VAR_LABELS[s]))

    ## (2) "変数:水準"
    if (grepl(":", s, fixed = TRUE)) {
      p    <- strsplit(s, ":", fixed = TRUE)[[1]]
      base <- p[1]
      lev  <- paste(p[-1], collapse = ":")
      blab <- relabel_vars(base, strict = FALSE)
      llab <- relabel_levels(base, lev, strict = FALSE)
      if (identical(llab, lev)) llab <- relabel_vars(lev, strict = FALSE)
      return(paste0(blab, ": ", llab))
    }

    ## (3) splines::ns() の基底（Table 2 に出る形）
    m <- regmatches(s, regexec("^ns\\(\\s*([^,)]+).*\\)([0-9]+)$", s))[[1]]
    if (length(m) == 3L) {
      return(paste0(relabel_vars(trimws(m[2]), strict = FALSE), ", spline ", m[3]))
    }

    ## (4) clm / lm の因子係数名「変数名＋水準名」
    ##     残りが既知の水準である場合にだけ採用する（誤った前方一致を避ける）
    for (v in cand) {
      if (startsWith(s, v) && nchar(s) > nchar(v)) {
        lev   <- substring(s, nchar(v) + 1L)
        known <- c(names(LEVELS_BY_VAR[[v]]), names(LEVEL_LABELS))
        if (lev %in% known) {
          return(paste0(unname(VAR_LABELS[v]), ": ",
                        relabel_levels(v, lev, strict = FALSE)))
        }
      }
    }

    ## (5) 接尾辞
    if (grepl("_measurable$", s)) {
      b <- sub("_measurable$", "", s)
      if (b %in% names(VAR_LABELS)) {
        return(paste0(unname(VAR_LABELS[b]), " (measurable)"))
      }
    }
    if (grepl("_c$", s)) {
      b <- sub("_c$", "", s)
      if (b %in% names(VAR_LABELS)) return(unname(VAR_LABELS[b]))
    }

    unknown <<- c(unknown, s)
    s
  }

  out <- vapply(x, one, character(1), USE.NAMES = FALSE)
  if (strict) .label_warn("変数名", unknown)
  out
}

#' 選抜の段階を英語ラベルに変換する
relabel_stage <- function(x) {
  x <- as.character(x)
  if (!length(x)) return(character(0))
  hit <- x %in% names(STAGE_LABELS)
  out <- x
  out[hit] <- unname(STAGE_LABELS[x[hit]])
  out
}

#' パイプライン名を英語ラベルに変換する
relabel_pipeline <- function(label) {
  label <- as.character(label)
  if (!length(label)) return(character(0))
  hit <- label %in% names(PIPE_LABELS)
  out <- label
  out[hit] <- unname(PIPE_LABELS[label[hit]])
  out
}


# -----------------------------------------------------------------------------
# 5. 被覆の点検
#    Alldata の各列が「ラベル済み」「意図して使わない」「未定義」のどれかを示す。
#    未定義の列があれば、ラベルの付け忘れか、計画に無い列である。
# -----------------------------------------------------------------------------

check_alldata_labels <- function(data, quiet = FALSE) {
  nm <- names(data)

  status <- ifelse(nm %in% names(VAR_LABELS), "labelled",
            ifelse(nm %in% names(VARS_EXCLUDED), "not used", "UNDEFINED"))
  note <- ifelse(status == "labelled",  unname(VAR_LABELS[nm]),
          ifelse(status == "not used",  unname(VARS_EXCLUDED[nm]), ""))
  var_report <- data.frame(variable = nm, status = status, note = note,
                           stringsAsFactors = FALSE)

  ## 水準の被覆（LEVELS_BY_VAR に定義のある変数のみ）
  lev_rows <- list()
  for (v in intersect(nm, names(LEVELS_BY_VAR))) {
    obs <- unique(as.character(data[[v]]))
    obs <- obs[!is.na(obs)]
    und <- setdiff(obs, names(LEVELS_BY_VAR[[v]]))
    lev_rows[[v]] <- data.frame(
      variable  = v,
      n_levels  = length(obs),
      undefined = if (length(und)) paste(und, collapse = " | ") else "",
      stringsAsFactors = FALSE
    )
  }
  lev_report <- if (length(lev_rows)) do.call(rbind, lev_rows) else
    data.frame(variable = character(0), n_levels = integer(0),
               undefined = character(0), stringsAsFactors = FALSE)

  ## 計画で使う変数のうち、データに無いもの（派生変数は 02 以降で作るため除く）
  derived <- c("age_group", "period", "fy_in", "los", "mFIM_in_r", "cFIM_in_r")
  missing_vars <- setdiff(setdiff(names(VAR_LABELS), derived), nm)

  if (!quiet) {
    cat("\n-- check_alldata_labels ------------------------------------------\n")
    print(var_report, row.names = FALSE)
    cat("\n-- levels --------------------------------------------------------\n")
    print(lev_report, row.names = FALSE)
    cat("\nUNDEFINED columns : ",
        if (any(status == "UNDEFINED")) paste(nm[status == "UNDEFINED"], collapse = ", ")
        else "(none)", "\n", sep = "")
    cat("Levels not labelled: ",
        if (any(nzchar(lev_report$undefined)))
          paste(lev_report$variable[nzchar(lev_report$undefined)], collapse = ", ")
        else "(none)", "\n", sep = "")
    cat("Planned but absent : ",
        if (length(missing_vars)) paste(missing_vars, collapse = ", ") else "(none)",
        "\n", sep = "")
    cat("------------------------------------------------------------------\n")
  }

  invisible(list(variables = var_report, levels = lev_report,
                 missing = missing_vars))
}
