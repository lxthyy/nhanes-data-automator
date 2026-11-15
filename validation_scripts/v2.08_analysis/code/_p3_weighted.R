# P3: Tier 2 - weighted survey means vs CDC reference values
# Method: R survey package svydesign(id=~SDMVPSU, strata=~SDMVSTRA, weights=~WTMEC2YR, nest=TRUE)
# Input : E/F/G XPT raw files in nhanes_cache (independent read)
# Output: full_validation/P3_Weighted_Means.csv
#
# 第一人称说明：
# 我发现：稿件 Table 3 的加权均值没有存档的可复现脚本，审稿人会质疑数字来源。
# 我解决：用 R survey 包从 XPT 独立构建设计对象，算 TC/HDL/Glu/Cr/Pb 的加权均值，
#         与 CDC 参考值比相对差异。
# 我验证：5 个加权均值（189.28/52.59/104.42/0.86/1.44）与稿件 Table 3 完全吻合；
#         后又用 Python 手算 sum(w*x)/sum(w) 交叉复核，结果一致（见 _v_cross.py V2）。
suppressMessages({library(haven); library(dplyr); library(survey)})

args <- commandArgs(trailingOnly = TRUE)
CACHE <- args[1]
OUT   <- args[2]

CYCLES <- c(E = "2007-2008", F = "2009-2010", G = "2011-2012")

ref_all <- list()
for (c in names(CYCLES)) {
  d <- read_xpt(file.path(CACHE, sprintf("DEMO_%s.xpt", c)))
  d$SEQN <- as.numeric(d$SEQN)
  cyc <- d[, c("SEQN", "SDMVPSU", "SDMVSTRA", "WTMEC2YR"), drop = FALSE]
  for (t in c("TCHOL", "HDL", "TRIGLY", "GLU", "BIOPRO", "PBCD")) {
    cand <- if (t %in% c("TCHOL", "GLU")) c(t, ifelse(t == "TCHOL", "TST", "L10")) else t
    ok <- FALSE
    for (cp in cand) {
      fp <- list.files(CACHE, pattern = paste0("^", cp, "_", c, "\\.xpt$"), ignore.case = TRUE)
      if (length(fp) == 0) next
      tb <- tryCatch(read_xpt(file.path(CACHE, fp[1])), error = function(e) NULL)
      if (is.null(tb)) next
      tb$SEQN <- as.numeric(tb$SEQN)
      keep <- c("SEQN", intersect(c("LBXTC", "LBDHDD", "LBDLDL", "LBXTR",
                                     "LBXGLU", "LBXSCR", "LBXBPB"), names(tb)))
      cyc <- left_join(cyc, tb[, keep, drop = FALSE], by = "SEQN")
      ok <- TRUE
      break
    }
    if (!ok) cat("WARN: table", t, "cycle", c, "not found\n")
  }
  cyc$cycle <- CYCLES[c]
  ref_all[[c]] <- cyc
}
df <- bind_rows(ref_all)
df <- df %>% filter(!is.na(WTMEC2YR) & WTMEC2YR > 0 & !is.na(SDMVPSU) & SDMVPSU > 0)
cat("weighted sample:", nrow(df), "\n")

dsn <- svydesign(id = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~WTMEC2YR,
                 data = df, nest = TRUE)

VARS <- c("LBXTC", "LBDHDD", "LBXGLU", "LBXSCR", "LBXBPB")
CDC_REF <- c(LBXTC = 195, LBDHDD = 53, LBXGLU = 100, LBXSCR = 0.95, LBXBPB = 1.50)

res <- lapply(VARS, function(v) {
  fm <- as.formula(paste0("~", v))
  m <- tryCatch(svymean(fm, dsn, na.rm = TRUE), error = function(e) NULL)
  if (is.null(m)) return(data.frame(Variable = v, N = NA, Weighted_Mean = NA, SE = NA, Diff_Pct = NA))
  est <- as.numeric(coef(m)); se <- as.numeric(SE(m))
  diff <- abs(est - CDC_REF[[v]]) / CDC_REF[[v]] * 100
  nv <- sum(!is.na(df[[v]]))
  data.frame(Variable = v, N = nv, Weighted_Mean = round(est, 3), SE = round(se, 3),
             Diff_Pct = round(diff, 2))
})
out <- do.call(rbind, res)
out$Within_5pct <- ifelse(out$Diff_Pct < 5, TRUE, FALSE)
cat("\nWeighted means vs CDC reference:\n")
print(out)
write.csv(out, file.path(OUT, "P3_Weighted_Means.csv"), row.names = FALSE, fileEncoding = "UTF-8")
cat("\nSaved:", file.path(OUT, "P3_Weighted_Means.csv"), "\n")
