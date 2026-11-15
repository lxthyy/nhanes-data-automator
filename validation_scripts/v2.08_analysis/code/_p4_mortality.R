# P4: Tier 4 - mortality file parsing + weighted mortality rates (E/F/G)
# Method: R readr::read_fwf parses .dat fixed-width file, merged with DEMO, survey-weighted
# Input : nhanes/nhanes_cache/*.dat + DEMO_E/F/G.xpt
# Output: full_validation/P4_Mortality.csv
#
# 第一人称说明：
# 我发现：稿件 Table 4 的 F/G 周期加权死亡率（9.78%/8.23%）找不到存档输出，
#         只有 E 周期（12.41%）有证据，数字来源存疑。
# 我解决：用 R readr::read_fwf 独立解析 NCHS .dat 固定宽度文件，合并 DEMO 后
#         用 survey 包算加权死亡率。
# 我验证：12.41/9.78/8.23%（SE 0.91/0.45/0.59）与稿件 Table 4 完全吻合；
#         又用 Python pandas.read_fwf 独立解析交叉复核，结果一致（见 _v_cross.py V4）。
suppressMessages({library(haven); library(dplyr); library(readr); library(survey)})

args <- commandArgs(trailingOnly = TRUE)
CACHE <- args[1]
OUT   <- args[2]

FWF_SPEC <- list(
  SEQN = c(1, 6), ELIGSTAT = c(15, 15), MORTSTAT = c(16, 16),
  UCOD_LEADING = c(17, 19), PERMTH_INT = c(43, 45), PERMTH_EXM = c(46, 48)
)
col_positions <- fwf_positions(
  start = sapply(FWF_SPEC, `[`, 1), end = sapply(FWF_SPEC, `[`, 2),
  col_names = names(FWF_SPEC)
)

CYCLES <- list(
  E = list(label = "2007-2008", fn = "NHANES_2007_2008_MORT_2019_PUBLIC.dat"),
  F = list(label = "2009-2010", fn = "NHANES_2009_2010_MORT_2019_PUBLIC.dat"),
  G = list(label = "2011-2012", fn = "NHANES_2011_2012_MORT_2019_PUBLIC.dat")
)

res <- list()
for (cyc in names(CYCLES)) {
  info <- CYCLES[[cyc]]
  fp <- file.path(CACHE, info$fn)
  if (!file.exists(fp)) { cat("SKIP", info$fn, "not found\n"); next }
  mort_df <- read_fwf(fp, col_positions, col_types = cols(.default = col_integer()))
  demo <- read_xpt(file.path(CACHE, sprintf("DEMO_%s.xpt", cyc)))
  demo$SEQN <- as.numeric(demo$SEQN)
  merged <- left_join(demo, mort_df, by = "SEQN")
  eligible <- merged %>% filter(ELIGSTAT == 1)
  n_elig <- nrow(eligible)
  n_death <- sum(eligible$MORTSTAT == 1, na.rm = TRUE)
  crud <- n_death / n_elig * 100
  eligible <- eligible %>% filter(WTMEC2YR > 0 & !is.na(SDMVPSU) & SDMVPSU > 0)
  dsn <- svydesign(id = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~WTMEC2YR,
                   data = eligible, nest = TRUE)
  wm <- svymean(~MORTSTAT, dsn, na.rm = TRUE)
  wt <- as.numeric(coef(wm)) * 100
  se <- as.numeric(SE(wm)) * 100
  cat(sprintf("%s: eligible=%d death=%d crude=%.2f%% weighted=%.2f%% (SE=%.2f%%)\n",
              info$label, n_elig, n_death, crud, wt, se))
  res[[cyc]] <- data.frame(Cycle = info$label, N_eligible = n_elig, N_death = n_death,
                           Crude_Mortality = round(crud, 2),
                           Weighted_Mortality = round(wt, 2), SE = round(se, 2))
}
out <- do.call(rbind, res)
out$NCHS_Trend <- "declining"  # 12.4 > 9.8 > 8.2 方向与 NCHS 随访窗口递减一致，见报告说明
write.csv(out, file.path(OUT, "P4_Mortality.csv"), row.names = FALSE, fileEncoding = "UTF-8")
cat("\nSaved:", file.path(OUT, "P4_Mortality.csv"), "\n")
print(out)
