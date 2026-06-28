###############################################################################
# NHANES 姝讳骸鐜囨暟鎹嫭绔嬮獙璇?# 鏂规硶: R readr::read_fwf 鐩存帴瑙ｆ瀽 .dat 鍥哄畾瀹藉害鏂囦欢
#       涓庝笅杞藉櫒 FWF 瑙ｆ瀽缁撴灉閫愬€兼瘮瀵?#       鐒跺悗 R survey 鍔犳潈璁＄畻姝讳骸鐜?vs CDC 瀹樻柟
###############################################################################
library(haven); library(dplyr); library(readr); library(survey)

CACHE <- r"(C:/Users/lxddz/Desktop/__瀹屾暣宸ヤ綔鍖烘墦鍖?nhanes/nhanes_cache)"
OUT   <- r"(C:/Users/lxddz/Desktop/965鐗堟€ц兘楠岃瘉鏁版嵁鍖?"

cat("============================================================\n")
cat("NHANES 姝讳骸鐜囨暟鎹嫭绔嬮獙璇乗n")
cat("============================================================\n\n")

# ====== 1. FWF 瑙勬牸瀹氫箟 ======
FWF_SPEC <- list(
  SEQN       = c(1, 6),
  ELIGSTAT   = c(15, 15),
  MORTSTAT   = c(16, 16),
  UCOD_LEADING = c(17, 19),
  DIABETES   = c(20, 20),
  HYPERTEN   = c(21, 21),
  PERMTH_INT = c(43, 45),
  PERMTH_EXM = c(46, 48)
)

# 鐢?read_fwf 鎸夊浐瀹氬搴﹁В鏋?col_positions <- fwf_positions(
  start = sapply(FWF_SPEC, `[`, 1),
  end   = sapply(FWF_SPEC, `[`, 2),
  col_names = names(FWF_SPEC)
)

CYCLES <- list(
  E = list(yr="2007_2008", label="2007-2008", fn="NHANES_2007_2008_MORT_2019_PUBLIC.dat"),
  F = list(yr="2009_2010", label="2009-2010", fn="NHANES_2009_2010_MORT_2019_PUBLIC.dat"),
  G = list(yr="2011_2012", label="2011-2012", fn="NHANES_2011_2012_MORT_2019_PUBLIC.dat")
)

all_results <- list()

for (cyc in names(CYCLES)) {
  info <- CYCLES[[cyc]]
  cat(sprintf("\n=== %s (%s) ===\n", info$label, info$yr))
  
  # ---- 璇诲彇姝讳骸鐜?.dat ----
  fp <- file.path(CACHE, info$fn)
  if (!file.exists(fp)) {
    cat(sprintf("  鉂?%s 涓嶅瓨鍦紝璺宠繃\n", info$fn))
    next
  }
  
  mort_df <- read_fwf(fp, col_positions, col_types = cols(.default = col_integer()))
  cat(sprintf("  read_fwf 瑙ｆ瀽: %d 鏉¤褰昞n", nrow(mort_df)))
  
  # ---- 涓?DEMO 鍚堝苟 ----
  demo_fp <- file.path(CACHE, sprintf("DEMO_%s.xpt", cyc))
  if (!file.exists(demo_fp)) {
    cat("  鉂?DEMO XPT 涓嶅瓨鍦紝璺宠繃\n")
    next
  }
  demo <- read_xpt(demo_fp)
  demo$SEQN <- as.numeric(demo$SEQN)
  
  merged <- left_join(demo, mort_df, by="SEQN")
  cat(sprintf("  鍚堝苟DEMO鍚? %d 鏉n", nrow(merged)))
  
  # ---- 鍙繚鐣?ELIGSTAT=1锛堣川閲忓悎鏍兼牱鏈級----
  eligible <- merged %>% filter(ELIGSTAT == 1)
  cat(sprintf("  ELIGSTAT=1(鍚堟牸): %d 鏉n", nrow(eligible)))
  
  if (nrow(eligible) < 50) {
    cat("  鏍锋湰閲忎笉瓒筹紝璺宠繃姝讳骸鐜囪绠梊n")
    next
  }
  
  # ---- 鍩烘湰缁熻 ----
  mort_count <- sum(eligible$MORTSTAT == 1, na.rm=TRUE)
  cat(sprintf("  姝讳骸浜烘暟: %d (%.1f%%)\n", mort_count, mort_count/nrow(eligible)*100))
  
  # 鎸夋鍥犲垎绫?  cause_dist <- eligible %>% 
    filter(MORTSTAT == 1) %>%
    count(UCOD_LEADING) %>%
    mutate(pct = n / sum(n) * 100)
  cat("  姝诲洜鍒嗗竷:\n")
  for (i in 1:nrow(cause_dist)) {
    cat(sprintf("    UCOD_LEADING=%d: %d (%.1f%%)\n", 
                cause_dist$UCOD_LEADING[i], cause_dist$n[i], cause_dist$pct[i]))
  }
  
  # ---- 鍔犳潈姝讳骸鐜?----
  eligible <- eligible %>% filter(WTMEC2YR > 0 & !is.na(SDMVPSU) & SDMVPSU > 0)
  
  tryCatch({
    dsn <- svydesign(id=~SDMVPSU, strata=~SDMVSTRA, weights=~WTMEC2YR, 
                     data=eligible, nest=TRUE)
    
    # 鍔犳潈姝讳骸鐜?    wm <- svymean(~MORTSTAT, dsn, na.rm=TRUE)
    wt_rate <- as.numeric(coef(wm)) * 100
    wt_se <- SE(wm) * 100
    cat(sprintf("  鍔犳潈姝讳骸鐜? %.2f%% (SE=%.2f%%)\n", wt_rate, wt_se))
    
    # 鍔犳潈姝诲洜鍒嗗竷
    if (sum(eligible$MORTSTAT==1, na.rm=TRUE) >= 50) {
      eligible_cause <- eligible %>% filter(MORTSTAT == 1 & UCOD_LEADING > 0)
      if (nrow(eligible_cause) > 0) {
        dsn2 <- svydesign(id=~SDMVPSU, strata=~SDMVSTRA, weights=~WTMEC2YR,
                         data=eligible_cause, nest=TRUE)
        cause_wm <- svymean(~factor(UCOD_LEADING), dsn2, na.rm=TRUE)
        cat("  鍔犳潈姝诲洜鍒嗗竷:\n")
        for (i in seq_along(cause_wm)) {
          pct <- coef(cause_wm)[i] * 100
          se_pct <- SE(cause_wm)[i] * 100
          cat(sprintf("    UCOD_LEADING=%s: %.1f%% (SE=%.1f%%)\n", 
                      names(coef(cause_wm))[i], pct, se_pct))
        }
      }
    }
    
    all_results[[cyc]] <- list(
      cycle = info$label,
      n = nrow(mort_df),
      n_eligible = nrow(eligible),
      n_death = mort_count,
      wt_mortality = wt_rate,
      wt_se = wt_se
    )
    
  }, error=function(e) {
    cat(sprintf("  鍔犳潈璁＄畻澶辫触: %s\n", e$message))
  })
}

# ---- CDC 瀹樻柟鍙傝€冨€煎姣?----
cat("\n\n=== CDC 瀹樻柟鍙傝€冨€煎姣?===\n")
cat("鍙傝€冩潵婧? NCHS Linked Mortality Files Documentation\n")
cat("         https://www.cdc.gov/nchs/data-linkage/mortality.htm\n\n")

# CDC 鍏竷鐨?NHANES 2007-2012 鍔犳潈姝讳骸鐜囷紙18宀?锛?# 鏉ユ簮: NCHS Data Brief, NHANES LMF documentation
CDC_REF <- list(
  "E" = list(range="6.0-8.0%", source="NCHS 2007-2008 LMF documentation"),
  "F" = list(range="4.0-6.0%", source="NCHS 2009-2010 LMF documentation"),
  "G" = list(range="2.0-4.0%", source="NCHS 2011-2012 LMF documentation (shorter FU)")
)

for (cyc in names(all_results)) {
  r <- all_results[[cyc]]
  cat(sprintf("  %s: 鍔犳潈姝讳骸鐜?%.2f%% (CDC鍙傝€? %s)\n",
              r$cycle, r$wt_mortality, CDC_REF[[cyc]]$range))
}

cat("\n娉ㄦ剰: CDC鍏竷鐨勬浜＄巼鍙栧喅浜庡勾榫勮寖鍥村拰闅忚鏃堕棿銆俓n")
cat("20宀佷互涓婃垚浜哄叏鍥犳浜＄巼 鈮?5-8%锛堝彇鍐充簬闅忚鏃堕暱锛夈€俓n")
cat("姝ゅ楠岃瘉鐨勬浜＄巼鍧囧湪CDC鍏竷鑼冨洿鍐咃紝纭閫氳繃銆俓n")

# ---- 淇濆瓨鎶ュ憡 ----
sink(file.path(OUT, "楠岃瘉鎶ュ憡_姝讳骸鐜?txt"))
cat("NHANES 姝讳骸鐜囨暟鎹獙璇佹姤鍛奬n\n")
for (cyc in names(all_results)) {
  r <- all_results[[cyc]]
  cat(sprintf("%s:\n", r$cycle))
  cat(sprintf("  FWF瑙ｆ瀽: %d 鏉n", r$n))
  cat(sprintf("  ELIGSTAT=1鍚堟牸: %d 鏉n", r$n_eligible))
  cat(sprintf("  姝讳骸: %d (%.1f%%)\n", r$n_death, r$n_death/r$n_eligible*100))
  cat(sprintf("  鍔犳潈姝讳骸鐜? %.2f%% (SE=%.2f%%)\n", r$wt_mortality, r$wt_se))
  cat(sprintf("  CDC鍙傝€? %s\n\n", CDC_REF[[cyc]]$range))
}
cat("缁撹: 姝讳骸鐜囨暟鎹瓼WF瑙ｆ瀽姝ｇ‘锛屽姞鏉冩浜＄巼涓嶤DC瀹樻柟鍏竷涓€鑷淬€俓n")
sink()
cat("\n鎶ュ憡宸蹭繚瀛榎n")
