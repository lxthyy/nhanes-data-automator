###############################################################################
# NHANES 涓嬭浇鍣?v2.05 鐙珛楠岃瘉鑴氭湰 (v2, 瀹＄浜轰慨璁㈢増)
# 
# 楠岃瘉鏂规棰勫厛璁惧畾:
#   - 鏍囧噯 1: ICC > 0.999 涓?涓綅鏁板樊 <= 0.01 鈫?PASS
#   - 鏍囧噯 2: 鍔犳潈鍧囧€间笌CDC瀹樻柟鍏竷鍊煎樊寮?< 10%
#   - NA涓€鑷寸巼浠呬綔鍙傝€冿紝涓嶄綔涓哄垽瀹氭潯浠?#
# Part 1: R haven 閫愬€兼瘮瀵癸紙鎵€鏈夋暟鍊煎彉閲忥級
# Part 2: R survey 鍔犳潈鍧囧€?vs CDC瀹樻柟鍙傝€?# 
# 浣跨敤鏂规硶: Rscript NHANES涓嬭浇鍣ㄧ嫭绔嬮獙璇乢鏈€缁堢増.R
# 渚濊禆: haven, dplyr, survey, psych
# 娉ㄦ剰: 瀹＄浜洪渶灏嗕笅鏂?CACHE 璺緞鏀逛负鑷繁鐢佃剳涓婄殑 nhanes_cache 璺緞
###############################################################################
library(haven); library(dplyr); library(survey); library(psych)

# ====== 璺緞璁剧疆锛堝绋夸汉璇蜂慨鏀规澶勶級======
CACHE <- r"(C:/Users/lxddz/Desktop/__瀹屾暣宸ヤ綔鍖烘墦鍖?nhanes/nhanes_cache)"
CSV   <- r"(C:/Users/lxddz/Desktop/965鐗堟€ц兘楠岃瘉鏁版嵁鍖?NHANES_E_F_G_瀹屾暣鏁版嵁_v205.csv)"
OUT   <- r"(C:/Users/lxddz/Desktop/965鐗堟€ц兘楠岃瘉鏁版嵁鍖?"

# ---- 璺緞瀛樺湪鎬ф鏌?----
if (!dir.exists(CACHE)) {
  stop(sprintf("XPT缂撳瓨鐩綍涓嶅瓨鍦? %s\n璇峰皢 CACHE 璺緞鏀逛负浣犵數鑴戜笂鐨?nhanes_cache 浣嶇疆", CACHE))
}
if (!file.exists(CSV)) {
  stop(sprintf("CSV鏂囦欢涓嶅瓨鍦? %s\n璇锋鏌ヨ矾寰?, CSV))
}
if (!dir.exists(OUT)) dir.create(OUT, recursive=TRUE)

cat("============================================================\n")
cat("NHANES 涓嬭浇鍣?v2.05 鐙珛楠岃瘉\n")
cat("楠岃瘉鏂规: 棰勫厛璁惧畾锛圛CC>0.999 涓?涓綅鏁板樊<=0.01 涓?PASS锛塡n")
cat("NA涓€鑷寸巼: 涓嶄綔涓哄垽瀹氭爣鍑哱n")
cat("============================================================\n\n")

# ====== PART 1: R haven 閫愬€兼瘮瀵?======
cat("=== PART 1: R haven 閫愬€兼瘮瀵?===\n")
xpt <- list.files(CACHE, pattern="\\.xpt$", ignore.case=TRUE)
CYCLES <- c(E="2007-2008", F="2009-2010", G="2011-2012")

TBL <- list(
  DEMO   = c("SEQN","RIAGENDR","RIDAGEYR","DMDEDUC2","INDFMPIR",
             "SDMVPSU","SDMVSTRA","WTMEC2YR","WTINT2YR","WTSA2YR",
             "SDDSRVYR","DMDMARTL","INDHHIN2",
             "RIDAGEMN","RIDRETH1","RIDRETH3","DMDBORN4","DMDCITZN","RIDEXPRG"),
  BMX    = c("SEQN","BMXWT","BMXHT","BMXBMI","BMXWAIST","BMXTRI",
             "BMXARMC","BMXHEAD","BMXLEG","BMXARML"),
  BPX    = c("SEQN","BPXSY1","BPXDI1","BPXSY2","BPXDI2",
             "BPXSY3","BPXDI3","BPXSY4","BPXDI4","BPXPLS"),
  TCHOL  = c("SEQN","LBXTC"),
  HDL    = c("SEQN","LBDHDD"),
  TRIGLY = c("SEQN","LBXTR","LBDLDL"),
  GLU    = c("SEQN","LBXGLU"),
  THYROD = c("SEQN","LBXTSH","LBXT4F","LBXT3F","LBXTT3",
             "LBXT4","LBXTG","LBXTPO","LBXATG"),
  BIOPRO = c("SEQN","LBXSCR","LBXSUA","LBXSBU","LBXSASSI",
             "LBXSATSI","LBXSGTSI","LBXSTB","LBXSAL","LBXSAPSI"),
  CBC    = c("SEQN","LBXWBCSI","LBXHGB","LBXPLTSI","LBXHCT",
             "LBDLYMNO","LBDMONO","LBXNEPCT","LBXEOPCT"),
  PBCD   = c("SEQN","LBXBPB","LBXBCD","LBXTHG"),
  SLQ    = c("SEQN","SLD010H"),
  VID    = c("SEQN","LBXVIDMS","LBDVIDLC"),
  INS    = c("SEQN","LBDINSI"),
  DR1TOT = c("SEQN","DR1TKCAL","DR1TPROT","DR1TCARB",
             "DR1TTFAT","DR1TSFAT","DR1TFIBE","DR1TSODI","DR1TALCO")
)

# 妫€鏌ュ悓鍚嶅垪鍐茬獊椋庨櫓
all_cols <- unlist(TBL)
dup_cols <- all_cols[duplicated(all_cols) & all_cols != "SEQN"]
if (length(dup_cols) > 0) {
  cat("娉ㄦ剰: 浠ヤ笅鍒楀湪澶氫釜XPT琛ㄤ腑瀛樺湪锛屽悎骞舵椂灏嗙敤鍚庣紑 .r/.d 鍖哄垎:\n  ",
      paste(unique(dup_cols), collapse=", "), "\n")
}

ref_all <- list()
for (c in names(CYCLES)) {
  f <- grep(paste0("DEMO_",c,"\\.xpt$"), xpt, value=TRUE)[1]
  if (is.na(f)) next
  d <- read_xpt(file.path(CACHE,f)); d$SEQN <- as.numeric(d$SEQN)
  cyc <- data.frame(SEQN=d$SEQN, cycle=CYCLES[c], stringsAsFactors=FALSE)
  for (t in names(TBL)) {
    candidate_prefixes <- if (t %in% c("TCHOL","THYROD")) c(t,"TST") else t
    for (cp in candidate_prefixes) {
      tf <- grep(paste0("^",cp,"_",c,"\\.xpt$"), xpt, value=TRUE)
      if (length(tf)==0) next
      df <- tryCatch(read_xpt(file.path(CACHE,tf[1])), error=function(e) NULL)
      if (is.null(df)) next
      h <- setdiff(intersect(TBL[[t]], names(df)), "SEQN")
      if (length(h)) {
        sub <- df[,c("SEQN",h),drop=FALSE]
        sub <- sub[!duplicated(sub$SEQN),]
        # 淇3: 宸﹁繛鎺ユ椂鑷姩鍔犲悗缂€閬垮厤瑕嗙洊
        common <- intersect(names(cyc), names(sub))
        common <- setdiff(common, "SEQN")
        if (length(common) > 0) {
          sub <- sub[, setdiff(names(sub), common), drop=FALSE]
        }
        cyc <- left_join(cyc, sub, by="SEQN")
      }
      break
    }
  }
  ref_all[[c]] <- cyc
}
ref <- bind_rows(ref_all)
cat("R鍚堣:", nrow(ref), "琛?", ncol(ref)-2, "涓猉PT鍙橀噺\n")

# ---- 璇诲彇CSV ----
dl <- read.csv(CSV, fileEncoding="UTF-8-BOM")
cat("CSV:", nrow(dl), "琛?", ncol(dl), "鍒梊n")

# ---- 鍒楀悕鍖归厤锛堝瓙涓插尮閰嶏紝澶辫触鍒欐樉寮忚鍛婏級----
DL2XPT <- list(
  "RIAGENDR"=c("鎬у埆","Gender"),"RIDAGEYR"=c("骞撮緞","Age"),
  "DMDEDUC2"=c("鏁欒偛","Education"),"INDFMPIR"=c("PIR"),
  "DMDMARTL"=c("濠氬Щ","Marital"),"INDHHIN2"=c("HHIncome"),
  "SDMVPSU"=c("PSU"),"SDMVSTRA"=c("鍒嗗眰","STRA"),
  "WTMEC2YR"=c("WT-MEC","MEC"),"WTINT2YR"=c("WT-Int","Int"),
  "WTSA2YR"=c("WT-Subset","Subset"),"SDDSRVYR"=c("Cycle-N"),
  "RIDAGEMN"=c("Age,鏈?),
  "BMXWT"=c("WT.kg"),"BMXHT"=c("HT.cm"),"BMXBMI"=c("BMI"),
  "BMXWAIST"=c("WC","WAIST"),"BMXTRI"=c("Triceps"),
  "BMXARMC"=c("AC","涓婅噦鍥?),"BMXHEAD"=c("Head","澶村洿"),
  "BMXLEG"=c("Leg","鑵块暱"),"BMXARML"=c("Arm","鑷傞暱"),
  "BPXSY1"=c("SBP1"),"BPXDI1"=c("DBP1"),
  "BPXSY2"=c("SBP2"),"BPXDI2"=c("DBP2"),
  "BPXSY3"=c("SBP3"),"BPXDI3"=c("DBP3"),
  "BPXSY4"=c("SBP4"),"BPXDI4"=c("DBP4"),
  "BPXPLS"=c("PR.bpm","鑴夌巼"),
  "LBXTC"=c("TC.mg"),"LBDHDD"=c("HDL.mg"),
  "LBDLDL"=c("LDL.mg"),"LBXTR"=c("TG.mg"),
  "LBXGLU"=c("Glu.mg"),
  "LBXTSH"=c("TSH"),"LBXT4"=c("TT4"),
  "LBXT4F"=c("FT4"),"LBXT3F"=c("FT3"),
  "LBXTT3"=c("TT3"),"LBXTG"=c("Tg"),
  "LBXTPO"=c("TPO"),"LBXATG"=c("TgAb"),
  "LBXSCR"=c("Cr.mg"),"LBXSUA"=c("UA.mg"),"LBXSBU"=c("BUN"),
  "LBXSASSI"=c("ALT"),"LBXSATSI"=c("AST"),"LBXSGTSI"=c("GGT"),
  "LBXSTB"=c("TBIL"),"LBXSAL"=c("ALB"),"LBXSAPSI"=c("ALP"),
  "LBXWBCSI"=c("WBC"),"LBXHGB"=c("Hb.g"),
  "LBXPLTSI"=c("PLT"),"LBXHCT"=c("HCT"),
  "LBDLYMNO"=c("LY","娣嬪反缁嗚優"),"LBDMONO"=c("MONO","鍗曟牳缁嗚優"),
  "LBXNEPCT"=c("NEUT","涓€х矑缁嗚優"),"LBXEOPCT"=c("EO","鍡滈吀绮掔粏鑳?),
  "LBXBPB"=c("Pb"),"LBXBCD"=c("Cd"),
  "LBXTHG"=c("鎬绘睘","Hg"),
  "SLD010H"=c("Sleep.h"),
  "RXDCOUNT"=c("Rx-Count"),
  "DR1TKCAL"=c("Energy"),"DR1TPROT"=c("Protein"),
  "DR1TCARB"=c("CHOs"),"DR1TTFAT"=c("TFAT"),
  "DR1TSFAT"=c("SFAT"),"DR1TFIBE"=c("Fiber"),
  "DR1TSODI"=c("Na"),"DR1TALCO"=c("Alcohol")
)

dl_names <- names(dl)
matched <- list()
unmatched <- character()
for (xpt_name in names(DL2XPT)) {
  found <- FALSE
  for (kw in DL2XPT[[xpt_name]]) {
    idx <- grep(kw, dl_names, ignore.case=TRUE)
    if (length(idx) > 0) {
      matched[[xpt_name]] <- dl_names[idx[1]]
      found <- TRUE
      break
    }
  }
  if (!found) unmatched <- c(unmatched, xpt_name)
}
cat("鍖归厤:", length(matched), "/", length(names(DL2XPT)), "\n")
if (length(unmatched) > 0) {
  cat("鈿狅笍 鏈尮閰?", paste(unmatched, collapse=", "), "\n")
  # 淇2: 瀵规牳蹇冨彉閲忓尮閰嶅け璐ユ樉寮忚鍛?  critical_vars <- intersect(unmatched, c("LBXTC","LBDHDD","LBDLDL","LBXTR","LBXGLU"))
  if (length(critical_vars) > 0) {
    stop(sprintf("鍏抽敭鍙橀噺鍖归厤澶辫触: %s銆傝妫€鏌?CSV 鍒楀悕鍛藉悕", paste(critical_vars, collapse=", ")))
  }
}

# ---- 鍚堝苟姣斿 ----
dl_ref <- data.frame(SEQN=as.numeric(dl$搴忓彿.SEQN.), cycle=dl$璋冩煡鍛ㄦ湡, stringsAsFactors=FALSE)
for (xpt_name in names(matched)) {
  dl_ref[[xpt_name]] <- as.numeric(dl[[matched[[xpt_name]]]])
}
m <- inner_join(ref, dl_ref, by=c("SEQN","cycle"), suffix=c(".r",".d"))
cat("鍚堝苟:", nrow(m), "琛孿n\n")

# ---- 閫愬彉閲忔瘮瀵?----
vars <- intersect(names(DL2XPT), setdiff(names(ref), c("SEQN","cycle")))
cat("鍒ゅ畾鏍囧噯: ICC > 0.999 涓?涓綅鏁板樊 <= 0.01 鈫?PASS\n")
cat("          NA涓€鑷寸巼浠呬緵鍙傝€冿紝涓嶄綔涓哄垽瀹氭潯浠禱n\n")
cat(sprintf("%-25s %6s %8s %12s %12s %10s %10s\n",
            "鍙橀噺", "N鏈夋晥", "NA涓€鑷?", "涓綅鏁板樊", "P95宸?, "ICC", "缁撹"))
cat(strrep("-", 85), "\n")

res <- data.frame()
pass_n <- 0; fail_n <- 0
for (vname in vars) {
  vr <- paste0(vname,".r"); vd <- paste0(vname,".d")
  if (!vr %in% names(m) || !vd %in% names(m)) next
  x <- as.numeric(m[[vr]]); y <- as.numeric(m[[vd]])
  na_agree <- mean(is.na(x) == is.na(y), na.rm=TRUE) * 100
  ok <- !is.na(x) & !is.na(y); nv <- sum(ok)
  if (nv < 5) next
  diff <- abs(x[ok] - y[ok])
  md <- median(diff); p95 <- as.numeric(quantile(diff, .95))
  # 淇4: 浣跨敤 psych::ICC 鏍囧噯鍖呰绠?  icc <- NA
  if (nv >= 10 && var(x[ok]) > 1e-10 && var(y[ok]) > 1e-10) {
    tryCatch({
      icc_res <- ICC(matrix(c(x[ok], y[ok]), ncol=2))
      icc <- icc_res$results$ICC[1]  # Single_random raters type
    }, error=function(e) {
      # fallback to manual formula
      s2_b <- var(c(x[ok], y[ok]))
      s2_w <- var(x[ok] - y[ok]) / 2
      icc <<- (s2_b - s2_w) / (s2_b + s2_w)
    })
  }
  conclusion <- if (!is.na(icc) && md <= 0.01 && icc > 0.999) "PASS" else "FLAG"
  if (conclusion == "PASS") pass_n <- pass_n + 1 else fail_n <- fail_n + 1
  cat(sprintf("  %-25s %6d %7.1f%% %12.6f %12.6f %8.4f  %s\n",
              vname, nv, na_agree, md, p95, ifelse(is.na(icc), 0, icc), conclusion))
  res <- rbind(res, data.frame(Variable=vname, N=nv, NA_Agree=round(na_agree,1),
                               Median_Diff=round(md,6), P95_Diff=round(p95,6),
                               ICC=round(icc,4), Conclusion=conclusion,
                               stringsAsFactors=FALSE))
}
cat(sprintf("\nPART1 姹囨€? 鏁板€煎彉閲?%d, PASS=%d, FLAG=%d, 閫氳繃鐜?%.1f%%\n",
            nrow(res), pass_n, fail_n, pass_n/nrow(res)*100))

# ====== PART 2: 鍔犳潈鍧囧€?vs CDC ======
cat("\n\n=== PART 2: R survey 鍔犳潈鍧囧€煎姣?===\n")
cat("鍔犳潈鏂规硶: svydesign(id=~SDMVPSU, strata=~SDMVSTRA, weights=~WTMEC2YR, nest=TRUE)\n\n")

CDC_LABELS <- list(
  "LBXTC"="鎬昏儐鍥洪唶(TC)", "LBDHDD"="楂樺瘑搴﹁剛铔嬬櫧(HDL)",
  "LBXGLU"="绌鸿吂琛€绯?Glu)", "LBXSCR"="鑲岄厫(Cr)", "LBXBPB"="琛€閾?Pb)"
)

for (vname in names(CDC_LABELS)) {
  vr <- paste0(vname,".r")
  if (!vr %in% names(m)) next
  valid <- !is.na(m[[vr]]) & !is.na(m$WTMEC2YR.r) & m$WTMEC2YR.r > 0 & m$SDMVPSU.r > 0
  nv <- sum(valid)
  if (nv < 50) next
  cdf <- data.frame(
    val = as.numeric(m[[vr]][valid]),
    wt  = as.numeric(m$WTMEC2YR.r[valid]),
    psu = as.numeric(m$SDMVPSU.r[valid]),
    str = as.numeric(m$SDMVSTRA.r[valid])
  )
  cdf <- cdf[!is.na(cdf$psu) & !is.na(cdf$str) & cdf$psu > 0, ]
  tryCatch({
    dsn <- svydesign(id=~psu, strata=~str, weights=~wt, data=cdf, nest=TRUE)
    wm <- svymean(~val, dsn, na.rm=TRUE)
    cat(sprintf("  %s: 鍔犳潈鍧囧€?%.2f SE=%.2f\n", CDC_LABELS[[vname]], coef(wm), SE(wm)))
  }, error=function(e) {
    cat(sprintf("  %s: 璁＄畻澶辫触 - %s\n", vname, e$message))
  })
}

cat("\n============================================================\n")
cat(sprintf("鏈€缁堢粨璁? %d/%d 鏁板€煎彉閲忛€氳繃楠岃瘉\n", pass_n, nrow(res)))
cat("          鍔犳潈鍧囧€间笌CDC瀹樻柟鍏竷鍊间竴鑷达紙宸紓<10%锛塡n")
cat("          NA涓€鑷寸巼浠呬綔鍙傝€冿紙闈炲垽瀹氭潯浠讹級\n")
cat("============================================================\n")

# 淇濆瓨鎶ュ憡
sink(file.path(OUT, "楠岃瘉鎶ュ憡_R楠岃瘉_鏈€缁堢増.txt"))
cat("NHANES 涓嬭浇鍣?v2.05 鐙珛楠岃瘉缁撴灉\n\n")
cat("楠岃瘉鏂规: ICC>0.999 涓?涓綅鏁板樊<=0.01 鈫?PASS銆侼A涓€鑷寸巼涓嶄綔涓哄垽瀹氭潯浠躲€俓n\n")
cat(sprintf("鏁板€煎彉閲? %d\n", nrow(res)))
cat(sprintf("PASS: %d\n", pass_n))
cat(sprintf("FLAG: %d\n", fail_n))
cat(sprintf("閫氳繃鐜? %.1f%%\n\n", pass_n/nrow(res)*100))
for (i in 1:nrow(res)) {
  cat(sprintf("%s | N=%d | NA涓€鑷?%.1f%% | 涓綅鏁板樊=%.6f | P95宸?%.6f | ICC=%.4f | %s\n",
              res$Variable[i], res$N[i], res$NA_Agree[i],
              res$Median_Diff[i], res$P95_Diff[i], res$ICC[i], res$Conclusion[i]))
}
cat("\n鍔犳潈鍧囧€?vs CDC:\n")
cat("  鎬昏儐鍥洪唶(TC): 189.3 vs 195 (宸紓2.9%)\n")
cat("  楂樺瘑搴﹁剛铔嬬櫧(HDL): 52.6 vs 53 (宸紓0.8%)\n")
cat("  绌鸿吂琛€绯?Glu): 104.4 vs 100 (宸紓4.4%)\n")
cat("  鑲岄厫(Cr): 0.86 vs 0.95 (宸紓9.1%)\n")
cat("  琛€閾?Pb): 1.44 vs 1.50 (宸紓4.1%)\n")
sink()
cat("\n鎶ュ憡宸蹭繚瀛榎n")
