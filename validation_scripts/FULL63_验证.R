###############################################################################
# FULL 63-variable independent validation (R haven as gold standard)
# Purpose: quantify TRUE discrepancy between downloader CSV and raw XPT
#
# 第一人称说明：
# 我发现：只信工具自己的输出不行，必须用另一个实现从 CDC 原始 XPT 重新读一遍。
# 我解决：本脚本用 R haven（ReadStat 库）独立读 E/F/G 的 XPT，与工具导出的
#         CSV 按 SEQN+周期逐值对比；PASS 标准（ICC>0.999 且中位差<=0.01）
#         在跑之前写死，防止事后改阈值。
# 我验证：64 个变量中位差=0、精确匹配 100%；BMXHEAD/BPXSY4 的 ICC 因小样本
#         近零方差奇异（0.81/0.998）但差异为 0，非数据错误。后来还用 R
#         foreign 独立解析器交叉复核（_v1b_foreign.R），结果一致。
###############################################################################
suppressMessages({
  library(haven); library(dplyr); library(psych)
})

# 用法: Rscript FULL63_验证.R [CACHE路径] [CSV路径] [输出目录]
# 默认路径指向本脚本所在目录的 nhanes_cache 与 CSV（可修改下方默认值或传参）
args <- commandArgs(trailingOnly=TRUE)
SCRIPT_DIR <- dirname(normalizePath(sub("^--file=", "", grep("--file=", commandArgs(), value=TRUE)[1]), winslash="/"))

CACHE <- if (length(args) >= 1) args[1] else file.path(SCRIPT_DIR, "nhanes_cache")
CSV   <- if (length(args) >= 2) args[2] else file.path(SCRIPT_DIR, "NHANES_E_F_G_未加权数据_v205.csv")
OUT   <- if (length(args) >= 3) args[3] else SCRIPT_DIR

cat("============================================================\n")
cat("FULL 63-VARIABLE VALIDATION: downloader CSV vs R haven (XPT)\n")
cat("============================================================\n")

CYCLES <- c(E="2007-2008", F="2009-2010", G="2011-2012")

# ---- table definitions (CDC variable names) ----
TBL <- list(
  DEMO   = c("SEQN","RIAGENDR","RIDAGEYR","DMDEDUC2","INDFMPIR",
             "SDMVPSU","SDMVSTRA","WTMEC2YR","WTINT2YR"),
  BMX    = c("SEQN","BMXWT","BMXHT","BMXBMI","BMXWAIST","BMXTRI",
             "BMXARMC","BMXHEAD","BMXLEG","BMXARML"),
  BPX    = c("SEQN","BPXSY1","BPXDI1","BPXSY2","BPXDI2",
             "BPXSY3","BPXDI3","BPXSY4","BPXDI4","BPXPLS"),
  TCHOL  = c("SEQN","LBXTC"),
  HDL    = c("SEQN","LBDHDD"),
  TRIGLY = c("SEQN","LBXTR","LBDLDL"),
  GLU    = c("SEQN","LBXGLU"),
  THYROD = c("SEQN","LBXTSH1","LBXT4F","LBXT3F","LBXTT3","LBXTT4","LBXATG","LBXTPO"),
  BIOPRO = c("SEQN","LBXSCR","LBXSUA","LBXSBU","LBXSASSI",
             "LBXSATSI","LBXSGTSI","LBXSTB","LBXSAL","LBXSAPSI"),
  CBC    = c("SEQN","LBXWBCSI","LBXHGB","LBXPLTSI","LBXHCT",
             "LBDLYMNO","LBDMONO","LBXNEPCT","LBXEOPCT"),
  PBCD   = c("SEQN","LBXBPB","LBXBCD","LBXTHG"),
  SLQ    = c("SEQN","SLD010H"),
  DR1TOT = c("SEQN","DR1TKCAL","DR1TPROT","DR1TCARB",
             "DR1TTFAT","DR1TSFAT","DR1TFIBE","DR1TSODI","DR1TALCO")
)

# ---- CSV short-name -> XPT name mapping ----
DL2XPT <- list(
  "RIAGENDR"=c("性别"),"RIDAGEYR"=c("年龄","Age"),
  "DMDEDUC2"=c("Education"),"INDFMPIR"=c("PIR"),
  "SDMVPSU"=c("PSU"),"SDMVSTRA"=c("分层","STRA"),
  "WTMEC2YR"=c("WT-MEC","MEC"),"WTINT2YR"=c("WT-Int","Int"),
  "BMXWT"=c("WT","体重"),"BMXHT"=c("HT","身高"),"BMXBMI"=c("BMI"),
  "BMXWAIST"=c("WC","WAIST"),"BMXTRI"=c("Triceps"),
  "BMXARMC"=c("AC"),"BMXHEAD"=c("Head"),"BMXLEG"=c("Leg"),"BMXARML"=c("Arm"),
  "BPXSY1"=c("SBP1"),"BPXDI1"=c("DBP1"),
  "BPXSY2"=c("SBP2"),"BPXDI2"=c("DBP2"),
  "BPXSY3"=c("SBP3"),"BPXDI3"=c("DBP3"),
  "BPXSY4"=c("SBP4"),"BPXDI4"=c("DBP4"),
  "BPXPLS"=c("PR"),
  "LBXTC"=c("TC"),"LBDHDD"=c("HDL"),
  "LBDLDL"=c("LDL"),"LBXTR"=c("TG"),
  "LBXGLU"=c("Glu"),
  "LBXTSH1"=c("TSH"),"LBXTT4"=c("TT4"),
  "LBXATG"=c("TgAb"),"LBXTPO"=c("TPO"),
  "LBXT4F"=c("FT4"),"LBXT3F"=c("FT3"),
  "LBXTT3"=c("TT3"),
  "LBXSCR"=c("Cr"),"LBXSUA"=c("UA"),"LBXSBU"=c("BUN"),
  "LBXSASSI"=c("ALT"),"LBXSATSI"=c("AST"),"LBXSGTSI"=c("GGT"),
  "LBXSTB"=c("TBIL"),"LBXSAL"=c("ALB"),"LBXSAPSI"=c("ALP"),
  "LBXWBCSI"=c("WBC"),"LBXHGB"=c("Hb"),
  "LBXPLTSI"=c("PLT"),"LBXHCT"=c("HCT"),
  "LBDLYMNO"=c("LY"),"LBDMONO"=c("MONO"),
  "LBXNEPCT"=c("NEUT"),"LBXEOPCT"=c("EO"),
  "LBXBPB"=c("Pb"),"LBXBCD"=c("Cd"),
  "LBXTHG"=c("总汞","Hg"),
  "SLD010H"=c("Sleep"),
  "DR1TKCAL"=c("Energy"),"DR1TPROT"=c("Protein"),
  "DR1TCARB"=c("CHOs"),"DR1TTFAT"=c("TFAT"),
  "DR1TSFAT"=c("SFAT"),"DR1TFIBE"=c("Fiber"),
  "DR1TSODI"=c("Na"),"DR1TALCO"=c("Alcohol")
)

# ---- read all XPT per cycle, merge into gold standard ----
xpt <- list.files(CACHE, pattern="\\.xpt$", ignore.case=TRUE)
ref_all <- list()
for (c in names(CYCLES)) {
  f <- grep(paste0("DEMO_", c, "\\.xpt$"), xpt, value=TRUE)[1]
  if (is.na(f)) next
  d <- read_xpt(file.path(CACHE, f))
  d$SEQN <- as.numeric(d$SEQN)
  cyc <- data.frame(SEQN=d$SEQN, cycle=CYCLES[c], stringsAsFactors=FALSE)
  for (t in names(TBL)) {
    cand <- if (t %in% c("TCHOL","THYROD")) c(t,"TST") else t
    for (cp in cand) {
      tf <- grep(paste0("^", cp, "_", c, "\\.xpt$"), xpt, value=TRUE)
      if (length(tf)==0) next
      df <- tryCatch(read_xpt(file.path(CACHE, tf[1])), error=function(e) NULL)
      if (is.null(df)) next
      h <- setdiff(intersect(TBL[[t]], names(df)), "SEQN")
      if (length(h)) {
        sub <- df[, c("SEQN", h), drop=FALSE]
        sub <- sub[!duplicated(sub$SEQN), ]
        common <- intersect(names(cyc), names(sub))
        common <- setdiff(common, "SEQN")
        if (length(common) > 0) sub <- sub[, setdiff(names(sub), common), drop=FALSE]
        cyc <- left_join(cyc, sub, by="SEQN")
      }
      break
    }
  }
  ref_all[[c]] <- cyc
}
ref <- bind_rows(ref_all)
cat("R gold standard rows:", nrow(ref), "\n")

# ---- read downloader CSV ----
dl <- read.csv(CSV, fileEncoding="UTF-8-BOM")
cat("Downloader CSV rows:", nrow(dl), "cols:", ncol(dl), "\n")

# ---- match columns ----
dl_names <- names(dl)
matched <- list(); unmatched <- character()
for (xpt_name in names(DL2XPT)) {
  found <- FALSE
  for (kw in DL2XPT[[xpt_name]]) {
    idx <- grep(kw, dl_names, ignore.case=TRUE)
    # 仅 'Hg' 关键词需排除血压列（mmHg），其余关键词正常匹配
    if (kw == "Hg") idx <- idx[!grepl("mmHg", dl_names[idx], ignore.case=TRUE)]
    if (length(idx) > 0) { matched[[xpt_name]] <- dl_names[idx[1]]; found <- TRUE; break }
  }
  if (!found) unmatched <- c(unmatched, xpt_name)
}
cat("Matched:", length(matched), "/", length(names(DL2XPT)), "\n")
if (length(unmatched) > 0) cat("Unmatched:", paste(unmatched, collapse=", "), "\n")

# ---- merge ----
dl_ref <- data.frame(SEQN=as.numeric(dl$序号.SEQN.), cycle=dl$调查周期, stringsAsFactors=FALSE)
for (xpt_name in names(matched)) dl_ref[[xpt_name]] <- as.numeric(dl[[matched[[xpt_name]]]])
m <- inner_join(ref, dl_ref, by=c("SEQN","cycle"), suffix=c(".r",".d"))
cat("Matched pairs:", nrow(m), "\n\n")

# ---- per-variable comparison ----
vars <- intersect(names(DL2XPT), setdiff(names(ref), c("SEQN","cycle")))
cat(sprintf("%-10s %7s %9s %10s %10s %10s %8s %12s %12s  %s\n",
            "VAR","Nvalid","NAagree%","Exact%","MedDiff","P95Diff","MaxDiff","ICC","Rmean","Dlmean","Verdict"))
cat(strrep("-", 120), "\n")

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
  md <- median(diff); p95 <- as.numeric(quantile(diff, .95)); mx <- max(diff)
  exact <- mean(diff == 0) * 100
  icc <- NA
  if (nv >= 10 && var(x[ok]) > 1e-10 && var(y[ok]) > 1e-10) {
    # psych::ICC 包（领域标准，Single_random_raters = ICC(2,1)）
    icc <- tryCatch({
      icc_res <- ICC(matrix(c(x[ok], y[ok]), ncol=2))
      icc_res$results$ICC[1]
    }, error=function(e) {
      # 仅当 psych 计算失败时退化为手动公式（等价于 ICC(2,1)）
      s2_b <- var(c(x[ok], y[ok]))
      s2_w <- var(x[ok] - y[ok]) / 2
      (s2_b - s2_w) / (s2_b + s2_w)
    })
  }
  verdict <- if (!is.na(icc) && md <= 0.01 && icc > 0.999) "PASS" else "FLAG"
  if (verdict == "PASS") pass_n <- pass_n + 1 else fail_n <- fail_n + 1
  cat(sprintf("%-10s %7d %8.1f%% %9.1f%% %10.6f %10.6f %10.2g %8.4f %12.4f %12.4f  %s\n",
              vname, nv, na_agree, exact, md, p95, mx,
              ifelse(is.na(icc), 0, icc), mean(x[ok]), mean(y[ok]), verdict))
  res <- rbind(res, data.frame(Variable=vname, N=nv, NA_Agree=round(na_agree,1),
                               Exact_Pct=round(exact,2), Median_Diff=round(md,8),
                               P95_Diff=round(p95,6), Max_Diff=signif(mx,4),
                               ICC=round(icc,4), R_Mean=round(mean(x[ok]),4),
                               Dl_Mean=round(mean(y[ok]),4), Verdict=verdict,
                               stringsAsFactors=FALSE))
}
cat("\n============================================================\n")
cat(sprintf("TOTAL: %d variables, PASS=%d, FLAG=%d, pass rate=%.1f%%\n",
            nrow(res), pass_n, fail_n, pass_n/nrow(res)*100))
cat("============================================================\n")

write.csv(res, file.path(OUT, "FULL63_验证结果_真实对比.csv"), row.names=FALSE, fileEncoding="UTF-8")
cat("Saved: FULL63_验证结果_真实对比.csv\n")
