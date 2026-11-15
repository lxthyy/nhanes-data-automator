# V3a: export E/F/G wide table read by R foreign::read.xport (for Python independent checks)
#
# 第一人称说明：
# 我发现：Python 侧要独立验证，但 pyreadstat 读 CDC XPT 会报字符集错误。
# 我解决：先用 R foreign::read.xport（独立解析器）把 E/F/G 读成一张宽表
#         （ref_foreign_wide.csv），作为 Python 侧"不经工具代码"的数据源。
# 我验证：该宽表与 haven 读取结果逐值一致（见 V1b），后续 V2-V7 的 Python
#         交叉验证都用它当输入。
suppressMessages({library(dplyr)})
args <- commandArgs(trailingOnly = TRUE)
CACHE <- args[1]; OUT <- args[2]

CYCLES <- c(E = "2007-2008", F = "2009-2010", G = "2011-2012")
TBL <- list(
  DEMO = c("SEQN","RIAGENDR","RIDAGEYR","DMDEDUC2","INDFMPIR","SDMVPSU","SDMVSTRA","WTMEC2YR","WTINT2YR"),
  BMX = c("SEQN","BMXWT","BMXHT","BMXBMI","BMXWAIST","BMXTRI","BMXARMC","BMXHEAD","BMXLEG","BMXARML"),
  BPX = c("SEQN","BPXSY1","BPXDI1","BPXSY2","BPXDI2","BPXSY3","BPXDI3","BPXSY4","BPXDI4","BPXPLS"),
  TCHOL = c("SEQN","LBXTC"), HDL = c("SEQN","LBDHDD"),
  TRIGLY = c("SEQN","LBXTR","LBDLDL"), GLU = c("SEQN","LBXGLU"),
  THYROD = c("SEQN","LBXTSH1","LBXT4F","LBXT3F","LBXTT3","LBXTT4","LBXATG","LBXTPO"),
  BIOPRO = c("SEQN","LBXSCR","LBXSUA","LBXSBU","LBXSASSI","LBXSATSI","LBXSGTSI","LBXSTB","LBXSAL","LBXSAPSI"),
  CBC = c("SEQN","LBXWBCSI","LBXHGB","LBXPLTSI","LBXHCT","LBDLYMNO","LBDMONO","LBXNEPCT","LBXEOPCT"),
  PBCD = c("SEQN","LBXBPB","LBXBCD","LBXTHG"), SLQ = c("SEQN","SLD010H"),
  DR1TOT = c("SEQN","DR1TKCAL","DR1TPROT","DR1TCARB","DR1TTFAT","DR1TSFAT","DR1TFIBE","DR1TSODI","DR1TALCO")
)
read_fx <- function(fp) {
  z <- tryCatch(foreign::read.xport(fp), error = function(e) NULL)
  if (is.null(z)) return(NULL)
  z$SEQN <- as.numeric(z$SEQN); z
}
ref_all <- list()
for (c in names(CYCLES)) {
  d <- read_fx(file.path(CACHE, sprintf("DEMO_%s.xpt", c)))
  if (is.null(d)) next
  cyc <- data.frame(SEQN = d$SEQN, cycle = CYCLES[c], stringsAsFactors = FALSE)
  for (t in names(TBL)) {
    cand <- if (t %in% c("TCHOL","THYROD")) c(t, "TST") else t
    for (cp in cand) {
      tf <- list.files(CACHE, pattern = paste0("^", cp, "_", c, "\\.xpt$"), ignore.case = TRUE)
      if (length(tf) == 0) next
      tb <- read_fx(file.path(CACHE, tf[1]))
      if (is.null(tb)) next
      keep <- setdiff(intersect(TBL[[t]], names(tb)), "SEQN")
      cyc <- left_join(cyc, tb[, c("SEQN", keep), drop = FALSE], by = "SEQN")
      break
    }
  }
  ref_all[[c]] <- cyc
}
ref <- bind_rows(ref_all)
out_csv <- file.path(OUT, "ref_foreign_wide.csv")
write.csv(ref, out_csv, row.names = FALSE, fileEncoding = "UTF-8")
cat("foreign 宽表已导出:", out_csv, "|", nrow(ref), "行 x", ncol(ref), "列\n")
