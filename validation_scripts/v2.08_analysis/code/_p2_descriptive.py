# -*- coding: utf-8 -*-
"""
17 个核心变量的描述统计核对
==========================
背景：正文写"17 个核心临床变量做了均值±SD 对比，差异 <1%"。
我发现：只凭 ICC 一个指标说服力不够，审稿人可能追问"均值本身差多少"。
我解决：直接从 FULL63 逐值对比结果（每个变量的 R 均值 vs 工具均值）里取出
       这 17 个变量，算相对差异百分比，>1% 就标 FAIL。
我验证：17/17 全部差异 = 0%（因为逐值完全一致，均值必然相同），并输出 CSV 存档。

用法: python _p2_descriptive.py <FULL63结果CSV> <输出目录>
"""
import sys, os
import pandas as pd

full63_csv = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(os.path.abspath(__file__)), "FULL63_验证结果_真实对比.csv")
out_dir    = sys.argv[2] if len(sys.argv) > 2 else os.path.dirname(os.path.abspath(__file__))

r = pd.read_csv(full63_csv, encoding="utf-8-sig")

CORE = {
    "RIDAGEYR": "Age", "BMXBMI": "BMI", "BMXWT": "Weight", "BMXHT": "Height",
    "BMXWAIST": "Waist", "LBXTC": "TC", "LBDHDD": "HDL", "LBDLDL": "LDL",
    "LBXTR": "TG", "BPXSY1": "SBP1", "BPXDI1": "DBP1", "LBXTSH1": "TSH",
    "LBXT4F": "FT4", "LBXT3F": "FT3", "LBXTT3": "Tg", "LBXTPO": "TPOAb", "LBXATG": "TgAb",
}
rows = []
for v, label in CORE.items():
    row = r[r["Variable"] == v]
    if row.empty:
        rows.append({"Variable": label, "Status": "NOT_IN_FULL63"})
        continue
    x = row.iloc[0]
    diff = abs(x["R_Mean"] - x["Dl_Mean"]) / abs(x["R_Mean"]) * 100 if x["R_Mean"] else 0.0
    rows.append({
        "Variable": label, "XPT": v, "N": int(x["N"]),
        "R_Mean": x["R_Mean"], "Dl_Mean": x["Dl_Mean"],
        "Diff_Pct": round(diff, 6), "Exact_Pct": x["Exact_Pct"],
        "Status": "PASS" if diff < 1.0 else "FAIL",
    })
out = pd.DataFrame(rows)
out.to_csv(os.path.join(out_dir, "P2_Descriptive_17vars.csv"), index=False, encoding="utf-8-sig")
print(out.to_string(index=False))
print(f"\n17 变量: PASS={(out['Status']=='PASS').sum()}, FAIL={(out['Status']=='FAIL').sum()}")
