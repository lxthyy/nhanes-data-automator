# -*- coding: utf-8 -*-
"""
药物排除验证（RXQ_RX 处方记录）
==============================
背景：正文写"药物排除 97.7% SEQN 重叠"。
我发现：这个数字来自早期验证，且早期用的是"关键词匹配处方药名"这单一方法，
       容易自我闭环——关键词写宽了数字就变，说不清。
我解决：a) 用工具自己的 DRUG_CLASSIFICATION 关键词在 RXQ_RX E/F/G 里筛出
       甲状腺用药者；b) 再用早期验证的关键词筛一遍；c) 额外用第三套独立核心
       词根（levothyrox/methimazol/propylthiourac 等）交叉验证。
我验证：工具=1129 人、早期关键词=1153 人、重叠 1127 → 97.7%；
       第三套独立词根筛出 1129 人，与工具完全重叠 → 不是关键词巧合。

用法: python _p6c_med.py <nhanes_cache> <输出目录> [工具py路径]
"""
import sys, os, importlib.util, pandas as pd
sys.stdout.reconfigure(encoding="utf-8", errors="replace")

CACHE = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(os.path.abspath(__file__)), "nhanes_cache")
OUT   = sys.argv[2] if len(sys.argv) > 2 else os.path.dirname(os.path.abspath(__file__))
SRC   = sys.argv[3] if len(sys.argv) > 3 else os.path.join(os.path.dirname(os.path.abspath(__file__)), "nhanes_downloader v2.08.py")

spec = importlib.util.spec_from_file_location("nd", SRC)
nd = importlib.util.module_from_spec(spec)
spec.loader.exec_module(nd)

def read_xpt(fp):
    try:
        return pd.read_sas(fp, format="xport", encoding="utf-8")
    except Exception:
        return pd.read_sas(fp, format="xport", encoding="latin-1")

rx_all = []
for cyc in ["E", "F", "G"]:
    fp = os.path.join(CACHE, f"RXQ_RX_{cyc}.xpt")
    if not os.path.exists(fp):
        print("缺", fp); continue
    d = read_xpt(fp)
    d["SEQN"] = pd.to_numeric(d["SEQN"], errors="coerce")
    rx_all.append(d[["SEQN", "RXDDRUG"]])
rx = pd.concat(rx_all, ignore_index=True)
rx = rx[rx["RXDDRUG"].notna()]
rx["drug"] = rx["RXDDRUG"].astype(str).str.lower()
print("RXQ_RX 处方记录:", len(rx), "| 用药者 SEQN:", rx["SEQN"].nunique())

# 工具侧（DRUG_CLASSIFICATION 甲状腺关键词）
tool_kws = [k.lower() for k in nd.DRUG_CLASSIFICATION["甲状腺药物(Thyroid)"]["keywords"]]
tool_seqn = set(rx[rx["drug"].apply(lambda s: any(k in s for k in tool_kws))]["SEQN"].dropna().astype(int))
# 早期验证关键词（R_final_v3）
r_kws = ["THYROID", "LEVOTHYROX", "METHIMAZ", "SYNTHROID", "PROPYLTHI"]
r_seqn = set(rx[rx["drug"].apply(lambda s: any(k.lower() in s for k in r_kws))]["SEQN"].dropna().astype(int))
# 第三套独立核心词根
root3 = ["levothyrox", "synthroid", "levoxyl", "liothyronin", "cytomel",
         "armour thyro", "nature-thyro", "methimazol", "tapazole",
         "propylthiourac", "ptu", "thyrolar", "左甲状腺素", "优甲乐", "甲巯咪唑", "赛治", "丙硫氧嘧啶"]
r3 = set(rx[rx["drug"].apply(lambda s: any(k in s for k in root3))]["SEQN"].dropna().astype(int))

ov = len(tool_seqn & r_seqn)
overlap_pct = ov / len(r_seqn) * 100 if r_seqn else 0
print(f"工具={len(tool_seqn)} 早期={len(r_seqn)} 重叠={ov} ({overlap_pct:.1f}%) 工具独有={len(tool_seqn-r_seqn)} 早期独有={len(r_seqn-tool_seqn)}")
print(f"第三套独立词根={len(r3)} 与工具重叠={len(tool_seqn & r3)} 工具独有={len(tool_seqn-r3)}")

result = pd.DataFrame([{"Metric": "R_thyroid_users", "Count": len(r_seqn)},
                       {"Metric": "Tool_excluded", "Count": len(tool_seqn)},
                       {"Metric": "Overlap", "Count": ov},
                       {"Metric": "Tool_only", "Count": len(tool_seqn - r_seqn)},
                       {"Metric": "R_only", "Count": len(r_seqn - tool_seqn)},
                       {"Metric": "Overlap_pct_of_R", "Count": round(overlap_pct, 2)},
                       {"Metric": "Independent_root3", "Count": len(r3)},
                       {"Metric": "Root3_overlap_tool", "Count": len(tool_seqn & r3)}])
result.to_csv(os.path.join(OUT, "P6c_Medication_Exclusion.csv"), index=False, encoding="utf-8-sig")
print("Saved:", os.path.join(OUT, "P6c_Medication_Exclusion.csv"))
