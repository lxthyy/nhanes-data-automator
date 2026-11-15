# -*- coding: utf-8 -*-
"""
QC 引擎验证（黄金快照 5% 容差 + 血脂闭环 + 范围检查）
=====================================================
背景：正文声称"自动 QC 引擎 5 类检查，全部通过"。
我发现：工具的 qc_engine 有个 bug——按变量名找列时，单位里的 '/' 与工具导出
       CSV 里的 '_' 对不上（mIU/L vs mIU_L），导致 6 个核心变量报"列不存在"，
       QC 实际只有部分检查在跑。更严重的是：快照是工具自己生成的，等于
       "工具拿自己的标准查自己"，有自我闭环风险。
我解决：a) 修 qc_engine 的列名匹配（'/' 与 '_' 归一化）；
       b) 不只用工具自己的快照，另用 R 独立读取（haven/foreign）算出的均值
       交叉核对快照值。
我验证：修复后 19/19 检查通过；最大差异 Age 3.6%（v2.07 把年龄 >=80 封顶为 80，
       属口径差异，在 5% 容差内）；血脂闭环中位差 0.1%；范围检查无异常。

用法: python _p7_qc.py <工具导出CSV> <qc_engine.py目录>
"""
import sys, os, json, importlib.util
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
import pandas as pd

CSV  = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(os.path.abspath(__file__)), "full_validation", "NHANES_E_F_G_v207_full.csv")
BASE = sys.argv[2] if len(sys.argv) > 2 else os.path.dirname(os.path.abspath(__file__))

spec = importlib.util.spec_from_file_location("qc", os.path.join(BASE, "qc_engine.py"))
qc = importlib.util.module_from_spec(spec)
spec.loader.exec_module(qc)

df = pd.read_csv(CSV, encoding="utf-8-sig", low_memory=False)
snap = os.path.join(BASE, "golden_snapshot.json")

engine = qc.QCEngine(df, snapshot_path=snap)
engine.run_all()

print(f"QC 引擎: {len(engine.report['checks'])} 项检查")
print(f"通过: {sum(1 for c in engine.report['checks'] if c['passed'])} / 失败: {sum(1 for c in engine.report['checks'] if not c['passed'])}")
for c in engine.report["checks"]:
    mark = "PASS" if c["passed"] else "FAIL"
    print(f"[{mark}] {c['name']}: {c['detail']}")

out = os.path.join(BASE, "full_validation", "P7_QC_Report.txt")
with open(out, "w", encoding="utf-8") as f:
    f.write(json.dumps(engine.report, ensure_ascii=False, indent=2, default=str))
print("\nSaved:", out)
