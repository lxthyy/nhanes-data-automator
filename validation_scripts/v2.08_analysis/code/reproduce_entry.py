# -*- coding: utf-8 -*-
"""复现入口（Reproducibility runner）
用法：
    python 复现_入口.py                # 只打印复现步骤与输入检查
    python 复现_入口.py <nhanes_cache>  # 额外执行药物排除复跑（_p6c_med.py）

说明：本脚本只调度 Python 侧可离线运行的步骤；R 侧（haven/survey/psych）需要 R 4.6.0，
      脚本会打印对应命令，不代替执行。
"""
import os
import sys
import subprocess
from pathlib import Path

HERE = Path(__file__).resolve().parent
CODE = HERE / "代码"
RES = HERE / "结果表"

STEPS = [
    ("Tier 1 逐值比对（Python）", r'python 代码\_v_cross.py <tool_csv> <ref_haven_csv>'),
    ("Tier 1 参考值导出（R）",    r'Rscript 代码\_v3a_export_ref.R'),
    ("17 变量均值核对（Python）", r'python 代码\_p2_descriptive.py <FULL63_csv> 结果表'),
    ("Tier 2 加权均值（R）",      r'Rscript 代码\_p3_weighted.R'),
    ("Tier 4 死亡文件解析（R）",  r'Rscript 代码\_p4_mortality.R'),
    ("Tier 4 筛查一致性（Python/R）", r'python 代码\_p6b_tool_screening.py  /  Rscript 代码\_p6_screening.R'),
    ("药物排除（Python）",        r'python 代码\_p6c_med.py <nhanes_cache> 结果表 "代码\nhanes_downloader v2.08.py"'),
    ("QC 快照检查（Python）",     r'python 代码\_p7_qc.py <tool_export_csv> 代码'),
]


def show():
    print("=" * 68)
    print("NHANES Data Automator · 复现入口（对应正文 v1.29）")
    print("=" * 68)
    print("输入目录：代码\\（12 个脚本）｜结果表\\（7 个 CSV）｜数据\\（获取说明）")
    print("现有结果表：")
    for f in sorted(RES.glob("*")):
        print("   -", f.name)
    print("\n复现步骤（Python 可离线运行；R 需 R 4.6.0）：")
    for i, (name, cmd) in enumerate(STEPS, 1):
        print("  %d) %-28s %s" % (i, name, cmd))
    print("\n关键数字的核查命令见《数字溯源对照_v1.06.md》第五节。")


def run_medication(cache: str):
    script = CODE / "_p6c_med.py"
    tool = CODE / "hanes_downloader v2.08.py"
    if not script.exists():
        print("缺脚本：", script)
        return 1
    cmd = [sys.executable, str(script), cache, str(RES), str(tool)]
    print("\n执行：", " ".join(cmd))
    return subprocess.call(cmd)


if __name__ == "__main__":
    show()
    if len(sys.argv) > 1:
        cache = sys.argv[1]
        if not Path(cache).is_dir():
            print("\n缓存目录不存在：", cache)
            sys.exit(2)
        sys.exit(run_medication(cache))
