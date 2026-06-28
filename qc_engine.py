# -*- coding: utf-8 -*-
"""生成黄金快照 + QCEngine"""

import json, os, sys
import pandas as pd
import numpy as np

sys.stdout.reconfigure(encoding='utf-8')

SNAPSHOT_PATH = r'C:\Users\lxddz\Desktop\NHANES_Data_Automator_v2.05_Validation\NHANES_E_F_G_golden_snapshot.json'

# 10个核心变量（覆盖不同检测类型）
CORE_VARS = [
    ('年龄(Age,岁)', 'numeric'),
    ('BMI(kg/m²)', 'numeric'),
    ('收缩压-1(SBP1,mmHg)', 'numeric'),
    ('促甲状腺激素(TSH,mIU/L)', 'numeric'),
    ('总胆固醇(TC,mg/dL)', 'numeric'),
    ('高密度脂蛋白(HDL,mg/dL)', 'numeric'),
    ('空腹血糖(Glu,mg/dL)', 'numeric'),
    ('白细胞(WBC,×10⁹/L)', 'numeric'),
    ('性别(Gender)', 'categorical'),
    ('教育程度(Education)', 'categorical'),
]


def generate_snapshot(csv_path):
    """从已验证的CSV生成黄金快照"""
    print(f'读取: {csv_path}')
    df = pd.read_csv(csv_path, encoding='utf-8-sig', low_memory=False)
    print(f'行数: {len(df)}, 列数: {len(df.columns)}')

    snapshot = {
            'n_total': int(len(df)),
            'n_cycles': {str(k): int(v) for k, v in df['调查周期'].value_counts().items()},
            'generated_at': pd.Timestamp.now().strftime('%Y-%m-%d %H:%M'),
            'source': 'NHANES E/F/G (2007-2012), validated v2.05',
            'core_vars': {}
        }

    for var_name, var_type in CORE_VARS:
        found = [c for c in df.columns if var_name in c]
        if not found:
            print(f'  ⚠️ 未找到: {var_name}')
            continue
        col = found[0]
        if var_type == 'numeric':
            vals = pd.to_numeric(df[col], errors='coerce')
            snapshot['core_vars'][var_name] = {
                'type': 'numeric',
                'mean': round(float(vals.mean()), 4),
                'std': round(float(vals.std()), 4),
                'missing_pct': round(float(vals.isna().mean() * 100), 2),
                'n_valid': int(vals.notna().sum()),
            }
        elif var_type == 'categorical':
            vc = df[col].value_counts()
            snapshot['core_vars'][var_name] = {
                'type': 'categorical',
                'categories': {str(k): int(v) for k, v in vc.items()},
                'missing_pct': round(float(df[col].isna().mean() * 100), 2),
            }
        print(f'  ✅ {col}')

    with open(SNAPSHOT_PATH, 'w', encoding='utf-8') as f:
        json.dump(snapshot, f, ensure_ascii=False, indent=2)
    print(f'\n黄金快照已保存: {SNAPSHOT_PATH}')
    print(f'  N={snapshot["n_total"]}, 核心变量={len(snapshot["core_vars"])}')
    return snapshot


class QCEngine:
    """轻量数据健康检查引擎"""

    def __init__(self, df, snapshot_path=None):
        self.df = df
        self.snapshot_path = snapshot_path or SNAPSHOT_PATH
        self.report = {
            'timestamp': pd.Timestamp.now().strftime('%Y-%m-%d %H:%M'),
            'n_current': len(df),
            'checks': [],
            'passed': True,
        }

        if os.path.exists(self.snapshot_path):
            with open(self.snapshot_path, 'r', encoding='utf-8') as f:
                self.snapshot = json.load(f)
        else:
            self.snapshot = None

    def add_check(self, name, passed, detail, value=None, threshold=None):
        self.report['checks'].append({
            'name': name,
            'passed': passed,
            'detail': detail,
            'value': value,
            'threshold': threshold,
        })
        if not passed:
            self.report['passed'] = False

    def check_mean(self, var_name, threshold=5.0):
        """关键变量均值 vs 黄金快照"""
        if not self.snapshot or var_name not in self.snapshot['core_vars']:
            self.add_check(var_name, False, '无黄金快照参考', None, None)
            return
        ref = self.snapshot['core_vars'][var_name]
        if ref['type'] != 'numeric':
            return
        found = [c for c in self.df.columns if var_name in c]
        if not found:
            self.add_check(var_name, False, '列不存在', None, None)
            return
        vals = pd.to_numeric(self.df[found[0]], errors='coerce')
        curr_mean = float(vals.mean())
        diff_pct = abs(curr_mean - ref['mean']) / ref['mean'] * 100 if ref['mean'] != 0 else 0
        passed = diff_pct < threshold
        self.add_check(
            f'均值: {var_name}',
            passed,
            f'当前={curr_mean:.2f}, 快照={ref["mean"]:.2f}, 差异={diff_pct:.2f}%',
            round(diff_pct, 2), f'<{threshold}%'
        )

    def check_missing(self, var_name, threshold=5.0):
        """缺失率 vs 黄金快照"""
        if not self.snapshot or var_name not in self.snapshot['core_vars']:
            return
        ref = self.snapshot['core_vars'][var_name]
        found = [c for c in self.df.columns if var_name in c]
        if not found:
            return
        curr_miss = round(float(self.df[found[0]].isna().mean() * 100), 2)
        diff = abs(curr_miss - ref['missing_pct'])
        passed = diff < threshold
        self.add_check(
            f'缺失率: {var_name}',
            passed,
            f'当前={curr_miss:.1f}%, 快照={ref["missing_pct"]:.1f}%, 差异={diff:.1f}%',
            round(diff, 1), f'<{threshold}%'
        )

    def check_cleaning_conservation(self, original_n, steps):
        """清洗存续守恒：原始 - ∑排除 = 保留"""
        accounted = sum(r for _, r in steps)
        remaining = original_n - accounted
        passed = remaining == self.report['n_current']
        self.add_check(
            '清洗存续守恒',
            passed,
            f'{original_n} - {accounted} = {remaining}, 实际保留={self.report["n_current"]}',
            remaining - self.report['n_current'], 0
        )

    def check_lipid_closure(self):
        """血脂闭环：LDL + HDL + TG/5 ≈ TC"""
        cols = {k: [c for c in self.df.columns if k in c] for k in
                ['LDL.mg', 'HDL.mg', 'TG.mg', 'TC.mg']}
        if not all(cols.values()):
            return
        ldl = pd.to_numeric(self.df[cols['LDL.mg'][0]], errors='coerce')
        hdl = pd.to_numeric(self.df[cols['HDL.mg'][0]], errors='coerce')
        tg = pd.to_numeric(self.df[cols['TG.mg'][0]], errors='coerce')
        tc = pd.to_numeric(self.df[cols['TC.mg'][0]], errors='coerce')
        valid = ldl.notna() & hdl.notna() & tg.notna() & tc.notna()
        if valid.sum() < 10:
            return
        estimated = ldl + hdl + tg / 5
        diff = (abs(estimated - tc) / tc * 100).median()
        passed = diff < 10
        self.add_check(
            '血脂闭环(LDL+HDL+TG/5≈TC)',
            passed,
            f'中位数差异={diff:.1f}%, N={int(valid.sum())}',
            round(diff, 1), '<10%'
        )

    def check_range_bounds(self):
        """所有数值变量是否在临床合理范围内"""
        bounds = {'SBP': (50, 250), 'DBP': (20, 150), 'BMI': (10, 60),
                  'TSH': (0.001, 100), 'Glu': (10, 500), 'WBC': (0.1, 50)}
        for kw, (lo, hi) in bounds.items():
            found = [c for c in self.df.columns if kw in c.upper()]
            if not found:
                continue
            vals = pd.to_numeric(self.df[found[0]], errors='coerce')
            out = vals[(vals < lo) | (vals > hi)].dropna()
            passed = len(out) / len(vals.notna()) * 100 < 1.0
            if not passed:
                self.add_check(
                    f'范围: {found[0]}',
                    passed,
                    f'{len(out)}条超出[{lo},{hi}], 范围=[{vals.min():.1f},{vals.max():.1f}]',
                    int(len(out)), 0
                )

    def run_all(self, df_before_cleaning=None, cleaning_steps=None):
        """运行全部健康检查"""
        # 均值检查（10个核心变量）
        for var_name, _ in CORE_VARS:
            self.check_mean(var_name)

        # 缺失率检查（10个核心变量）
        for var_name, _ in CORE_VARS:
            self.check_missing(var_name)

        # 清洗存续守恒
        if df_before_cleaning is not None and cleaning_steps:
            self.check_cleaning_conservation(len(df_before_cleaning), cleaning_steps)

        # 血脂闭环
        self.check_lipid_closure()

        # 范围边界
        self.check_range_bounds()

        # 摘要
        total = len(self.report['checks'])
        passed_n = sum(1 for c in self.report['checks'] if c['passed'])
        self.report['summary'] = f'{passed_n}/{total} 检查通过'
        self.report['passed'] = passed_n == total
        return self.report


if __name__ == '__main__':
    # 生成黄金快照
    csv = r'C:\Users\lxddz\Desktop\NHANES_Data_Automator_v2.05_Validation\NHANES_E_F_G_完整数据_v205.csv'
    generate_snapshot(csv)
