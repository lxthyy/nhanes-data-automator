# NHANES Data Automator

**One-click NHANES data extraction with built-in R validation.**

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

---

## Overview

NHANES Data Automator is a Python desktop application (Tkinter GUI) that automates the downloading, merging, unit conversion, cleaning, and export of NHANES data across **13 survey cycles (1999–2024)**. It extracts **63 numeric variables** and **over 55 categorical variables** across 11 clinical domains.

**Key features:**
- **No-code GUI** — select cycles and variables, click export
- **Automated cross-cycle merge** — variable names harmonized across cycles
- **Built-in unit conversion** — mg/dL ↔ mmol/L, μmol/L ↔ mg/dL, etc.
- **Independent R validation** — all 63 numeric variables validated against R `haven` (100% PASS, ICC=1.0)
- **Automated QC** — 8 checks per export, golden snapshot drift detection
- **Local SQLite database** — optional pre-download of all cycles for fast queries

---

## Quick Start

```bash
# Install dependencies
pip install pandas numpy scipy statsmodels matplotlib openpyxl python-docx

# Run the application
python nhanes_data_automator.py
```

The GUI will guide you through:
1. Select survey cycles (check boxes)
2. Choose variable groups (5×5 grid)
3. Click ▶ **Start** → export CSV
4. Click ✅ **0.1N Strict** or 🔶 **0.1N Lenient** for instant validation

---

## Project Structure

```
nhanes-data-automator/
├── LICENSE                     # MIT License
├── README.md                   # This file
├── requirements.txt            # Python dependencies
├── setup.py                    # Package setup
├── nhanes_data_automator.py    # Main program (GUI + engine)
├── qc_engine.py                # QC engine + golden snapshot
├── validation_engine.py        # 0.1N validation engine
├── tests/
│   ├── __init__.py
│   └── test_core.py            # Unit tests (pytest)
└── validation_scripts/
    ├── NHANES独立验证_最终版.R    # R haven validation
    ├── 死亡率验证.R              # Mortality FWF validation
    ├── 分类变量验证.py            # Categorical mapping check
    └── 数据完整性验证.py          # Completeness check
```

---

## Validation Summary

| Tier | What | Result |
|:----:|:-----|:------|
| 1 | R `haven` cell-by-cell comparison (63 variables) | **63/63 PASS** (ICC=1.0) |
| 2 | R `survey` weighted means vs CDC reference | **All <5% diff** (except Cr 9.1%) |
| 3 | Categorical variable mapping (57+ variables) | **100% correct** |
| 4 | Mortality FWF parsing vs R `read_fwf` | **3/3 cycles PASS** |

Full validation report: [manuscript_draft_v1.docx](https://github.com/lxthyy/nhanes-data-automator)

---

## Citation

If you use this tool in your research, please cite:

> [Paper title — submitted to BMC Medical Research Methodology, 2026]

---

## License

MIT License. See [LICENSE](LICENSE) for details.
