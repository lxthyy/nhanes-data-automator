# NHANES Data Automator

**One-click NHANES data extraction with built-in R validation.**

## Overview

NHANES Data Automator is a Python desktop application (Tkinter GUI) that automates the downloading, merging, unit conversion, cleaning, and export of NHANES data across **13 survey cycles (1999–2024)**. It extracts **64 numeric variables** and **55 categorical variables** across 11 clinical domains.

**Key features:**

- **No-code GUI** — select cycles and variables, click export
- **Automated cross-cycle merge** — variable names harmonized across cycles
- **Built-in unit conversion** — mg/dL ↔ mmol/L, μmol/L ↔ mg/dL, etc.
- **Cross-software consistency verification** — all 64 numeric variables compared element-wise against R `haven` (zero median difference; 62/64 variables ICC ≥0.999, with the remaining two showing exact match but reduced ICC due to near-zero variance in small subsamples)
- **Automated QC** — snapshot-based drift detection with 5% tolerance
- **Local SQLite database** — optional pre-download of all cycles for fast queries

## Quick Start

```
# Install dependencies
pip install pandas numpy scipy statsmodels matplotlib openpyxl python-docx

# Run the application
python "nhanes_downloader v2.08.py"
```
The GUI will guide you through:

1. Select survey cycles (check boxes)

2. Choose variable groups (5×5 grid)

3. Click ▶ **Start** → export CSV

4. Click ✅ **0.1N Strict** or 🔶 **0.1N Lenient** for instant validation

## Changelog

### v2.08 (2026-07-27)

**QC column-matching fix:** the automated QC module now resolves column names with unit slash/underscore variants (e.g., `mIU/L` vs `mIU_L`) through a dedicated `_find_col` lookup, applied to mean consistency, missing-rate consistency, and lipid-panel internal consistency checks. Full re-validation with 64 numeric variables (including TSH and TT4) against R `haven` confirmed zero median difference across all variables; all 19 QC checks passed (max deviation 3.6%).

### v2.06 (2026-07-27)

#### Bug Fixes

| Bug | Root Cause | Fix |
|---|---|---|
| **Pregnancy exclusion silently skipped** despite checkbox checked | ① CORE_PATTERNS missing "怀孕" → pregnancy column dropped at extraction; ② `_make_short_names` converts Chinese to English names, exclusion module then fails to find Chinese column names | ① Added `"怀孕"/"Pregnant"/"RIDEXPRG"` to CORE_PATTERNS; ② Added short-name fallback search; ③ Log level changed from skip to warning |
| **Statin checkbox clipped** at window bottom | Drug exclusion area had only 2 columns, statins sorted to bottom of second column | Changed to 4-column layout |

#### New Features

| Feature | Description |
|---|---|
| **Edit All Parameters dialog** | Age/BMI/9 QC range limits + derived indices — all editable in one popup; toggle on/off, modify min/max, instant apply |
| **2 new templates** | "China-US TSH Main Analysis (replicate Chinese paper)" — thyroid meds only + QC + age/BMI; "China-US TSH Sensitivity Analysis" — main + pregnancy exclusion |
| **BMI inclusion criteria** | BMI input field added (in popup dialog) |

#### v2.06 Updates (2026-07-27)

| Update | Description |
|---|---|
| **Version number alignment** | Code internal version (`__version__`) now matches file name |
| **Repository URL fix** | Fixed incorrect repo URL in docstring |
| **Package version** | Updated `setup.py` version to 2.06 |
| **Documentation** | Updated README.md to reference v2.06 |

## Project Structure

```
nhanes-data-automator/
├── LICENSE # MIT License
├── README.md # This file
├── requirements.txt # Python dependencies
├── setup.py # Package setup
├── nhanes_downloader v2.08.py # Main program (GUI + engine)
├── qc_engine.py # QC engine + golden snapshot
├── validation_engine.py # 0.1N validation engine
├── tests/
│ ├── __init__.py
│ └── test_core.py # Unit tests (pytest)
└── validation_scripts/
 ├── NHANES_R_validation.R # R haven validation
 ├── mortality_validation.R # Mortality FWF validation
 ├── categorical_mapping_check.py # Categorical mapping check
 └── data_completeness_check.py # Completeness check
```

## Validation Summary

The manuscript reports a four-tier cross-software consistency verification strategy (30,442 participants, NHANES cycles 2007–2012):

| Tier | What | Result |
|---|---|---|
| 1 | R `haven` element-wise comparison (64 variables) | **64/64 exact match, zero median difference**; ICC ≥0.999 for 62 variables; BMXHEAD (0.81) and BPXSY4 (0.998) showed exact match with reduced ICC due to near-zero variance in small subsamples |
| 2 | R `survey` weighted means vs CDC reference | **4/5 within 5% diff** (TC 2.9%, HDL 0.8%, glucose 4.4%, lead 4.1%; Cr 9.1%) |
| 3 | Categorical variable mapping (55 variables) | **100% correct** |
| 4 | Mortality FWF parsing vs R `read_fwf` | **3/3 cycles PASS**, expected NCHS trends |

Screening integrity: **6/6 scenarios at 100% participant-level agreement with zero false positives**. All validation scripts are available in `validation_scripts/` for independent reproduction.

## Citation

If you use this tool in your research, please cite:

Li X. Bridging the Programming Gap in Population Health Research: Development and Validation of a No-Code, Open-Source Tool for Automated NHANES Data Extraction. *Computer Methods and Programs in Biomedicine*. 2026. (in review)

## License

MIT License. See [LICENSE](https://github.com/lxthyy/nhanes-data-automator/blob/main/LICENSE) for details.
