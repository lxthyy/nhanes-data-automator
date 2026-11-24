# validation_scripts / v2.08_analysis

Validation archive for:

> **From SAS Transport Files to Analysis-Ready Data Without Code: Development and Four-Tier
> Cross-Software Validation of the NHANES Data Automator**

Tool source: https://github.com/lxthyy/nhanes-data-automator (tag **v2.08**, MIT License).

## Layout

```
code/                 12 validation scripts (Python + R) used for the four-tier comparison
code/reproduce_entry.py  runnable entry point (prints the pipeline; optional medication re-run)
results/              8 derived result files referenced by the manuscript
environment.txt       software versions and rerun record
```

## Quick start

```bash
python code/reproduce_entry.py                 # prints the pipeline and checks inputs
python code/reproduce_entry.py <nhanes_cache>  # also re-runs the medication-exclusion check
```

Inputs: NHANES 2007-2012 (cycles E/F/G) XPT files, and the NCHS linked mortality file
(2019 public-use release). Download URLs and cache layout are in `environment.txt`
and in the manuscript's Data availability statement.

R steps require R 4.6.0 with `haven`, `survey` and `psych`; `reproduce_entry.py` prints the
exact commands rather than running them.

## Result files

| File | Contents |
|---|---|
| `FULL63_逐值对比_64vars.csv` | element-wise comparison of all 64 numeric variables |
| `P2_Descriptive_17vars.csv` | 17 core variables, mean difference |
| `P3_Weighted_Means.csv` | survey-weighted estimates vs published references |
| `P4_Mortality.csv` | three-cycle mortality, weighted |
| `P6_Screening_Consistency.csv` | six screening scenarios, participant-level |
| `P6c_Medication_Exclusion.csv` | medication exclusion overlap (1,127 / 1,153) |
| `Supplementary_Table_S2_ICC_Detail.csv` | ICC detail table (Supplementary S2) |
