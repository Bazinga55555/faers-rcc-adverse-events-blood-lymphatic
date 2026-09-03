# Blood and lymphatic system disorders across five mechanistic classes of renal cell carcinoma therapy

A mechanism-stratified pharmacovigilance study of haematologic adverse events associated with
renal cell carcinoma (RCC) drugs, using the FDA Adverse Event Reporting System (FAERS/AEMS,
2004 Q1 – 2026 Q2). Signals are reproduced in two independent databases (JADER, Canada Vigilance).

## Study design

- **15 RCC drugs** grouped into **five mechanistic classes**: HIF-2α inhibitor (belzutifan),
  VEGFR-TKI (7 drugs), VEGF mAb (bevacizumab), mTOR inhibitor (everolimus/temsirolimus),
  immune checkpoint inhibitor (4 drugs).
- **Exposure / role**: primary or secondary suspect (PS/SS) only.
- **Outcomes**: 251 core + 75 laboratory BLSD Preferred Terms (MedDRA SOC 10005329).
- **Disproportionality**: ROR, PRR, BCPNN/IC, MGPS/EBGM; Haldane–Anscombe +0.5.
  Primary signal = ROR 95% CI lower bound > 1 **and** IC025 > 0, n ≥ 3.
- **Comparator**: *restricted* — each index drug/class is contrasted against the other 14 RCC drugs
  (controls for indication and protopathic bias).
- **Time-to-onset**: Weibull β shape parameter; only complete 8-digit start/event dates accepted.
- **External validation**: identical 15-drug list and restricted comparator in JADER (PMDA, 2026 Q2)
  and Canada Vigilance (to 2026-04-30). Reproduction = external ROR lower CI > 1 with a ≥ 3.

## Cohort

522,813 unique FAERS reports (557,950 report–mechanism pairs).

## Pipeline (run in order)

All analysis scripts live in the `scripts/` directory of this repository.
They read raw/cleaned data and write results through the path configuration in
`scripts/paths.py`, which by default expects the on-disk layout
`01_原始数据/`, `02_解压数据/`, `03_清洗数据/`, `05_结果/`, `00_外部库/`
under the project root. Override the root with the `FAERS_RCC_ROOT` environment
variable or a local, non-committed `paths_local.py` (see `paths.py`).

Run scripts with the Python interpreter recorded in `scripts/99_run_pipeline.py`
(or your own Python 3.9+ with pandas / NumPy / SciPy installed).

| Script | Purpose |
|---|---|
| `02_parse_merge_v2.py` | Parse 90 FAERS quarters, CRC-verify, merge to parquet |
| `04_drug_map.py` | Map drug names → 15 drugs + mechanism (regex over INN/brand/salt) |
| `05_cohort.py` | Build cohort, apply PS/SS filter |
| `07_blsd_pt.py` | Normalise PTs → MedDRA SOC 10005329; 251 core + 75 lab |
| `06_signal.py` | Mechanism- and drug-level disproportionality (4 algorithms) |
| `08_weibull.py` | Time-to-onset Weibull (strict 8-digit date convention) |
| `09_key_pt_table.py` | Key-PT summary tables |
| `10_sensitivity.py` | Five sensitivity scenarios |
| `11_figures.py` | Figures 1–4 |
| `12_summary_report.py` | Results material (44 key-PT mechanism signals) |
| `14_drug_gradient.py` | Drug-level gradient + reliability grading |
| `15_temporal_validation.py` | FAERS internal time-split validation (2004–2018 / 2019–2026) |
| `16_baseline_table.py` | Table 1 baseline characteristics (PS/SS basis) |
| `17_case_series.py` | Case-by-case review of 4 mechanistically informative signals |
| `18_diag_tto.py` | Diagnostic for TTO date-parsing bugs |
| `19_flow_diagram.py` | STROBE-style study flow diagram (Fig 0) |
| `13_ext_validation.py` | **External validation** in JADER + Canada Vigilance (restricted comparator) |

## Key results

- 44 primary signals among 38 key BLSD PTs (83 across the full vocabulary); 196 drug-level.
- belzutifan→anaemia ROR 30.0 (FAERS), reproduced in JADER at 52.2.
- External validation: JADER strictly reproduced 19/44 key-PT mechanism-level and 29/196 drug-level
  signals; Canada Vigilance directionally consistent but underpowered for rare/new agents.

## Data availability & registration

- FAERS/AEMS: US FDA (public). JADER: PMDA (public). Canada Vigilance: Health Canada (public).
- Protocol, scripts, and outputs deposited on the Open Science Framework:
  **https://doi.org/10.17605/OSF.IO/5TR96**
- Code and result tables mirrored on GitHub:
  **https://github.com/Bazinga55555/faers-rcc-adverse-events-blood-lymphatic**

## Environment

Python 3.9.13; pandas, NumPy, SciPy. Streaming parsers handle the ~1.3 GB external databases
(JADER cp932; Canada Vigilance `$`-delimited) in ~1 minute.

## Note on large data

Raw and cleaned data (`01_原始数据/`, `02_解压数据/`, `03_清洗数据/`, `00_外部库/`) are excluded
from version control (see `.gitignore`); only code and small result tables are committed.

## Licence

Released under the MIT Licence — see [`LICENSE`](LICENSE).
