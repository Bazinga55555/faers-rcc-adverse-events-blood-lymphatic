# -*- coding: utf-8 -*-
"""
05_cohort.py -- cohort construction + baseline-characteristics table

Join drug_mapped with demo to produce the five mechanism cohorts plus baseline
statistics (with missing rates). This step is the prerequisite for signal
detection (06) and mechanism stratification (07).

Outputs:
  03_clean_data/cohort.parquet    cohort master table (one row = one de-duplicated report of a target drug)
  05_results/tables/baseline.csv  baseline features of the five cohorts (sex/age/reporter/outcome/source country)
"""
import os
import sys
import io
import time

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import pandas as pd

from paths import ROOT, DIR_CLEAN, DIR_RESULTS, SUB_TABLES
OUT  = os.path.join(ROOT, DIR_CLEAN)
RES  = os.path.join(ROOT, DIR_RESULTS, SUB_TABLES)

MECH_LABEL = {
    "A_HIF2a": "HIF-2α inhibitor (belzutifan)",
    "B_VEGFR_TKI": "VEGFR-TKI",
    "C_VEGF_mAb": "VEGF mAb (bevacizumab)",
    "D_mTOR": "mTOR inhibitor",
    "E_ICI": "Immune checkpoint inhibitor",
}


def main():
    os.makedirs(RES, exist_ok=True)
    print("[%s] Loading data..." % time.strftime("%H:%M:%S"))
    drug = pd.read_parquet(os.path.join(OUT, "drug_mapped.parquet"))
    demo = pd.read_parquet(os.path.join(OUT, "demo_dedup.parquet"))
    reac = pd.read_parquet(os.path.join(OUT, "reac.parquet"))
    outc = pd.read_parquet(os.path.join(OUT, "outc.parquet")) if os.path.exists(os.path.join(OUT, "outc.parquet")) else pd.DataFrame()

    # de-duplicated reports
    keep_ids = set(demo["primaryid"])
    drug = drug[drug["primaryid"].isin(keep_ids)]
    role_ok = drug["role_cod"].str.upper().isin(["PS", "SS"])
    drug = drug[role_ok]

    # one report may contain several target drugs; keep the report-drug pairs and
    # de-duplicate counts by primaryid when building the cohort
    cohort = drug.merge(demo, on="primaryid", how="left")

    # outcome join
    if not outc.empty and "outc_cod" in outc.columns:
        outc_agg = outc.groupby("primaryid")["outc_cod"].agg(
            lambda s: "|".join(sorted(set(s.dropna()))))
        cohort = cohort.merge(outc_agg.rename("outc"), on="primaryid", how="left")

    # age/sex cleaning
    cohort["age"] = pd.to_numeric(cohort["age"], errors="coerce")
    cohort["gndr_cod"] = cohort["gndr_cod"].str.upper()

    # output cohort
    cohort.to_parquet(os.path.join(OUT, "cohort.parquet"), index=False)

    # ---- baseline-characteristics table ----
    rows = []
    for mech, label in MECH_LABEL.items():
        sub = cohort[cohort["mechanism"] == mech]
        n_reports = sub["primaryid"].nunique()
        n_drug_rows = len(sub)
        rows.append({
            "mechanism": label,
            "n_reports": n_reports,
            "n_drug_rows": n_drug_rows,
            "age_median": round(sub["age"].median(), 1),
            "age_missing_%": round(sub["age"].isna().mean() * 100, 1),
            "male_%": round((sub["gndr_cod"] == "M").mean() * 100, 1),
            "female_%": round((sub["gndr_cod"] == "F").mean() * 100, 1),
            "sex_missing_%": round(sub["gndr_cod"].isna().mean() * 100, 1),
        })

    baseline = pd.DataFrame(rows)
    baseline.to_csv(os.path.join(RES, "baseline.csv"), index=False, encoding="utf-8-sig")

    print("\n=== Five-cohort baseline ===")
    print(baseline.to_string(index=False))
    print("\n[%s] Done; cohort.parquet has %d rows" % (time.strftime("%H:%M:%S"), len(cohort)))


if __name__ == "__main__":
    main()
