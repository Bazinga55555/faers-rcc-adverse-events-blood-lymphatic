# -*- coding: utf-8 -*-
"""
10_sensitivity.py -- sensitivity analysis (test robustness of the core findings)

Scenarios:
  S0 base        main analysis (PS/SS + BLSD core PT + all reports)
  S1 lab         additionally include lab-investigation PTs (broaden outcome definition)
  S2 rcc_ind     keep only reports with renal cell carcinoma as the indication (control indication confounding)
  S3 ps_only     keep only primary-suspect drugs (PS) (exclude SS, reduce polypharmacy noise)
  S4 single_drug keep only reports involving a single target drug (exclude drug-drug interactions)

Run mechanism-level disproportionality for each scenario and output the ROR of "headline signals"
under each scenario, to judge whether the core findings (H1/H2) are robust.

Outputs:
  05_results/tables/sensitivity_headline.csv  headline signal x scenario ROR comparison
  05_results/figures/                          (figures drawn by later scripts)
"""
import os
import sys
import io
import time
import importlib.util

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import numpy as np
import pandas as pd

# dynamically import 06_signal (filename starts with a digit)
_spec = importlib.util.spec_from_file_location("sig", os.path.join(HERE, "06_signal.py"))
sig = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(sig)

from paths import ROOT, DIR_CLEAN, DIR_RESULTS, SUB_TABLES
CLEAN = os.path.join(ROOT, DIR_CLEAN)
RES = os.path.join(ROOT, DIR_RESULTS, SUB_TABLES)

# renal cell carcinoma indication regex (avoid non-oncologic terms like renal failure / renal impairment)
RCC_PAT = r"renal cell|renal cancer|kidney cancer|renal carcinoma|kidney carcinoma|" \
          r"renal neoplasm|kidney neoplasm|renal tumour|kidney tumour|" \
          r"renal tumor|kidney tumor|renal adenocarcinoma"

# headline signals (key combinations driving the H1/H2 conclusions)
HEADLINE = [
    ("A_HIF2a", "Anaemia"),
    ("A_HIF2a", "Neutropenia"),                      # negative control
    ("A_HIF2a", "Thrombocytopenia"),                 # negative control
    ("B_VEGFR_TKI", "Polycythaemia"),
    ("B_VEGFR_TKI", "Thrombocytopenia"),
    ("B_VEGFR_TKI", "Disseminated intravascular coagulation"),
    ("C_VEGF_mAb", "Myelosuppression"),
    ("C_VEGF_mAb", "Thrombotic microangiopathy"),
    ("C_VEGF_mAb", "Neutropenia"),
    ("D_mTOR", "Lymphocele"),
    ("D_mTOR", "Iron deficiency anaemia"),
    ("D_mTOR", "Thrombotic microangiopathy"),
    ("E_ICI", "Aplasia pure red cell"),
    ("E_ICI", "Autoimmune haemolytic anaemia"),
    ("E_ICI", "Immune thrombocytopenia"),
    ("E_ICI", "Agranulocytosis"),
    ("E_ICI", "Aplastic anaemia"),
]


def build(scenario):
    """Build (drug, reac, demo) for a given scenario."""
    drug = pd.read_parquet(os.path.join(CLEAN, "drug_mapped.parquet"))
    reac = pd.read_parquet(os.path.join(CLEAN, "reac_norm.parquet"))
    demo = pd.read_parquet(os.path.join(CLEAN, "demo_dedup.parquet"))

    # ---- outcome scope
    blsd = pd.read_csv(os.path.join(CLEAN, "blsd_pt.csv"))
    keep_pt = set(blsd["pt_std"])
    if scenario == "lab":
        keep_pt |= set(pd.read_csv(os.path.join(CLEAN, "blsd_pt_lab.csv"))["pt_std"])
    reac = reac[reac["pt_std"].isin(keep_pt)]

    # ---- report scope
    keep_ids = set(demo["primaryid"])
    if scenario == "rcc_ind":
        indi = pd.read_parquet(os.path.join(CLEAN, "indi.parquet"))
        rcc = indi[indi["indi_pt"].astype(str).str.contains(RCC_PAT, case=False, na=False)]
        keep_ids &= set(rcc["primaryid"])
    if scenario == "single_drug":
        cnt = drug.groupby("primaryid")["drug_std"].nunique()
        keep_ids &= set(cnt[cnt == 1].index)

    demo = demo[demo["primaryid"].isin(keep_ids)]
    drug = drug[drug["primaryid"].isin(keep_ids)]
    reac = reac[reac["primaryid"].isin(keep_ids)]

    # ---- drug role
    if scenario == "ps_only":
        drug = drug[drug["role_cod"].str.upper() == "PS"]
    else:
        drug = drug[drug["role_cod"].str.upper().isin(["PS", "SS"])]
    return drug, reac, demo


def main():
    os.makedirs(RES, exist_ok=True)
    scenarios = [
        ("S0_base", "base"),
        ("S1_lab", "lab"),
        ("S2_rccInd", "rcc_ind"),
        ("S3_psOnly", "ps_only"),
        ("S4_singleDrug", "single_drug"),
    ]
    frames = {}
    for name, sc in scenarios:
        t0 = time.time()
        drug, reac, demo = build(sc)
        res = sig.run(drug, reac, by_mechanism=True, verbose=False)
        res = res.set_index(["mechanism", "pt_disp"])
        frames[name] = res
        print("[%s] %-14s reports %6d, drug rows %7d, outcome rows %6d, combos %5d, dual-positive %4d (%.0fs)"
              % (time.strftime("%H:%M:%S"), name, demo["primaryid"].nunique(), len(drug),
                 len(reac), len(res), int(res["primary_signal"].sum()), time.time() - t0))

    rows = []
    for mech, pt in HEADLINE:
        row = {"mechanism": mech, "pt": pt}
        for name, _ in scenarios:
            res = frames[name]
            if (mech, pt) in res.index:
                r = res.loc[(mech, pt)]
                row[name] = "%.2f (%.2f-%.2f)%s" % (r["ror"], r["ror_low"], r["ror_high"],
                                                    "*" if r["primary_signal"] else "")
                row[name + "_n"] = int(r["n"])
            else:
                row[name] = "—"
                row[name + "_n"] = np.nan
        rows.append(row)
    out = pd.DataFrame(rows)
    cols = ["mechanism", "pt"] + [n for n, _ in scenarios] + \
           ["%s_n" % n for n, _ in scenarios]
    out = out[cols]
    out.to_csv(os.path.join(RES, "sensitivity_headline.csv"), index=False, encoding="utf-8-sig")

    print("\n=== Headline-signal ROR by scenario (95% CI), * = ROR+IC dual-positive ===")
    print(out[[c for c in out.columns if not c.endswith("_n")]].to_string(index=False))
    print("\n=== Effective sample size n by scenario ===")
    print(out[["mechanism", "pt"] + ["%s_n" % n for n, _ in scenarios]].to_string(index=False))
    print("\nWrote %s" % os.path.join(RES, "sensitivity_headline.csv"))


if __name__ == "__main__":
    main()
