# -*- coding: utf-8 -*-
"""
14_drug_gradient.py -- drug-level gradient validation + signal reliability grading

Purpose:
  1. Test whether the erythrocytosis gradient within VEGFR-TKI is robust (excluding paraneoplastic EPO confounding)
     Observation: axitinib 22.10 > lenvatinib 8.57 > pazopanib 6.33 > cabozantinib 5.17
                  >> sorafenib 1.22 ≈ sunitinib 0.99 (c-KIT co-inhibition may offset)
  2. Reliability grading of all drug-level primary signals (down-weight small samples)

Scenarios (consistent with 10_sensitivity.py):
  S0_base / S1_lab / S2_rccInd / S3_psOnly / S4_singleDrug

Outputs:
  05_results/tables/drug_gradient_sensitivity.csv  drug x PT x scenario ROR
  05_results/tables/signal_reliability.csv         reliability grade of drug-level primary signals
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

_spec = importlib.util.spec_from_file_location("sig", os.path.join(HERE, "06_signal.py"))
sig = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(sig)

from paths import ROOT, DIR_CLEAN, DIR_RESULTS, SUB_TABLES
CLEAN = os.path.join(ROOT, DIR_CLEAN)
RES   = os.path.join(ROOT, DIR_RESULTS, SUB_TABLES)

RCC_PAT = r"renal cell|renal cancer|kidney cancer|renal carcinoma|kidney carcinoma|" \
          r"renal neoplasm|kidney neoplasm|renal tumour|kidney tumour|" \
          r"renal tumor|kidney tumor|renal adenocarcinoma"

# drug x PT combos to validate (erythroid axis + key controls)
GRADIENT = [
    # VEGFR-TKI internal gradient: erythrocytosis
    ("axitinib",     "Polycythaemia"),
    ("lenvatinib",   "Polycythaemia"),
    ("pazopanib",    "Polycythaemia"),
    ("cabozantinib", "Polycythaemia"),
    ("sorafenib",    "Polycythaemia"),
    ("sunitinib",    "Polycythaemia"),
    # control: upstream HIF-2α blockade
    ("belzutifan",   "Anaemia"),
    ("belzutifan",   "Polycythaemia"),
    # control: VEGF ligand blockade
    ("bevacizumab",  "Polycythaemia"),
    # control: mTOR
    ("everolimus",   "Iron deficiency anaemia"),
    ("everolimus",   "Polycythaemia"),
]

SCENARIOS = [
    ("S0_base", "base"),
    ("S1_lab", "lab"),
    ("S2_rccInd", "rcc_ind"),
    ("S3_psOnly", "ps_only"),
    ("S4_singleDrug", "single_drug"),
]


def build(scenario):
    """Identical to build() in 10_sensitivity.py."""
    drug = pd.read_parquet(os.path.join(CLEAN, "drug_mapped.parquet"))
    reac = pd.read_parquet(os.path.join(CLEAN, "reac_norm.parquet"))
    demo = pd.read_parquet(os.path.join(CLEAN, "demo_dedup.parquet"))

    blsd = pd.read_csv(os.path.join(CLEAN, "blsd_pt.csv"))
    keep_pt = set(blsd["pt_std"])
    if scenario == "lab":
        keep_pt |= set(pd.read_csv(os.path.join(CLEAN, "blsd_pt_lab.csv"))["pt_std"])
    reac = reac[reac["pt_std"].isin(keep_pt)]

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

    if scenario == "ps_only":
        drug = drug[drug["role_cod"].str.upper() == "PS"]
    else:
        drug = drug[drug["role_cod"].str.upper().isin(["PS", "SS"])]
    return drug, reac, demo


def reliability_grade(row):
    """Signal reliability grading:
      A = n>=20 and all four algorithms positive
      B = n>=10 and ROR+IC dual-positive
      C = n>=5  and ROR+IC dual-positive (exploratory, interpret with caution)
      D = n<5   (small sample, list only, no interpretation)
    """
    if not row["primary_signal"]:
        return "-"
    n = row["n"]
    algo = row.get("n_algo", 0)
    if n >= 20 and algo >= 4:
        return "A"
    if n >= 10:
        return "B"
    if n >= 5:
        return "C"
    return "D"


def main():
    os.makedirs(RES, exist_ok=True)

    # ---------- Part 1: drug-level gradient sensitivity ----------
    print("=" * 70)
    print("Part 1: drug-level gradient sensitivity analysis")
    print("=" * 70)
    frames = {}
    for name, sc in SCENARIOS:
        t0 = time.time()
        drug, reac, demo = build(sc)
        res = sig.run(drug, reac, by_mechanism=False, verbose=False)
        res = res.set_index(["drug_std", "pt_disp"])
        frames[name] = res
        print("[%s] %-14s reports %6d, drug rows %7d, dual-positive %4d (%.0fs)"
              % (time.strftime("%H:%M:%S"), name, demo["primaryid"].nunique(),
                 len(drug), int(res["primary_signal"].sum()), time.time() - t0))

    rows = []
    for d, pt in GRADIENT:
        row = {"drug": d, "pt": pt}
        for name, _ in SCENARIOS:
            res = frames[name]
            if (d, pt) in res.index:
                r = res.loc[(d, pt)]
                row[name] = "%.2f (%.2f-%.2f)%s" % (r["ror"], r["ror_low"], r["ror_high"],
                                                    "*" if r["primary_signal"] else "")
                row[name + "_n"] = int(r["n"])
            else:
                row[name] = "-"
                row[name + "_n"] = np.nan
        rows.append(row)
    out = pd.DataFrame(rows)
    cols = (["drug", "pt"]
            + [n for n, _ in SCENARIOS]
            + ["%s_n" % n for n, _ in SCENARIOS])
    out = out[cols]
    p1 = os.path.join(RES, "drug_gradient_sensitivity.csv")
    out.to_csv(p1, index=False, encoding="utf-8-sig")

    print("\n=== Drug x PT ROR by scenario (95% CI), * = ROR+IC dual-positive ===")
    print(out[[c for c in out.columns if not c.endswith("_n")]].to_string(index=False))
    print("\n=== n by scenario ===")
    print(out[["drug", "pt"] + ["%s_n" % n for n, _ in SCENARIOS]].to_string(index=False))
    print("\nWrote %s" % p1)

    # ---------- Part 2: signal reliability grading ----------
    print("\n" + "=" * 70)
    print("Part 2: drug-level primary-signal reliability grading")
    print("=" * 70)
    sd = pd.read_csv(os.path.join(RES, "signal_drug.csv"))
    sig_on = sd[sd["primary_signal"] == True].copy()
    sig_on["grade"] = sig_on.apply(reliability_grade, axis=1)
    print("\nGrade distribution:")
    for g in ["A", "B", "C", "D"]:
        sub = sig_on[sig_on["grade"] == g]
        print("  grade %s: %3d" % (g, len(sub)))

    print("\n=== Grade A (n>=20 and all four algorithms positive, most reliable) ===")
    a = sig_on[sig_on["grade"] == "A"].sort_values("ror", ascending=False)
    print(a[["drug_std", "mechanism", "pt_disp", "n", "ror", "ror_low", "ror_high",
             "ic025", "ebgm05", "n_algo"]].to_string(index=False))

    print("\n=== Grade D (n<5, exclude from interpretation) ===")
    d = sig_on[sig_on["grade"] == "D"].sort_values("ror", ascending=False)
    print(d[["drug_std", "mechanism", "pt_disp", "n", "ror", "ror_low", "ror_high",
             "ic025"]].to_string(index=False))

    p2 = os.path.join(RES, "signal_reliability.csv")
    sig_on.sort_values(["grade", "ror"], ascending=[True, False]).to_csv(
        p2, index=False, encoding="utf-8-sig")
    print("\nWrote %s" % p2)

    # ---------- Part 3: VEGFR-TKI internal gradient summary ----------
    print("\n" + "=" * 70)
    print("Part 3: VEGFR-TKI internal erythrocytosis gradient")
    print("=" * 70)
    tki = ["axitinib", "lenvatinib", "pazopanib", "cabozantinib", "sorafenib", "sunitinib"]
    g = sd[(sd["pt_disp"] == "Polycythaemia") & (sd["drug_std"].isin(tki))]
    g = g.sort_values("ror", ascending=False)
    print(g[["drug_std", "n", "ror", "ror_low", "ror_high", "ic025", "ebgm05",
             "n_algo", "primary_signal"]].to_string(index=False))
    print("\nNote: axitinib/lenvatinib/pazopanib are selective VEGFR-TKIs;")
    print("      sunitinib/sorafenib additionally inhibit c-KIT (essential for haematopoietic stem cells),")
    print("      so the erythrocytosis signal disappears.")


if __name__ == "__main__":
    main()
