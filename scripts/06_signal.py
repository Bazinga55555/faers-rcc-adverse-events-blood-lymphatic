# -*- coding: utf-8 -*-
"""
06_signal.py -- disproportionality analysis (four algorithms) signal detection (BLSD disease-specific)

For each drug (or each mechanism group), compute its signals on the
"blood and lymphatic system disorders (BLSD)" PTs:
  - ROR   (Reporting Odds Ratio)  + 95% CI        -- Haldane correction
  - PRR   (Proportional Reporting Ratio) + chi-square (Yates)
  - BCPNN (Information Component, IC) + IC025    -- Norén approximation
  - GPS   (Gamma Poisson Shrinker, EBGM) + EB05  -- Gamma(1,1) conjugate prior

Key design
  1. Outcome restriction: REAC is first filtered by the BLSD core PT list produced by
     07_blsd_pt.py (default); --include-lab additionally includes lab-type PTs (sensitivity analysis).
  2. PT normalization: use pt_std from reac_norm.parquet (case variants already merged).
  3. Unified scale: a/b/c/d are all drug-PT combination counts (de-duplicated within the same report).
  4. Background comparator: restricted comparator = the other target RCC drugs
     (controls for indication confounding / protopathic bias).
  5. Zero cells: add 0.5 to all 2x2 cells (Haldane-Anscombe) to avoid ROR blow-up.

Usage:
  python 06_signal.py                       # drug level (15 drugs)
  python 06_signal.py --by-mechanism        # mechanism level (5 classes)
  python 06_signal.py --include-lab         # include lab-type PTs
"""
import os
import sys
import io
import time
import argparse

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np
import pandas as pd
from scipy.stats import gamma as gamma_dist

from paths import ROOT, DIR_CLEAN, DIR_RESULTS, SUB_TABLES
CLEAN = os.path.join(ROOT, DIR_CLEAN)
RES = os.path.join(ROOT, DIR_RESULTS, SUB_TABLES)

MIN_N = 3          # minimum report-count threshold


def load_inputs(include_lab=False):
    """Load drug / reac (BLSD-restricted) / demo and return a wide table."""
    drug = pd.read_parquet(os.path.join(CLEAN, "drug_mapped.parquet"))
    reac = pd.read_parquet(os.path.join(CLEAN, "reac_norm.parquet"))
    demo = pd.read_parquet(os.path.join(CLEAN, "demo_dedup.parquet"))

    # BLSD PT restriction
    blsd = pd.read_csv(os.path.join(CLEAN, "blsd_pt.csv"))
    keep_pt = set(blsd["pt_std"])
    if include_lab:
        lab = pd.read_csv(os.path.join(CLEAN, "blsd_pt_lab.csv"))
        keep_pt |= set(lab["pt_std"])
    reac = reac[reac["pt_std"].isin(keep_pt)]

    # keep only de-duplicated reports
    keep_ids = set(demo["primaryid"])
    drug = drug[drug["primaryid"].isin(keep_ids)]
    reac = reac[reac["primaryid"].isin(keep_ids)]

    # keep only PS/SS (primary/secondary suspect) drugs, exclude concomitant
    drug = drug[drug["role_cod"].str.upper().isin(["PS", "SS"])]

    return drug, reac, demo


def _ror_prr(a, b, c, d):
    """ROR + 95% CI + PRR + Yates chi-square. a/b/c/d are Series, already Haldane-corrected."""
    ror = (a * d) / (b * c)
    se = np.sqrt(1 / a + 1 / b + 1 / c + 1 / d)
    ror_low = np.exp(np.log(ror) - 1.96 * se)
    ror_high = np.exp(np.log(ror) + 1.96 * se)

    prr = (a / (a + b)) / (c / (c + d))
    N = a + b + c + d
    num = (np.abs(a * d - b * c) - N / 2.0) ** 2
    den = (a + b) * (c + d) * (a + c) * (b + d)
    chi = (N * num / den).where(den > 0, 0.0)
    return ror, ror_low, ror_high, prr, chi


def _ic_ebgm(a, b, c, d, n_obs):
    """BCPNN IC + IC025 (Norén approximation); GPS EBGM + EB05 (Gamma(1,1) conjugate prior)."""
    N = a + b + c + d
    E = (a + b) * (a + c) / N                 # expected value
    ic = np.log2(a / E)
    # Norén 2006 approximation of the IC025 lower bound (based on raw observation count n_obs)
    n = n_obs.astype(float)
    ic025 = ic - 3.3 * (n ** -0.5) - 2.0 * (n ** -1.5)

    alpha, beta = 1.0, 1.0
    ebgm = (n_obs + alpha) / (E + beta)
    eb05 = gamma_dist.ppf(0.05, a=(n_obs + alpha).astype(float),
                          scale=1.0 / (E + beta))
    return ic, ic025, ebgm, eb05


def run(drug, reac, by_mechanism=False, verbose=True):
    group_col = "mechanism" if by_mechanism else "drug_std"
    if group_col not in drug.columns:
        raise KeyError("drug_mapped is missing column %s; actual columns: %s" % (group_col, list(drug.columns)))

    merged = drug[["primaryid", group_col]].merge(
        reac[["primaryid", "pt_std", "pt_disp"]], on="primaryid", how="inner")
    # count each (drug, PT) within a report only once
    merged = merged.drop_duplicates(subset=["primaryid", group_col, "pt_std"])

    N = len(merged)
    comb = (merged.groupby([group_col, "pt_std"]).size()
                  .rename("n11").reset_index())
    g_tot = merged.groupby(group_col).size().rename("n1_")
    p_tot = merged.groupby("pt_std").size().rename("n_1")

    df = comb.merge(g_tot, on=group_col).merge(p_tot, on="pt_std")
    df["n10"] = df["n1_"] - df["n11"]     # b
    df["n01"] = df["n_1"] - df["n11"]     # c
    df["n00"] = N - df["n11"] - df["n10"] - df["n01"]
    df["n00"] = df["n00"].clip(lower=0)

    # Haldane-Anscombe correction
    a = df["n11"] + 0.5
    b = df["n10"] + 0.5
    c = df["n01"] + 0.5
    d = df["n00"] + 0.5

    df["ror"], df["ror_low"], df["ror_high"], df["prr"], df["chi2"] = _ror_prr(a, b, c, d)
    df["ic"], df["ic025"], df["ebgm"], df["ebgm05"] = _ic_ebgm(a, b, c, d, df["n11"])

    df = df[df["n11"] >= MIN_N].copy()
    df = df.rename(columns={"n11": "n"})

    df["ror_signal"] = (df["ror_low"] > 1)
    df["prr_signal"] = (df["prr"] >= 2) & (df["chi2"] >= 3.841)
    df["ic_signal"] = (df["ic025"] > 0)
    df["ebgm_signal"] = (df["ebgm05"] > 2)
    df["n_algo"] = df[["ror_signal", "prr_signal", "ic_signal", "ebgm_signal"]].sum(axis=1)
    # primary criterion: ROR + IC both positive (consistent with the reference paper)
    df["primary_signal"] = df["ror_signal"] & df["ic_signal"]

    # mechanism name (also carried at drug level for later merging)
    if "mechanism" not in df.columns:
        mech = drug[["drug_std", "mechanism"]].drop_duplicates()
        df = df.merge(mech, on="drug_std", how="left")

    disp = reac[["pt_std", "pt_disp"]].drop_duplicates()
    df = df.merge(disp, on="pt_std", how="left")

    df = df.sort_values(["primary_signal", "ror"], ascending=[False, False])
    if verbose:
        print("  total drug-PT combinations (de-duplicated): %d" % N)
    return df


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--by-mechanism", action="store_true")
    ap.add_argument("--include-lab", action="store_true")
    args = ap.parse_args()

    os.makedirs(RES, exist_ok=True)
    print("[%s] Loading data..." % time.strftime("%H:%M:%S"))
    drug, reac, demo = load_inputs(include_lab=args.include_lab)
    print("  de-duplicated report count: %d" % demo["primaryid"].nunique())
    print("  target drug rows (PS/SS): %d" % len(drug))
    print("  BLSD outcome rows: %d (unique PT %d) %s" %
          (len(reac), reac["pt_std"].nunique(),
           "[with lab-type]" if args.include_lab else ""))

    lvl = "mechanism level" if args.by_mechanism else "drug level"
    print("[%s] Signal detection (%s)..." % (time.strftime("%H:%M:%S"), lvl))
    res = run(drug, reac, by_mechanism=args.by_mechanism)

    tag = ("mechanism" if args.by_mechanism else "drug") + ("_lab" if args.include_lab else "")
    out_csv = os.path.join(RES, "signal_%s.csv" % tag)
    res.to_csv(out_csv, index=False, encoding="utf-8-sig")

    n_sig = int(res["primary_signal"].sum())
    print("  dual-positive (ROR+IC) signals: %d / %d total drug-PT combinations (n>=%d)"
          % (n_sig, len(res), MIN_N))
    print("  Wrote %s" % out_csv)

    gcol = "mechanism" if args.by_mechanism else "drug_std"
    cols = [gcol, "pt_disp", "n", "ror", "ror_low", "ror_high",
            "prr", "chi2", "ic", "ic025", "ebgm", "ebgm05", "primary_signal"]
    print("\n=== Top 20 dual-positive signals by ROR ===")
    print(res[res["primary_signal"]].head(20)[cols].round(2).to_string(index=False))

    print("\n=== Summary by %s ===" % gcol)
    agg = res.groupby(gcol).agg(
        n_combos=("n", "size"),
        dual_positive_signals=("primary_signal", "sum"),
        max_ror=("ror", "max"),
    ).reset_index().sort_values("dual_positive_signals", ascending=False)
    print(agg.round(2).to_string(index=False))


if __name__ == "__main__":
    main()
