# -*- coding: utf-8 -*-
"""
08_weibull.py -- Weibull time-to-event (TTO) modelling (BLSD disease-specific, tests H3)

Fit a Weibull distribution to the time-to-onset (TTO) of "mechanism group x key BLSD PT",
using the shape parameter beta to distinguish mechanistic vs opportunistic events:
  - beta < 1 and 95% CI upper bound < 1 : Early failure  -- on-target acute effect
  - CI includes 1                        : Random failure -- background event (unrelated to time on drug)
  - beta > 1 and 95% CI lower bound > 1 : Wear-out failure -- cumulative toxicity / delayed event

beta's 95% CI uses the standard asymptotic formula for Weibull MLE (standard in pharmacovigilance literature):
    95% CI = beta * exp( ±1.96 * sqrt(0.6079 / n) )

TTO definition: event_dt - start_dt (days). start_dt comes from the THER table (earliest per primaryid).
  Filter 1 <= TTO <= 3650 (Weibull requires > 0; > 10 years treated as anomalous).

[v2 revision | date-granularity bias fix, 2026-09-02]
FAERS start_dt / event_dt have three granularities: YYYYMMDD (full), YYYYMM, YYYY.
The original implementation padded YYYYMM/YYYY with "01" and treated them as exact dates in the
computation, which interpreted "2023" as "2023-01-01", thereby creating the maximum possible TTO
within that year and systematically overestimating time-to-onset for coarse-granularity records.
Measured on the belzutifan x Anaemia stratum:
    8-digit full dates   n=41  median TTO 41 days
    6-digit (year-month) n=36  median TTO 491 days
    4-digit (year only)  n=8   median TTO 332 days
The median TTO of coarse-granularity records is 8-12x that of full-date records, and worsens
monotonically with coarser granularity, confirming this is an artefact of date imputation rather
than a real temporal difference.

The revision adopts a strict criterion (also the prevailing requirement in FAERS TTO literature):
  * both start_dt and event_dt must be full 8-digit dates simultaneously, otherwise the case's TTO is treated as missing
  * start_dt is taken as the earliest value only from THER rows with full dates
  * no date imputation is applied

Outputs:
  05_results/tables/weibull_mechanism.csv       mechanism group x key PT
  05_results/tables/weibull_drug.csv            drug x key PT
  05_results/tables/weibull_mechanism_all.csv   mechanism group x all BLSD events (combined)
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
from scipy.stats import weibull_min

from paths import ROOT, DIR_CLEAN, DIR_RESULTS, SUB_TABLES
CLEAN = os.path.join(ROOT, DIR_CLEAN)
RES = os.path.join(ROOT, DIR_RESULTS, SUB_TABLES)

# key PT list and mechanism order (shared definition, see pt_defs.py)
from pt_defs import KEY_PT_NAMES as KEY_PT, MECH_ORDER

MIN_N = 10          # minimum sample size for Weibull fitting
MAX_TTO = 3650      # days


def to_dt(s: pd.Series) -> pd.Series:
    """Vectorized parsing of FAERS full dates (YYYYMMDD).

    [v2] no longer pad YYYYMM / YYYY with imputation -- that would treat '2023' as
    '2023-01-01', creating the maximum possible TTO within that year and systematically
    overestimating time-to-onset. Coarse-granularity records are uniformly treated as missing.
    """
    s = s.astype(str).str.strip()
    return pd.to_datetime(s.where(s.str.match(r"^\d{8}$").fillna(False)),
                          format="%Y%m%d", errors="coerce")


def build_tto():
    drug = pd.read_parquet(os.path.join(CLEAN, "drug_mapped.parquet"))
    reac = pd.read_parquet(os.path.join(CLEAN, "reac_norm.parquet"))
    demo = pd.read_parquet(os.path.join(CLEAN, "demo_dedup.parquet"))
    ther = pd.read_parquet(os.path.join(CLEAN, "ther.parquet"))

    # BLSD restriction
    blsd = pd.read_csv(os.path.join(CLEAN, "blsd_pt.csv"))
    reac = reac[reac["pt_std"].isin(set(blsd["pt_std"]))]

    keep_ids = set(demo["primaryid"])
    drug = drug[drug["primaryid"].isin(keep_ids)]
    reac = reac[reac["primaryid"].isin(keep_ids)]
    drug = drug[drug["role_cod"].str.upper().isin(["PS", "SS"])]

    # start_dt (THER, earliest per report; [v2] only full 8-digit dates accepted)
    ther = ther[ther["primaryid"].isin(keep_ids)].copy()
    ther["_sd"] = to_dt(ther["start_dt"])
    start_map = ther.dropna(subset=["_sd"]).groupby("primaryid")["_sd"].min()

    demo = demo[["primaryid", "event_dt"]].drop_duplicates("primaryid").copy()
    demo["_ev"] = to_dt(demo["event_dt"])   # [v2] coarse-granularity event_dt treated as missing

    merged = (drug[["primaryid", "drug_std", "mechanism"]]
              .merge(reac[["primaryid", "pt_std", "pt_disp"]], on="primaryid", how="inner")
              .drop_duplicates(subset=["primaryid", "drug_std", "pt_std"]))
    merged = merged.merge(demo[["primaryid", "_ev"]], on="primaryid", how="left")
    merged = merged.merge(start_map.rename("_sd"), on="primaryid", how="left")

    merged["tto"] = (merged["_ev"] - merged["_sd"]).dt.days
    valid = merged["tto"].between(1, MAX_TTO)
    print("  Usable TTO: %d / %d (%.1f%%)" %
          (int(valid.sum()), len(merged), 100.0 * valid.sum() / max(len(merged), 1)))
    # [v2] report-date completeness retention rate, for the Methods statement in the paper
    has_ev = merged["_ev"].notna()
    has_sd = merged["_sd"].notna()
    print("  Date completeness: event full %d (%.1f%%); start full %d (%.1f%%); "
          "both full %d (%.1f%%)" %
          (int(has_ev.sum()), 100 * has_ev.mean(),
           int(has_sd.sum()), 100 * has_sd.mean(),
           int((has_ev & has_sd).sum()), 100 * (has_ev & has_sd).mean()))
    return merged[valid].copy()


def fit_weibull(ttos):
    """Return (beta, beta_low95, beta_high95, scale, median)."""
    x = np.asarray(ttos, dtype=float)
    x = x[x > 0]
    n = len(x)
    if n < MIN_N:
        return np.nan, np.nan, np.nan, np.nan, np.nan
    try:
        c, loc, scale = weibull_min.fit(x, floc=0)
    except Exception:
        return np.nan, np.nan, np.nan, np.nan, np.nan
    if not np.isfinite(c) or c <= 0:
        return np.nan, np.nan, np.nan, np.nan, np.nan
    se_log = np.sqrt(0.6079 / n)                    # standard asymptotic formula
    low = c * np.exp(-1.96 * se_log)
    high = c * np.exp(1.96 * se_log)
    median = scale * (np.log(2) ** (1.0 / c))
    return c, low, high, scale, median


def classify(beta, low, high):
    if not np.isfinite(beta):
        return "unknown"
    if np.isfinite(high) and high < 1:
        return "early_failure"
    if np.isfinite(low) and low > 1:
        return "wearout_failure"
    return "random_failure"


def analyze(merged, by, pt_filter=None, label=""):
    sub = merged
    if pt_filter is not None:
        sub = sub[sub["pt_disp"].isin(pt_filter)]
    rows = []
    for (g, pt), grp in sub.groupby([by, "pt_disp"]):
        if len(grp) < MIN_N:
            continue
        b, lo, hi, sc, med = fit_weibull(grp["tto"].values)
        rows.append({
            by: g, "pt": pt, "n": len(grp),
            "median_tto": float(grp["tto"].median()),
            "beta": b, "beta_low95": lo, "beta_high95": hi,
            "scale": sc, "weibull_median": med,
            "type": classify(b, lo, hi),
        })
    df = pd.DataFrame(rows)
    if not df.empty:
        df = df.sort_values([by, "n"], ascending=[True, False])
    print("  %s: %d groups meet n>=%d" % (label, len(df), MIN_N))
    return df


def main():
    os.makedirs(RES, exist_ok=True)
    print("[%s] Building TTO..." % time.strftime("%H:%M:%S"))
    merged = build_tto()

    print("[%s] Fitting Weibull..." % time.strftime("%H:%M:%S"))
    mres = analyze(merged, "mechanism", KEY_PT, "mechanism x key PT")
    dres = analyze(merged, "drug_std", KEY_PT, "drug x key PT")

    # mechanism group x all BLSD events (combined)
    allrows = []
    for m in MECH_ORDER:
        sub = merged[merged["mechanism"] == m]
        if len(sub) < MIN_N:
            continue
        b, lo, hi, sc, med = fit_weibull(sub["tto"].values)
        allrows.append({"mechanism": m, "pt": "[all BLSD events]", "n": len(sub),
                        "median_tto": float(sub["tto"].median()),
                        "beta": b, "beta_low95": lo, "beta_high95": hi,
                        "scale": sc, "weibull_median": med,
                        "type": classify(b, lo, hi)})
    ares = pd.DataFrame(allrows)

    mres.to_csv(os.path.join(RES, "weibull_mechanism.csv"), index=False, encoding="utf-8-sig")
    dres.to_csv(os.path.join(RES, "weibull_drug.csv"), index=False, encoding="utf-8-sig")
    ares.to_csv(os.path.join(RES, "weibull_mechanism_all.csv"), index=False, encoding="utf-8-sig")

    print("\n=== Weibull typing for mechanism group x all BLSD events (H3 core) ===")
    print(ares[["mechanism", "n", "median_tto", "beta", "beta_low95", "beta_high95", "type"]]
          .round(3).to_string(index=False))

    print("\n=== Mechanism group x key PT (by n descending, top 8 per group) ===")
    if not mres.empty:
        top = mres.sort_values("n", ascending=False).groupby("mechanism").head(8)
        print(top[["mechanism", "pt", "n", "median_tto", "beta",
                   "beta_low95", "beta_high95", "type"]].round(3).to_string(index=False))

    print("\n=== Type distribution (mechanism group x key PT) ===")
    print(mres["type"].value_counts().to_string() if not mres.empty else "(empty)")


if __name__ == "__main__":
    main()
