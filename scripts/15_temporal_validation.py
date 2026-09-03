# -*- coding: utf-8 -*-
"""
15_temporal_validation.py -- FAERS internal time-split validation (signal reproducibility)

Purpose: compensate for the lack of an in-house retrospective cohort and independent external data.
  Use FAERS itself to split into discovery / validation sets and test the reproducibility of
  primary signals.

Design:
    discovery set : 2004 Q1 - 2018 Q4
    validation set: 2019 Q1 - 2025 Q3

    For each mechanism x key-PT combination, compute the four algorithms in both periods;
    the decision rule is identical to 06_signal.py (primary = ROR lower bound > 1 and IC025 > 0).

Output metrics:
    - replication rate (share of discovery-positive signals still positive in validation)
    - direction consistency (whether the ROR point estimate is on the same side)
    - ROR point-estimate stability (ratio of the two periods)

Relation to "external validation":
    Time-split validation cannot replace independent-population validation (shared reporting
    culture and coding system), but it rules out the most common false-positive source -- a
    signal driven by a one-off event in a single period (e.g. a regulatory announcement, publication bias).
    For Drug Safety reviewers this is an acceptable minimum level of reproducibility evidence.

Usage:
    python 15_temporal_validation.py
"""
import io
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import gamma as gamma_dist

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from paths import ROOT, DIR_CLEAN, DIR_RESULTS, SUB_TABLES
CLEAN = os.path.join(ROOT, DIR_CLEAN)
RES = os.path.join(ROOT, DIR_RESULTS, SUB_TABLES)

# period split point (inclusive)
SPLIT_YEAR = 2019  # >= 2019 is the validation set

# pt_std -> pt_disp mapping (filled by main, for display names in signals_for)
pt_disp_map = pd.DataFrame(columns=["pt_std", "pt_disp"])


def qnum(q):
    """'2018q4' -> (2018, 4)."""
    y, n = q.lower().split("q")
    return int(y), int(n)


def period_of(q):
    y, _ = qnum(q)
    return "validation" if y >= SPLIT_YEAR else "discovery"


# ---------------------------------------------------------------- four algorithms
def _ror_prr(a, b, c, d):
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
    N = a + b + c + d
    E = (a + b) * (a + c) / N
    ic = np.log2(a / E)
    n = n_obs.astype(float)
    ic025 = ic - 3.3 * (n ** -0.5) - 2.0 * (n ** -1.5)   # Norén approximation
    alpha, beta = 1.0, 1.0
    ebgm = (n_obs + alpha) / (E + beta)
    eb05 = gamma_dist.ppf(0.05, a=(n_obs + alpha).astype(float), scale=1.0 / (E + beta))
    return ic, ic025, ebgm, eb05


def signals_for(df, group_col="mechanism"):
    """Compute the four algorithms once for a given subset. df: one row per (primaryid, mechanism/drug, PT) combo."""
    out = []
    for grp, sub_g in df.groupby(group_col):
        # numerator: combo count per PT within this group
        n11 = sub_g.groupby("pt_std").size()
        n1_ = len(sub_g)                       # total combos in this group
        # the **total** combo count of this PT (including this group) -- same convention as 06_signal.py,
        # so below n01 = n_1 - n11 is "not this group and this PT"
        n_1 = df.groupby("pt_std").size()
        N = len(df)                            # total combos

        tbl = pd.DataFrame({"n11": n11}).fillna(0)
        tbl["n_1"] = n_1.reindex(tbl.index).fillna(0)
        tbl["n10"] = n1_ - tbl["n11"]
        tbl["n01"] = tbl["n_1"] - tbl["n11"]
        tbl["n00"] = (N - tbl["n11"] - tbl["n10"] - tbl["n01"]).clip(lower=0)

        # Haldane-Anscombe correction
        a = tbl["n11"] + 0.5
        b = tbl["n10"] + 0.5
        c = tbl["n01"] + 0.5
        d = tbl["n00"] + 0.5

        ror, rlo, rhi, prr, chi = _ror_prr(a, b, c, d)
        ic, ic025, ebgm, eb05 = _ic_ebgm(a, b, c, d, tbl["n11"])

        t = pd.DataFrame({
            group_col: grp,
            "pt_std": tbl.index,
            "n": tbl["n11"].astype(int),
            "ror": ror, "ror_low": rlo, "ror_high": rhi,
            "prr": prr, "chi2": chi,
            "ic": ic, "ic025": ic025,
            "ebgm": ebgm, "ebgm05": eb05,
        })
        # keep only combos with actual observations (IC/EBGM undefined when n=0)
        t = t[t["n"] > 0]
        t["primary_signal"] = (t["ror_low"] > 1) & (t["ic025"] > 0)
        out.append(t)
    res = pd.concat(out, ignore_index=True)
    # attach PT display names
    res = res.merge(pt_disp_map, on="pt_std", how="left")
    res["pt_disp"] = res["pt_disp"].fillna(res["pt_std"])
    return res


def main():
    print("=" * 70)
    print("Time-split validation: discovery set vs validation set")
    print("=" * 70)

    drug = pd.read_parquet(os.path.join(CLEAN, "drug_mapped.parquet"))
    reac = pd.read_parquet(os.path.join(CLEAN, "reac_norm.parquet"))
    demo = pd.read_parquet(os.path.join(CLEAN, "demo_dedup.parquet"))

    # keep only key PTs (38, consistent with 09_key_pt_table)
    from pt_defs import KEY_PT_NAMES, MECH_ORDER

    key_std = set(p.lower().strip() for p in KEY_PT_NAMES)
    reac = reac[reac["pt_std"].isin(key_std)]

    # combine into (primaryid, mechanism, pt) combos
    merged = (
        drug[["primaryid", "mechanism"]]
        .merge(reac[["primaryid", "pt_std", "pt_disp"]], on="primaryid", how="inner")
        .drop_duplicates(subset=["primaryid", "mechanism", "pt_std"])
    )

    # attach quarter
    qmap = demo[["primaryid", "_q"]].drop_duplicates("primaryid")
    merged = merged.merge(qmap, on="primaryid", how="left")
    merged = merged[merged["_q"].notna()]
    merged["period"] = merged["_q"].apply(period_of)

    print(f"\nAvailable combos: {len(merged):,}")
    for p in ["discovery", "validation"]:
        sub = merged[merged["period"] == p]
        print(f"  {p:11s}: {len(sub):>8,} combos, {sub['primaryid'].nunique():>7,} reports, "
              f"{sub['_q'].nunique()} quarters")

    global pt_disp_map
    pt_disp_map = (
        reac[["pt_std", "pt_disp"]]
        .dropna(subset=["pt_std"])
        .drop_duplicates("pt_std")
    )

    # ---- compute signals for the two periods separately ----
    res = {}
    for p in ["discovery", "validation"]:
        sub = merged[merged["period"] == p]
        res[p] = signals_for(sub, "mechanism")
        n_sig = int(res[p]["primary_signal"].sum())
        print(f"\n{p}: {len(res[p])} combos, {n_sig} primary signals")

    # ---- merge and compare ----
    d = res["discovery"].rename(columns={
        "n": "n_disc", "ror": "ror_disc", "ror_low": "ror_low_disc",
        "ror_high": "ror_high_disc", "ic025": "ic025_disc",
        "primary_signal": "sig_disc"})
    v = res["validation"].rename(columns={
        "n": "n_val", "ror": "ror_val", "ror_low": "ror_low_val",
        "ror_high": "ror_high_val", "ic025": "ic025_val",
        "primary_signal": "sig_val"})

    cmp = d[["mechanism", "pt_std", "pt_disp", "n_disc", "ror_disc",
             "ror_low_disc", "ror_high_disc", "sig_disc"]].merge(
          v[["mechanism", "pt_std", "n_val", "ror_val",
             "ror_low_val", "ror_high_val", "sig_val"]],
          on=["mechanism", "pt_std"], how="outer")

    cmp["n_disc"] = cmp["n_disc"].fillna(0).astype(int)
    cmp["n_val"] = cmp["n_val"].fillna(0).astype(int)
    cmp["sig_disc"] = cmp["sig_disc"].fillna(False)
    cmp["sig_val"] = cmp["sig_val"].fillna(False)

    # replication flag
    cmp["replicated"] = cmp["sig_disc"] & cmp["sig_val"]
    # direction consistency (point estimates on the same side, both non-missing)
    def same_side(r):
        if pd.isna(r["ror_disc"]) or pd.isna(r["ror_val"]):
            return False
        return (r["ror_disc"] > 1) == (r["ror_val"] > 1)
    cmp["direction_consistent"] = cmp.apply(same_side, axis=1)
    # stability: ratio of the two-period RORs (log is more intuitive)
    cmp["ror_ratio"] = cmp["ror_val"] / cmp["ror_disc"]

    cmp = cmp.sort_values(["mechanism", "ror_disc"], ascending=[True, False])

    os.makedirs(RES, exist_ok=True)
    out_csv = os.path.join(RES, "temporal_validation.csv")
    cmp.to_csv(out_csv, index=False, encoding="utf-8-sig")
    print(f"\n[written] {out_csv}")

    # ---- summary ----
    print("\n" + "=" * 70)
    print("Reproducibility summary")
    print("=" * 70)

    n_disc = int(cmp["sig_disc"].sum())
    n_rep = int(cmp["replicated"].sum())
    print(f"\nDiscovery primary signals: {n_disc}")
    print(f"Still positive in validation: {n_rep}  replication rate = {n_rep/max(n_disc,1):.1%}")

    # direction consistency only over combos with enough n in both periods
    both = cmp[(cmp["n_disc"] >= 5) & (cmp["n_val"] >= 5)]
    print(f"\nCombos with n>=5 in both periods: {len(both)}")
    print(f"  direction consistent: {int(both['direction_consistent'].sum())} "
          f"({both['direction_consistent'].mean():.1%})")

    # replication rate by mechanism
    print("\nReplication by mechanism group:")
    print(f"{'mechanism':<14s} {'disc_pos':>8s} {'val_pos':>8s} {'rate':>8s}")
    for m in MECH_ORDER:
        sub = cmp[cmp["mechanism"] == m]
        nd = int(sub["sig_disc"].sum())
        nr = int(sub["replicated"].sum())
        rate = f"{nr/nd:.0%}" if nd else "—"
        print(f"{m:<14s} {nd:>8d} {nr:>8d} {rate:>8s}")

    # ---- per-signal performance ----
    print("\n" + "=" * 70)
    print("Per-signal validation of discovery-positive signals")
    print("=" * 70)
    sig = cmp[cmp["sig_disc"]].copy()
    sig["status"] = np.where(sig["sig_val"], "✓ replicated",
                     np.where(sig["n_val"] < 5, "? insufficient validation-sample size", "✗ not replicated"))
    show = sig[["mechanism", "pt_disp", "n_disc", "ror_disc",
                "n_val", "ror_val", "ror_ratio", "status"]].copy()
    show["ror_disc"] = show["ror_disc"].round(2)
    show["ror_val"] = show["ror_val"].round(2)
    show["ror_ratio"] = show["ror_ratio"].round(2)
    print(show.to_string(index=False))

    # status distribution
    print("\nStatus distribution:")
    for s, c in sig["status"].value_counts().items():
        print(f"  {s}: {c}")

    # ---- write summary table ----
    summary = pd.DataFrame([{
        "metric": "discovery primary signals", "value": n_disc},
        {"metric": "validation replicated", "value": n_rep},
        {"metric": "replication rate", "value": f"{n_rep/max(n_disc,1):.1%}"},
        {"metric": "combos with n>=5 both periods", "value": len(both)},
        {"metric": "direction-consistency rate", "value": f"{both['direction_consistent'].mean():.1%}"},
        {"metric": "discovery combos", "value": len(res['discovery'])},
        {"metric": "validation combos", "value": len(res['validation'])},
        {"metric": "discovery reports",
         "value": int(merged[merged['period']=='discovery']['primaryid'].nunique())},
        {"metric": "validation reports",
         "value": int(merged[merged['period']=='validation']['primaryid'].nunique())},
    ])
    out_sum = os.path.join(RES, "temporal_validation_summary.csv")
    summary.to_csv(out_sum, index=False, encoding="utf-8-sig")
    print(f"\n[written] {out_sum}")

    print("\nDone.")


if __name__ == "__main__":
    main()
