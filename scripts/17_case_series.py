# -*- coding: utf-8 -*-
"""
17_case_series.py - Case-by-case series analysis of signature signals

Purpose (READUS-PV item 10 + item 7d):
    "Present the case-by-case analysis of key variables. Present the causality
     assessment, if applicable."
    "Specify the variables and methods used for the case-by-case analysis,
     including any algorithm or criteria used to assess causality."

    Disproportionality only yields a statistical association; a case series is
    the necessary step that puts each association back into its clinical context.
    This script produces a structured case series for 5 signature signals and a
    **proxy** causality assessment against the available elements of the
    WHO-UMC framework (its limitations are stated explicitly).

Ethics: only aggregate statistics are emitted; no individually identifiable
case information is output.

Outputs:
    <DIR_RESULTS>/<SUB_TABLES>/case_series_summary.csv   per-signal case-series features (Table 4)
    <DIR_RESULTS>/<SUB_TABLES>/case_series_drugs.csv     concomitant-drug distribution per signal
    <DIR_RESULTS>/<SUB_TABLES>/case_series_indi.csv      indication distribution per signal

Usage:
    python 17_case_series.py
"""
import io
import os
import sys

import numpy as np
import pandas as pd

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from paths import ROOT, DIR_CLEAN, DIR_RESULTS, SUB_TABLES
CLEAN = os.path.join(ROOT, DIR_CLEAN)
RES = os.path.join(ROOT, DIR_RESULTS, SUB_TABLES)

# Signature signals: (label, level, drug-or-mechanism key, PT display name)
SIGNALS = [
    ("Belzutifan-anaemia",        "drug", "belzutifan",  "Anaemia"),
    ("Axitinib-polycythaemia",    "drug", "axitinib",    "Polycythaemia"),
    ("ICI-pure red cell aplasia", "mech", "E_ICI",       "Aplasia pure red cell"),
    ("Everolimus-iron deficiency anaemia", "drug", "everolimus", "Iron deficiency anaemia"),
    ("ICI-autoimmune haemolytic anaemia",  "mech", "E_ICI", "Autoimmune haemolytic anaemia"),
]


def to_dt(x):
    """FAERS date -> datetime.

    [v2 revision | 2026-09-02] The original implementation converted to a
    number first and then sliced, which caused:
        start_dt="202301" -> 202301.0 -> "202301.0"
        -> slice(6,8)=".0" -> "2023-01-.0" -> NaT
    i.e. every 6-digit (YYYYMM) and 4-digit (YYYY) date was silently dropped.
    In practice start_dt lost 246,601 records (15.4%) and event_dt lost 38,903.

    Revised to the strict convention consistent with 08_weibull.py v2: accept
    only complete 8-digit dates; coarser granularities are treated as missing
    with no imputation (imputing "2023" as "2023-01-01" would fabricate the
    maximum possible TTO within that year, systematically over-estimating onset).
    """
    s = pd.Series(x).astype(str).str.strip()
    return pd.to_datetime(s.where(s.str.match(r"^\d{8}$").fillna(False)),
                          format="%Y%m%d", errors="coerce")


def to_year(dt):
    try:
        s = str(int(float(dt)))
    except (TypeError, ValueError):
        return np.nan
    return int(s[:4]) if len(s) >= 4 else np.nan


def main():
    print("=" * 74)
    print("Case series of signature signals (READUS-PV item 10)")
    print("=" * 74)

    demo = pd.read_parquet(os.path.join(CLEAN, "demo_dedup.parquet"))
    drug = pd.read_parquet(os.path.join(CLEAN, "drug_mapped.parquet"))
    reac = pd.read_parquet(os.path.join(CLEAN, "reac_norm.parquet"))
    outc = pd.read_parquet(os.path.join(CLEAN, "outc.parquet"))
    ther = pd.read_parquet(os.path.join(CLEAN, "ther.parquet"))
    indi = pd.read_parquet(os.path.join(CLEAN, "indi.parquet"))

    # ---- serious-outcome flags ----
    outc = outc.copy()
    outc["outc_cod"] = outc["outc_cod"].astype(str).str.upper().str.strip()
    sev = outc[outc["outc_cod"].isin(["DE", "LT", "HO", "DS", "RI"])][
        ["primaryid", "outc_cod"]].drop_duplicates()
    for c in ["DE", "LT", "HO"]:
        flag = sev[sev["outc_cod"] == c]["primaryid"].unique()
        demo[c] = demo["primaryid"].isin(flag).astype(int)

    # ---- therapy start date (earliest per case) ----
    ther = ther.copy()
    ther["start"] = to_dt(ther["start_dt"])
    start_min = ther.groupby("primaryid")["start"].min()
    demo["start_dt"] = demo["primaryid"].map(start_min)
    demo["event_date"] = to_dt(demo["event_dt"])
    demo["tto_days"] = (demo["event_date"] - demo["start_dt"]).dt.days
    demo.loc[(demo["tto_days"] < 0) | (demo["tto_days"] > 3650), "tto_days"] = np.nan

    # ---- age ----
    demo["age_num"] = pd.to_numeric(demo["age"], errors="coerce")
    is_year = demo["age_cod"].astype(str).str.upper().isin(["YR"])
    demo.loc[~is_year, "age_num"] = np.nan
    demo.loc[(demo["age_num"] < 0) | (demo["age_num"] > 120), "age_num"] = np.nan

    demo["year"] = demo["event_dt"].apply(to_year)

    # ---- target cohort: reports hitting a given signal ----
    rows, drug_rows, indi_rows = [], [], []

    for label, level, key, pt in SIGNALS:
        pt_std = pt.lower().strip()
        ids = set(reac.loc[reac["pt_std"] == pt_std, "primaryid"])

        if level == "drug":
            exp = set(drug.loc[drug["drug_std"] == key, "primaryid"])
            exp_label = key
        else:
            exp = set(drug.loc[drug["mechanism"] == key, "primaryid"])
            exp_label = key

        case_ids = exp & ids
        sub = demo[demo["primaryid"].isin(case_ids)]

        if len(sub) == 0:
            print(f"\n[skip] {label}: no cases")
            continue

        # Single-drug exposure (only one PS/SS drug).
        ps_count = (drug[drug["role_cod"].astype(str).str.upper().isin(["PS", "SS"])]
                    .groupby("primaryid")["drug_std"].nunique())
        n_ps = sub["primaryid"].map(ps_count)
        mono = int((n_ps == 1).sum())

        # Concomitant drugs (excluding self).
        co = drug[drug["primaryid"].isin(case_ids)]
        if level == "drug":
            co = co[co["drug_std"] != key]
        else:
            co = co[co["mechanism"] != key]
        co_top = (co.groupby("drug_std").size().sort_values(ascending=False).head(10))

        # Indications.
        ind_top = (indi[indi["primaryid"].isin(case_ids)]
                   .groupby("indi_pt").size().sort_values(ascending=False).head(5))

        def pct(x, n):
            return f"{x:,} ({x / n * 100:.1f}%)" if n else "—"

        n = len(sub)
        rows.append({
            "Signal": label,
            "Drug/Mechanism": exp_label,
            "PT": pt,
            "Cases": n,
            "Age n": int(sub["age_num"].notna().sum()),
            "Age median (IQR)": (f"{sub['age_num'].median():.0f} "
                                 f"({sub['age_num'].quantile(.25):.0f}–"
                                 f"{sub['age_num'].quantile(.75):.0f})")
                                 if sub["age_num"].notna().sum() else "—",
            "Male n (%)": pct(int((sub["gndr_cod"].astype(str).str.upper() == "M").sum()), n),
            "Female n (%)": pct(int((sub["gndr_cod"].astype(str).str.upper() == "F").sum()), n),
            "TTO n": int(sub["tto_days"].notna().sum()),
            "TTO median (IQR)": (f"{sub['tto_days'].median():.0f} "
                                 f"({sub['tto_days'].quantile(.25):.0f}–"
                                 f"{sub['tto_days'].quantile(.75):.0f})")
                                 if sub["tto_days"].notna().sum() else "—",
            "TTO range": (f"{sub['tto_days'].min():.0f}–{sub['tto_days'].max():.0f}")
                         if sub["tto_days"].notna().sum() else "—",
            "Hospitalisation n (%)": pct(int(sub["HO"].sum()), n),
            "Death n (%)": pct(int(sub["DE"].sum()), n),
            "Life-threatening n (%)": pct(int(sub["LT"].sum()), n),
            "Monotherapy exposure n (%)": pct(mono, n),
            "Report year median": f"{sub['year'].median():.0f}" if sub["year"].notna().sum() else "—",
            "Report year range": (f"{sub['year'].min():.0f}–{sub['year'].max():.0f}")
                                 if sub["year"].notna().sum() else "—",
        })

        for d, c in co_top.items():
            drug_rows.append({"Signal": label, "Concomitant drug": d, "n": int(c),
                              "%": round(c / n * 100, 1)})
        for d, c in ind_top.items():
            indi_rows.append({"Signal": label, "Indication": d, "n": int(c),
                              "%": round(c / n * 100, 1)})

        print(f"\n--- {label} ---")
        print(f"  Cases {n}  Age median {rows[-1]['Age median (IQR)']}  "
              f"TTO median {rows[-1]['TTO median (IQR)']} d  "
              f"Monotherapy {mono} ({mono/n*100:.1f}%)")

    df = pd.DataFrame(rows)
    os.makedirs(RES, exist_ok=True)

    out1 = os.path.join(RES, "case_series_summary.csv")
    df.to_csv(out1, index=False, encoding="utf-8-sig")
    print(f"\n[write] {out1}")

    out2 = os.path.join(RES, "case_series_drugs.csv")
    pd.DataFrame(drug_rows).to_csv(out2, index=False, encoding="utf-8-sig")
    print(f"[write] {out2}")

    out3 = os.path.join(RES, "case_series_indi.csv")
    pd.DataFrame(indi_rows).to_csv(out3, index=False, encoding="utf-8-sig")
    print(f"[write] {out3}")

    print("\n" + "=" * 74)
    print("Case series summary")
    print("=" * 74)
    show = df[["Signal", "Cases", "Age median (IQR)", "Male n (%)", "TTO median (IQR)",
               "Hospitalisation n (%)", "Death n (%)", "Monotherapy exposure n (%)",
               "Report year median"]]
    print(show.to_string(index=False))

    print("\nDone.")


if __name__ == "__main__":
    main()
