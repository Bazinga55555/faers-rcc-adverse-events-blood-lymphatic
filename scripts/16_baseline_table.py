# -*- coding: utf-8 -*-
"""
16_baseline_table.py - Table 1: baseline characteristics by mechanistic class

Purpose (READUS-PV item 8b):
    "Provide key demographic and clinical characteristics of cases,
     if possible comparing cases with any appropriate reference group."

    This table operationalises that item. Reference group = the other
    mechanistic classes (each class is compared against the rest).

Outputs:
    <DIR_RESULTS>/<SUB_TABLES>/table1_baseline.csv         long format (easy to reshape)
    <DIR_RESULTS>/<SUB_TABLES>/table1_baseline_wide.csv    wide format (ready for Table 1)
    <DIR_RESULTS>/<SUB_TABLES>/table1_baseline.md          Markdown version (paste-ready)

Usage:
    python 16_baseline_table.py
"""
import io
import os
import sys

import numpy as np
import pandas as pd

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from pt_defs import MECH_ORDER, MECH_LABEL

from paths import ROOT, DIR_CLEAN, DIR_RESULTS, SUB_TABLES
CLEAN = os.path.join(ROOT, DIR_CLEAN)
RES = os.path.join(ROOT, DIR_RESULTS, SUB_TABLES)

# ---------------------------------------------------------------- region map
REGION = {
    # North America
    "US": "North America", "CA": "North America", "MX": "North America",
    # Europe (Western + Northern + Southern + Eastern)
    "GB": "Europe", "DE": "Europe", "FR": "Europe", "IT": "Europe",
    "ES": "Europe", "NL": "Europe", "BE": "Europe", "SE": "Europe",
    "NO": "Europe", "DK": "Europe", "FI": "Europe", "IE": "Europe",
    "AT": "Europe", "CH": "Europe", "PT": "Europe", "GR": "Europe",
    "PL": "Europe", "CZ": "Europe", "HU": "Europe", "RO": "Europe",
    "SK": "Europe", "SI": "Europe", "HR": "Europe", "BG": "Europe",
    "RS": "Europe", "LT": "Europe", "LV": "Europe", "EE": "Europe",
    "LU": "Europe", "IS": "Europe", "MT": "Europe", "CY": "Europe",
    "RU": "Europe", "TR": "Europe", "UA": "Europe", "BY": "Europe",
    # Asia
    "JP": "Asia", "CN": "Asia", "IN": "Asia", "KR": "Asia",
    "TW": "Asia", "HK": "Asia", "SG": "Asia", "TH": "Asia",
    "MY": "Asia", "ID": "Asia", "PH": "Asia", "VN": "Asia",
    "PK": "Asia", "BD": "Asia", "LK": "Asia", "NP": "Asia",
    # Latin America
    "BR": "Latin America", "AR": "Latin America", "CL": "Latin America",
    "CO": "Latin America", "PE": "Latin America", "VE": "Latin America",
    "EC": "Latin America", "UY": "Latin America", "PY": "Latin America",
    "BO": "Latin America", "CR": "Latin America", "PA": "Latin America",
    "GT": "Latin America", "DO": "Latin America", "CU": "Latin America",
    # Middle East / Africa
    "SA": "Middle East / Africa", "AE": "Middle East / Africa",
    "IL": "Middle East / Africa", "EG": "Middle East / Africa",
    "ZA": "Middle East / Africa", "NG": "Middle East / Africa",
    "KE": "Middle East / Africa", "MA": "Middle East / Africa",
    "JO": "Middle East / Africa", "LB": "Middle East / Africa",
    "KW": "Middle East / Africa", "QA": "Middle East / Africa",
    # Oceania
    "AU": "Oceania", "NZ": "Oceania",
}

# Reporter occupation (FAERS occp_cod)
OCCP = {
    "MD": "Physician", "DO": "Physician",
    "PH": "Pharmacist",
    "CN": "Consumer/non-HCP", "LW": "Lawyer",
    "OT": "Other HCP", "RN": "Other HCP", "HP": "Other HCP",
    "VET": "Other HCP", "DEN": "Other HCP",
}

# Outcome (FAERS outc_cod)
OUTC_LABEL = {
    "DE": "Death", "LT": "Life-threatening", "HO": "Hospitalisation",
    "DS": "Disability", "CA": "Congenital anomaly",
    "RI": "Required intervention", "OT": "Other serious",
}


def to_region(c):
    if not isinstance(c, str):
        return "Not specified"
    c = c.strip().upper()
    if c in ("", "COUNTRY NOT SPECIFIED", "UNKNOWN", "NS"):
        return "Not specified"
    return REGION.get(c, "Other")


def to_year(dt):
    """event_dt may be YYYYMMDD / YYYYMM / YYYY"""
    try:
        s = str(int(float(dt)))
    except (TypeError, ValueError):
        return np.nan
    if len(s) >= 4:
        return int(s[:4])
    return np.nan


def main():
    print("=" * 70)
    print("Table 1: baseline characteristics by mechanistic class")
    print("=" * 70)

    demo = pd.read_parquet(os.path.join(CLEAN, "demo_dedup.parquet"))
    drug = pd.read_parquet(os.path.join(CLEAN, "drug_mapped.parquet"))
    outc = pd.read_parquet(os.path.join(CLEAN, "outc.parquet"))

    # Consistent with signal analysis: keep only primary/secondary suspects
    # (PS/SS); exclude concomitant (C) and interacting (I).
    drug = drug[drug["role_cod"].str.upper().isin(["PS", "SS"])]

    # Report -> mechanism (a report may span several mechanisms and is thus
    # counted once per mechanism; this must be footnoted).
    rep_mech = drug[["primaryid", "mechanism"]].drop_duplicates()
    demo = demo.merge(rep_mech, on="primaryid", how="inner")

    n_rep_total = demo["primaryid"].nunique()
    print(f"\nReport-mechanism pairs included: {len(demo):,}")
    print(f"Unique reports: {n_rep_total:,}")

    # ---------------- pre-process fields ----------------
    # Age: keep records expressed in years only; drop out-of-range values.
    demo["age_num"] = pd.to_numeric(demo["age"], errors="coerce")
    is_year = demo["age_cod"].astype(str).str.upper().isin(["YR", "YEAR", "YEARS"])
    demo.loc[~is_year, "age_num"] = np.nan
    demo.loc[(demo["age_num"] < 0) | (demo["age_num"] > 120), "age_num"] = np.nan

    demo["region"] = demo["reporter_country"].apply(to_region)
    demo["occupation"] = demo["occp_cod"].astype(str).str.upper().map(OCCP).fillna("Unknown")
    demo["year"] = demo["event_dt"].apply(to_year)

    # Serious outcomes: a report may carry several outc rows; expand + dedupe.
    outc = outc.copy()
    outc["outc_cod"] = outc["outc_cod"].astype(str).str.upper().str.strip()
    severe = outc[outc["outc_cod"].isin(OUTC_LABEL)][["primaryid", "outc_cod"]].drop_duplicates()
    sev_flags = pd.get_dummies(
        severe.assign(v=1).pivot_table(index="primaryid", columns="outc_cod",
                                       values="v", aggfunc="max", fill_value=0)
    ).reset_index() if len(severe) else pd.DataFrame(columns=["primaryid"])
    demo = demo.merge(sev_flags, on="primaryid", how="left")
    for c in ["DE", "LT", "HO", "DS", "CA", "RI", "OT"]:
        if c not in demo.columns:
            demo[c] = 0
        demo[c] = demo[c].fillna(0).astype(int)

    rows = []

    def add_row(variable, level, func):
        cells = {}
        for m in MECH_ORDER + ["ALL"]:
            sub = demo if m == "ALL" else demo[demo["mechanism"] == m]
            cells[m] = func(sub)
        rows.append({"Variable": variable, "Stratum": level, **cells})

    # ---- number of reports ----
    add_row("Reports", "n", lambda s: f"{s['primaryid'].nunique():,}")

    # ---- age ----
    add_row("Age (years)", "n available", lambda s: f"{s['age_num'].notna().sum():,}")
    add_row("Age (years)", "Median (IQR)",
            lambda s: (f"{s['age_num'].median():.0f} "
                       f"({s['age_num'].quantile(.25):.0f}–{s['age_num'].quantile(.75):.0f})")
            if s["age_num"].notna().sum() > 0 else "—")
    for lo, hi, lbl in [(0, 65, "<65"), (65, 75, "65–74"), (75, 200, "≥75")]:
        add_row("Age group (years)", lbl,
                lambda s, lo=lo, hi=hi: (
                    f"{((s['age_num'] >= lo) & (s['age_num'] < hi)).sum():,} "
                    f"({((s['age_num'] >= lo) & (s['age_num'] < hi)).mean() * 100:.1f}%)"
                    if s["age_num"].notna().sum() > 0 else "—"))

    # ---- sex ----
    g = demo["gndr_cod"].astype(str).str.upper().str.strip()
    for code, lbl in [("M", "Male"), ("F", "Female"), ("UNK", "Unknown"), ("NS", "Not stated")]:
        add_row("Sex", lbl,
                lambda s, code=code: f"{(s['gndr_cod'].astype(str).str.upper().str.strip() == code).sum():,} "
                                     f"({(s['gndr_cod'].astype(str).str.upper().str.strip() == code).mean() * 100:.1f}%)")

    # ---- region ----
    for r in ["North America", "Europe", "Asia", "Latin America",
              "Middle East / Africa", "Oceania", "Other", "Not specified"]:
        add_row("Report region", r,
                lambda s, r=r: f"{(s['region'] == r).sum():,} ({(s['region'] == r).mean() * 100:.1f}%)")

    # ---- reporter ----
    for o in ["Physician", "Pharmacist", "Other HCP", "Consumer/non-HCP", "Lawyer", "Unknown"]:
        add_row("Reporter", o,
                lambda s, o=o: f"{(s['occupation'] == o).sum():,} ({(s['occupation'] == o).mean() * 100:.1f}%)")

    # ---- serious outcomes ----
    for c in ["DE", "LT", "HO", "DS", "RI"]:
        if c not in demo.columns:
            continue
        add_row("Serious outcome", OUTC_LABEL[c],
                lambda s, c=c: f"{int(s[c].sum()):,} ({s[c].mean() * 100:.1f}%)")

    # ---- report year ----
    add_row("Report year", "Median (IQR)",
            lambda s: (f"{s['year'].median():.0f} "
                       f"({s['year'].quantile(.25):.0f}–{s['year'].quantile(.75):.0f})")
            if s["year"].notna().sum() > 0 else "—")
    for lo, hi, lbl in [(2004, 2015, "2004–2014"), (2015, 2021, "2015–2020"), (2021, 2030, "2021–2026")]:
        add_row("Report year", lbl,
                lambda s, lo=lo, hi=hi: (
                    f"{((s['year'] >= lo) & (s['year'] < hi)).sum():,} "
                    f"({((s['year'] >= lo) & (s['year'] < hi)).mean() * 100:.1f}%)"
                    if s["year"].notna().sum() > 0 else "—"))

    # ---- assemble ----
    df = pd.DataFrame(rows)
    df = df.rename(columns={m: f"{m}" for m in MECH_ORDER + ["ALL"]})

    os.makedirs(RES, exist_ok=True)
    out_csv = os.path.join(RES, "table1_baseline.csv")
    df.to_csv(out_csv, index=False, encoding="utf-8-sig")
    print(f"\n[write] {out_csv}")

    # ---- wide table (paper Table 1) ----
    wide = df.copy()
    wide.columns = ["Variable", "Stratum"] + [MECH_LABEL.get(c, c) for c in MECH_ORDER] + ["All"]
    out_wide = os.path.join(RES, "table1_baseline_wide.csv")
    wide.to_csv(out_wide, index=False, encoding="utf-8-sig")
    print(f"[write] {out_wide}")

    # ---- Markdown ----
    lines = ["**Table 1** Baseline characteristics of included reports, by mechanistic class",
             "",
             "| Characteristic | " + " | ".join(
                 [MECH_LABEL.get(m, m) for m in MECH_ORDER] + ["All reports"]) + " |",
             "|---|" + "---|" * (len(MECH_ORDER) + 1)]
    for _, r in df.iterrows():
        cells = [r[m] for m in MECH_ORDER] + [r["ALL"]]
        lines.append(f"| {r['Variable']} — {r['Stratum']} | " + " | ".join(cells) + " |")
    lines += ["",
              "Data are n (%) unless stated otherwise. Percentages are calculated within column. "
              "A report listing drugs from more than one mechanistic class is counted once in each "
              "relevant column; column totals therefore exceed the number of unique reports.",
              "IQR = interquartile range; HCP = health care professional."]
    out_md = os.path.join(RES, "table1_baseline.md")
    with open(out_md, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))
    print(f"[write] {out_md}")

    print("\n" + "=" * 70)
    print(wide.to_string(index=False))
    print("\nDone.")


if __name__ == "__main__":
    main()
