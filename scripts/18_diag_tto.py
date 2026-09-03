# -*- coding: utf-8 -*-
"""
18_diag_tto.py - Diagnostic for the median_tto discrepancy (08_weibull vs 17_case_series)

Issue: the two scripts implement to_dt differently
  - 08: string parsing, YYYYMMDD -> exact; YYYYMM -> pad 01; YYYY -> pad 0101
  - 17: convert to number first then slice, which yields strings like "202301.0"
        for 6-/4-digit dates -> NaT

If 17 drops a large share of coarse-grained start_dt values while 08 treats YYYY
as Jan 1 (i.e. the maximum possible TTO within that year), 08's median is
systematically over-estimated.

This script decomposes the problem layer by layer:
  A. parse-success rate of the two to_dt implementations across start_dt formats
  B. belzutifan x Anaemia under 4 conventions: n / median / mean
  C. same for each mechanism x Anaemia
Only aggregate statistics are emitted; no individual records are output.
"""
import io
import os
import sys

import numpy as np
import pandas as pd

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from paths import ROOT, DIR_CLEAN
CLEAN = os.path.join(ROOT, DIR_CLEAN)

MECH_ORDER = ["A_HIF2a", "B_VEGFR_TKI", "C_VEGF_mAb", "D_mTOR", "E_ICI"]


# ---------- two date parsers ----------
def to_dt_str(s: pd.Series) -> pd.Series:
    """08_weibull.py implementation: string parsing + padding"""
    s = s.astype(str).str.strip()
    out = pd.to_datetime(s, format="%Y%m%d", errors="coerce")
    m = out.isna() & s.str.match(r"^\d{6}$").fillna(False)
    if m.any():
        out.loc[m] = pd.to_datetime(s[m] + "01", format="%Y%m%d", errors="coerce")
    m2 = out.isna() & s.str.match(r"^\d{4}$").fillna(False)
    if m2.any():
        out.loc[m2] = pd.to_datetime(s[m2] + "0101", format="%Y%m%d", errors="coerce")
    return out


def to_dt_num(x):
    """17_case_series.py implementation: convert to number then slice"""
    s = pd.to_numeric(x, errors="coerce")
    s = s.astype("Float64").astype(str)
    y = s.str.slice(0, 4)
    m = s.str.slice(4, 6).replace({"": "01", "nan": None, "<NA>": None})
    d = s.str.slice(6, 8).replace({"": "01", "nan": None, "<NA>": None})
    return pd.to_datetime(
        y.astype(str) + "-" + m.fillna("01") + "-" + d.fillna("01"),
        errors="coerce",
    )


def fmt_bucket(s: pd.Series) -> pd.Series:
    """Bucket raw date strings by length to gauge granularity"""
    t = s.astype(str).str.strip()
    b = pd.Series("other", index=s.index)
    b[t.str.match(r"^\d{8}$").fillna(False)] = "8-digit (exact day)"
    b[t.str.match(r"^\d{6}$").fillna(False)] = "6-digit (year-month)"
    b[t.str.match(r"^\d{4}$").fillna(False)] = "4-digit (year)"
    b[t.isin(["nan", "None", "<NA>", ""])] = "Missing"
    return b


def main():
    print("=" * 78)
    print("TTO convention diagnostic: 08_weibull vs 17_case_series")
    print("=" * 78)

    drug = pd.read_parquet(os.path.join(CLEAN, "drug_mapped.parquet"))
    reac = pd.read_parquet(os.path.join(CLEAN, "reac_norm.parquet"))
    demo = pd.read_parquet(os.path.join(CLEAN, "demo_dedup.parquet"))
    ther = pd.read_parquet(os.path.join(CLEAN, "ther.parquet"))
    blsd = pd.read_csv(os.path.join(CLEAN, "blsd_pt.csv"))

    reac = reac[reac["pt_std"].isin(set(blsd["pt_std"]))]
    keep = set(demo["primaryid"])
    drug = drug[drug["primaryid"].isin(keep)]
    reac = reac[reac["primaryid"].isin(keep)]
    drug = drug[drug["role_cod"].astype(str).str.upper().isin(["PS", "SS"])]

    # ---------- A. parse-success rate ----------
    print("\n[A] format distribution of start_dt / event_dt and parse success of both methods")
    print("-" * 78)
    for col, df in [("start_dt", ther), ("event_dt", demo)]:
        raw = df[col]
        bk = fmt_bucket(raw)
        tab = pd.DataFrame({
            "n": bk.value_counts(),
            "08_str parsed": to_dt_str(raw).notna().groupby(bk).sum(),
            "17_num parsed": to_dt_num(raw).notna().groupby(bk).sum(),
        }).fillna(0).astype(int)
        tab["17/08 success ratio"] = (tab["17_num parsed"] /
                                     tab["08_str parsed"].replace(0, np.nan)).round(3)
        print(f"\n-- {col} --")
        print(tab.to_string())

    # ---------- build merged (identical to 08) ----------
    merged = (drug[["primaryid", "drug_std", "mechanism"]]
              .merge(reac[["primaryid", "pt_std", "pt_disp"]], on="primaryid", how="inner")
              .drop_duplicates(subset=["primaryid", "drug_std", "pt_std"]))

    demo2 = demo[["primaryid", "event_dt"]].copy()
    demo2["ev_str"] = to_dt_str(demo2["event_dt"])
    demo2["ev_num"] = to_dt_num(demo2["event_dt"])

    ther2 = ther[ther["primaryid"].isin(keep)].copy()
    ther2["sd_str"] = to_dt_str(ther2["start_dt"])
    ther2["sd_num"] = to_dt_num(ther2["start_dt"])
    # Record granularity of each case's start_dt (bucket of its earliest record).
    ther2["_bk"] = fmt_bucket(ther2["start_dt"])
    bk_rank = {"8-digit (exact day)": 3, "6-digit (year-month)": 2,
               "4-digit (year)": 1, "other": 0, "Missing": 0}
    ther2["_rk"] = ther2["_bk"].map(bk_rank)
    gran = ther2.groupby("primaryid")["_rk"].min()   # smaller = coarser

    # === strict convention: earliest start only from 8-digit THER rows ===
    ther8 = ther2[ther2["_rk"] == 3]
    sm_strict = ther8.dropna(subset=["sd_str"]).groupby("primaryid")["sd_str"].min()

    sm_str = ther2.dropna(subset=["sd_str"]).groupby("primaryid")["sd_str"].min()
    sm_num = ther2.dropna(subset=["sd_num"]).groupby("primaryid")["sd_num"].min()

    m = merged.merge(demo2[["primaryid", "ev_str", "ev_num"]], on="primaryid", how="left")
    m = m.merge(sm_str.rename("sd_str"), on="primaryid", how="left")
    m = m.merge(sm_num.rename("sd_num"), on="primaryid", how="left")
    m = m.merge(sm_strict.rename("sd_st8"), on="primaryid", how="left")
    m["gran"] = m["primaryid"].map(gran)

    m["tto_str"] = (m["ev_str"] - m["sd_str"]).dt.days
    m["tto_num"] = (m["ev_num"] - m["sd_num"]).dt.days
    # strict convention: start from 8-digit rows; event accepts 8-digit only
    ev8 = demo2[demo2["event_dt"].astype(str).str.match(r"^\d{8}$").fillna(False)]
    ev8map = ev8.drop_duplicates("primaryid").set_index("primaryid")["ev_str"]
    m["ev_st8"] = m["primaryid"].map(ev8map)
    m["tto_strict"] = (m["ev_st8"] - m["sd_st8"]).dt.days

    def summ(x, tag):
        x = x.dropna()
        if len(x) == 0:
            return f"{tag}: n=0"
        return (f"{tag}: n={len(x):4d}  median={x.median():7.1f}  "
                f"mean={x.mean():7.1f}  IQR=({x.quantile(.25):.0f}-{x.quantile(.75):.0f})")

    # ---------- B. belzutifan x Anaemia ----------
    print("\n" + "=" * 78)
    print("[B] belzutifan x Anaemia (08 reports median_tto=116.5, n=86; 17 reports 40.0, n=85)")
    print("=" * 78)
    b = m[(m["drug_std"] == "belzutifan") & (m["pt_std"] == "anaemia")]

    print("\n-- convention 1: 08 original (str parse, tto in [1,3650]) --")
    v = b["tto_str"][b["tto_str"].between(1, 3650)]
    print("  " + summ(v, "result"))

    print("\n-- convention 2: switch to 17's num parse (else identical to 08) --")
    v = b["tto_num"][b["tto_num"].between(1, 3650)]
    print("  " + summ(v, "result"))

    print("\n-- convention 3: 08 parse but drop tto=0 (>=0 filter, same lower bound as 17) --")
    v = b["tto_str"][b["tto_str"].between(0, 3650)]
    print("  " + summ(v, "result"))

    print("\n-- stratified by start_dt granularity (08 str parse) --")
    lab = {3: "8-digit (exact day)", 2: "6-digit (year-month)", 1: "4-digit (year)", 0: "other/Missing"}
    for g, sub in b.groupby("gran"):
        v = sub["tto_str"][sub["tto_str"].between(1, 3650)]
        print("  " + summ(v, f"{lab.get(g, g):<12}"))

    print("\n-- keep exact-day only (gran==3) --")
    v = b.loc[b["gran"] == 3, "tto_str"]
    v = v[v.between(1, 3650)]
    print("  " + summ(v, "result"))

    print("\n-- [recommended] strict convention: both start and event must be 8-digit --")
    v = b["tto_strict"][b["tto_strict"].between(1, 3650)]
    print("  " + summ(v, "result"))

    # ---------- C. each mechanism x Anaemia ----------
    print("\n" + "=" * 78)
    print("[C] each mechanism x Anaemia: 4 conventions compared")
    print("=" * 78)
    a = m[m["pt_std"] == "anaemia"]
    rows = []
    for mech in MECH_ORDER:
        sub = a[a["mechanism"] == mech]
        v1 = sub["tto_str"][sub["tto_str"].between(1, 3650)]
        v2 = sub["tto_num"][sub["tto_num"].between(1, 3650)]
        v3 = sub.loc[sub["gran"] == 3, "tto_str"]
        v3 = v3[v3.between(1, 3650)]
        v4 = sub["tto_strict"][sub["tto_strict"].between(1, 3650)]
        rows.append({
            "Mechanism": mech,
            "08 n": len(v1), "08 median": round(v1.median(), 1) if len(v1) else np.nan,
            "17 n": len(v2), "17 median": round(v2.median(), 1) if len(v2) else np.nan,
            "exact-day n": len(v3), "exact-day median": round(v3.median(), 1) if len(v3) else np.nan,
            "strict n": len(v4), "strict median": round(v4.median(), 1) if len(v4) else np.nan,
            "strict IQR": (f"{v4.quantile(.25):.0f}-{v4.quantile(.75):.0f}") if len(v4) else "—",
        })
    print(pd.DataFrame(rows).to_string(index=False))

    # ---------- D. mechanism x all BLSD ----------
    print("\n" + "=" * 78)
    print("[D] each mechanism x all BLSD events: 08 convention vs exact-day convention")
    print("=" * 78)
    rows = []
    for mech in MECH_ORDER:
        sub = m[m["mechanism"] == mech]
        v1 = sub["tto_str"][sub["tto_str"].between(1, 3650)]
        v3 = sub.loc[sub["gran"] == 3, "tto_str"]
        v3 = v3[v3.between(1, 3650)]
        v4 = sub["tto_strict"][sub["tto_strict"].between(1, 3650)]
        rows.append({
            "Mechanism": mech,
            "08 n": len(v1), "08 median": round(v1.median(), 1) if len(v1) else np.nan,
            "exact-day n": len(v3), "exact-day median": round(v3.median(), 1) if len(v3) else np.nan,
            "strict n": len(v4), "strict median": round(v4.median(), 1) if len(v4) else np.nan,
            "strict retention %": round(100 * len(v4) / max(len(v1), 1), 1),
        })
    print(pd.DataFrame(rows).to_string(index=False))

    print("\nDone.")


if __name__ == "__main__":
    main()
