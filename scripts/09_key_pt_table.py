# -*- coding: utf-8 -*-
"""
09_key_pt_table.py -- core results table for the paper: key BLSD PT x mechanism-group signal matrix

Outputs (can drop directly into Table 2 / Table 3 of the paper):
  05_results/tables/table_key_pt_mechanism.csv    long table (mechanism, pt, n, ROR (95% CI), IC (IC025), signal flag)
  05_results/tables/table_key_pt_matrix.csv       matrix table (rows = PT, cols = mechanism, cell = "ROR (lower-upper)")
  05_results/tables/table_key_pt_drug.csv         drug-level long table (supplementary material)
  05_results/tables/table_mechanism_overview.csv  BLSD event overview for each mechanism group

Also prints a readable grouped result per mechanism group.
"""
import os
import sys
import io
import importlib.util

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import numpy as np
import pandas as pd
import types

# dynamically import 06_signal (filename starts with a digit, cannot be imported normally)
_spec = importlib.util.spec_from_file_location("sig", os.path.join(HERE, "06_signal.py"))
sig = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(sig)
assert isinstance(sig, types.ModuleType)

# key PT list and mechanism order (shared definition, see pt_defs.py)
from pt_defs import KEY_PT, KEY_PT_NAMES, LINEAGE, MECH_ORDER, MECH_LABEL

from paths import ROOT, DIR_CLEAN, DIR_RESULTS, SUB_TABLES
CLEAN = os.path.join(ROOT, DIR_CLEAN)
RES = os.path.join(ROOT, DIR_RESULTS, SUB_TABLES)

def main():
    os.makedirs(RES, exist_ok=True)
    print("Loading data...")
    drug, reac, demo = sig.load_inputs(include_lab=False)
    print("  de-duplicated reports %d; target drug rows %d; BLSD outcome rows %d"
          % (demo["primaryid"].nunique(), len(drug), len(reac)))

    # ---------- mechanism level ----------
    print("\nMechanism-level signal detection...")
    mres = sig.run(drug, reac, by_mechanism=True, verbose=False)
    dres = sig.run(drug, reac, by_mechanism=False, verbose=False)

    key_pt = [p for p, _ in KEY_PT]
    lineage = {p: g for p, g in KEY_PT}

    long = mres[mres["pt_disp"].isin(key_pt)].copy()
    long["lineage"] = long["pt_disp"].map(lineage)
    long["ci"] = long.apply(lambda r: "%.2f (%.2f-%.2f)" % (r["ror"], r["ror_low"], r["ror_high"]), axis=1)
    long = long[["lineage", "mechanism", "pt_disp", "n", "ror", "ror_low", "ror_high",
                 "prr", "chi2", "ic", "ic025", "ebgm", "ebgm05", "primary_signal", "ci"]]
    long = long.sort_values(["lineage", "pt_disp", "mechanism"])
    long.to_csv(os.path.join(RES, "table_key_pt_mechanism.csv"), index=False, encoding="utf-8-sig")
    print("  Wrote table_key_pt_mechanism.csv (%d rows)" % len(long))

    # ---------- matrix table ----------
    mat = {}
    for _, r in long.iterrows():
        mat.setdefault(r["pt_disp"], {})[r["mechanism"]] = "%s%s" % (r["ci"], "*" if r["primary_signal"] else "")
    rows = []
    for p in key_pt:
        if p not in mat:
            continue
        row = {"lineage": lineage[p], "pt": p}
        for m in MECH_ORDER:
            row[m] = mat[p].get(m, "—")
        rows.append(row)
    matrix = pd.DataFrame(rows)
    matrix.to_csv(os.path.join(RES, "table_key_pt_matrix.csv"), index=False, encoding="utf-8-sig")
    print("  Wrote table_key_pt_matrix.csv (%d PTs)" % len(matrix))

    # ---------- drug level (supplementary material) ----------
    dlong = dres[dres["pt_disp"].isin(key_pt)].copy()
    dlong["lineage"] = dlong["pt_disp"].map(lineage)
    dlong["ci"] = dlong.apply(lambda r: "%.2f (%.2f-%.2f)" % (r["ror"], r["ror_low"], r["ror_high"]), axis=1)
    dlong = dlong[["lineage", "drug_std", "mechanism", "pt_disp", "n", "ror",
                   "ror_low", "ror_high", "ic", "ic025", "primary_signal", "ci"]]
    dlong = dlong.sort_values(["lineage", "pt_disp", "drug_std"])
    dlong.to_csv(os.path.join(RES, "table_key_pt_drug.csv"), index=False, encoding="utf-8-sig")
    print("  Wrote table_key_pt_drug.csv (%d rows)" % len(dlong))

    # ---------- BLSD overview per mechanism group ----------
    merged = drug[["primaryid", "mechanism"]].merge(
        reac[["primaryid", "pt_std"]], on="primaryid", how="inner").drop_duplicates()
    rep = drug[["primaryid", "mechanism"]].drop_duplicates()
    ov = []
    for m in MECH_ORDER:
        n_rep = int((rep["mechanism"] == m).sum())
        ids = set(rep.loc[rep["mechanism"] == m, "primaryid"])
        sub = merged[merged["mechanism"] == m]
        n_evt = len(sub)
        n_with = sub["primaryid"].nunique()
        ov.append({
            "mechanism": m,
            "n_reports": n_rep,
            "n_blsd_events": n_evt,
            "n_reports_with_blsd": n_with,
            "blsd_report_pct": round(100.0 * n_with / n_rep, 2) if n_rep else np.nan,
            "blsd_events_per_report": round(n_evt / n_rep, 3) if n_rep else np.nan,
            "unique_blsd_pt": sub["pt_std"].nunique(),
        })
    ovdf = pd.DataFrame(ov)
    ovdf.to_csv(os.path.join(RES, "table_mechanism_overview.csv"), index=False, encoding="utf-8-sig")
    print("\n=== BLSD event overview by mechanism group ===")
    print(ovdf.to_string(index=False))

    # ---------- print matrix ----------
    print("\n=== Key PT x mechanism group ROR (95% CI), * = ROR+IC dual-positive ===")
    for lin in ["Erythroid", "Leukocyte", "Platelet",
                "Multilineage/Bone marrow", "Coagulation/Microvascular", "Lymphatic/Spleen"]:
        sub = matrix[matrix["lineage"] == lin]
        if sub.empty:
            continue
        print("\n[%s]" % lin)
        print(sub[["pt"] + MECH_ORDER].to_string(index=False))


if __name__ == "__main__":
    main()
