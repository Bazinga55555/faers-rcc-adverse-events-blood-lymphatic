# -*- coding: utf-8 -*-
"""
04_drug_map.py -- FAERS drug-name standardization + mechanism grouping

Map the free-text drugname in drug.parquet to the 15 target drugs and tag each
with a mechanism label.

Matching strategy (four progressive layers):
  1. Exact uppercase match (generic / brand / code name)
  2. Regex boundary match (covers salt forms, hyphen variants)
  3. Normalized substring match (strip spaces/hyphens/salt suffixes)
  4. No match -> assign to "OTHER" (excluded from the cohort)

Outputs:
  03_clean_data/drug_mapped.parquet   only target-drug rows (with drug_std + mechanism)
  03_clean_data/match_audit.csv       matched row count audit per drug
"""
import os
import sys
import io
import re
import time

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import pandas as pd

from paths import ROOT, DIR_CLEAN
OUT  = os.path.join(ROOT, DIR_CLEAN)

# ---- drug regex table (generic + brand + code/salt-form variants) ----
# each mechanism -> dict(standard_name -> regex list)
MECHANISMS = {
    "A_HIF2a": {
        "belzutifan": [r"\bBELZUTIFAN\b", r"\bWELIREG\b", r"\bMK[ -]?6482\b", r"\bPT ?2385\b"],
    },
    "B_VEGFR_TKI": {
        "sunitinib":    [r"\bSUNITINIB\b", r"\bSUTENT\b", r"\bSU ?0?11248\b"],
        "sorafenib":    [r"\bSORAFENIB\b", r"\bNEXAVAR\b", r"\bBAY ?43[ -]?9006\b"],
        "pazopanib":    [r"\bPAZOPANIB\b", r"\bVOTRIENT\b", r"\bGW ?786034\b"],
        "axitinib":     [r"\bAXITINIB\b", r"\bINLYTA\b", r"\bAG ?0?13736\b"],
        "cabozantinib": [r"\bCABOZANTINIB\b", r"\bCABOMETYX\b", r"\bCOMETRIQ\b",
                         r"\bXL ?184\b", r"\bBMS[ -]?907351\b"],
        "lenvatinib":   [r"\bLENVATINIB\b", r"\bLENVIMA\b", r"\bLENVIPA\b", r"\bE ?7080\b"],
        "tivozanib":    [r"\bTIVOZANIB\b", r"\bFOTIVDA\b", r"\bAV[ -]?951\b", r"\bKRN ?951\b"],
    },
    "C_VEGF_mAb": {
        "bevacizumab":  [r"\bBEVACIZUMAB\b", r"\bAVASTIN\b", r"\bMVASI\b", r"\bZIRABEV\b",
                         r"\bALYMSYS\b", r"\bVEGZELMA\b", r"\bONBEVZI\b"],
    },
    "D_mTOR": {
        "everolimus":   [r"\bEVEROLIMUS\b", r"\bAFINITOR\b", r"\bRAD ?001\b", r"\bSDZ[ -]?RAD\b"],
        "temsirolimus": [r"\bTEMSIROLIMUS\b", r"\bTORISEL\b", r"\bCCI[ -]?779\b"],
    },
    "E_ICI": {
        "nivolumab":     [r"\bNIVOLUMAB\b", r"\bOPDIVO\b", r"\bBMS[ -]?936558\b", r"\bMDX[ -]?1106\b"],
        "ipilimumab":    [r"\bIPILIMUMAB\b", r"\bYERVOY\b", r"\bBMS[ -]?734016\b", r"\bMDX[ -]?010\b"],
        "pembrolizumab": [r"\bPEMBROLIZUMAB\b", r"\bKEYTRUDA\b", r"\bMK[ -]?3475\b", r"\bSCH ?900475\b"],
        "avelumab":      [r"\bAVELUMAB\b", r"\bBAVENCIO\b", r"\bMSB ?0010718C\b"],
    },
}

# precompile regexes
COMPILED = {}
for mech, drugs in MECHANISMS.items():
    for std, pats in drugs.items():
        COMPILED[std] = (mech, [re.compile(p) for p in pats])

# salt-form suffixes (used for normalized fallback matching)
SALT_SUFFIX = re.compile(r"\s*(MALATE|TOSYLATE|TOSILATE|HYDROCHLORIDE|HCL|BESYLATE|MESYLATE|SUCCINATE)\s*$")


def normalize(s):
    if not isinstance(s, str):
        return ""
    s = s.upper().strip()
    s = SALT_SUFFIX.sub("", s)
    s = re.sub(r"[\s\-/]+", "", s)  # strip spaces/hyphens/slashes
    return s


def map_drugname(val):
    """Return (standard_name, mechanism) or (None, None)."""
    if not isinstance(val, str):
        return None, None
    up = val.upper().strip()

    # 1) regex boundary match
    for std, (mech, pats) in COMPILED.items():
        for p in pats:
            if p.search(up):
                return std, mech

    # 2) normalized exact match (fallback)
    norm = normalize(up)
    if not norm:
        return None, None
    for std, (mech, pats) in COMPILED.items():
        for p in pats:
            # take the regex core word (strip \b and groups), compare after normalization
            core = p.pattern.replace("\\b", "").replace(" ", "").replace("?", "")
            core = re.sub(r"[\[\]\(\)\-\s]+", "", core)
            if core and core == norm:
                return std, mech
    return None, None


def main():
    os.makedirs(OUT, exist_ok=True)
    drug_path = os.path.join(OUT, "drug.parquet")
    if not os.path.exists(drug_path):
        print("!! Cannot find %s; run 02_parse_merge.py first" % drug_path)
        return

    print("[%s] Reading drug.parquet..." % time.strftime("%H:%M:%S"))
    drug = pd.read_parquet(drug_path)
    print("  Total drug rows: %d" % len(drug))

    # mapping (vectorized: build a lookup from unique drugname first to avoid row-by-row O(n))
    uniq = drug["drugname"].dropna().unique()
    print("  Unique drugname count: %d" % len(uniq))
    mapdict = {}
    for u in uniq:
        std, mech = map_drugname(u)
        if std:
            mapdict[u] = (std, mech)

    print("  Unique names matching target drugs: %d" % len(mapdict))

    # apply mapping
    mapped = drug[drug["drugname"].isin(mapdict.keys())].copy()
    mapped["drug_std"] = mapped["drugname"].map(lambda v: mapdict[v][0])
    mapped["mechanism"] = mapped["drugname"].map(lambda v: mapdict[v][1])

    # audit
    audit = (
        mapped.groupby(["mechanism", "drug_std"])
        .size()
        .reset_index(name="n_drug_rows")
        .sort_values(["mechanism", "n_drug_rows"], ascending=[True, False])
    )
    audit.to_csv(os.path.join(OUT, "match_audit.csv"), index=False, encoding="utf-8-sig")

    mapped.to_parquet(os.path.join(OUT, "drug_mapped.parquet"), index=False)

    print("\n=== Matching audit ===")
    print(audit.to_string(index=False))
    print("\n[%s] Done; drug_mapped.parquet has %d rows" % (time.strftime("%H:%M:%S"), len(mapped)))


if __name__ == "__main__":
    main()
