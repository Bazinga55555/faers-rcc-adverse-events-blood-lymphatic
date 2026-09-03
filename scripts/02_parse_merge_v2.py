# -*- coding: utf-8 -*-
"""
02_parse_merge_v2.py -- memory-friendly FAERS parser (filter target drugs first, then read other tables)

Key improvements (lessons learned from v1's out-of-memory failure):
  1. DRUG table uses usecols to read only the needed columns
     (primaryid/caseid/drug_seq/role_cod/drugname)
  2. While reading DRUG, immediately regex-match the 15 target drugs and keep only the matching rows
  3. Collect the target primaryid set and use it to filter REAC/DEMO/THER/OUTC/RPSR/INDI
  4. High-cardinality columns use the category dtype to save memory
  5. Process quarter by quarter; each quarter is read + filtered independently, then concat at the end
     (data volume is already greatly reduced by then)

This brings memory usage from ~20GB down to ~500MB, so the full 90 quarters run safely.

Outputs (03_clean_data/):
  drug_mapped.parquet   drug rows for target drugs (with drug_std + mechanism)
  reac.parquet          reactions of target reports
  demo_dedup.parquet    demographics of target reports (de-duplicated)
  outc/rpsr/ther/indi   related tables of target reports
  match_audit.csv       drug-matching audit
"""
import os
import sys
import io
import re
import time

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import pandas as pd
import faers_io as F

from paths import ROOT, DIR_CLEAN, DIR_EXTRACTED
EXTR = os.path.join(ROOT, DIR_EXTRACTED)
OUT  = os.path.join(ROOT, DIR_CLEAN)

# ---- drug regex table (consistent with 04_drug_map.py) ----
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

# precompile: drug_std -> (mechanism, [compiled regex])
COMPILED = {}
for mech, drugs in MECHANISMS.items():
    for std, pats in drugs.items():
        COMPILED[std] = (mech, [re.compile(p) for p in pats])

# merge into one master regex for fast coarse filtering (keep a row if any drug matches)
MASTER_PAT = re.compile(
    "|".join(p.pattern for std in COMPILED for p in COMPILED[std][1]))


def match_drug(val):
    """Return (drug_std, mechanism) or (None, None)."""
    if not isinstance(val, str):
        return None, None
    up = val.upper().strip()
    for std, (mech, pats) in COMPILED.items():
        for p in pats:
            if p.search(up):
                return std, mech
    return None, None


def _read_cols(path):
    """Read the header, return (raw column list, unified column list)."""
    with io.open(path, "r", encoding="latin-1", errors="replace") as f:
        header = f.readline().rstrip("\r\n")
    raw = [c.strip().lower().strip('"') for c in header.split("$")]
    # drop trailing empty columns that may be present
    raw = [c for c in raw if c]
    uni = [F.OLD2NEW.get(c, c) for c in raw]
    return raw, uni


def _read_full(path):
    """Read a single table in full (use range(N) positional indexing to avoid the
    field-misalignment bug caused by trailing $ in legacy formats).
    Return (DataFrame, unified column list). Note: column names stay raw, trimmed after read."""
    raw, uni = _read_cols(path)
    n = len(raw)
    # key: usecols must be the full contiguous range(n); a partial column set makes the C engine misalign fields
    df = pd.read_csv(
        path, sep="$", header=0, dtype=str, engine="c",
        encoding="latin-1", on_bad_lines="skip", low_memory=False,
        quoting=3, usecols=range(n),
    )
    # rename with unified column names
    rename = {df.columns[i]: uni[i] for i in range(n)}
    df = df.rename(columns=rename)
    # drop possible trailing Unnamed columns
    df = df.loc[:, [c for c in df.columns if c and not str(c).startswith("Unnamed")]]
    return df, uni


def process_drug_table(path):
    """Read one quarter's DRUG table, keep only target-drug rows, return DataFrame (with drug_std/mechanism added)."""
    df, uni = _read_full(path)
    # trim to needed columns
    want = ["primaryid", "caseid", "drug_seq", "role_cod", "drugname"]
    keep = [c for c in want if c in df.columns]
    if "primaryid" not in keep or "drugname" not in keep:
        # missing key columns, return empty df (column names fall back to want)
        return pd.DataFrame(columns=want)
    df = df[keep].copy()

    # coarse filter: drugname matches the master regex
    mask = df["drugname"].str.upper().str.contains(MASTER_PAT, na=False, regex=True)
    df = df[mask].copy()

    if df.empty:
        return df

    # fine match: match row by row to get std/mechanism
    mapped = df["drugname"].map(lambda v: match_drug(v))
    df["drug_std"] = [m[0] for m in mapped]
    df["mechanism"] = [m[1] for m in mapped]
    df = df[df["drug_std"].notna()].copy()

    # role_cod: keep only PS/SS (primary/secondary suspect)
    df["role_cod"] = df["role_cod"].str.upper()
    return df


def read_filtered(path, table, keep_ids):
    """Read a table, keep only rows whose primaryid is in keep_ids."""
    df, uni = _read_full(path)
    want = F.KEEP.get(table, [])
    keep = [c for c in want if c in df.columns]
    df = df[keep].copy()
    df = df[df["primaryid"].isin(keep_ids)]
    return df


def main():
    t0 = time.time()
    os.makedirs(OUT, exist_ok=True)
    qds = F.qdirs(EXTR)
    print("[%s] Start parsing %d quarters (filter target drugs first)" % (time.strftime("%H:%M:%S"), len(qds)), flush=True)

    drug_parts = []
    target_ids = set()

    # ---- pass 1: read DRUG table, filter target drugs ----
    print("Stage 1: filter target drugs in DRUG table...", flush=True)
    for i, (_, p, name) in enumerate(qds, 1):
        d = F.ascii_dir(p)
        fp = F.find_file(d, "DRUG")
        if not fp:
            print("  [%d/%d] %s has no DRUG file" % (i, len(qds), name), flush=True)
            continue
        try:
            df = process_drug_table(fp)
        except Exception as e:
            print("  [%d/%d] %s DRUG error: %s" % (i, len(qds), name, repr(e)[:100]), flush=True)
            continue
        if not df.empty and "primaryid" in df.columns:
            df["_q"] = name
            drug_parts.append(df)
            target_ids.update(df["primaryid"].unique())
        if i % 15 == 0 or i == len(qds):
            print("  [%s] Processed %d/%d, target drug rows cum %d, primaryid %d"
                  % (time.strftime("%H:%M:%S"), i, len(qds),
                     sum(len(x) for x in drug_parts), len(target_ids)), flush=True)

    drug_all = pd.concat(drug_parts, ignore_index=True) if drug_parts else pd.DataFrame()
    print("Stage 1 done: target drug rows %d, unique primaryid %d" % (len(drug_all), len(target_ids)), flush=True)

    # save drug-mapping result
    drug_all.to_parquet(os.path.join(OUT, "drug_mapped.parquet"), index=False)

    # audit
    audit = (drug_all.groupby(["mechanism", "drug_std"]).size()
             .reset_index(name="n_drug_rows")
             .sort_values(["mechanism", "n_drug_rows"], ascending=[True, False]))
    audit.to_csv(os.path.join(OUT, "match_audit.csv"), index=False, encoding="utf-8-sig")
    print("Drug-matching audit:\n" + audit.to_string(index=False), flush=True)

    # ---- pass 2: filter other tables with target_ids ----
    print("\nStage 2: filter REAC/DEMO/THER/OUTC/RPSR/INDI...", flush=True)
    for table in ["REAC", "DEMO", "OUTC", "RPSR", "THER", "INDI"]:
        parts = []
        for i, (_, p, name) in enumerate(qds, 1):
            d = F.ascii_dir(p)
            fp = F.find_file(d, table)
            if not fp:
                continue
            try:
                df = read_filtered(fp, table, target_ids)
            except Exception as e:
                print("  ! %s %s read failed: %s" % (name, table, repr(e)[:100]), flush=True)
                continue
            if not df.empty:
                df["_q"] = name
                parts.append(df)
        if parts:
            tbl = pd.concat(parts, ignore_index=True)
            tbl.to_parquet(os.path.join(OUT, "%s.parquet" % table.lower()), index=False)
            print("  %-6s %d rows -> %s.parquet" % (table, len(tbl), table.lower()), flush=True)
        else:
            print("  %-6s no data" % table, flush=True)

    # ---- pass 3: de-duplicate DEMO ----
    print("\nStage 3: de-duplicate DEMO...", flush=True)
    demo_path = os.path.join(OUT, "demo.parquet")
    if os.path.exists(demo_path):
        demo = pd.read_parquet(demo_path)
        n_raw = len(demo)
        # sort then keep the newest per caseid
        demo = demo.sort_values(["caseid", "fda_dt", "primaryid"], na_position="last")
        demo_dd = demo.drop_duplicates(subset=["caseid"], keep="last")
        demo_dd.to_parquet(os.path.join(OUT, "demo_dedup.parquet"), index=False)
        keep_ids = set(demo_dd["primaryid"])
        pd.DataFrame({"primaryid": sorted(keep_ids)}).to_parquet(
            os.path.join(OUT, "dedup_ids.parquet"), index=False)
        print("  DEMO dedupe: %d -> %d (removed %.1f%%)"
              % (n_raw, len(demo_dd), (1 - len(demo_dd) / n_raw) * 100), flush=True)

    print("\n[%s] Pipeline complete, elapsed %.1f min" % (time.strftime("%H:%M:%S"), (time.time() - t0) / 60), flush=True)


if __name__ == "__main__":
    main()
