# -*- coding: utf-8 -*-
"""
02_parse_merge.py -- FAERS full-quarter parsing + merge + de-duplication

Steps:
  1. Iterate 02_extracted/{quarter}/ascii/*.TXT, read the 7 tables per quarter
  2. Unify column names (legacy ISR/CASE -> primaryid/caseid)
  3. Vertically concatenate across quarters into 7 full tables
  4. De-duplicate (FDA three-step rule: within CASEID keep max FDA_DT ->
     for the same FDA_DT keep max PRIMARYID)
  5. Write parquet to 03_clean_data/

Outputs:
  03_clean_data/demo.parquet        demographics (de-duplicated)
  03_clean_data/drug.parquet        drugs (NOT de-duplicated; at dedup time
                                    filtered by demo's primaryid)
  03_clean_data/reac.parquet        reactions
  03_clean_data/outc.parquet        outcomes
  03_clean_data/rpsr.parquet        report sources
  03_clean_data/ther.parquet        therapy dates
  03_clean_data/indi.parquet        indications
  03_clean_data/demo_dedup.parquet  de-duplicated demographics (with dedup flag)
"""
import os
import sys
import io
import time
import glob

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import pandas as pd
import faers_io as F

from paths import ROOT, DIR_CLEAN, DIR_EXTRACTED
EXTR = os.path.join(ROOT, DIR_EXTRACTED)
OUT  = os.path.join(ROOT, DIR_CLEAN)

TABLES = ["DEMO", "DRUG", "REAC", "OUTC", "RPSR", "THER", "INDI"]


def parse_all():
    os.makedirs(OUT, exist_ok=True)
    qds = F.qdirs(EXTR)
    print("[%s] Start parsing %d quarters" % (time.strftime("%H:%M:%S"), len(qds)))

    merged = {t: [] for t in TABLES}
    for i, (_, p, name) in enumerate(qds, 1):
        d = F.ascii_dir(p)
        for t in TABLES:
            fp = F.find_file(d, t)
            if not fp:
                continue
            try:
                df = F.read_table(fp, t)
            except Exception as e:
                print("  ! %s %s read failed: %s" % (name, t, repr(e)[:100]))
                continue
            # record source quarter (handy for troubleshooting)
            df["_q"] = name
            merged[t].append(df)
        if i % 10 == 0 or i == len(qds):
            print("[%s] Parsed %d/%d quarters" % (time.strftime("%H:%M:%S"), i, len(qds)))

    print("[%s] Vertical concatenation..." % time.strftime("%H:%M:%S"))
    out = {}
    for t in TABLES:
        if not merged[t]:
            print("  ! %s has no data" % t)
            continue
        out[t] = pd.concat(merged[t], ignore_index=True)
        print("  %-6s %d rows after merge" % (t, len(out[t])))

    return out


def dedup(demo):
    """FDA three-step de-duplication: within CASEID keep max FDA_DT ->
    for the same FDA_DT keep max PRIMARYID."""
    d = demo.copy()
    # 1) sort
    d = d.sort_values(["caseid", "fda_dt", "primaryid"],
                      na_position="last")
    # 2) within the same CASEID keep the row with max FDA_DT
    d = d.drop_duplicates(subset=["caseid"], keep="last")
    # 3) for the same (caseid, fda_dt) keep max primaryid
    #    (already guaranteed by the previous step; fallback here)
    d = d.drop_duplicates(subset=["caseid", "fda_dt"], keep="last")
    return d


def main():
    t0 = time.time()
    tables = parse_all()

    # write the raw (non-deduplicated) 7 tables
    for t, df in tables.items():
        df.to_parquet(os.path.join(OUT, "%s.parquet" % t.lower()), index=False)
        print("  Wrote %s.parquet (%d rows)" % (t.lower(), len(df)))

    # de-duplicate
    print("[%s] De-duplicating..." % time.strftime("%H:%M:%S"))
    demo_raw = tables.get("DEMO")
    if demo_raw is not None:
        n_raw = len(demo_raw)
        demo_dd = dedup(demo_raw)
        # keep the key columns used by dedup + all remaining columns
        demo_dd.to_parquet(os.path.join(OUT, "demo_dedup.parquet"), index=False)
        print("  DEMO dedupe: %d -> %d rows (kept %d, removed %.1f%%)"
              % (n_raw, len(demo_dd), len(demo_dd), (1 - len(demo_dd) / n_raw) * 100))
        # whitelist of primaryid after de-duplication
        keep_ids = set(demo_dd["primaryid"])
        keep_ids.to_parquet  # noop
        pd.DataFrame({"primaryid": sorted(keep_ids)}).to_parquet(
            os.path.join(OUT, "dedup_ids.parquet"), index=False)
        print("  Dedup whitelist primaryid count: %d" % len(keep_ids))

    print("[%s] Done, elapsed %.1f min" % (time.strftime("%H:%M:%S"), (time.time() - t0) / 60))


if __name__ == "__main__":
    main()
