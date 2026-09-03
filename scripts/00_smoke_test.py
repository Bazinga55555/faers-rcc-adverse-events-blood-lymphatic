# -*- coding: utf-8 -*-
"""Smoke test: verify that the reader in 02_parse_merge_v2 handles both the
legacy (pre-2012Q4) and the current FAERS ASCII layout correctly.
"""
import importlib.util
import io
import os
import sys

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

import faers_io as F
import pandas as pd
from paths import P, DIR_EXTRACTED

spec = importlib.util.spec_from_file_location(
    "v2", os.path.join(HERE, "02_parse_merge_v2.py"))
v2 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(v2)

EXTR = P(DIR_EXTRACTED)

# Check the DRUG filter on one legacy, one boundary, and one recent quarter
for name in ["2004q1", "2012q3", "2026q2"]:
    d = F.ascii_dir(os.path.join(EXTR, name))
    fp = F.find_file(d, "DRUG")
    if not fp:
        print(name, "NOT FOUND")
        continue
    df = v2.process_drug_table(fp)
    print("%s DRUG shape after filtering = %s" % (name, df.shape))
    if not df.empty:
        head = df[["drug_std", "mechanism", "drugname"]].head(2)
        print("   sample:", head.to_string(index=False).replace("\n", " | "))

# Check the REAC filter against a small set of primary ids
print("\n=== REAC filter test ===")
d = F.ascii_dir(os.path.join(EXTR, "2004q1"))
fp = F.find_file(d, "REAC")
if fp:
    # pull a few primary ids out of the DRUG table first
    dfp = F.find_file(d, "DRUG")
    drug = v2.process_drug_table(dfp)
    ids = set(drug["primaryid"].head(50))
    reac = v2.read_filtered(fp, "REAC", ids)
    print("REAC shape after filtering = %s, unique primaryid = %d"
          % (reac.shape, reac["primaryid"].nunique()))
