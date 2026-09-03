# -*- coding: utf-8 -*-
"""Generic reader for FAERS/AEMS ASCII files
------------------------------------------
Handles the three format differences across 2004Q1-2026Q2:
  1. delimiter '$', no quoting (some drug names contain '$', needs tolerant rows)
  2. encoding latin-1 (includes cp1252 characters)
  3. field-name drift:
        legacy format (2004Q1-2012Q3) : ISR / CASE
        new format    (2012Q4-)       : primaryid / caseid
     both are unified to primaryid / caseid
"""

import os
import io
import glob
import pandas as pd

# legacy field names -> unified names
OLD2NEW = {
    "isr":      "primaryid",
    "case":     "caseid",
    "i_f_cod":  "i_f_code",
    "foll_seq": "foll_seq",
    "image":    "image",
    "event_dt": "event_dt",
    "mfr_dt":   "mfr_dt",
    "fda_dt":   "fda_dt",
    "rept_cod": "rept_cod",
    "mfr_num":  "mfr_num",
    "mfr_sndr": "mfr_sndr",
    "age":      "age",
    "age_cod":  "age_cod",
    "gndr_cod": "gndr_cod",
    "sex":      "gndr_cod",   # new format (post-2015Q1) uses 'sex' for sex
    "e_sub":    "e_sub",
    "wt":       "wt",
    "wt_cod":   "wt_cod",
    "rept_dt":  "rept_dt",
    "occp_cod": "occp_cod",
    "death_dt": "death_dt",
    "to_mfr":   "to_mfr",
    "confid":   "confid",
    "drug_seq": "drug_seq",
    "role_cod": "role_cod",
    "drugname": "drugname",
    "val_vbm":  "val_vbm",
    "route":    "route",
    "dose_vbm": "dose_vbm",
    "dechal":   "dechal",
    "rechal":   "rechal",
    "lot_num":  "lot_num",
    "exp_dt":   "exp_dt",
    "nda_num":  "nda_num",
    "pt":       "pt",
    "outc_cod": "outc_cod",
    "rpsr_cod": "rpsr_cod",
    "start_dt": "start_dt",
    "end_dt":   "end_dt",
    "dur":      "dur",
    "dur_cod":  "dur_cod",
    "indi_drug_seq": "indi_drug_seq",
    "indi_pt":  "indi_pt",
}

# columns to keep (per table); drop the rest to save memory
KEEP = {
    "DEMO": ["primaryid", "caseid", "event_dt", "fda_dt", "age", "age_cod",
             "gndr_cod", "wt", "wt_cod", "rept_dt", "occp_cod", "mfr_sndr",
             "rept_cod", "i_f_code", "foll_seq", "reporter_country", "occr_country"],
    "DRUG": ["primaryid", "caseid", "drug_seq", "role_cod", "drugname",
             "route", "dose_amt", "dose_unit", "val_vbm"],
    "REAC": ["primaryid", "caseid", "pt"],
    "OUTC": ["primaryid", "caseid", "outc_cod"],
    "RPSR": ["primaryid", "caseid", "rpsr_cod"],
    "THER": ["primaryid", "caseid", "start_dt", "end_dt", "dur", "dur_cod"],
    "INDI": ["primaryid", "caseid", "indi_drug_seq", "indi_pt"],
}


def find_file(qdir, table):
    """Locate a table file inside a quarter directory (case-insensitive)."""
    for pat in ("%s*.TXT" % table, "%s*.txt" % table):
        hits = glob.glob(os.path.join(qdir, pat))
        if hits:
            return hits[0]
    # some quarters place files directly under the quarter root
    for pat in ("%s*.TXT" % table, "%s*.txt" % table):
        hits = glob.glob(os.path.join(os.path.dirname(qdir), pat))
        if hits:
            return hits[0]
    return None


def read_table(path, table, usecols=None):
    """Read a single FAERS ASCII table, unify column names, return DataFrame (all str)."""
    with io.open(path, "r", encoding="latin-1", errors="replace") as f:
        header = f.readline().rstrip("\r\n")
    cols = [c.strip().lower().strip('"') for c in header.split("$")]
    cols = [OLD2NEW.get(c, c) for c in cols]

    df = pd.read_csv(
        path,
        sep="$",
        header=0,
        names=cols,
        dtype=str,
        engine="c",
        encoding="latin-1",
        on_bad_lines="skip",
        low_memory=False,
        quoting=3,          # QUOTE_NONE
    )
    # drop trailing all-empty columns (some files have a trailing '$')
    df = df.loc[:, [c for c in df.columns if c and not c.startswith("unnamed")]]

    want = usecols or KEEP.get(table)
    if want:
        keep = [c for c in want if c in df.columns]
        df = df[keep]
    return df


def qdirs(extr_root):
    """Return quarter directories sorted chronologically."""
    ds = []
    for name in os.listdir(extr_root):
        p = os.path.join(extr_root, name)
        if not os.path.isdir(p):
            continue
        # directory names look like 2012q4
        try:
            y, q = name.lower().split("q")
            ds.append(((int(y), int(q)), p, name.lower()))
        except Exception:
            continue
    ds.sort()
    return ds


def ascii_dir(qdir):
    """Quarter directory -> actual sub-directory holding the .TXT files (case-insensitive)."""
    for sub in (os.path.join(qdir, "ascii"), os.path.join(qdir, "ASCII")):
        if os.path.isdir(sub):
            return sub
    # fallback: find any sub-directory containing .txt files
    for entry in os.listdir(qdir):
        full = os.path.join(qdir, entry)
        if os.path.isdir(full):
            if glob.glob(os.path.join(full, "*.txt")) or glob.glob(os.path.join(full, "*.TXT")):
                return full
    return qdir
