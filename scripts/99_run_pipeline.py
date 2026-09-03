# -*- coding: utf-8 -*-
"""
99_run_pipeline.py - run the full analysis pipeline in one step (v2)

Prerequisite: 01_download_faers.py has downloaded and extracted all 89+ quarters
into the extracted-data folder (DIR_EXTRACTED). If some quarters are missing, run
01b_fix_missing.py first to backfill them.

Executed in order:
  02_parse_merge_v2.py  -> parse, merge, dedupe (writes DIR_CLEAN/*.parquet)
  07_blsd_pt.py         -> BLSD-restricted PT list + PT normalisation (reac_norm.parquet)
  06_signal.py          -> four-algorithm signal detection (drug + mechanism level)
  09_key_pt_table.py    -> key PT x mechanism core table
  08_weibull.py         -> Weibull TTO classification
  10_sensitivity.py     -> five-scenario sensitivity analysis
  11_figures.py         -> four paper figures

Usage:
  python 99_run_pipeline.py
"""
import os
import sys
import io
import time
import subprocess

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")

PY = r"D:\Python\Python39\python.exe"
SCR = os.path.dirname(os.path.abspath(__file__))

STEPS = [
    ("02_parse_merge_v2.py",   "parse + merge + dedupe (v2 memory-friendly)"),
    ("07_blsd_pt.py",          "BLSD PT restriction + PT normalisation"),
    ("06_signal.py",           "disproportionality (drug level)"),
    ("06_signal.py --by-mechanism", "disproportionality (mechanism level)"),
    ("09_key_pt_table.py",     "key PT x mechanism table"),
    ("08_weibull.py",          "Weibull TTO classification"),
    ("10_sensitivity.py",      "five-scenario sensitivity analysis"),
    ("11_figures.py",          "four paper figures"),
]


def main():
    t0 = time.time()
    for i, (cmd, label) in enumerate(STEPS, 1):
        print("\n" + "=" * 60)
        print("[%s] step %d/%d: %s" % (time.strftime("%H:%M:%S"), i, len(STEPS), label))
        print("=" * 60)
        full = [PY] + cmd.split()
        rc = subprocess.run(full, cwd=SCR)
        if rc.returncode != 0:
            print("!! step failed: %s (exit=%d)" % (label, rc.returncode))
            sys.exit(rc.returncode)
    print("\n[%s] pipeline complete, elapsed %.1f min" % (time.strftime("%H:%M:%S"), (time.time() - t0) / 60))


if __name__ == "__main__":
    main()
