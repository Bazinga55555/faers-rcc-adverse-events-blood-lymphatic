# -*- coding: utf-8 -*-
"""Central path configuration for the FAERS RCC pharmacovigilance pipeline.

Every analysis script imports ``ROOT`` and the directory-name constants from
this module instead of hard-coding absolute paths, so the same code runs
unchanged on any machine once the project root is known.

Supply the project root in either of two ways:

1. the environment variable ``FAERS_RCC_ROOT``; or
2. a ``paths_local.py`` file next to this one, defining ``PROJECT_ROOT``.
   That file is deliberately not committed (see ``.gitignore``), so no
   machine-specific path ever enters version control.

Nothing else in this repository contains a machine-specific path.
"""

import os
import sys

# Directory names inside the project root. Change them here if you lay the
# project out under different folder names on your own machine.
DIR_RAW = "01_原始数据"            # downloaded quarterly ZIP files
DIR_EXTRACTED = "02_解压数据"       # unzipped quarterly ASCII files
DIR_CLEAN = "03_清洗数据"           # parsed, de-duplicated parquet files
DIR_RESULTS = "05_结果"             # analysis output
DIR_EXTERNAL = "00_外部库"          # JADER / Canada Vigilance extracts
DIR_SCRIPTS = "04_分析脚本"          # analysis scripts (log files land here)
SUB_TABLES = "表"                   # results/tables
SUB_FIGURES = "图"                  # results/figures


def _find_root():
    """Resolve the project root from the environment or a local file."""
    env = os.environ.get("FAERS_RCC_ROOT")
    if env:
        return env
    here = os.path.dirname(os.path.abspath(__file__))
    if here not in sys.path:
        sys.path.insert(0, here)
    try:
        from paths_local import PROJECT_ROOT  # optional, not committed
        return PROJECT_ROOT
    except ImportError:
        return None


ROOT = _find_root()

if not ROOT:
    sys.stderr.write(
        "FAERS_RCC_ROOT is not set, so the project root is unknown.\n"
        "Point it at the directory that holds the 0*/ folders, e.g.\n"
        "    export FAERS_RCC_ROOT=/path/to/project     (Linux/macOS)\n"
        "    set FAERS_RCC_ROOT=C:\\path\\to\\project      (Windows)\n"
        "Alternatively create paths_local.py beside paths.py containing\n"
        "    PROJECT_ROOT = r'/path/to/project'\n"
    )
    raise SystemExit(1)


def P(*parts):
    """Join one or more path fragments onto the project root."""
    return os.path.join(ROOT, *parts)


def clean(*parts):
    """Path inside the cleaned-data directory."""
    return P(DIR_CLEAN, *parts)


def results(*parts):
    """Path inside the results directory."""
    return P(DIR_RESULTS, *parts)


def tables(*parts):
    """Path inside the results/tables directory."""
    return P(DIR_RESULTS, SUB_TABLES, *parts)


def figures(*parts):
    """Path inside the results/figures directory."""
    return P(DIR_RESULTS, SUB_FIGURES, *parts)


def external(*parts):
    """Path inside the external-database directory."""
    return P(DIR_EXTERNAL, *parts)


def ensure(*dirs):
    """Create the given directories if they do not already exist."""
    for d in dirs:
        if d and not os.path.isdir(d):
            os.makedirs(d)
