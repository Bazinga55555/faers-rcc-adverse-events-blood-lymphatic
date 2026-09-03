# -*- coding: utf-8 -*-
"""
19_flow_diagram.py - Study flow diagram (READUS-PV item 8a / STROBE-style)

Outputs:
    <DIR_RESULTS>/<SUB_FIGURES>/Fig0_flow_diagram.png   (300 dpi submission figure)
    <DIR_RESULTS>/<SUB_FIGURES>/Fig0_flow_diagram.svg    vector version (editable in Inkscape)
    <DIR_RESULTS>/<SUB_FIGURES>/Fig0_flow_diagram.txt    text version (for Supplementary)

All counts come from the actual pipeline run (rerun_full.log and Table 1), so
they stay strictly consistent with the main text, Table 1, and the signal results.

Usage:
    python 19_flow_diagram.py
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

# Keep text as editable <text> nodes in SVG instead of glyph outlines,
# so the diagram can be fine-tuned in Inkscape.
plt.rcParams["svg.fonttype"] = "none"
plt.rcParams["font.sans-serif"] = ["DejaVu Sans", "Arial", "Helvetica"]

sys.stdout = __import__("io").TextIOWrapper(sys.stdout.buffer, encoding="utf-8")

from paths import ROOT, DIR_RESULTS, SUB_FIGURES
FIGDIR = os.path.join(ROOT, DIR_RESULTS, SUB_FIGURES)
os.makedirs(FIGDIR, exist_ok=True)

# ------------------------------------------------------------------ colours
C_BOX = "#1f4e79"      # dark blue (included)
C_EXC = "#b23b3b"      # red (excluded)
C_EDGE = "#444444"
C_MAIN_TEXT = "#1a1a1a"
C_SUB = "#5a5a5a"

# ------------------------------------------------------------------ canvas
fig, ax = plt.subplots(figsize=(10.5, 12.5), dpi=300)
ax.set_xlim(0, 100)
ax.set_ylim(0, 130)
ax.axis("off")


def box(cx, cy, w, h, title, sub, fc=C_BOX, tc="white", ts_title=12.5, ts_sub=9.5):
    """Main box: rounded rectangle + title + optional sub-line(s)."""
    ax.add_patch(FancyBboxPatch(
        (cx - w / 2, cy - h / 2), w, h,
        boxstyle="round,pad=0.8,rounding_size=1.6",
        linewidth=1.4, edgecolor=C_EDGE, facecolor=fc, zorder=3))
    ax.text(cx, cy + h * 0.16, title, ha="center", va="center",
            fontsize=ts_title, fontweight="bold", color=tc, zorder=4)
    if isinstance(sub, str):
        sub = [sub]
    yy = cy - h * 0.10
    for i, s in enumerate(sub):
        ax.text(cx, yy - i * (h * 0.20), s, ha="center", va="center",
                fontsize=ts_sub, color=tc, zorder=4)


def excl_box(cx, cy, w, h, title, sub):
    """Exclusion box: right side, red border, white fill, red text."""
    ax.add_patch(FancyBboxPatch(
        (cx - w / 2, cy - h / 2), w, h,
        boxstyle="round,pad=0.6,rounding_size=1.2",
        linewidth=1.3, edgecolor=C_EXC, facecolor="white", zorder=3,
        linestyle="--"))
    ax.text(cx, cy + h * 0.18, title, ha="center", va="center",
            fontsize=10.5, fontweight="bold", color=C_EXC, zorder=4)
    if isinstance(sub, str):
        sub = [sub]
    yy = cy - h * 0.10
    for i, s in enumerate(sub):
        ax.text(cx, yy - i * (h * 0.24), s, ha="center", va="center",
                fontsize=8.6, color=C_MAIN_TEXT, zorder=4)


def arrow(x1, y1, x2, y2, color="#333333"):
    ax.add_patch(FancyArrowPatch(
        (x1, y1), (x2, y2), arrowstyle="-|>", mutation_scale=18,
        linewidth=1.5, color=color, zorder=2))


def excl_arrow(x1, y1, x2, y2):
    ax.add_patch(FancyArrowPatch(
        (x1, y1), (x2, y2), arrowstyle="-|>", mutation_scale=14,
        linewidth=1.1, color=C_EXC, zorder=2, linestyle=":"))


# ============================================================ layout
MAIN_X = 46.0          # main-flow axis
EXC_X = 85.0           # exclusion-box axis
W_MAIN = 60.0
H = 15.0

# ---- Box 1: data source ----
box(MAIN_X, 120, W_MAIN, H + 4,
    "FAERS / AEMS public data, 2004Q1 - 2026Q2",
    ["90 quarterly files, full download with CRC verification"])

# ---- Box 2: drug mapping ----
box(MAIN_X, 100, W_MAIN, H,
    "Drug records mapped to 15 RCC drugs",
    ["1,130,974 drug-event records   .   710,384 unique reports (primaryid)"])

# ---- Box 3: deduplication ----
box(MAIN_X, 81, W_MAIN, H - 2,
    "Unique reports after case deduplication",
    ["549,021 reports"])

# ---- Box 4: PS/SS ----
box(MAIN_X, 62, W_MAIN, H,
    "Reports with >=1 drug as primary / secondary suspect",
    ["522,813 unique reports   (557,950 report-mechanism pairs)"])
excl_box(EXC_X, 71.5, 24, 9,
    "Excluded (n = 26,208)",
    ["drug not PS/SS\n(concomitant /\ninteracting / unclear)"])

# ---- Box 5: BLSD ----
box(MAIN_X, 44, W_MAIN, H,
    "Outcome restricted to 251 core BLSD Preferred Terms",
    ["65,751 BLSD outcome records   .   49,254 reports with >=1 BLSD event"])

# ---- Box 6: signals ----
box(MAIN_X, 26, W_MAIN, H,
    "Disproportionality analysis (ROR / PRR / IC / EBGM)",
    ["primary signal = ROR 95% CI lower bound > 1  and  IC025 > 0, n >= 3"])

# ---- Box 7: results ----
box(MAIN_X, 8, W_MAIN, H + 2,
    "Signals detected",
    ["Mechanism level (38 key PTs): 44 primary signals",
     "Drug level: 196 primary signals"])

# ============================================================ arrows (main flow)
arrow(MAIN_X, 117.4, MAIN_X, 109.0)
arrow(MAIN_X, 91.0, MAIN_X, 89.0)
arrow(MAIN_X, 72.0, MAIN_X, 70.0)
arrow(MAIN_X, 53.0, MAIN_X, 52.0)
arrow(MAIN_X, 35.0, MAIN_X, 34.0)
arrow(MAIN_X, 17.0, MAIN_X, 15.5)

# ============================================================ exclusion arrows
excl_arrow(71.0, 62.0, EXC_X - 11, 71.5)

# ============================================================ title
ax.text(MAIN_X, 128.5,
        "Study flow diagram - READUS-PV item 8a",
        ha="center", va="center", fontsize=14.5, fontweight="bold",
        color=C_MAIN_TEXT)
ax.text(MAIN_X, 125.8,
        "Mechanism-stratified disproportionality analysis of blood and lymphatic system disorders",
        ha="center", va="center", fontsize=10, color=C_SUB, style="italic")

# ============================================================ save
out_png = os.path.join(FIGDIR, "Fig0_flow_diagram.png")
out_svg = os.path.join(FIGDIR, "Fig0_flow_diagram.svg")
plt.tight_layout()
plt.savefig(out_png, dpi=300, bbox_inches="tight", facecolor="white")
plt.savefig(out_svg, format="svg", bbox_inches="tight", facecolor="white")
plt.close()
print(f"[write] {out_png}")
print(f"[write] {out_svg}")

# ============================================================ text version
txt = """Study flow diagram (READUS-PV item 8a)

FAERS / AEMS public data, 2004Q1 - 2026Q2 (90 quarterly files, CRC-verified)
  |
  v
Drug records mapped to 15 RCC drugs
  n = 1,130,974 drug-event records; 710,384 unique reports (primaryid)
  |
  v
Unique reports after case deduplication
  n = 549,021
  |
  v
Reports with >=1 drug as primary/secondary suspect (PS/SS)
  n = 522,813 unique reports (557,950 report-mechanism pairs)
    \\--> Excluded n = 26,208: drug not PS/SS (concomitant / interacting / unclear)
  |
  v
Outcome restricted to 251 core BLSD Preferred Terms
  n = 65,751 BLSD outcome records; 49,254 reports with >=1 BLSD event
  |
  v
Disproportionality analysis (ROR / PRR / BCPNN-IC / MGPS-EBGM)
  primary signal = ROR 95% CI lower bound > 1 and IC025 > 0, n >= 3
  |
  v
Signals detected
  Mechanism level (38 key PTs): 44 primary signals
  Drug level: 196 primary signals
"""
out_txt = os.path.join(FIGDIR, "Fig0_flow_diagram.txt")
with open(out_txt, "w", encoding="utf-8") as f:
    f.write(txt)
print(f"[write] {out_txt}")
print("Done.")
