# -*- coding: utf-8 -*-
"""
11_figures.py -- Journal figures (Drug Safety style: clean, greyscale-friendly, 300 dpi)

Outputs (both a 300-dpi PNG and an editable SVG for each):
    Fig1_heatmap_ror.png / .svg      ROR heatmap, key BLSD PT x mechanism class (log2 scale)
    Fig2_forest_headline.png / .svg  Forest plot of headline signals (ROR + 95% CI, log axis)
    Fig3_tto_cumulative.png / .svg   Cumulative time-to-onset curves of BLSD events by class
    Fig4_mechanism_profile.png / .svg BLSD report share + lineage composition by class

The SVG is written with rcParams["svg.fonttype"] = "none", so every label remains a
real <text> element and can be edited directly in Inkscape.

Depends on: 05_results/tables/table_key_pt_mechanism.csv, sensitivity_headline.csv,
            weibull_mechanism_all.csv, table_mechanism_overview.csv
"""
import os
import sys
import io

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

from paths import ROOT, DIR_CLEAN, DIR_RESULTS, SUB_FIGURES, SUB_TABLES
RES = os.path.join(ROOT, DIR_RESULTS, SUB_TABLES)
FIG = os.path.join(ROOT, DIR_RESULTS, SUB_FIGURES)
CLEAN = os.path.join(ROOT, DIR_CLEAN)

# English-only rendering. No CJK font is registered, so any residual Chinese
# would render as tofu boxes -- an intentional, visible failure signal.
plt.rcParams["font.sans-serif"] = ["DejaVu Sans", "Arial", "Helvetica"]
plt.rcParams["axes.unicode_minus"] = False
plt.rcParams["font.size"] = 9
plt.rcParams["figure.dpi"] = 300
# Keep text as editable <text> nodes in SVG instead of glyph outlines,
# so the figures can be fine-tuned in Inkscape.
plt.rcParams["svg.fonttype"] = "none"

MECH_ORDER = ["A_HIF2a", "B_VEGFR_TKI", "C_VEGF_mAb", "D_mTOR", "E_ICI"]
MECH_LABEL = ["HIF-2\u03b1i\n(belzutifan)", "VEGFR-TKI\n(7 drugs)", "VEGF mAb\n(bevacizumab)",
              "mTORi\n(2 drugs)", "ICI\n(4 drugs)"]

# The lineage column in the result tables is already English; map it to a compact display label.
LINEAGE_EN = {
    "Erythroid": "Erythroid",
    "Leukocyte": "Leukocyte",
    "Platelet": "Platelet",
    "Multilineage/Bone marrow": "Multilineage/marrow",
    "Coagulation/Microvascular": "Coagulation/microvascular",
    "Lymphatic/Spleen": "Lymphatic/splenic",
}

COLORS = ["#B2182B", "#D6604D", "#F4A582", "#92C5DE", "#4393C3"]   # red -> blue

# Lineage ordering (keys match the lineage values stored in the CSV).
LIN_RANK = {"Erythroid": 0, "Leukocyte": 1, "Platelet": 2,
            "Multilineage/Bone marrow": 3, "Coagulation/Microvascular": 4, "Lymphatic/Spleen": 5}


def save_both(fig, out_png):
    """Write a 300-dpi PNG plus an editable SVG beside it (text stays editable)."""
    fig.savefig(out_png, dpi=300, bbox_inches="tight", facecolor="white")
    out_svg = os.path.splitext(out_png)[0] + ".svg"
    fig.savefig(out_svg, format="svg", bbox_inches="tight", facecolor="white")
    print("wrote", out_png)
    print("wrote", out_svg)


# ------------------------------------------------------------------ Fig1 heatmap
def fig_heatmap():
    df = pd.read_csv(os.path.join(RES, "table_key_pt_mechanism.csv"))
    df = df[df["n"] >= 3].copy()
    df["log2ror"] = np.log2(df["ror"].clip(lower=1e-3))

    order = [p for p in df.drop_duplicates("pt_disp")["pt_disp"]]
    # sort by haematopoietic lineage
    df["_lr"] = df["lineage"].map(LIN_RANK)
    order = (df.sort_values(["_lr", "pt_disp"]).drop_duplicates("pt_disp")["pt_disp"].tolist())

    piv = df.pivot_table(index="pt_disp", columns="mechanism", values="log2ror")
    piv = piv.reindex(index=order, columns=MECH_ORDER)
    sig = df.pivot_table(index="pt_disp", columns="mechanism", values="primary_signal")
    sig = sig.reindex(index=order, columns=MECH_ORDER).fillna(False).astype(bool)

    vmax = 5.0
    norm = TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)
    fig, ax = plt.subplots(figsize=(6.4, 9.4))
    im = ax.imshow(piv.values, cmap="RdBu_r", norm=norm, aspect="auto")

    ax.set_xticks(range(len(MECH_ORDER)))
    ax.set_xticklabels(MECH_LABEL, fontsize=7.8)
    ax.set_yticks(range(len(order)))
    ax.set_yticklabels(order, fontsize=7.6)
    ax.set_title("ROR of key haematological adverse events\n(log2 scale; red = signal enrichment)",
                 fontsize=9.5, pad=10)

    # cell annotations
    for i in range(piv.shape[0]):
        for j in range(piv.shape[1]):
            v = piv.values[i, j]
            if not np.isfinite(v):
                ax.text(j, i, "\u2014", ha="center", va="center", fontsize=6.2, color="#999999")
                continue
            ror = 2 ** v
            txt = ("%.1f" % ror) if ror >= 10 else ("%.2f" % ror)
            if sig.values[i, j]:
                txt = txt + "*"
            ax.text(j, i, txt, ha="center", va="center", fontsize=6.4,
                    color="white" if abs(v) > 2.6 else "black")

    # lineage separator lines
    lin_of = df.drop_duplicates("pt_disp").set_index("pt_disp")["lineage"].reindex(order)
    prev = None
    for i, (pt, ln) in enumerate(lin_of.items()):
        if prev is not None and ln != prev:
            ax.axhline(i - 0.5, color="white", lw=1.4)
        prev = ln
    # lineage labels
    for ln, rk in LIN_RANK.items():
        idx = [i for i, x in enumerate(lin_of.values) if x == ln]
        if idx:
            ax.text(-0.62, np.mean(idx), LINEAGE_EN.get(ln, ln), ha="center", va="center",
                    rotation=90, fontsize=8, color="#333333")

    cbar = fig.colorbar(im, ax=ax, fraction=0.030, pad=0.02)
    cbar.set_label("log2(ROR)", fontsize=8)
    cbar.ax.tick_params(labelsize=7)
    for s in ("top", "right", "left", "bottom"):
        ax.spines[s].set_visible(False)
    ax.tick_params(length=0)

    fig.tight_layout()
    save_both(fig, os.path.join(FIG, "Fig1_heatmap_ror.png"))
    plt.close(fig)


# ------------------------------------------------------------------ Fig2 forest
def fig_forest():
    df = pd.read_csv(os.path.join(RES, "sensitivity_headline.csv"))
    rows = []
    for _, r in df.iterrows():
        v = str(r["S0_base"])
        if v in ("\u2014", "nan"):
            continue
        est = float(v.split(" (")[0])
        ci = v.split("(")[1].rstrip(")*").replace("\u2013", "-").replace("\u2014", "-")
        lo, hi = [float(t) for t in ci.split("-")]
        rows.append({"label": "%s  %s" % (r["mechanism"], r["pt"]),
                     "est": est, "lo": lo, "hi": hi,
                     "sig": r["S0_base"].endswith("*")})
    d = pd.DataFrame(rows).sort_values("est")
    d = d[d["est"] > 0]

    fig, ax = plt.subplots(figsize=(6.6, 6.4))
    y = np.arange(len(d))
    ax.errorbar(d["est"], y,
                xerr=[d["est"] - d["lo"], d["hi"] - d["est"]],
                fmt="none", ecolor="#555555", elinewidth=1.0, capsize=2.5)
    colors = ["#B2182B" if s else "#999999" for s in d["sig"]]
    ax.scatter(d["est"], y, s=26, c=colors, zorder=3)

    ax.axvline(1.0, color="black", ls="--", lw=0.9)
    ax.set_xscale("log")
    ax.set_yticks(y)
    ax.set_yticklabels(d["label"], fontsize=7.6)
    ax.set_xlabel("ROR (95% CI), log scale", fontsize=9)
    ax.set_title("Headline haematological signals\n(mechanism level, restricted comparator)",
                 fontsize=9.5)
    ax.set_xlim(0.05, max(d["hi"]) * 1.6)
    ax.grid(axis="x", ls=":", color="#DDDDDD", lw=0.6)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    fig.tight_layout()
    save_both(fig, os.path.join(FIG, "Fig2_forest_headline.png"))
    plt.close(fig)


# ------------------------------------------------------------------ Fig3 TTO curves
def fig_tto():
    import importlib.util
    _s = importlib.util.spec_from_file_location("w", os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "08_weibull.py"))
    w = importlib.util.module_from_spec(_s)
    _s.loader.exec_module(w)
    merged = w.build_tto()

    fig, ax = plt.subplots(figsize=(6.4, 4.4))
    for i, m in enumerate(MECH_ORDER):
        x = np.sort(merged.loc[merged["mechanism"] == m, "tto"].values)
        if len(x) < 10:
            continue
        y = np.arange(1, len(x) + 1) / len(x)
        ax.step(x, y * 100, where="post", lw=1.6, color=COLORS[i],
                label="%s (n=%d, median %.0f d)" % (
                    MECH_LABEL[i].replace("\n", " "), len(x), np.median(x)))

    ax.set_xlim(0, 730)
    ax.set_ylim(0, 100)
    ax.set_xlabel("Time to onset (days)", fontsize=9)
    ax.set_ylabel("Cumulative BLSD events (%)", fontsize=9)
    ax.set_title("Cumulative onset of haematological adverse events\nby mechanism class",
                 fontsize=9.5)
    ax.legend(fontsize=7.4, loc="lower right", frameon=False)
    ax.grid(ls=":", color="#DDDDDD", lw=0.6)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    fig.tight_layout()
    save_both(fig, os.path.join(FIG, "Fig3_tto_cumulative.png"))
    plt.close(fig)


# ------------------------------------------------------------------ Fig4 mechanism profile
def fig_profile():
    ov = pd.read_csv(os.path.join(RES, "table_mechanism_overview.csv"))
    ov = ov.set_index("mechanism").reindex(MECH_ORDER)

    fig, axes = plt.subplots(1, 2, figsize=(8.2, 3.5))
    ax = axes[0]
    b = ax.bar(MECH_LABEL, ov["blsd_report_pct"], color=COLORS)
    for r, v in zip(b, ov["blsd_report_pct"]):
        ax.text(r.get_x() + r.get_width() / 2, v + 0.25, "%.1f%%" % v,
                ha="center", fontsize=7.6)
    ax.set_ylabel("Reports with >=1 BLSD event (%)", fontsize=8.5)
    ax.set_title("A  Share of reports containing BLSD events", fontsize=9.5)
    ax.tick_params(axis="x", labelsize=7)

    # lineage composition (count of double-positive signals)
    df = pd.read_csv(os.path.join(RES, "table_key_pt_mechanism.csv"))
    df = df[df["primary_signal"] == True]
    df["_lr"] = df["lineage"].map(LIN_RANK)
    ax2 = axes[1]
    bottom = np.zeros(len(MECH_ORDER))
    for ln, rk in sorted(LIN_RANK.items(), key=lambda kv: kv[1]):
        vals = []
        for m in MECH_ORDER:
            vals.append(int(((df["mechanism"] == m) & (df["lineage"] == ln)).sum()))
        ax2.bar(MECH_LABEL, vals, bottom=bottom, label=LINEAGE_EN.get(ln, ln),
                color=plt.cm.Set2(rk / 5.0))
        bottom += np.array(vals, dtype=float)
    ax2.set_ylabel("Number of double-positive signals", fontsize=8.5)
    ax2.set_title("B  Lineage composition of\ndouble-positive signals", fontsize=9.5)
    ax2.tick_params(axis="x", labelsize=7)
    ax2.legend(fontsize=6.8, frameon=False, ncol=2)

    for a in axes:
        for s in ("top", "right"):
            a.spines[s].set_visible(False)
    fig.tight_layout()
    save_both(fig, os.path.join(FIG, "Fig4_mechanism_profile.png"))
    plt.close(fig)


if __name__ == "__main__":
    os.makedirs(FIG, exist_ok=True)
    fig_heatmap()
    fig_forest()
    fig_profile()
    fig_tto()
    print("all figures written (PNG 300 dpi + editable SVG)")
