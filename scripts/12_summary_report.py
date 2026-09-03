# -*- coding: utf-8 -*-
"""
12_summary_report.py -- aggregate this analysis's "Results section materials" into paste-ready markdown

Outputs:
  05_results/Results_materials.md
  05_results/Results_materials.csv  (flat table of all primary signals)
"""
import io
import os
import sys
import pandas as pd

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
from paths import ROOT, DIR_RESULTS, SUB_TABLES
RES  = os.path.join(ROOT, DIR_RESULTS, SUB_TABLES)
OUT  = os.path.join(ROOT, DIR_RESULTS)

MECH_ORDER = ["A_HIF2a", "B_VEGFR_TKI", "C_VEGF_mAb", "D_mTOR", "E_ICI"]
MECH_LABEL = {
    "A_HIF2a":   "HIF-2α inhibitor (belzutifan, 1 drug)",
    "B_VEGFR_TKI":"VEGFR-TKI (7 drugs)",
    "C_VEGF_mAb": "VEGF mAb (bevacizumab, 1 drug)",
    "D_mTOR":    "mTOR inhibitor (2 drugs)",
    "E_ICI":     "Immune checkpoint inhibitor (4 drugs)",
}


def fmt_ci(lo, hi):
    return f"({lo:.2f}\u2013{hi:.2f})"


def main():
    ov    = pd.read_csv(os.path.join(RES, "table_mechanism_overview.csv"))
    mt    = pd.read_csv(os.path.join(RES, "table_key_pt_mechanism.csv"))
    sh    = pd.read_csv(os.path.join(RES, "sensitivity_headline.csv"))
    w_all = pd.read_csv(os.path.join(RES, "weibull_mechanism_all.csv"))
    w_pt  = pd.read_csv(os.path.join(RES, "weibull_mechanism.csv"))

    sig = mt[mt["primary_signal"] == True].copy()
    sig = sig.sort_values(["mechanism", "ror"], ascending=[True, False])

    # flat table of all primary signals (CSV)
    flat = sig[["mechanism", "lineage", "pt_disp", "n",
                "ror", "ror_low", "ror_high", "ic", "ic025",
                "ebgm", "ebgm05"]].copy()
    flat.to_csv(os.path.join(OUT, "Results_materials.csv"), index=False, encoding="utf-8-sig")
    print("Wrote Results_materials.csv:", len(flat), "rows")

    # write markdown report
    lines = []
    lines.append("# Results materials: five mechanism classes x BLSD signals / TTO / sensitivity\n")
    lines.append("> Data: FAERS/AEMS 2004 Q1 \u2013 2026 Q2, 5 mechanism classes, 15 RCC drugs, n=522,813 reports.\n")
    lines.append("> Outcome scope: core PTs of MedDRA BLSD SOC (10005329) (251).\n")
    lines.append("> Primary criterion: ROR 95% CI lower bound > 1 **and** IC025 > 0 dual-positive (marked \"*\").\n")
    lines.append("> restricted comparator: the other 14 target RCC drugs (controls for indication confounding and protopathic bias).\n")
    lines.append("\n---\n")

    # ---- 3.1 cohort overview ----
    col_map = {"n_reports": "reports", "n_blsd_events": "events",
               "n_reports_with_blsd": "blsd_reports",
               "blsd_report_pct": "blsd_pct",
               "unique_blsd_pt": "unique_pt"}
    ov2 = ov.rename(columns=col_map)
    lines.append("## 3.1 Cohort overview\n")
    lines.append("| Mechanism | Reports | BLSD events | Reports w/ BLSD | BLSD report % | Unique BLSD PT |")
    lines.append("|---|---:|---:|---:|---:|---:|")
    for m in MECH_ORDER:
        r = ov2[ov2["mechanism"] == m].iloc[0]
        lines.append("| %s | %d | %d | %d | %.2f | %d |"
                     % (MECH_LABEL[m], int(r['reports']), int(r['events']),
                        int(r['blsd_reports']), r['blsd_pct'], int(r['unique_pt'])))
    lines.append("")
    lines.append("**HIF-2αi has the highest report share (19.78%), suggesting belzutifan haematologic events are highly focused reporting.**\n")

    # ---- 3.2 key PT signals ----
    lines.append("## 3.2 Key PT signals (by mechanism)\n")
    for m in MECH_ORDER:
        sub = sig[sig["mechanism"] == m]
        lines.append(f"### {MECH_LABEL[m]} ({len(sub)} primary signals)\n")
        lines.append("| PT | n | ROR (95% CI) | IC (IC025) | EBGM (EB05) |")
        lines.append("|---|---:|---|---|---|")
        for _, r in sub.iterrows():
            lines.append(f"| {r['pt_disp']} | {int(r['n']):,} | "
                         f"{r['ror']:.2f} {fmt_ci(r['ror_low'], r['ror_high'])} | "
                         f"{r['ic']:.2f} ({r['ic025']:.2f}) | "
                         f"{r['ebgm']:.2f} ({r['ebgm05']:.2f}) |")
        lines.append("")

    # ---- 3.3 H1 ----
    lines.append("## 3.3 H1 test: upstream vs downstream opposite directions on the erythroid axis\n")
    ana = sig[sig["pt_disp"] == "Anaemia"]
    poly = sig[sig["pt_disp"] == "Polycythaemia"]
    lines.append("| Mechanism | Anaemia ROR (95% CI) | Polycythaemia ROR (95% CI) | Interpretation |")
    lines.append("|---|---|---|---|")
    for m in MECH_ORDER:
        a = ana[ana["mechanism"] == m]
        p = poly[poly["mechanism"] == m]
        astr = (f"{a['ror'].iloc[0]:.2f} {fmt_ci(a['ror_low'].iloc[0], a['ror_high'].iloc[0])}" if len(a) else "\u2014")
        pstr = (f"{p['ror'].iloc[0]:.2f} {fmt_ci(p['ror_low'].iloc[0], p['ror_high'].iloc[0])}" if len(p) else "\u2014")
        interp = ""
        if m == "A_HIF2a":
            interp = "**upstream EPO transcription shutdown \u2192 very strong anaemia signal**"
        elif m == "B_VEGFR_TKI":
            interp = "downstream hypoxia blockade \u2192 HIF stabilization \u2192 compensatory EPO rise \u2192 erythrocytosis"
        lines.append(f"| {MECH_LABEL[m]} | {astr} | {pstr} | {interp} |")
    lines.append("")
    lines.append("**Group A Anaemia ROR=30.00 (21.23\u201342.39), Group B Polycythaemia ROR=7.50 (5.43\u201310.35), opposite directions, H1 holds.**\n")

    # ---- 3.4 H3 Weibull ----
    lines.append("## 3.4 H3 test: Weibull TTO typing\n")
    lines.append("| Mechanism | All BLSD n | Median TTO (d) | beta (95% CI) | Type |")
    lines.append("|---|---:|---:|---|---|")
    for _, r in w_all.iterrows():
        lines.append(f"| {MECH_LABEL[r['mechanism']]} | {int(r['n']):,} | {r['median_tto']:.0f} | "
                     f"{r['beta']:.3f} ({r['beta_low95']:.3f}\u2013{r['beta_high95']:.3f}) | {r['type']} |")
    lines.append("")
    lines.append("**All 5 groups are early_failure (beta<1, CI upper<1), suggesting haematologic events occur early on therapy (indication confounding / first-dose effect). Median TTO: mTOR 56 d \u2192 belzutifan 41 d \u2192 ICI 40 d \u2192 VEGF mAb 37 d \u2192 VEGFR-TKI 36 d.**\n")

    wearout = w_pt[w_pt["type"] == "wearout_failure"].sort_values("median_tto", ascending=False)
    if len(wearout):
        lines.append("**Wearout_failure (2 items, cumulative / delayed):**\n")
        lines.append("| Mechanism x PT | n | Median TTO (d) | beta (95% CI) |")
        lines.append("|---|---:|---:|---|")
        for _, r in wearout.iterrows():
            lines.append(f"| {MECH_LABEL[r['mechanism']]} x {r['pt']} | {int(r['n'])} | {r['median_tto']:.0f} | "
                         f"{r['beta']:.3f} ({r['beta_low95']:.3f}\u2013{r['beta_high95']:.3f}) |")
        lines.append("")

    # ---- 3.5 sensitivity ----
    lines.append("## 3.5 Sensitivity analysis (five scenarios)\n")
    lines.append("| Mechanism | PT | S0 base | S1 lab | S2 RCC indication | S3 PS only | S4 single-drug |")
    lines.append("|---|---|---|---|---|---|")
    for _, r in sh.iterrows():
        lines.append(f"| {MECH_LABEL[r['mechanism']]} | {r['pt']} | "
                     f"{r['S0_base']} | {r['S1_lab']} | {r['S2_rccInd']} | {r['S3_psOnly']} | {r['S4_singleDrug']} |")
    lines.append("")
    lines.append("**Key robustness:**\n")
    lines.append("- **belzutifan x Anaemia**: robust across all 5 scenarios (10.48 \u2013 33.54)")
    lines.append("- **ICI immune cytopenia spectrum** (PRCA / AIHA / ITP): robust across all 5 scenarios")
    lines.append("- **VEGFR-TKI x Polycythaemia**: attenuates to 1.18 (0.74\u20131.86) under S2 RCC indication, **downgrade to exploratory conclusion** (confounded by RCC's paraneoplastic erythrocytosis)")
    lines.append("- **bevacizumab x TMA**: amplifies to 12.42 (8.29\u201318.61) under S2 indication, suggesting TMA is highly concentrated in the RCC-indication population")
    lines.append("")

    out_md = os.path.join(OUT, "Results_materials.md")
    with io.open(out_md, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))
    print("Wrote", out_md, len(lines), "lines")


if __name__ == "__main__":
    main()
