# -*- coding: utf-8 -*-
"""
13_ext_validation.py -- cross-database external validation (corrected version v2)
=====================================================================
Reproduce the FAERS primary signals in two independent databases -- JADER (PMDA, Japan) and
Canada Vigilance (Health Canada) -- to test signal-direction consistency (reproducibility).

Methodological notes (v2 correction)
------------------------------------
* Restricted comparator (consistent with FAERS internal): the external databases also use the
  other 14 RCC drugs as the comparator (not the whole database), ensuring apple-to-apple comparability.
  The early version used a whole-database comparator; because the whole database is full of
  chemo etc. with extremely high haematologic toxicity, it systematically depressed the RCC-drug
  signals (JADER SOC-level B/C/D/E classes ROR<1 is this artefact), now abandoned.
* 2x2 (only reports containing the 15 target RCC drugs):
      a = contains this entity (suspected) AND contains this PT
      b = contains this entity AND does not contain this PT
      c = contains the other 14 drugs AND contains this PT
      d = contains the other 14 drugs AND does not contain this PT
  Haldane-Anscombe +0.5; ROR 95% CI by Woolf method.
* Reproduced: FAERS primary signal positive (ror_low>1) AND external ROR lower bound > 1 AND a >= 3.

Event vocabulary
----------------
* Full BLSD PT = blsd_pt.csv (251 core) + blsd_pt_lab.csv (75 lab), lowercase pt_std.
* JADER: Japanese PT names -> lowercase pt_std (covers 38 core PTs; the remaining PTs lack a
  systematically built Japanese mapping and are marked N/A, only SOC-level any-BLSD is used for class-level testing).
* Canada: English PT names (token-order normalization) -> lowercase pt_std, covers the full vocabulary.

Outputs
-------
05_results/tables/ext_validation_mechanism.csv
05_results/tables/ext_validation_drug.csv
05_results/tables/ext_validation_soc.csv
05_results/tables/ext_validation_summary.csv
"""
import os
import io
import sys
import csv
import time
import re

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import pandas as pd
import numpy as np

from paths import ROOT, DIR_CLEAN, DIR_EXTERNAL, DIR_RESULTS, SUB_TABLES
EXT  = os.path.join(ROOT, DIR_EXTERNAL)
TAB  = os.path.join(ROOT, DIR_RESULTS, SUB_TABLES)
CLEAN= os.path.join(ROOT, DIR_CLEAN)

MECHANISMS = {
    "A_HIF2a":     ["belzutifan"],
    "B_VEGFR_TKI": ["sunitinib", "sorafenib", "pazopanib", "axitinib",
                    "cabozantinib", "lenvatinib", "tivozanib"],
    "C_VEGF_mAb":  ["bevacizumab"],
    "D_mTOR":      ["everolimus", "temsirolimus"],
    "E_ICI":       ["nivolumab", "ipilimumab", "pembrolizumab", "avelumab"],
}
ALL_DRUGS = [d for ds in MECHANISMS.values() for d in ds]

# ---- full BLSD PT vocabulary (lowercase) ----
def load_vocab():
    vocab = set()
    for fn in ["blsd_pt.csv", "blsd_pt_lab.csv"]:
        p = os.path.join(CLEAN, fn)
        if os.path.exists(p):
            df = pd.read_csv(p)
            if "pt_std" in df.columns:
                vocab |= set(df["pt_std"].astype(str).str.lower())
    return vocab

VOCAB = load_vocab()

# ---- JADER Japanese drug names -> INN ----
JP_DRUG = {
    "ベルズチファン":"belzutifan","スニチニブ":"sunitinib","ソラフェニブ":"sorafenib",
    "パゾパニブ":"pazopanib","アキシチニブ":"axitinib","カボザンチニブ":"cabozantinib",
    "レンバチニブ":"lenvatinib","チボザニブ":"tivozanib","ベバシズマブ":"bevacizumab",
    "エベロリムス":"everolimus","テムシロリムス":"temsirolimus","ニボルマブ":"nivolumab",
    "イピリムマブ":"ipilimumab","ペムブロリズマブ":"pembrolizumab","アベルマブ":"avelumab",
}
# ---- JADER Japanese PT -> lowercase pt_std (38 core PTs) ----
JP_PT = {
    "貧血":"anaemia","大球性貧血":"anaemia macrocytic","鉄欠乏性貧血":"iron deficiency anaemia",
    "腎性貧血":"nephrogenic anaemia","多血症":"polycythaemia","溶血性貧血":"haemolytic anaemia",
    "自己免疫性溶血性貧血":"autoimmune haemolytic anaemia",
    "微小血管障害性溶血性貧血":"microangiopathic haemolytic anaemia",
    "純赤芽球癆":"aplasia pure red cell","再生不良性貧血":"aplastic anaemia","溶血":"haemolysis",
    "好中球減少症":"neutropenia","発熱性好中球減少症":"febrile neutropenia",
    "白血球減少症":"leukopenia","無顆粒球症":"agranulocytosis","リンパ球減少症":"lymphopenia",
    "白血球増多症":"leukocytosis","好中球増多症":"neutrophilia","好酸球増多症":"eosinophilia",
    "血小板減少症":"thrombocytopenia","免疫性血小板減少症":"immune thrombocytopenia",
    "血小板増多症":"thrombocytosis","血小板障害":"platelet disorder","汎血球減少症":"pancytopenia",
    "骨髄抑制":"myelosuppression","骨髄不全":"bone marrow failure","血球減少症":"cytopenia",
    "二血球減少症":"bicytopenia","発熱性骨髄無形成":"febrile bone marrow aplasia",
    "播種性血管内凝固":"disseminated intravascular coagulation","凝固異常":"coagulopathy",
    "血栓性微小血管症":"thrombotic microangiopathy",
    "血栓性血小板減少性紫斑病":"thrombotic thrombocytopenic purpura",
    "溶血性尿毒症症候群":"haemolytic uraemic syndrome","リンパ節腫脹":"lymphadenopathy",
    "リンパ浮腫":"lymphoedema","脾腫":"splenomegaly","リンパ瘤":"lymphocele",
}

# ---- Canada English drug-name matching ----
EN_DRUG = {
    "belzutifan":["belzutifan","welireg","mk6482","pt2385"],
    "sunitinib":["sunitinib","sutent","su011248"],
    "sorafenib":["sorafenib","nexavar","bay439006"],
    "pazopanib":["pazopanib","votrient","gw786034"],
    "axitinib":["axitinib","inlyta","ag013736"],
    "cabozantinib":["cabozantinib","cabometyx","cometriq","xl184","bms907351"],
    "lenvatinib":["lenvatinib","lenvima","lenvipa","e7080"],
    "tivozanib":["tivozanib","fotivda","av951","krn951"],
    "bevacizumab":["bevacizumab","avastin","mvsi","zirabev","alymsys","vegzelma","onbevzi"],
    "everolimus":["everolimus","afinitor","rad001","sdzrad"],
    "temsirolimus":["temsirolimus","torisel","cci779"],
    "nivolumab":["nivolumab","opdivo","bms936558","mdx1106"],
    "ipilimumab":["ipilimumab","yervoy","bms734016","mdx010"],
    "pembrolizumab":["pembrolizumab","keytruda","mk3475","sch900475"],
    "avelumab":["avelumab","bavencio","msb0010718c"],
}
EN_DRUG_NORM = {d:[re.sub(r"[^a-z0-9]","",t.lower()) for t in toks]
                for d,toks in EN_DRUG.items()}
CANADA_BLSD_SOC = "Blood and lymphatic system disorders"

def norm_tokens(s):
    return " ".join(sorted(re.findall(r"[a-z0-9]+", s.lower())))

NORM2PT = {norm_tokens(p): p for p in VOCAB}

def map_en_drug(name):
    if not isinstance(name, str):
        return []
    n = re.sub(r"[^a-z0-9]","", name.lower())
    if not n:
        return []
    out = []
    for d, toks in EN_DRUG_NORM.items():
        for t in toks:
            if t and (t in n or n in t):
                out.append(d); break
    return out

def map_jp_drug(generic, brand):
    out = []
    for jp, inn in JP_DRUG.items():
        if jp in (generic or "") or jp in (brand or ""):
            out.append(inn)
    return out


# =====================================================================
def parse_jader():
    print("[%s] Parsing JADER ..." % time.strftime("%H:%M:%S"))
    jader = os.path.join(EXT, "JADER")
    drug_case, pt_case, blsd_cases = {}, {}, set()
    with open(os.path.join(jader, "drug202608.csv"), encoding="cp932", newline="") as f:
        for row in csv.reader(f):
            if len(row) < 6: continue
            if row[3] != "被疑薬": continue
            for d in map_jp_drug(row[4], row[5]):
                drug_case.setdefault(d, set()).add(row[0])
    with open(os.path.join(jader, "reac202608.csv"), encoding="cp932", newline="") as f:
        for row in csv.reader(f):
            if len(row) < 4: continue
            p = JP_PT.get(row[3])
            if p:
                pt_case.setdefault(p, set()).add(row[0]); blsd_cases.add(row[0])
    N = set()
    with open(os.path.join(jader, "demo202608.csv"), encoding="cp932", newline="") as f:
        for row in csv.reader(f):
            if row: N.add(row[0])
    print("  JADER: target-drug reports=%d, BLSD PT(38) reports=%d, N=%d"
          % (sum(len(v) for v in drug_case.values()), len(blsd_cases), len(N)))
    return drug_case, pt_case, blsd_cases, len(N)


def parse_canada():
    print("[%s] Parsing Canada Vigilance ..." % time.strftime("%H:%M:%S"))
    cv = os.path.join(EXT, "CanadaVigilance")
    drug_case, pt_case, blsd_cases = {}, {}, set()
    with open(os.path.join(cv, "report_drug.txt"), encoding="utf-8", newline="") as f:
        for row in csv.reader(f, delimiter="$"):
            if len(row) < 5: continue
            if "Suspect" not in row[4]: continue
            for d in map_en_drug(row[3]):
                drug_case.setdefault(d, set()).add(row[0])
    with open(os.path.join(cv, "reactions.txt"), encoding="utf-8", newline="") as f:
        for row in csv.reader(f, delimiter="$"):
            if len(row) < 8: continue
            if row[7] == CANADA_BLSD_SOC:
                blsd_cases.add(row[1])
            p = NORM2PT.get(norm_tokens(row[5]))
            if p:
                pt_case.setdefault(p, set()).add(row[1])
    N = set()
    with open(os.path.join(cv, "reports.txt"), encoding="utf-8", newline="") as f:
        for row in csv.reader(f, delimiter="$"):
            if len(row) >= 2 and row[1]: N.add(row[1])
    print("  Canada: target-drug reports=%d, BLSD(any) reports=%d, BLSD PT reports=%d, N=%d"
          % (sum(len(v) for v in drug_case.values()), len(blsd_cases),
             sum(len(v) for v in pt_case.values()), len(N)))
    return drug_case, pt_case, blsd_cases, len(N)


def ror_restricted(ent, other, pt):
    """restricted comparator: ent vs other (the other 14 drugs); pt = case set for this PT."""
    if not ent:
        return (np.nan, np.nan, np.nan, 0)
    a = len(ent & pt); b = len(ent) - a
    c = len(other & pt); d = len(other) - c
    aa, bb, cc, dd = a+0.5, b+0.5, c+0.5, d+0.5
    ror = (aa*dd)/(bb*cc)
    se = np.sqrt(1/aa+1/bb+1/cc+1/dd)
    return (ror, np.exp(np.log(ror)-1.96*se), np.exp(np.log(ror)+1.96*se), a)


def mech_union(drug_case, mech):
    s = set()
    for d in MECHANISMS[mech]:
        s |= drug_case.get(d, set())
    return s


def main():
    t0 = time.time()
    mech_ref = pd.read_csv(os.path.join(TAB, "signal_mechanism.csv"))
    drug_ref = pd.read_csv(os.path.join(TAB, "signal_drug.csv"))

    jd, jp, jb, jN = parse_jader()
    cd, cp, cb, cN = parse_canada()
    j_all = set().union(*jd.values()) if jd else set()
    c_all = set().union(*cd.values()) if cd else set()

    # ---- mechanism-level primary signals ----
    rows_m = []
    for _, r in mech_ref[mech_ref["primary_signal"] == True].iterrows():
        m, pt = r["mechanism"], str(r["pt_std"]).lower()
        ent_j, ent_c = mech_union(jd, m), mech_union(cd, m)
        j_ror, j_lo, _, j_n = ror_restricted(ent_j, j_all-ent_j, jp.get(pt, set()))
        c_ror, c_lo, _, c_n = ror_restricted(ent_c, c_all-ent_c, cp.get(pt, set()))
        fs = bool(r["ror_low"] > 1)
        j_appl = pt in jp
        j_rep = bool(fs and j_appl and j_lo > 1 and j_n >= 3)
        c_rep = bool(fs and c_lo > 1 and c_n >= 3)
        rows_m.append({
            "mechanism": m, "pt_std": pt, "pt_disp": r["pt_disp"],
            "faers_ror": round(r["ror"], 2), "faers_low": round(r["ror_low"], 2), "faers_n": int(r["n"]),
            "jader_applicable": j_appl,
            "jader_ror": round(j_ror, 2) if j_appl and j_n else "",
            "jader_low": round(j_lo, 2) if j_appl and j_n else "",
            "jader_n": j_n, "jader_reproduced": j_rep,
            "canada_ror": round(c_ror, 2) if c_n else "",
            "canada_low": round(c_lo, 2) if c_n else "",
            "canada_n": c_n, "canada_reproduced": c_rep,
            "db_reproduced": int(j_rep)+int(c_rep),
        })
    out_m = pd.DataFrame(rows_m)
    out_m.to_csv(os.path.join(TAB, "ext_validation_mechanism.csv"), index=False, encoding="utf-8-sig")

    # ---- drug-level primary signals ----
    rows_d = []
    for _, r in drug_ref[drug_ref["primary_signal"] == True].iterrows():
        d, pt = r["drug_std"], str(r["pt_std"]).lower()
        ent_j, ent_c = jd.get(d, set()), cd.get(d, set())
        j_ror, j_lo, _, j_n = ror_restricted(ent_j, j_all-ent_j, jp.get(pt, set()))
        c_ror, c_lo, _, c_n = ror_restricted(ent_c, c_all-ent_c, cp.get(pt, set()))
        fs = bool(r["ror_low"] > 1)
        j_appl = pt in jp
        j_rep = bool(fs and j_appl and j_lo > 1 and j_n >= 3)
        c_rep = bool(fs and c_lo > 1 and c_n >= 3)
        rows_d.append({
            "drug_std": d, "mechanism": r["mechanism"], "pt_std": pt, "pt_disp": r["pt_disp"],
            "faers_ror": round(r["ror"], 2), "faers_low": round(r["ror_low"], 2), "faers_n": int(r["n"]),
            "jader_applicable": j_appl,
            "jader_ror": round(j_ror, 2) if j_appl and j_n else "",
            "jader_low": round(j_lo, 2) if j_appl and j_n else "",
            "jader_n": j_n, "jader_reproduced": j_rep,
            "canada_ror": round(c_ror, 2) if c_n else "",
            "canada_low": round(c_lo, 2) if c_n else "",
            "canada_n": c_n, "canada_reproduced": c_rep,
            "db_reproduced": int(j_rep)+int(c_rep),
        })
    out_d = pd.DataFrame(rows_d)
    out_d.to_csv(os.path.join(TAB, "ext_validation_drug.csv"), index=False, encoding="utf-8-sig")

    # ---- mechanism-level any-BLSD (restricted comparator, class-level signal) ----
    rows_s = []
    for m in ["A_HIF2a","B_VEGFR_TKI","C_VEGF_mAb","D_mTOR","E_ICI"]:
        ent_j, ent_c = mech_union(jd, m), mech_union(cd, m)
        j_ror, j_lo, _, j_n = ror_restricted(ent_j, j_all-ent_j, jb)
        c_ror, c_lo, _, c_n = ror_restricted(ent_c, c_all-ent_c, cb)
        rows_s.append({
            "mechanism": m,
            "jader_blsd_ror": round(j_ror, 2) if j_n else "", "jader_blsd_low": round(j_lo, 2) if j_n else "", "jader_n": j_n,
            "canada_blsd_ror": round(c_ror, 2) if c_n else "", "canada_blsd_low": round(c_lo, 2) if c_n else "", "canada_n": c_n,
        })
    pd.DataFrame(rows_s).to_csv(os.path.join(TAB, "ext_validation_soc.csv"), index=False, encoding="utf-8-sig")

    summary = {
        "mechanism_primary_signals": len(out_m),
        "mechanism_reproduced_canada": int(out_m["canada_reproduced"].sum()),
        "mechanism_reproduced_jader": int(out_m["jader_reproduced"].sum()),
        "mechanism_reproduced_both": int((out_m["db_reproduced"] == 2).sum()),
        "drug_primary_signals": len(out_d),
        "drug_reproduced_canada": int(out_d["canada_reproduced"].sum()),
        "drug_reproduced_jader": int(out_d["jader_reproduced"].sum()),
        "drug_reproduced_both": int((out_d["db_reproduced"] == 2).sum()),
    }
    pd.DataFrame([summary]).to_csv(os.path.join(TAB, "ext_validation_summary.csv"), index=False, encoding="utf-8-sig")

    print("\n===== External validation summary (restricted comparator) =====")
    for k, v in summary.items():
        print(f"  {k}: {v}")
    print("Elapsed %.1f s" % (time.time() - t0))


if __name__ == "__main__":
    main()
