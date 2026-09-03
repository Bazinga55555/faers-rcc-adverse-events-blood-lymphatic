# -*- coding: utf-8 -*-
"""
07_blsd_pt.py -- build the "blood and lymphatic system disorders (BLSD)" PT list + PT name normalization

Background: the public FAERS ASCII REAC table has only PT (Preferred Term), no SOC field,
  so PTs must be mapped to MedDRA SOC 'Blood and lymphatic system disorders' (10005329) on our own.
  This script produces an auditable list in two steps:
    Step A  PT name normalization (case/whitespace) to eliminate one PT being split across multiple rows
    Step B  regex coarse screening -> rule-based classification (core / lab / exclude)

Outputs:
  03_clean_data/reac_norm.parquet       normalized REAC (with pt_std / pt_disp)
  03_clean_data/blsd_pt.csv            core PT list (for the main analysis)
  03_clean_data/blsd_pt_lab.csv        lab PT list (for sensitivity analysis)
  05_results/tables/blsd_pt_audit.csv  all candidate PTs and their classification reasons
                                      (for review; can be used directly as a supplementary table)

Classification:
  core    = PTs clearly within BLSD SOC, diagnostic type (anaemia / cytopenia / coagulopathy /
           haemolysis / haemoglobinopathy / bone-marrow failure / spleen / lymphatic)
  lab     = laboratory-investigation PTs (strictly within Investigations SOC), used only in sensitivity analysis
  exclude = removed (neoplasms, site-specific bleeding, thrombosis, skin purpura, infection,
           organ-specific eosinophilic disorders, etc.)
"""
import os
import sys
import io
import re

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")

import pandas as pd

from paths import ROOT, DIR_CLEAN, DIR_RESULTS, SUB_TABLES
CLEAN = os.path.join(ROOT, DIR_CLEAN)
RES = os.path.join(ROOT, DIR_RESULTS, SUB_TABLES)

# =============================================================== Step A: normalization
def norm_pt(s: pd.Series) -> pd.Series:
    return (s.astype(str).str.strip()
             .str.replace(r"\s+", " ", regex=True)
             .str.lower())


# =============================================================== Step B: coarse screen (broad)
INCLUDE = re.compile(r"|".join([
    # ---- anaemia / erythrocyte
    r"anaemi", r"anemi",
    r"hypochrom", r"macrocytic", r"microcytic", r"normocytic", r"normochrom",
    r"sideroblast", r"spherocyt", r"elliptocyt", r"reticulocyt",
    r"erythrocyt", r"erythropen", r"erythroblastopen",
    r"red cell", r"red blood cell",
    r"iron deficien", r"iron-deficien", r"b12 deficien", r"vitamin b12",
    r"folate deficien", r"folic acid deficien",
    r"megaloblast", r"pernicious", r"aplastic", r"aplasia",
    r"hypoplasia", r"pure red cell",
    # ---- haemoglobin / haemolysis / haemoglobinopathy
    r"haemoglobin", r"hemoglobin",
    r"haemolys", r"hemolys", r"haemolyt", r"hemolyt",
    r"methaemoglobin", r"methemoglobin", r"sulphaemoglobin", r"sulfhemoglobin",
    r"paroxysmal nocturnal",
    r"haemosideros", r"hemosideros", r"haemochromatos", r"hemochromatos",
    r"porphyr", r"thalassaem", r"thalassem", r"sickle",
    r"haemoglobinopath", r"hemoglobinopath",
    # ---- cytopenia / bone marrow
    r"neutropen", r"agranulocyt", r"granulocytopen",
    r"leukopen", r"leucopen", r"lymphopen", r"thrombocytopen",
    r"pancytopen", r"cytopen", r"monocytopen", r"eosinopen", r"basopen",
    r"myelosuppress", r"myelosuppression", r"bone marrow",
    r"marrow failure", r"marrow suppression", r"marrow aplasia",
    r"haematopoie", r"hematopoie", r"haemopoie",
    # ---- cytosis
    r"leukocytos", r"leucocytos", r"lymphocytos", r"monocytos",
    r"eosinophili", r"basophili", r"neutrophili",
    r"thrombocytos", r"thrombocythaemi", r"thrombocythemi",
    r"polycythaemi", r"polycythemi",
    # ---- coagulation / platelet
    r"coagul", r"disseminated intravascular",
    r"antithrombin", r"von willebrand", r"haemophili", r"hemophili",
    r"factor (?:i{1,3}|iv|v|vi|vii|viii|ix|x|xi|xii|xiii) deficien",
    r"protein c deficien", r"protein s deficien",
    r"fibrinogen", r"afibrinogen", r"hypofibrinogen", r"dysfibrinogen",
    r"antiphospholipid", r"lupus anticoagulant",
    r"thrombotic microangiopath", r"thrombotic thrombocytopenic",
    r"haemolytic uraemic", r"hemolytic uremic",
    r"platelet", r"thrombocytopath", r"thrombastheni",
    r"purpura", r"petechi", r"ecchymos", r"bruising", r"bruise",
    # ---- spleen / lymphatic
    r"splenomegal", r"splenic", r"hypersplen", r"asplen", r"spleen",
    r"lymphadenopath", r"lymphoedema", r"lymphedema", r"lymphocele",
    r"lymph node", r"lymphang", r"lymphatic", r"lymph gland",
]), re.I)

# ---------------------------------------------- Exclude 1: neoplasms / proliferative disease (Neoplasms SOC)
EXCL_NEOPLASM = re.compile(r"|".join([
    r"leukaemi", r"leukemi", r"leukaemoid", r"leukemoid",
    r"lymphoma", r"myeloma", r"myelodysplas", r"myeloproliferat",
    r"myelofibros", r"waldenstrom", r"castleman",
    r"histiocytos", r"mastocytos", r"lymphangioleiomyomatosis",
    r"polycythaemia vera", r"polycythemia vera",
    r"essential thrombocythaemi", r"essential thrombocythemi",
    r"chronic myeloid", r"acute myeloid", r"acute lymphoblastic",
    r"acute lymphocytic", r"chronic lymphocytic", r"chronic myelomonocytic",
    r"blast crisis", r"chloroma", r"granulocytic sarcoma",
    r"refractory anaemia", r"refractory anemia", r"excess of blasts",
    r"neoplasm", r"tumour", r"tumor", r"carcinoma", r"sarcoma",
    r"metastas", r"malignan", r"cancer", r"adenocarcinoma", r"blastoma",
    r"carcinomatosa",
]), re.I)

# ---------------------------------------------- Exclude 2: site-specific bleeding / haematoma (non-BLSD SOC)
EXCL_BLEED = re.compile(r"|".join([
    r"haemorrhage", r"hemorrhage", r"haemorrhagic", r"hemorrhagic",
    r"bleed", r"haematoma", r"hematoma", r"blood loss",
    r"haemoptys", r"hemoptys", r"haematemes", r"hematemes",
    r"melena", r"melaena", r"epistaxis", r"contusion",
    r"haematuria", r"hematuria", r"petechi",
]), re.I)

# ---------------------------------------------- Exclude 3: thrombosis / embolism / ischaemia (Vascular SOC)
EXCL_THROMBOSIS = re.compile(r"|".join([
    r"thrombos", r"thromboembol", r"embolism", r"embolus",
    r"deep vein", r"phlebitis", r"thrombophlebitis",
    r"stroke", r"infarct", r"ischaemi", r"ischemi",
]), re.I)

# ---------------------------------------------- Exclude 4: infection / procedure / non-AE
EXCL_PROC_INF = re.compile(r"|".join([
    r"transfusion", r"transplant",
    r"bone marrow (?:biopsy|aspirat|examin|harvest|culture)",
    r"splenectom",
    r"splenic (?:injur|ruptur|infarct|laceration|cyst|mass|lesion|necrosis|abscess|haemorrhage)",
    r"blood (?:group|type|donor|sample|test|culture|transfusion)",
    r"anticoagulant therapy", r"anticoagulation", r"drug monitoring",
    r"^blood$", r"^bone marrow$",
    r"neutropenic (?:sepsis|infection)", r"neutropenic colitis",
    r"sepsis", r"septic", r"\binfection\b", r"abscess", r"pneumonia",
    r"tuberculosis", r"infected lymphocele",
    r"acute febrile neutrophilic dermatosis", r"drug reaction with eosinophilia",
    r"eosinophilic (?:fasciitis|pneumonia|myocarditis|oesophagitis|esophagitis|"
    r"gastritis|gastroenteritis|colitis|dermatitis|granuloma)",
    r"eosinophilic$",                    # e.g. "Gastroenteritis eosinophilic"
    r"^lymphangitis$", r"necrotic lymphadenopathy", r"necrotising lymphadenitis",
    r"granulomatosis with polyangiitis",
    r"metaplasia", r"bone marrow oedema", r"bone marrow infiltration",
    r"red blood cells urine", r"^iron deficiency$",
    r"^vitamin b12 (?:decreased|abnormal|increased)$",
    r"^pseudoporphyria$", r"congenital",
    r"candidiasis", r"aspergillosis", r"mycosis",
    r"fungaemia", r"fungemia", r"septicaemia", r"septicemia",
    r"stem cell mobilisation", r"stem cell mobilization",
    r"apheresis", r"harvest", r"mobilisation", r"mobilization",
    r"lymphangiectasia intestinal", r"intestinal lymphangiectasia",
    r"exposure", r"off label use", r"product used",
    r"condition aggravated", r"disease progression",
]), re.I)

# ---------------------------------------------- Exclude 5: skin-type purpura / ecchymosis (Skin SOC)
#   exception: thrombocytopenic purpura / ITP / TTP belong to BLSD
SKIN_PURPURA = re.compile(r"purpura|petechi|ecchymos|bruising|bruise", re.I)
KEEP_IF_TCP = re.compile(r"thrombocytopenic purpura", re.I)

# ---------------------------------------------- lab-investigation PTs (Investigations SOC)
LAB_PAT = re.compile(r"|".join([
    r"(?:count|level|index|haemoglobin|hemoglobin|haematocrit|hematocrit|"
    r"platelet|white blood|red blood|neutrophil|lymphocyte|monocyte|"
    r"eosinophil|basophil|reticulocyte|erythrocyte|leukocyte|leucocyte|"
    r"blood|full blood|fbc|wbc|rbc|hb|hgb)"
    r".*(?:increased|decreased|abnormal|elevated|reduced|low|high|fluctuat|"
    r"present|absent|normal|raised|depressed|below|above|prolonged|shortened)",
    r"(?:increased|decreased|abnormal|elevated|reduced|low|high|raised|depressed|prolonged)"
    r".*(?:count|level|haemoglobin|hemoglobin|haematocrit|hematocrit|"
    r"platelet|white blood|red blood|neutrophil|lymphocyte|monocyte|"
    r"eosinophil|basophil|reticulocyte|erythrocyte|leukocyte|leucocyte|blood)",
    r"^blood (?:count|test|film|smear|picture)", r"^full blood count",
    r"bone marrow (?:biopsy|aspirat|examin|culture|finding)",
    r"^laboratory", r"^investigation", r"differential",
    r"^haemoglobin$", r"^haematocrit$", r"^platelet count$",
    r"red cell distribution width", r"mean cell haemoglobin",
    r"mean platelet volume", r"red blood cell sedimentation rate",
    r"sedimentation rate", r"lymph nodes scan", r"coagulation time",
    r"coagulation test", r"glycosylated haemoglobin",
]), re.I)

# ---------------------------------------------- whitelist (highest priority, overrides the excludes above)
#   these are the most important haematologic endpoints of chemo/targeted therapy, clearly within BLSD SOC
WHITELIST = re.compile(r"|".join([
    r"^febrile neutropenia$",
    r"^febrile bone marrow aplasia$",
    r"^anaemia$", r"^neutropenia$", r"^thrombocytopenia$",
    r"^leukopenia$", r"^pancytopenia$", r"^lymphopenia$",
    r"^agranulocytosis$", r"^aplastic anaemia$", r"^aplastic anemia$",
    r"^bone marrow failure$", r"^myelosuppression$", r"^cytopenia$",
    r"^disseminated intravascular coagulation$",
    r"^thrombotic microangiopathy$", r"^thrombotic thrombocytopenic purpura$",
    r"^idiopathic thrombocytopenic purpura$",
    r"^immune thrombocytopenic purpura$", r"^immune thrombocytopenia$",
    r"^haemolytic anaemia$", r"^hemolytic anemia$",
    r"^autoimmune haemolytic anaemia$",
    r"^aplastic anaemia$", r"^pure red cell aplasia$",
    r"^hypereosinophilic syndrome$",
]), re.I)


def classify(pt: str):
    """Classify by priority, return (category, reason)."""
    s = str(pt).strip().lower()

    if WHITELIST.search(s):
        return "core", "whitelist: definite BLSD SOC term"
    if EXCL_NEOPLASM.search(s):
        return "exclude", "neoplasm/proliferative (Neoplasms SOC)"
    if EXCL_THROMBOSIS.search(s):
        return "exclude", "thrombosis/embolism/ischaemia (Vascular SOC)"
    if EXCL_BLEED.search(s):
        return "exclude", "site-specific haemorrhage/haematoma (non-BLSD SOC)"
    if EXCL_PROC_INF.search(s):
        return "exclude", "infection/procedure/non-AE/organ-specific term"
    if SKIN_PURPURA.search(s) and not KEEP_IF_TCP.search(s):
        return "exclude", "skin purpura/ecchymosis/petechiae (Skin SOC)"
    if LAB_PAT.search(s):
        return "lab", "laboratory/investigation term (Investigations SOC)"
    return "core", "BLSD SOC diagnosis term"


def main():
    reac = pd.read_parquet(os.path.join(CLEAN, "reac.parquet"))
    print("Total REAC rows: %d, original unique PT count: %d" % (len(reac), reac["pt"].nunique()))

    # ---- Step A: normalization
    reac["pt_std"] = norm_pt(reac["pt"])
    print("Unique PT count after normalization: %d (merged %d case/whitespace variants)"
          % (reac["pt_std"].nunique(), reac["pt"].nunique() - reac["pt_std"].nunique()))

    # display name: the most frequent original spelling under this pt_std
    disp = (reac.groupby("pt_std")["pt"]
                .agg(lambda s: s.value_counts().index[0])
                .rename("pt_disp"))
    reac = reac.merge(disp, on="pt_std", how="left")
    reac[["primaryid", "pt_std", "pt_disp"]].to_parquet(
        os.path.join(CLEAN, "reac_norm.parquet"), index=False)
    print("Wrote reac_norm.parquet")

    # ---- Step B: classification
    vc = reac["pt_std"].value_counts().rename_axis("pt_std").reset_index(name="n")
    vc["hit"] = vc["pt_std"].str.contains(INCLUDE)
    cand = vc[vc["hit"]].merge(disp.rename_axis("pt_std"), on="pt_std", how="left")
    print("Coarse-screen candidate PTs: %d (covering %d report rows)" % (len(cand), int(cand["n"].sum())))

    cls = cand["pt_std"].apply(classify)
    cand["category"] = [c[0] for c in cls]
    cand["reason"] = [c[1] for c in cls]

    os.makedirs(RES, exist_ok=True)
    cand.sort_values(["category", "n"], ascending=[True, False]).to_csv(
        os.path.join(RES, "blsd_pt_audit.csv"), index=False, encoding="utf-8-sig")

    core = cand[cand["category"] == "core"].copy()
    lab = cand[cand["category"] == "lab"].copy()
    core.to_csv(os.path.join(CLEAN, "blsd_pt.csv"), index=False, encoding="utf-8-sig")
    lab.to_csv(os.path.join(CLEAN, "blsd_pt_lab.csv"), index=False, encoding="utf-8-sig")

    print("\n=== Classification summary ===")
    print(cand["category"].value_counts().to_string())
    print("  core BLSD PTs : %d (covering %d report rows)" % (len(core), int(core["n"].sum())))
    print("  lab-type PTs  : %d (covering %d report rows)" % (len(lab), int(lab["n"].sum())))

    print("\n=== Core BLSD PTs (reports >= 30, total %d) ===" % int((core["n"] >= 30).sum()))
    print(core[core["n"] >= 30][["pt_disp", "n"]].to_string(index=False))

    print("\n=== Excluded high-frequency PTs (top 50, for review) ===")
    exc = cand[cand["category"] == "exclude"].sort_values("n", ascending=False)
    print(exc.head(50)[["pt_disp", "n", "reason"]].to_string(index=False))

    print("\n=== Lab-type PTs (top 30) ===")
    print(lab.head(30)[["pt_disp", "n"]].to_string(index=False))


if __name__ == "__main__":
    main()
