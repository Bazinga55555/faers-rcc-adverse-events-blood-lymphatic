# -*- coding: utf-8 -*-
"""
pt_defs.py -- key BLSD PT list and mechanism-group order (shared by 06/08/09 and others)

Note: this module must NOT do any sys.stdout redirection and must have no side effects,
  so it can be safely imported by other scripts.
"""

# Key PTs: (PT name, haematologic lineage). Names match the MedDRA PT text exactly (pt_disp)
KEY_PT = [
    # Erythroid
    ("Anaemia", "Erythroid"),
    ("Anaemia macrocytic", "Erythroid"),
    ("Iron deficiency anaemia", "Erythroid"),
    ("Nephrogenic anaemia", "Erythroid"),
    ("Polycythaemia", "Erythroid"),
    ("Haemolytic anaemia", "Erythroid"),
    ("Autoimmune haemolytic anaemia", "Erythroid"),
    ("Microangiopathic haemolytic anaemia", "Erythroid"),
    ("Aplasia pure red cell", "Erythroid"),
    ("Aplastic anaemia", "Erythroid"),
    ("Haemolysis", "Erythroid"),
    # Leukocyte
    ("Neutropenia", "Leukocyte"),
    ("Febrile neutropenia", "Leukocyte"),
    ("Leukopenia", "Leukocyte"),
    ("Agranulocytosis", "Leukocyte"),
    ("Lymphopenia", "Leukocyte"),
    ("Leukocytosis", "Leukocyte"),
    ("Neutrophilia", "Leukocyte"),
    ("Eosinophilia", "Leukocyte"),
    # Platelet
    ("Thrombocytopenia", "Platelet"),
    ("Immune thrombocytopenia", "Platelet"),
    ("Thrombocytosis", "Platelet"),
    ("Platelet disorder", "Platelet"),
    # Multilineage/Bone marrow
    ("Pancytopenia", "Multilineage/Bone marrow"),
    ("Myelosuppression", "Multilineage/Bone marrow"),
    ("Bone marrow failure", "Multilineage/Bone marrow"),
    ("Cytopenia", "Multilineage/Bone marrow"),
    ("Bicytopenia", "Multilineage/Bone marrow"),
    ("Febrile bone marrow aplasia", "Multilineage/Bone marrow"),
    # Coagulation/Microvascular
    ("Disseminated intravascular coagulation", "Coagulation/Microvascular"),
    ("Coagulopathy", "Coagulation/Microvascular"),
    ("Thrombotic microangiopathy", "Coagulation/Microvascular"),
    ("Thrombotic thrombocytopenic purpura", "Coagulation/Microvascular"),
    ("Haemolytic uraemic syndrome", "Coagulation/Microvascular"),
    # Lymphatic/Spleen
    ("Lymphadenopathy", "Lymphatic/Spleen"),
    ("Lymphoedema", "Lymphatic/Spleen"),
    ("Splenomegaly", "Lymphatic/Spleen"),
    ("Lymphocele", "Lymphatic/Spleen"),
]

MECH_ORDER = ["A_HIF2a", "B_VEGFR_TKI", "C_VEGF_mAb", "D_mTOR", "E_ICI"]

KEY_PT_NAMES = [p for p, _ in KEY_PT]
LINEAGE = {p: g for p, g in KEY_PT}

MECH_LABEL = {
    "A_HIF2a": "HIF-2α inhibitor (belzutifan)",
    "B_VEGFR_TKI": "VEGFR-TKI (7 drugs)",
    "C_VEGF_mAb": "VEGF mAb (bevacizumab)",
    "D_mTOR": "mTOR inhibitor (everolimus/temsirolimus)",
    "E_ICI": "Immune checkpoint inhibitor (4 drugs)",
}
