MATURITY_MODULES = {
    "progenitor": ["SIX2", "CITED1", "EYA1", "PAX2", "LHX1"],
    "podocyte": ["NPHS1", "NPHS2", "PODXL", "WT1", "SYNPO"],
    "pt_mature": ["LRP2", "CUBN", "SLC34A1", "SLC5A2", "HNF4A"],
    "loh": ["UMOD", "SLC12A1", "CLDN16"],
    "distal": ["SLC12A3", "PVALB", "TRPM6"],
    "collecting_duct": ["AQP2", "GATA3", "FOXI1", "ATP6V1B1"],
    "glom_endo": ["EHD3", "PLVAP", "PECAM1", "KDR"],
    "offtarget": ["MAP2", "STMN2", "MYOG", "TTN"],
    "fibrosis": ["COL1A1", "COL3A1", "POSTN", "ACTA2", "FAP"],
}
REPAIR_MODULES = {
    "healthy_pt": ["LRP2", "CUBN", "SLC34A1", "HNF4A"],
    "acute_injury": ["HAVCR1", "LCN2", "SOX9"],
    "failed_repair": ["VCAM1", "PROM1", "DCDC2", "IL32"],
    "podocyte_health": ["NPHS1", "NPHS2", "WT1", "SYNPO"],
    "inflammation": ["TNF", "MIF", "CD74", "CXCL1", "NFKBIA"],
    "fibrosis": ["ACTA2", "COL1A1", "COL3A1", "POSTN", "FAP"],
}
ADULT_TARGET = {
    "progenitor": 0.00, "podocyte": 0.08, "pt_mature": 0.40, "loh": 0.18,
    "distal": 0.10, "collecting_duct": 0.12, "glom_endo": 0.07,
    "offtarget": 0.00, "fibrosis": 0.05,
}
