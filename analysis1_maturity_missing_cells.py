#!/usr/bin/env python3
"""Analysis 1 — Kidney organoid maturity and missing cell types.

Score how close a batch is to a working nephron versus progenitors,
off-target cells, and scar stroma. Default: simulated public-atlas-like
organoid (days 7, 14, 26). Swap simulate_organoid() for GEO/Zenodo matrices.
Public: GSE213152, GSE269904 + Zenodo 18732338, KPMP atlas.kpmp.org
"""
from __future__ import annotations
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

OUT = Path(__file__).resolve().parent / "analysis_outputs"
OUT.mkdir(parents=True, exist_ok=True)
RNG = np.random.default_rng(26)

MODULES = {
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
ADULT_TARGET = {
    "progenitor": 0.00, "podocyte": 0.08, "pt_mature": 0.40, "loh": 0.18,
    "distal": 0.10, "collecting_duct": 0.12, "glom_endo": 0.07,
    "offtarget": 0.00, "fibrosis": 0.05,
}

def unique_genes():
    seen, out = set(), []
    for genes in MODULES.values():
        for g in genes:
            if g not in seen:
                seen.add(g); out.append(g)
    return out

def simulate_organoid(n_per_day=400):
    genes = unique_genes()
    days = [7, 14, 26]
    type_by_day = {
        7: {"progenitor": 0.55, "podocyte": 0.08, "pt_mature": 0.07, "loh": 0.02,
            "distal": 0.02, "collecting_duct": 0.01, "glom_endo": 0.02,
            "offtarget": 0.13, "fibrosis": 0.10},
        14: {"progenitor": 0.22, "podocyte": 0.16, "pt_mature": 0.22, "loh": 0.06,
             "distal": 0.06, "collecting_duct": 0.04, "glom_endo": 0.03,
             "offtarget": 0.09, "fibrosis": 0.12},
        26: {"progenitor": 0.08, "podocyte": 0.20, "pt_mature": 0.28, "loh": 0.08,
             "distal": 0.08, "collecting_duct": 0.05, "glom_endo": 0.04,
             "offtarget": 0.06, "fibrosis": 0.13},
    }
    rows, meta, cell_i = [], [], 0
    for day in days:
        fracs = type_by_day[day]
        types = list(fracs)
        counts = RNG.multinomial(n_per_day, [fracs[t] for t in types])
        for t, n in zip(types, counts):
            block = np.clip(RNG.normal(0.3, 0.25, size=(n, len(genes))), 0, None)
            for g in MODULES[t]:
                block[:, genes.index(g)] += RNG.normal(3.2 if day >= 14 or t == "progenitor" else 2.4, 0.4, size=n)
            if day == 7 and t != "progenitor":
                for g in MODULES["progenitor"]:
                    block[:, genes.index(g)] += RNG.normal(0.8, 0.3, size=n)
            rows.append(np.clip(block, 0, None))
            for _ in range(n):
                meta.append({"cell": f"c{cell_i}", "day": day, "true_type": t})
                cell_i += 1
    return pd.DataFrame(np.vstack(rows), columns=genes), pd.DataFrame(meta)

def zscore_modules(expr):
    scores = pd.DataFrame(index=expr.index)
    for name, genes in MODULES.items():
        use = [g for g in genes if g in expr.columns]
        sub = expr[use]
        z = (sub - sub.mean()) / sub.std(ddof=0).replace(0, 1)
        scores[name] = z.mean(axis=1)
    return scores

def maturity_score(scores):
    adult = scores[["podocyte", "pt_mature", "loh", "distal", "collecting_duct", "glom_endo"]].mean(axis=1)
    return adult - scores["progenitor"]

def radar(fracs, path):
    labels = list(ADULT_TARGET)
    adult = np.array([ADULT_TARGET[k] for k in labels], float); adult = adult / adult.sum()
    angles = np.linspace(0, 2 * np.pi, len(labels), endpoint=False)
    angles = np.concatenate([angles, angles[:1]])
    fig, ax = plt.subplots(figsize=(7.2, 7.2), subplot_kw=dict(polar=True))
    ax.plot(angles, np.concatenate([adult, adult[:1]]), color="#1B4F72", lw=2, label="working-kidney target")
    colors = {7: "#AF7AC5", 14: "#148F77", 26: "#B9770E"}
    for day, row in fracs.iterrows():
        vals = np.array([row[k] for k in labels], float)
        vals = vals / vals.sum() if vals.sum() else vals
        ax.plot(angles, np.concatenate([vals, vals[:1]]), color=colors.get(int(day), "gray"), lw=2, label=f"organoid d{int(day)}")
        ax.fill(angles, np.concatenate([vals, vals[:1]]), color=colors.get(int(day), "gray"), alpha=0.08)
    ax.set_xticks(angles[:-1]); ax.set_xticklabels(labels, size=8)
    ax.set_title("Distance to a working nephron census", pad=16)
    ax.legend(loc="upper right", bbox_to_anchor=(1.35, 1.1), fontsize=8)
    fig.tight_layout(); fig.savefig(path, dpi=160); plt.close()

def stacked_bars(fracs, path):
    fig, ax = plt.subplots(figsize=(8, 4.6))
    bottom = np.zeros(len(fracs))
    palette = {"progenitor": "#9B59B6", "podocyte": "#1ABC9C", "pt_mature": "#2980B9",
               "loh": "#E67E22", "distal": "#27AE60", "collecting_duct": "#8E44AD",
               "glom_endo": "#E74C3C", "offtarget": "#7F8C8D", "fibrosis": "#C0392B"}
    x = np.arange(len(fracs))
    for col in fracs.columns:
        ax.bar(x, fracs[col], bottom=bottom, color=palette.get(col, "#333"), label=col, width=0.7)
        bottom += fracs[col].to_numpy()
    ax.set_xticks(x); ax.set_xticklabels([f"day {i}" for i in fracs.index])
    ax.set_ylabel("Cell fraction"); ax.set_title("What is present vs what a working kidney needs")
    ax.legend(fontsize=7, ncol=3, frameon=False); ax.set_ylim(0, 1)
    fig.tight_layout(); fig.savefig(path, dpi=160); plt.close()

def main():
    expr, meta = simulate_organoid()
    scores = zscore_modules(expr)
    scores["maturity"] = maturity_score(scores)
    scores["assigned"] = scores.drop(columns=["maturity"]).idxmax(axis=1)
    scores["day"] = meta["day"].to_numpy()
    scores["true_type"] = meta["true_type"].to_numpy()
    fracs = scores.groupby(["day", "assigned"]).size().unstack(fill_value=0).apply(lambda r: r / r.sum(), axis=1)
    for col in MODULES:
        if col not in fracs.columns:
            fracs[col] = 0.0
    fracs = fracs[list(MODULES)]
    gaps = []
    for day, row in fracs.iterrows():
        for lineage, target in ADULT_TARGET.items():
            have = float(row[lineage])
            gaps.append({
                "day": int(day), "lineage": lineage,
                "organoid_fraction": round(have, 4), "working_kidney_target": target,
                "gap_target_minus_have": round(target - have, 4),
                "status": ("excess_unwanted" if lineage in ("progenitor", "offtarget", "fibrosis") and have > target + 0.03
                            else "missing_for_function" if lineage not in ("progenitor", "offtarget", "fibrosis") and have < target - 0.05
                            else "closer"),
            })
    gap_df = pd.DataFrame(gaps)
    mat = scores.groupby("day")["maturity"].mean().rename("mean_maturity_score")
    scores.to_csv(OUT / "analysis1_cell_scores.csv", index=False)
    fracs.to_csv(OUT / "analysis1_lineage_fractions.csv")
    gap_df.to_csv(OUT / "analysis1_missing_cell_report.csv", index=False)
    mat.to_csv(OUT / "analysis1_maturity_by_day.csv")
    radar(fracs, OUT / "analysis1_working_kidney_radar.png")
    stacked_bars(fracs, OUT / "analysis1_lineage_stacked.png")
    print("Analysis 1 written to", OUT)
    print(fracs.round(3).to_string())
    print(mat.round(3).to_string())

if __name__ == "__main__":
    main()
