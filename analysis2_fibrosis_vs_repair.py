#!/usr/bin/env python3
"""Analysis 2 — Fibrosis vs repair after metabolic / ischaemic stress.

Recovery Index = (healthy PT + podocyte) - (failed-repair + fibrosis)
Default: simulated GSE342984-like 11 mM vs 33 mM glucose vs inhibitor.
Public: GSE342984, Combes ischaemic organoid AKI, KPMP failed-repair labels.
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
RNG = np.random.default_rng(33)

MODULES = {
    "healthy_pt": ["LRP2", "CUBN", "SLC34A1", "HNF4A"],
    "acute_injury": ["HAVCR1", "LCN2", "SOX9"],
    "failed_repair": ["VCAM1", "PROM1", "DCDC2", "IL32"],
    "podocyte_health": ["NPHS1", "NPHS2", "WT1", "SYNPO"],
    "inflammation": ["TNF", "MIF", "CD74", "CXCL1", "NFKBIA"],
    "fibrosis": ["ACTA2", "COL1A1", "COL3A1", "POSTN", "FAP"],
}

def unique_genes():
    seen, out = set(), []
    for gs in MODULES.values():
        for g in gs:
            if g not in seen:
                seen.add(g); out.append(g)
    return out

def simulate():
    genes = unique_genes()
    arms = {
        "glucose_11mM": {"healthy_pt": 0.34, "acute_injury": 0.08, "failed_repair": 0.06,
            "podocyte_health": 0.22, "inflammation": 0.10, "fibrosis": 0.20},
        "glucose_33mM": {"healthy_pt": 0.16, "acute_injury": 0.16, "failed_repair": 0.18,
            "podocyte_health": 0.10, "inflammation": 0.18, "fibrosis": 0.22},
        "glucose_33mM_inhibitor": {"healthy_pt": 0.26, "acute_injury": 0.10, "failed_repair": 0.10,
            "podocyte_health": 0.18, "inflammation": 0.12, "fibrosis": 0.24},
    }
    blocks, meta, i, n = [], [], 0, 500
    for arm, fracs in arms.items():
        types = list(fracs)
        counts = RNG.multinomial(n, [fracs[t] for t in types])
        for t, k in zip(types, counts):
            X = np.clip(RNG.normal(0.25, 0.2, size=(k, len(genes))), 0, None)
            for g in MODULES[t]:
                X[:, genes.index(g)] += RNG.normal(3.0, 0.35, size=k)
            if arm == "glucose_33mM":
                for g in MODULES["inflammation"]:
                    X[:, genes.index(g)] += RNG.normal(0.9, 0.25, size=k)
                for g in MODULES["healthy_pt"] + MODULES["podocyte_health"]:
                    X[:, genes.index(g)] *= 0.75
            if arm == "glucose_33mM_inhibitor":
                for g in MODULES["inflammation"]:
                    X[:, genes.index(g)] += RNG.normal(0.25, 0.2, size=k)
            blocks.append(np.clip(X, 0, None))
            for _ in range(k):
                meta.append({"cell": f"c{i}", "arm": arm, "true_state": t}); i += 1
    return pd.DataFrame(np.vstack(blocks), columns=genes), pd.DataFrame(meta)

def zscore_modules(expr):
    scores = pd.DataFrame(index=expr.index)
    for name, genes in MODULES.items():
        use = [g for g in genes if g in expr.columns]
        sub = expr[use]
        z = (sub - sub.mean()) / sub.std(ddof=0).replace(0, 1)
        scores[name] = z.mean(axis=1)
    return scores

def recovery_index(s):
    return (s["healthy_pt"] + s["podocyte_health"]) - (s["failed_repair"] + s["fibrosis"])

def violins(scores, path):
    arms = ["glucose_11mM", "glucose_33mM", "glucose_33mM_inhibitor"]
    keys = ["healthy_pt", "failed_repair", "fibrosis", "recovery_index"]
    fig, axes = plt.subplots(1, 4, figsize=(11.5, 3.8))
    colors = ["#148F77", "#C0392B", "#B9770E"]
    for ax, key in zip(axes, keys):
        data = [scores.loc[scores.arm == a, key].to_numpy() for a in arms]
        parts = ax.violinplot(data, showmeans=True, showextrema=False)
        for i, body in enumerate(parts["bodies"]):
            body.set_facecolor(colors[i]); body.set_alpha(0.7)
        ax.set_xticks([1, 2, 3])
        ax.set_xticklabels(["11 mM", "33 mM", "33 mM+inh"], rotation=20, ha="right", fontsize=8)
        ax.set_title(key.replace("_", " "), fontsize=10)
        ax.axhline(0, color="#bbb", lw=0.8)
    fig.suptitle("Does the tissue return toward working programs?", fontsize=12)
    fig.tight_layout(); fig.savefig(path, dpi=160); plt.close()

def bars_index(summary, path):
    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    xs = np.arange(len(summary))
    ax.bar(xs, summary["recovery_index"], color=["#148F77", "#C0392B", "#B9770E"], width=0.65)
    ax.axhline(0, color="#333", lw=0.8)
    ax.set_ylabel("Recovery Index")
    ax.set_title("Working programs minus failed-repair and scar")
    ax.set_xticks(xs)
    ax.set_xticklabels(["11 mM\ncontrol", "33 mM\nhigh glucose", "33 mM +\nMIF/TNF inh"], fontsize=8)
    fig.tight_layout(); fig.savefig(path, dpi=160); plt.close()

def main():
    expr, meta = simulate()
    scores = zscore_modules(expr)
    scores["recovery_index"] = recovery_index(scores)
    scores["arm"] = meta["arm"].to_numpy()
    scores["true_state"] = meta["true_state"].to_numpy()
    scores["assigned"] = scores[list(MODULES)].idxmax(axis=1)
    summary = scores.groupby("arm")[
        ["healthy_pt", "acute_injury", "failed_repair", "podocyte_health", "inflammation", "fibrosis", "recovery_index"]
    ].mean().reset_index()
    order = ["glucose_11mM", "glucose_33mM", "glucose_33mM_inhibitor"]
    summary["arm"] = pd.Categorical(summary["arm"], order, ordered=True)
    summary = summary.sort_values("arm")
    fracs = scores.groupby(["arm", "assigned"]).size().unstack(fill_value=0).apply(lambda r: r / r.sum(), axis=1).reindex(order)
    scores.to_csv(OUT / "analysis2_cell_scores.csv", index=False)
    summary.to_csv(OUT / "analysis2_arm_summary.csv", index=False)
    fracs.to_csv(OUT / "analysis2_state_fractions.csv")
    violins(scores, OUT / "analysis2_module_violins.png")
    bars_index(summary, OUT / "analysis2_recovery_index.png")
    print("Analysis 2 written to", OUT)
    print(summary.round(3).to_string(index=False))

if __name__ == "__main__":
    main()
