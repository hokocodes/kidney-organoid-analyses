#!/usr/bin/env python3
"""Score real organoid count matrices. See DATA_SOURCES.md and download_public_data.sh."""
from __future__ import annotations
import argparse
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from kidney_modules import ADULT_TARGET, MATURITY_MODULES, REPAIR_MODULES
from load_counts import all_module_genes, load_any, load_h5ad, log1p_norm, qc_filter

OUT = Path(__file__).resolve().parent / "analysis_outputs"
OUT.mkdir(parents=True, exist_ok=True)

def module_scores(expr, modules):
    present = set(expr.columns)
    scores = pd.DataFrame(index=expr.index)
    used = {}
    for name, genes in modules.items():
        have = [g for g in genes if g in present]
        used[name] = have
        if not have:
            scores[name] = 0.0
            continue
        sub = expr[have]
        z = (sub - sub.mean()) / sub.std(ddof=0).replace(0, 1)
        scores[name] = z.mean(axis=1)
    return scores, used

def load_labeled(pairs, keep_genes):
    blocks, labels = [], []
    for path, label in pairs:
        p = Path(path)
        if str(p).endswith(".h5ad"):
            expr, obs = load_h5ad(p, keep_genes=keep_genes)
            expr = log1p_norm(qc_filter(expr))
            blocks.append(expr)
            labels.append(obs.loc[expr.index, label].astype(str) if label in obs.columns else pd.Series(label, index=expr.index))
            continue
        expr = log1p_norm(qc_filter(load_any(p, keep_genes=keep_genes)))
        blocks.append(expr)
        labels.append(pd.Series(label, index=expr.index))
    expr = pd.concat(blocks, axis=0)
    lab = pd.concat(labels).reindex(expr.index)
    return expr, lab

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--mode", choices=["maturity", "repair", "both"], default="both")
    p.add_argument("--matrix", action="append", default=[])
    p.add_argument("--label", action="append", default=[])
    args = p.parse_args()
    if not args.matrix:
        raise SystemExit("Pass --matrix paths. Run: bash download_public_data.sh")
    labels = args.label if args.label else [Path(m).stem for m in args.matrix]
    expr, lab = load_labeled(list(zip(args.matrix, labels)), all_module_genes(MATURITY_MODULES, REPAIR_MODULES))
    print(f"Loaded {expr.shape[0]} cells x {expr.shape[1]} genes")
    if args.mode in ("maturity", "both"):
        ms, used = module_scores(expr, MATURITY_MODULES)
        print("maturity genes", {k: len(v) for k, v in used.items()})
        ms = ms.copy(); ms["group"] = lab.to_numpy()
        ms["assigned"] = ms[list(MATURITY_MODULES)].idxmax(axis=1)
        fracs = ms.groupby(["group", "assigned"]).size().unstack(fill_value=0).apply(lambda r: r / r.sum(), axis=1)
        fracs.to_csv(OUT / "real_maturity_fractions.csv")
        print(fracs.round(3))
    if args.mode in ("repair", "both"):
        rs, used = module_scores(expr, REPAIR_MODULES)
        print("repair genes", {k: len(v) for k, v in used.items()})
        rs = rs.copy(); rs["group"] = lab.to_numpy()
        rs["recovery_index"] = (rs["healthy_pt"] + rs["podocyte_health"]) - (rs["failed_repair"] + rs["fibrosis"])
        summary = rs.groupby("group")[list(REPAIR_MODULES) + ["recovery_index"]].mean()
        summary.to_csv(OUT / "real_repair_summary.csv")
        print(summary.round(3))

if __name__ == "__main__":
    main()
