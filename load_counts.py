from __future__ import annotations
from pathlib import Path
import numpy as np
import pandas as pd

def _decode(xs):
    return [x.decode() if isinstance(x, (bytes, np.bytes_)) else str(x) for x in xs]

def _norm_genes(names):
    return [str(n).split("\t")[-1] if str(n).startswith("ENSG") and "\t" in str(n) else str(n) for n in names]

def _subset(mat, bcs, genes, keep_genes=None):
    genes = list(genes)
    if keep_genes is not None:
        want = set(keep_genes)
        idx = [i for i, g in enumerate(genes) if g in want]
        mat = mat[:, idx]
        genes = [genes[i] for i in idx]
    expr = pd.DataFrame(mat.toarray(), index=bcs, columns=genes)
    if expr.columns.duplicated().any():
        expr = expr.T.groupby(level=0).sum().T
    return expr

def load_10x_h5(path, keep_genes=None):
    import h5py
    from scipy import sparse
    path = Path(path)
    with h5py.File(path, "r") as f:
        grp = f["matrix"] if "matrix" in f else f[list(f.keys())[0]]
        shape = tuple(int(x) for x in grp["shape"][:])
        mat = sparse.csc_matrix((grp["data"][:], grp["indices"][:], grp["indptr"][:]), shape=shape)
        feats = grp["features"] if "features" in grp else grp["genes"]
        raw = feats["name"][:] if "name" in feats else feats["id"][:]
        genes = _norm_genes(_decode(raw))
        bcs = _decode(grp["barcodes"][:]) if "barcodes" in grp else [f"c{i}" for i in range(shape[1])]
    return _subset(mat.T.tocsr(), bcs, genes, keep_genes)

def load_10x_mtx(folder, keep_genes=None):
    from scipy.io import mmread
    folder = Path(folder)
    mtx = next(folder.glob("*matrix.mtx*"))
    genes_f = next((p for p in list(folder.glob("*features.tsv*")) + list(folder.glob("*genes.tsv*"))), None)
    bc_f = next(folder.glob("*barcodes.tsv*"), None)
    mat = mmread(mtx).tocsr()
    if genes_f:
        gdf = pd.read_csv(genes_f, sep="\t", header=None)
        genes = _norm_genes((gdf.iloc[:, 1] if gdf.shape[1] > 1 else gdf.iloc[:, 0]).astype(str).tolist())
    else:
        genes = [f"g{i}" for i in range(mat.shape[0])]
    bcs = pd.read_csv(bc_f, header=None)[0].astype(str).tolist() if bc_f else [f"c{i}" for i in range(mat.shape[1])]
    if mat.shape[0] == len(genes):
        mat = mat.T.tocsr()
    return _subset(mat, bcs, genes, keep_genes)

def load_h5ad(path, keep_genes=None):
    import anndata as ad
    a = ad.read_h5ad(path, backed="r")
    if keep_genes is not None:
        have = [g for g in keep_genes if g in a.var_names]
        X = a[:, have].X
        genes = have
    else:
        X = a.X
        genes = a.var_names.astype(str).tolist()
    if hasattr(X, "toarray"):
        X = X.toarray()
    expr = pd.DataFrame(np.asarray(X), index=a.obs_names.astype(str), columns=genes)
    obs = a.obs.copy()
    a.file.close()
    return expr, obs

def all_module_genes(*module_dicts):
    seen, out = set(), []
    for d in module_dicts:
        for genes in d.values():
            for g in genes:
                if g not in seen:
                    seen.add(g); out.append(g)
    return out

def load_any(path, keep_genes=None):
    path = Path(path)
    if path.is_dir():
        return load_10x_mtx(path, keep_genes=keep_genes)
    suf = "".join(path.suffixes).lower()
    if suf.endswith(".h5ad"):
        expr, _ = load_h5ad(path, keep_genes=keep_genes); return expr
    if suf.endswith(".h5"):
        return load_10x_h5(path, keep_genes=keep_genes)
    if suf.endswith(".csv") or suf.endswith(".csv.gz"):
        df = pd.read_csv(path, index_col=0)
        if keep_genes is not None:
            df = df[[c for c in keep_genes if c in df.columns]]
        return df
    raise ValueError(path)

def qc_filter(expr, min_counts=1.0):
    return expr.loc[expr.sum(axis=1) >= min_counts]

def log1p_norm(expr):
    lib = expr.sum(axis=1).replace(0, np.nan)
    return np.log1p(expr.div(lib, axis=0).fillna(0) * 1e4)
