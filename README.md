# Kidney organoid analyses

Score how close lab-grown kidney tissue is to a **working kidney** (filter, reabsorbing tubule, low scar) and whether injured tissue can come back.

Owner: [hokocodes](https://github.com/hokocodes)

## Analyses

1. **Maturity and missing cell types** (`analysis1_maturity_missing_cells.py`)
   - Module scores: progenitor, podocyte, mature PT, loop of Henle, distal, collecting duct, glomerular endothelium, off-target, fibrosis
   - Outputs a missing-cell report vs a working-nephron census and a maturity score (adult programs minus progenitors)
   - Public targets: GEO GSE213152, GSE269904 + [Zenodo 18732338](https://zenodo.org/records/18732338), [KPMP](https://atlas.kpmp.org/)

2. **Fibrosis vs repair** (`analysis2_fibrosis_vs_repair.py`)
   - Recovery Index = (healthy PT + podocyte) − (failed-repair + fibrosis)
   - Simulated 11 mM vs 33 mM glucose vs inhibitor arm (GSE342984-style biology)
   - Public targets: GEO GSE342984, Combes ischaemic organoid AKI, KPMP failed-repair labels

Default mode uses simulated atlas-like cells so the repo runs without multi-GB downloads. Swap the `simulate_*` functions for a real gene-by-cell matrix (columns = HGNC symbols).

## Run

```bash
python3 -m pip install numpy pandas matplotlib
python3 analysis1_maturity_missing_cells.py
python3 analysis2_fibrosis_vs_repair.py
```

Writes CSVs and PNGs to `analysis_outputs/`.

## Honest limit

Marker scores are not filtration or concentrating ability. They tell you which lineages to add next (assembloid, transplant, better patterning) so tissue can work again.
