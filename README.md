# Kidney organoid analyses

Score how close lab-grown kidney tissue is to a **working kidney** and whether injured tissue can come back.

## Real data (what you actually need)

The demo scripts invent cells. Real scoring uses `analyze_real.py` on count matrices you download. Full URL table: [DATA_SOURCES.md](DATA_SOURCES.md).

Start here — four filtered 10x H5 files from **GSE342984** (~115 MB total):

```bash
pip install numpy pandas matplotlib h5py scipy
bash download_public_data.sh
python3 analyze_real.py --mode both \
  --matrix data/GSM9944363_Sample1_filtered_feature_bc_matrix.h5 --label glucose_11mM \
  --matrix data/GSM9944364_Sample2_filtered_feature_bc_matrix.h5 --label glucose_11mM \
  --matrix data/GSM9944365_Sample3_filtered_feature_bc_matrix.h5 --label glucose_33mM \
  --matrix data/GSM9944366_Sample4_filtered_feature_bc_matrix.h5 --label glucose_33mM
```

| Sample | Condition |
|---|---|
| GSM9944363 Sample1 | 11 mM control |
| GSM9944364 Sample2 | 11 mM control |
| GSM9944365 Sample3 | 33 mM high glucose |
| GSM9944366 Sample4 | 33 mM high glucose |

Larger maturity atlas: [Zenodo 18732338](https://zenodo.org/records/18732338) Stage4 `h5ad` (1.3 GB) or GEO [GSE269904](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE269904). Adult / spatial reference: [KPMP](https://atlas.kpmp.org/).

## Demo scripts (no download)

```bash
python3 analysis1_maturity_missing_cells.py
python3 analysis2_fibrosis_vs_repair.py
```
