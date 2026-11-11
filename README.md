# PBMC scRNA-seq integration with Harmony

A short, runnable walk-through of integrating several 10x single-cell RNA-seq
samples: quality control, doublet removal, Harmony batch correction,
clustering, marker genes, provisional cell-type labels and composition.

Data: mouse blood (PBMC), five 10x Cell Ranger `filtered_feature_bc_matrix.h5`
files (mm10), three young (9 weeks) and two old (19 months) female SPF mice.
The whole notebook runs in about a minute on a laptop.

## Run it

```bash
git clone <this repo's URL>
cd pbmc-harmony-integration
conda env create -f environment.yml
conda activate pbmc-integration
python -m ipykernel install --user --name pbmc-integration --display-name "Python (pbmc-integration)"
jupyter lab integration.ipynb
```

Then pick the **Python (pbmc-integration)** kernel and run all cells. Or run it
headless:

```bash
jupyter nbconvert --to notebook --execute --inplace integration.ipynb
```

The notebook in the repo is saved **with its outputs**, so you can compare
your run against it cell by cell.

## Layout

| Path | What |
|---|---|
| `integration.ipynb` | the analysis |
| `data/metadata.txt` | one row per sample: tissue, age, sex, mouse id |
| `data/pbmc/<sample>/filtered_feature_bc_matrix.h5` | Cell Ranger output per sample |
| `environment.yml` | conda environment (versions pinned to a tested set) |
| `scripts/make_notebook.py` | builds `integration.ipynb` from source (edit here, then rebuild) |
| `scripts/cluster_copy/` | how the `.h5` files were copied off the cluster |
| `results/` | written by the notebook (not tracked by git) |

## What to know before you start

- **The two old samples (A170, A172) are almost all red blood cells**: about
  90 % of their reads are haemoglobin and a typical cell has under 100 genes.
  QC removes nearly all of them, so only the three young samples are
  integrated and **this data cannot compare old and young blood**.
- Cluster labels are assigned automatically from marker-gene scores. They are
  a starting point; check each one against the marker table and dot plot.
- Doublets are called with a fixed scrublet threshold (0.25). Check the
  histograms in section 3.
