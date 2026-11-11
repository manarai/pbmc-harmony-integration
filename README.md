# Mouse scRNA-seq integration with Harmony

Runnable walk-throughs of integrating 10x single-cell RNA-seq samples:
quality control, doublet removal, Harmony batch correction, clustering,
marker genes, provisional cell-type labels and composition.

Data: 28 Cell Ranger `filtered_feature_bc_matrix.h5` files (mm10) from female
SPF mice, young (9 weeks) and old (19 months), across five tissues. Several
mice were sampled in more than one tissue (see `data/metadata.txt`).

| Tissue | Young | Old |
|---|---|---|
| Bone marrow | 3 | 3 |
| Large intestine | 3 | 3 |
| Lung | 3 | 3 |
| PBMC | 3 | 2 |
| Spleen | 3 | 2 |

## Two notebooks

| Notebook | Data | Runtime | What it covers |
|---|---|---|---|
| `integration.ipynb` | PBMC only (5 samples) | ~1 min | the basic workflow, start here |
| `integration_all_tissues.ipynb` | all 28 samples | ~5 min | the same workflow on all five tissues, plus a per-tissue check of sample mixing and a pseudobulk old-vs-young test (DESeq2) per tissue × cell type |

## Run it

```bash
git clone https://github.com/manarai/pbmc-harmony-integration.git
cd pbmc-harmony-integration
conda env create -f environment.yml
conda activate pbmc-integration
python -m ipykernel install --user --name pbmc-integration --display-name "Python (pbmc-integration)"
jupyter lab
```

Open a notebook, pick the **Python (pbmc-integration)** kernel and run all
cells. Or run one headless:

```bash
jupyter nbconvert --to notebook --execute --inplace integration_all_tissues.ipynb
```

Both notebooks are saved **with their outputs**, so you can compare your run
against them cell by cell. The all-tissue notebook writes a ~425 MB
`results/all_tissues/integrated.h5ad`.

If you already created the environment before the all-tissue notebook was
added, update it (adds `pydeseq2`):

```bash
conda env update -n pbmc-integration -f environment.yml
```

## Layout

| Path | What |
|---|---|
| `integration.ipynb` | PBMC analysis |
| `integration_all_tissues.ipynb` | all-tissue analysis |
| `data/metadata.txt` | one row per sample: tissue, age, sex, mouse id |
| `data/<tissue>/<sample>/filtered_feature_bc_matrix.h5` | Cell Ranger output per sample |
| `environment.yml` | conda environment (versions pinned to a tested set) |
| `scripts/make_notebook.py`, `scripts/make_notebook_all_tissues.py` | build the notebooks from source (edit there, then rebuild) |
| `scripts/cluster_copy/` | how the `.h5` files were copied off the cluster |
| `results/` | written by the notebooks (not tracked by git) |

## What to know before you start

- **The two old PBMC samples (A170, A172) are almost all red blood cells**:
  about 90 % of their reads are haemoglobin and a typical cell has under 100
  genes. QC removes nearly all of them, so PBMC has young samples only and
  **cannot be used to compare old and young blood**.
- Samples are corrected for `sample`, and each sample is one mouse at one age,
  so the integrated space can absorb part of a real age effect. It is used for
  clustering and labelling only. The age test uses raw counts summed per
  mouse (pseudobulk). With 2–3 mice per group, read the hit lists as
  candidates.
- Large intestine keeps cells with up to 50 % mitochondrial reads (25 % in the
  other tissues), because intestinal epithelium normally runs high.
- Cluster labels are assigned automatically from marker-gene scores. They are
  a starting point. Check each one against the marker table and dot plot
  (for example, in the all-tissue run one bone-marrow/spleen cluster is
  weakly labelled "Epithelial").
- Doublets are called with a fixed scrublet threshold (0.25). Check the
  histograms. In lung and large intestine fewer doublets are called than the
  10x rule of thumb expects.
