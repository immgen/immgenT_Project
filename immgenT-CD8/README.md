# immgenT-CD8 analysis companion

This repository is the complete, ready-to-run companion code and data package for the **immgenT-CD8** CITE-seq analysis (CD8ab T cell atlas). It includes the curated MuData object, supporting matrices and gene lists, and three Jupyter notebooks that reproduce the main figures and Extended Data figures.

## Data availability

The complete data package is available here on GitHub and on Zenodo:  
**https://zenodo.org/records/21839963/files/immgenT-CD8.zip**

## Contents

```
immgenT-CD8/
├── data/                     # All required input files
│   ├── immgenT-CD8.h5mu      # Curated MuData (RNA + ADT)
│   ├── immgenT-CD8.xlsx      # Sample-level metadata / proportions
│   ├── cell_factor_matrix.txt.gz
│   ├── gene_factor_matrix.txt
│   ├── ttlist_OneVsAll.txt
│   ├── CD5hi_Fulton.txt
│   ├── CD5lo_Fulton.txt
│   ├── TRM_Crowl.txt
│   ├── TRM_Mackay.txt
│   └── TRM_Milner.txt
├── immgenT-CD8_1.ipynb       # Analysis notebook (part 1)
├── immgenT-CD8_2.ipynb       # Analysis notebook (part 2)
├── immgenT-CD8_3.ipynb       # Analysis notebook (part 3)
├── environment.yml           # Conda environment specification
├── README.md
└── .gitignore
```

All required input files live in the `data/` folder. The notebook uses relative paths of the form `data/...` and should be run from this project root.

## Environment setup

Create and activate the conda environment from the provided specification:

```bash
conda env create -f environment.yml
conda activate immgenT-CD8
```

If you prefer an explicit install path with your own conda installation:

```bash
conda env create -f environment.yml -n immgenT-CD8
conda activate immgenT-CD8
```

Key packages include `muon`, `scanpy`, `gseapy`, `pandas`, `numpy`, `scipy`, `matplotlib`, `seaborn`, and `adjustText`, plus supporting libraries such as `h5py` and `openpyxl` (see `environment.yml` for the full pinned list).

## How to run

1. Activate the environment (`conda activate immgenT-CD8`).
2. From this directory (project root), open the notebooks with any Jupyter-compatible frontend that uses the active conda environment (for example VS Code / Cursor, or JupyterLab if installed separately).
3. Select the `immgenT-CD8` kernel and run the notebooks in order: `immgenT-CD8_1.ipynb` → `immgenT-CD8_2.ipynb` → `immgenT-CD8_3.ipynb`. Run cells top to bottom within each notebook. The notebooks themselves contain detailed instructions, figure annotations, and analysis notes.

To register the kernel for Jupyter frontends:

```bash
python -m ipykernel install --user --name immgenT-CD8 --display-name "Python (immgenT-CD8)"
```

The environment provides `ipykernel` / `ipython` (see `environment.yml`). Install `jupyterlab` or `nbconvert` into the env if you want classic `jupyter lab` / `jupyter nbconvert` workflows.

**Note:** Loading `data/immgenT-CD8.h5mu` and the large factor matrices is memory-intensive. Plan for tens of GB of RAM for a full run.

## Data notes

- `immgenT-CD8.h5mu` holds sparse raw counts, metadata, and the main embedding (`MDE_INCREMENTAL`).
- Gene-program / factor matrices and DEG tables support figure panels that do not recompute the full NMF/DEG pipeline.
- Signature lists (`CD5hi_Fulton.txt`, `CD5lo_Fulton.txt`, `TRM_*`) are used for module scoring and GSEA-style analyses.

## Citation

If you use this code or data, please cite the immgenT-CD8 manuscript / resource associated with this analysis.
