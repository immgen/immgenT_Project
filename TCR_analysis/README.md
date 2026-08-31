# abTCR Repertoire Analysis in the immgenT Dataset

## Description

This repository contains the analysis scripts and Jupyter notebooks used to analyze T-cell receptor (TCR) repertoire data from the **immgenT dataset**.

The analysis includes:

* Extraction and processing of IMGT/HighV-QUEST outputs: ExtractionData_HighVQUEST.ipynb
* Analysis of TCR Vα-Jα (**VaJa**) gene-pair frequencies and visualization of VaJa frequencies using heatmaps: Heatmap_vaja.ipynb
* Analysis and visualization of TCR mutations: Plot_Mutations.ipynb
* Additional data processing performed using R scripts.

The jupyter notebooks are implemented in **Python 3.12.3**.

---

## Repository Structure

The analysis files are located in the `T_analysis` directory:

```text
immgenT_Project/
└── TCR_analysis/
    ├── ExtractionData_HighVQUEST.ipynb
    ├── Heatmap_VaJa.ipynb
    ├── Plot_Mutations.ipynb
    ├── *.R
    └── requirements.txt
```

The input data required to run the analyses are provided separately through Zenodo:

**Zenodo DOI:** `10.5281/zenodo.21839963`

The Zenodo repository contains the data directory:

```text
immgenT-TCR/
└── data/
    ├── Frequency/
    ├── HighV_QUEST/
    └── Mutations/
```

The data are not stored directly in this GitHub repository.

---

## Requirements

### Python

The notebooks were developed and tested with:

* Python 3.12.3
* Jupyter Notebook or JupyterLab

Python dependencies are listed in `requirements.txt`.

Install them with:

```bash
pip install -r requirements.txt
```

For reproducibility, using Python 3.12.3 is recommended.

### R

Some analysis steps require R. The required R packages and R version should be specified here if applicable.

For example:

```text
R version: X.X.X

Required R packages:
- package1
- package2
- package3
```

If the R scripts are not required to run the notebooks, this should be stated explicitly.

---

## Input Data

The input data required for the analyses are available from Zenodo:

[Zenodo record — DOI 10.5281/zenodo.21839963](https://zenodo.org/records/21839963?utm_source=chatgpt.com)

Download the dataset and place the `data/` directory in the folder coantaining the notebooks.

The data directory contains, among other files:

* `Frequency/` — data used for Vα-Jα frequency analysis and heatmap generation.
* `HighV_QUEST/` — IMGT/HighV-QUEST results.
* `Mutations/` — data used for TCR mutation analysis.

**Important:** The notebooks may contain relative paths to these directories. Before running the analyses, verify that the paths in the notebooks correspond to the location of the downloaded `data/` directory.

---

## How to Run

### 1. Clone the repository

Clone the GitHub repository and move to the analysis directory:

```bash
git clone <repository-url>
cd immgenT_Project/TCR_analysis
```

Replace `<repository-url>` with the URL of this GitHub repository.

### 2. Download the data

Download the dataset from Zenodo:

[Zenodo — 10.5281/zenodo.21839963](https://zenodo.org/records/21839963?utm_source=chatgpt.com)

Place the downloaded `data/` directory in the location expected by the analysis scripts (in the folder with the scripts).

### 3. Make the outputs folder

Create an `outputs/` folder where the outputs figure and tables will be saved.


### 4. Set up the Python environment

It is recommended to create a dedicated virtual environment:

```bash
python3.12 -m venv .venv
```

Activate it:

**Linux/macOS:**

```bash
source .venv/bin/activate
```

**Windows:**

```bash
.venv\Scripts\activate
```

Then install the dependencies:

```bash
pip install -r requirements.txt
```

### 5. Launch Jupyter

From the `TCR_analysis` directory, run:

```bash
jupyter notebook
```

or:

```bash
jupyter lab
```

### 6. Run the notebooks

The main notebooks are:

| Notebook                          | Description                                        |
| --------------------------------- | -------------------------------------------------- |
| `ExtractionData_HighVQUEST.ipynb` | Extraction and processing of IMGT/HighV-QUEST data |
| `Heatmap_VaJa.ipynb`              | Analysis and visualization of Vα-Jα frequencies    |
| `Plot_Mutations.ipynb`            | Analysis and visualization of TCR mutations        |


---

## Outputs

The analysis produces figures, tables, and/or processed datasets used to characterize the immgenT TCR repertoire.

The main outputs include:

* Vα-Jα frequency heatmaps.
* Processed IMGT/HighV-QUEST data.
* TCR mutation analysis plots.
* Processed tables generated during the analysis.

The exact output files will be stored in the `outputs/` folder and their locations are described within the corresponding notebooks.

---

## Reproducibility

To reproduce the analyses:

1. Clone this repository.
2. Download the corresponding dataset from Zenodo.
3. Use Python 3.12.3 and install the dependencies listed in `requirements.txt`.
4. Install the required R version and packages, if applicable.
5. Place the Zenodo `data/` directory in the expected location.
6. Run the notebooks

The GitHub repository contains the analysis code, while the Zenodo record contains the input data required to reproduce the analyses.

---


