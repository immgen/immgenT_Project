# immgenT Treg

This repository contains the analysis and reproducibility code accompanying the **immgenT Treg** manuscript. The rendered workflowr analysis can be browsed here: **https://immgen.github.io/immgenT_Project/Treg_paper/** This website provides a figure-by-figure view of the analysis, with the corresponding R code and workflowr reproducibility information.

## Data availability

The complete data package is available on Zenodo:

**https://zenodo.org/records/21839963**

The Zenodo archive contains the curated Treg Seurat object and supporting input files required to reproduce the analyses in this directory. Large single-cell data objects are distributed through Zenodo rather than stored directly on GitHub.

## Repository structure

The Treg analysis is contained within the `Treg_paper/` directory of the immgenT project repository:

```text
immgenT_Project/
└── Treg_paper/
    ├── analysis/
    │   ├── Treg_Workflow.Rmd       # Main Treg analysis workflow
    │   ├── index.Rmd               # workflowr website homepage
    │   ├── about.Rmd
    │   └── license.Rmd
    ├── code/                       # Supporting R scripts and functions
    ├── data/                       # Input data downloaded from Zenodo
    ├── docs/                       # Rendered workflowr website
    ├── output/                     # Analysis outputs
    ├── _workflowr.yml              # workflowr configuration
    ├── treg_github.Rproj           # RStudio project
    ├── README.md
    └── .gitignore
```

The main analysis is contained in `analysis/Treg_Workflow.Rmd`. Supporting data files are stored in `data/`, and analysis outputs are written to `output/` or workflowr-generated figure directories under `docs/`.

## Environment setup

The analysis is written in R and uses Seurat for single-cell analysis together with workflowr for organization and reproducible rendering.

Clone or download the `immgenT_Project` repository and navigate to the `Treg_paper/` directory. The Treg analysis can then be opened using the included RStudio project:

```text
treg_github.Rproj
```

Install `workflowr` if it is not already available:

```r
install.packages("workflowr")
```

The analysis additionally uses packages including `Seurat`, `tidyverse`, `dplyr`, `ggplot2`, `ggrastr`, `pheatmap`, `RColorBrewer`, `scales`, `cowplot`, and `ggrepel`, together with additional packages specified within the workflow.

## Data setup

Download the Treg data package from Zenodo and place the required input files in the `data/` directory.

After downloading the data, the relevant portion of the project should have the general structure:

```text
Treg_paper/
├── data/
│   ├── [Treg Seurat object]
│   └── [supporting analysis files]
├── analysis/
├── code/
├── output/
└── ...
```
Large intermediate objects and other files that are impractical to distribute through GitHub are provided through the associated Zenodo archive.

## Running the analysis

From the `Treg_paper/` project root:

1. Download the associated data package from Zenodo and place the required files in `data/`.
2. Open `treg_github.Rproj` in RStudio.
3. Install the required R packages if necessary.
4. Build the workflow:

```r
library(workflowr)

wflow_build("analysis/Treg_Workflow.Rmd")
```
The rendered analysis is written to:

```text
docs/Treg_Workflow.html
```
The R Markdown workflow can also be run interactively in RStudio. Code chunks should generally be run in order because later analyses may depend on objects generated in earlier sections.

## Data and analysis notes

- The curated Treg Seurat object contains the single-cell RNA and CITE-seq information used throughout the analysis.
- The primary Treg embedding used for visualization is `mde_incremental`.
- Treg populations are annotated using the `annotation_level2` metadata field, including the major Treg states `Treg.A`–`Treg.F` and proliferative Tregs.
- Supporting files distributed with the Zenodo data package provide precomputed results or metadata for analyses that do not require recomputation of the full upstream single-cell pipeline.
- Large intermediate objects and computationally intensive upstream processing steps are not regenerated when a curated or precomputed input is sufficient to reproduce the corresponding figure panel.

## Reproducibility

This analysis uses [workflowr](https://github.com/workflowr/workflowr) to organize the R Markdown source, rendered results, and version information associated with the analysis.

To rebuild the Treg workflow after making changes:

```r
library(workflowr)

wflow_build("analysis/Treg_Workflow.Rmd")
```
The generated workflowr website and figure outputs are stored under `docs/`.

The published analysis can be viewed at:

**https://immgen.github.io/immgenT_Project/Treg_paper/**

## Citation

If you use this code or data, please cite the associated immgenT Treg manuscript: Freuchet et al. **A Reference Landscape of Regulatory T Cell States in Mice**. https://www.biorxiv.org/content/10.64898/2026.01.30.702856

