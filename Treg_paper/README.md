# treg_github

A [workflowr][] project.

[workflowr]: https://github.com/workflowr/workflowr

# immgenT-Treg analysis companion

This repository is the complete, ready-to-run companion code for the **immgenT-Treg** CITE-seq analysis (regulatory T cell atlas). It contains the R/workflowr analysis used to reproduce the main and supplementary figure panels, together with supporting scripts and links to the required input data.

## Data availability

The complete data package is available on Zenodo:

**https://zenodo.org/records/21839963**

The Zenodo archive contains the curated Treg Seurat object and supporting input files required to reproduce the analyses in this repository.

## Contents

```text
immgenT-Treg/
├── analysis/
│   ├── Treg_Workflow.Rmd       # Main Treg analysis workflow
│   ├── index.Rmd               # workflowr project homepage
│   ├── about.Rmd
│   └── license.Rmd
├── code/                       # Supporting R scripts and functions
├── data/                       # Input data downloaded from Zenodo
├── docs/                       # Rendered workflowr site
├── output/                     # Analysis outputs
├── _workflowr.yml              # workflowr configuration
├── immgenT-Treg.Rproj          # RStudio project
├── README.md
└── .gitignore
```

The main analysis is contained in `analysis/Treg_Workflow.Rmd`. Supporting data files are stored in the `data/` directory and analysis outputs are written to `output/` or the workflowr-generated figure directories.

## Environment setup

The analysis is written in R and uses Seurat for single-cell analysis together with workflowr for organization and reproducible rendering.

Clone or download this repository and open the RStudio project:

```text
immgenT-Treg.Rproj
```

Install `workflowr` if it is not already available:

```r
install.packages("workflowr")
```

The analysis additionally uses packages including `Seurat`, `tidyverse`, `dplyr`, `ggplot2`, `ggrastr`, `pheatmap`, `RColorBrewer`, `scales`, `cowplot`, `ggrepel`, and other packages specified within the workflow.

## Data setup

Download the Treg data package from Zenodo and place the required files in the `data/` directory.

The repository is designed so that large single-cell data objects are distributed through Zenodo rather than stored directly in GitHub.

After downloading the data, the project should have the general structure:

```text
immgenT-Treg/
├── data/
│   ├── [Treg Seurat object]
│   └── [supporting analysis files]
├── analysis/
├── code/
└── ...
```

## How to run

1. Download or clone this repository.
2. Download the associated data package from Zenodo and place the required files in `data/`.
3. Open `immgenT-Treg.Rproj` in RStudio.
4. Install the required R packages if necessary.
5. Build the workflow from the project root:

```r
library(workflowr)

wflow_build("analysis/Treg_Workflow.Rmd")
```

The rendered analysis will be written to:

```text
docs/Treg_Workflow.html
```

The R Markdown workflow can also be run interactively in RStudio. Code chunks should be run in order because later analyses may depend on objects generated in earlier sections.

## Data notes

- The curated Treg Seurat object contains the single-cell RNA and CITE-seq information used throughout the analysis.
- The primary Treg embedding used for visualization is `mde_incremental`.
- Treg populations are annotated using the `annotation_level2` metadata field, including the major Treg states `Treg.A`-`Treg.F` and proliferative Tregs.
- Supporting files distributed with the Zenodo data package provide precomputed results or metadata for analyses that do not require recomputation of the full upstream single-cell pipeline.
- Large intermediate objects and computationally intensive upstream processing steps are not regenerated where a curated or precomputed input is sufficient to reproduce the corresponding figure panel.

## Reproducibility

This repository uses workflowr to organize the analysis and track the relationship between source code and rendered results.

To rebuild the complete Treg workflow after making changes:

```r
library(workflowr)

wflow_build("analysis/Treg_Workflow.Rmd")
```

The generated workflowr site and figure outputs are stored under `docs/`.

## Citation

If you use this code or data, please cite the immgenT-Treg manuscript / resource associated with this analysis.
