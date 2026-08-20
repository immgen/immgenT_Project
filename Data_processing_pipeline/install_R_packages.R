# install_R_packages.R
# Alternative to environment.yml for people who want to install
# into an existing R (>= 4.1) installation instead of using conda.
#
#   Rscript install_R_packages.R
#
# Seurat is pinned to 4.4.0 on purpose - see the note in environment.yml
# about v4 vs v5 slot syntax.

if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

remotes::install_version("Seurat", version = "4.4.0", upgrade = "never")

cran_pkgs = c(
  "ggplot2", "gridExtra", "dplyr", "RColorBrewer", "reshape2", "fitdistrplus",
  "ggExtra", "scales", "ggrastr", "gatepoints", "rafalib", "inflection",
  "seqinr", "Matrix"
)
new_cran = cran_pkgs[!cran_pkgs %in% installed.packages()[, "Package"]]
if (length(new_cran)) install.packages(new_cran)

bioc_pkgs = c("DropletUtils", "SingleR", "celldex", "UCell", "SingleCellExperiment", "BiocGenerics")
new_bioc = bioc_pkgs[!bioc_pkgs %in% installed.packages()[, "Package"]]
if (length(new_bioc)) BiocManager::install(new_bioc, update = FALSE, ask = FALSE)

message("Done. Remember: Cell Ranger, bcl2fastq, and IMGT/HighV-QUEST are NOT R packages - see README.md.")
