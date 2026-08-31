# ==============================================================================
# run_kappa_heatmap.R
#
# Standalone script: plots the kappa joint-tail heatmap for all pairwise
# TH program comparisons.
#
#   kappa(x, y) = P(X > x, Y > y)          <- numerator:   empirical joint tail
#                 ──────────────────
#                 P(X > x) * P(Y > y)       <- denominator: independence baseline
#
#   kappa < 1  =>  joint tail DEPLETION   (programs mutually exclusive at extremes)
#   kappa = 1  =>  independence
#   kappa > 1  =>  joint tail ENRICHMENT  (programs co-occur at extremes)
#
# Output: one PDF per pair in figure/kappa_heatmap/
# ==============================================================================

source("tail_dependence_functions.R")
library(tidyverse)


# ==============================================================================
# SETTINGS
# ==============================================================================

DATA_PATH <- "data/metadata_with_modulescore_annotation.csv"
FIG_DIR   <- "figure/kappa_heatmap"          # output folder
N_GRID    <- 80                               # grid resolution (n x n)
V_FIXED   <- 0.80                            # reference percentile shown on both axes

PROGRAMS  <- c("TH1_Score1", "TH2_Score1", "TH17_Score1", "TFH_TCHRON_Score1")


# ==============================================================================
# SETUP
# ==============================================================================

dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
metadata  <- read.csv(DATA_PATH)
short_label <- function(nm) sub("_Score1$", "", nm)

pretty_module_label <- function(score_name) {
  lab <- short_label(score_name)
  if (toupper(lab) == "TH17") return("Th17-like module score")
  if (toupper(lab) == "TH1")  return("Th1-like module score")
  if (toupper(lab) == "TH2")  return("Th2-like module score")
  if (toupper(lab) == "TFH_TCHRON") return("Tfh/Tchron-like module score")
  sprintf("%s-like module score", lab)
}


# ==============================================================================
# PLOTTING FUNCTION
# ==============================================================================

plot_kappa_heatmap <- function(x, y,
                                x_label,
                                y_label,
                                n_grid  = 80,
                                v_fixed = 0.80,
                                zlim_kappa = NULL) {
  plot_kappa_heatmap_square(
    x          = x,
    y          = y,
    x_label    = x_label,
    y_label    = y_label,
    n_grid     = n_grid,
    v_fixed    = v_fixed,
    zlim_kappa = zlim_kappa,
    main       = ""
  )
}


# ==============================================================================
# MAIN LOOP — all unique pairwise comparisons
# ==============================================================================

pairs <- combn(PROGRAMS, 2, simplify = FALSE)
cat(sprintf("Generating kappa heatmaps for %d pairs  ->  %s\n\n",
            length(pairs), FIG_DIR))

for (pair in pairs) {
  x_name <- pair[1]; y_name <- pair[2]
  x_lab  <- short_label(x_name); y_lab  <- short_label(y_name)
  x_axis <- pretty_module_label(x_name)
  y_axis <- pretty_module_label(y_name)
  tag    <- paste0(x_lab, "_vs_", y_lab)

  cat(sprintf("  %s  vs  %s ...\n", x_lab, y_lab))

  pdf_path <- file.path(FIG_DIR, paste0(tag, "_kappa_heatmap.pdf"))
  pdf(pdf_path, width = 7, height = 7)

  plot_kappa_heatmap(
    x       = metadata[[x_name]],
    y       = metadata[[y_name]],
    x_label = x_axis,
    y_label = y_axis,
    n_grid  = N_GRID,
    v_fixed = V_FIXED
  )

  dev.off()
  cat(sprintf("    saved: %s\n", basename(pdf_path)))
}

cat(sprintf("\nDone. All %d figures saved to:  %s/\n", length(pairs), FIG_DIR))
