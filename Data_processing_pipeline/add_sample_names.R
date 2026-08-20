# add_sample_names.R
# Add a human-readable sample_name column (from samples.csv) to the final
# dataset.Rds / dataset_clean.Rds objects, based on HTO_classification.simplified.
#
# This replaces naming_samples2.R. (naming_samples_RNA.R is retired entirely -
# it used to double as a stand-in for the skipped SingleR step; the wrapper now
# handles that with a plain file copy, so this script only has one job.)
#
# Works the same whether HTO=YES (real hashtags in samples.csv) or HTO=NO
# (a single row mapping "no_hashing" to your sample name) - see README.
#
# Usage:
#   Rscript add_sample_names.R <dataset.Rds> <dataset_clean.Rds> <path_to_sample_names> <EXP.txt>

args = commandArgs(TRUE)
dataset_path          = args[1]
dataset_clean_path    = args[2]
path_to_sample_names  = args[3]
EXP                   = args[4]

libs = c("Seurat", "ggplot2", "inflection", "grid", "gridExtra", "UCell")
sapply(libs, function(x) suppressMessages(suppressWarnings(library(x, character.only = TRUE, quietly = T, warn.conflicts = F))))

hashing = read.csv(file.path(path_to_sample_names, "samples.csv"), header = TRUE)
hashing = na.omit(hashing)

assign_sample_name = function(path_in, path_out) {
  sc = readRDS(path_in)
  sc$sample_name = hashing[match(sc$HTO_classification.simplified, hashing[, 1]), 2]
  saveRDS(sc, file = path_out)
}

assign_sample_name(dataset_path, "dataset.Rds")
assign_sample_name(dataset_clean_path, "dataset_clean.Rds")
