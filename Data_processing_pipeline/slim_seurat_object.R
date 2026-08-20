# slim_seurat_object.R
# Strip a Seurat object down to just counts/data + the assays and UMAPs that
# actually exist, and drop QC-only metadata columns before handing the object
# to end users.
#
# This replaces both dietSeurat.R (RNA+ADT) and dietSeurat_RNA.r (RNA only).
# Rather than taking an ADT=YES/NO flag, this script auto-detects which
# assays/reductions are present on the object (Assays()/Reductions()) and
# only keeps/drops what is actually there - so it works unchanged no matter
# which combination of ADT/HTO steps ran upstream.
#
# Usage:
#   Rscript slim_seurat_object.R <dataset_clean input.Rds> <dataset input.Rds> <path_to_sample_names> <EXP.txt>
# (path_to_sample_names and EXP are accepted for call-signature compatibility
#  with the rest of the pipeline but are not used by this script.)

args = commandArgs(TRUE)
dataset_clean_in      = args[1]
dataset_in            = args[2]
path_to_sample_names  = args[3]
EXP                   = args[4]

libs = c("Seurat", "ggplot2", "inflection", "grid", "gridExtra", "UCell")
sapply(libs, function(x) suppressMessages(suppressWarnings(library(x, character.only = TRUE, quietly = T, warn.conflicts = F))))

drop_if_present = function(sc, col) {
  if (col %in% colnames(sc@meta.data)) sc[[col]] = NULL
  sc
}

diet_and_clean = function(path_in, path_out) {

  sc = readRDS(path_in)

  assays_to_keep    = intersect(c('RNA', 'ADT'), Assays(sc))
  dimreducs_to_keep = intersect(c('umap_rna', 'umap_adt'), Reductions(sc))

  sc_slim = DietSeurat(
    sc,
    counts = TRUE,
    data = TRUE,
    scale.data = FALSE,
    features = NULL,
    assays = assays_to_keep,
    dimreducs = dimreducs_to_keep,
    graphs = NULL,
    misc = TRUE
  )

  qc_only_columns = c(
    'Cell_ID', 'HTO_maxID', 'HTO_secondID', 'HTO_margin', 'outliers_nGenes',
    'lowCountADT', 'percent_mito',
    'outliers_deadcells', 'autofluo',
    'signature_1T.enriched', 'signature_1mnp.enriched', 'signature_1B.enriched', 'signature_1ILC.enriched',
    'HTO_classification', 'HTO_classification.global', 'HTO_classification_orig', 'hash.ID',
    'RNA_snn_res.0.25', 'RNA_snn_res.0.5', 'RNA_snn_res.1.5', 'RNA_snn_res.2', 'RNA_snn_res.3', 'RNA_snn_res.4',
    'nFeature_HTO', 'seurat_clusters'
  )
  for (col in qc_only_columns) sc_slim = drop_if_present(sc_slim, col)

  if ('RNA_snn_res.1' %in% colnames(sc_slim@meta.data)) {
    sc_slim$RNA_clusters = sc_slim$RNA_snn_res.1
    sc_slim$RNA_snn_res.1 = NULL
  }
  if ('ADT_snn_res.1' %in% colnames(sc_slim@meta.data)) {
    sc_slim$Protein_clusters = sc_slim$ADT_snn_res.1
    sc_slim$ADT_snn_res.1 = NULL
  }

  saveRDS(sc_slim, path_out)
}

diet_and_clean(dataset_clean_in, 'dataset_clean.Rds')
diet_and_clean(dataset_in, 'dataset.Rds')
