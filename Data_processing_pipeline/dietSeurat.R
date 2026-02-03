args = commandArgs(TRUE)
dataset_clean = args[1]
dataset = args[2]
path_to_sample_names = args[3]
IGT = args[4]

libs = c("Seurat", "ggplot2","inflection","grid","gridExtra","UCell")
sapply(libs, function(x) suppressMessages(suppressWarnings(library(x, character.only = TRUE, quietly = T, warn.conflicts  = F))))

sc <- readRDS(dataset_clean)


sc_slim<- DietSeurat(
  sc,
  counts = TRUE,
  data = TRUE,
  scale.data = FALSE,
  features = NULL,
  assays = c('RNA','ADT'),
  dimreducs = c('umap_rna','umap_adt'),
  graphs = NULL,
  misc = TRUE
)


sc_slim$Cell_ID <- NULL
sc_slim$HTO_maxID <- NULL
sc_slim$HTO_secondID <- NULL
sc_slim$HTO_margin <- NULL
sc_slim$outliers_nGenes <- NULL
sc_slim$lowCountADT <- NULL
sc_slim$singler_immgen <- NULL
sc_slim$singler_mousernaseq <- NULL
sc_slim$percent_mito <- NULL
sc_slim$outliers_deadcells <- NULL
sc_slim$autofluo <- NULL
sc_slim$signature_1T.enriched <- NULL
sc_slim$signature_1mnp.enriched <- NULL
sc_slim$signature_1B.enriched <- NULL
sc_slim$signature_1ILC.enriched <- NULL
sc_slim$HTO_classification <- NULL
sc_slim$HTO_classification.global<- NULL
sc_slim$hash.ID<- NULL
sc_slim$RNA_snn_res.1.5 <- NULL
sc_slim$nFeature_HTO <- NULL
sc_slim$RNA_snn_res.0.5 <- NULL
sc_slim$RNA_snn_res.2 <- NULL
sc_slim$RNA_snn_res.3 <- NULL
sc_slim$RNA_snn_res.4 <- NULL
sc_slim$RNA_snn_res.0.25 <- NULL
sc_slim$HTO_classification_orig <- NULL
sc_slim$seurat_clusters <- NULL
sc_slim$RNA_clusters <- sc_slim$RNA_snn_res.1
sc_slim$Protein_clusters <- sc_slim$ADT_snn_res.1
sc_slim$RNA_snn_res.1 <- NULL
sc_slim$ADT_snn_res.1 <- NULL



saveRDS(sc_slim,'dataset_clean.Rds')

sc <- readRDS(dataset)

sc_slim<- DietSeurat(
  sc,
  counts = TRUE,
  data = TRUE,
  scale.data = FALSE,
  features = NULL,
  assays = c('RNA','ADT'),
  dimreducs = c('umap_rna','umap_adt'),
  graphs = NULL,
  misc = TRUE
)

sc_slim$Cell_ID <- NULL
sc_slim$HTO_maxID <- NULL
sc_slim$HTO_secondID <- NULL
sc_slim$HTO_margin <- NULL
sc_slim$outliers_nGenes <- NULL
sc_slim$lowCountADT <- NULL
sc_slim$singler_immgen <- NULL
sc_slim$singler_mousernaseq <- NULL
sc_slim$percent_mito <- NULL
sc_slim$outliers_deadcells <- NULL
sc_slim$autofluo <- NULL
sc_slim$signature_1T.enriched <- NULL
sc_slim$signature_1mnp.enriched <- NULL
sc_slim$signature_1B.enriched <- NULL
sc_slim$signature_1ILC.enriched <- NULL
sc_slim$HTO_classification <- NULL
sc_slim$HTO_classification.global<- NULL
sc_slim$hash.ID<- NULL
sc_slim$RNA_snn_res.1.5 <- NULL
sc_slim$nFeature_HTO <- NULL
sc_slim$RNA_snn_res.0.5 <- NULL
sc_slim$RNA_snn_res.2 <- NULL
sc_slim$RNA_snn_res.3 <- NULL
sc_slim$RNA_snn_res.4 <- NULL
sc_slim$RNA_snn_res.0.25 <- NULL
sc_slim$HTO_classification_orig <- NULL
sc_slim$seurat_clusters <- NULL
sc_slim$RNA_clusters <- sc_slim$RNA_snn_res.1
sc_slim$Protein_clusters <- sc_slim$ADT_snn_res.1
sc_slim$RNA_snn_res.1 <- NULL
sc_slim$ADT_snn_res.1 <- NULL

saveRDS(sc_slim,'dataset.Rds')
