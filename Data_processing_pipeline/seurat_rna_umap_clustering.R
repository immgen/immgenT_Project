# seurat_rna_umap_clustering.R
#
# Portability note: `rna_qc.R` (unchanged - see README.md) calls
# `source("/n/groups/cbdm_lab/immgen_t/old/seurat_rna_umap_clustering.R")`,
# which is a path on the original HMS O2 cluster. If you are not on that
# cluster (e.g. running the GSE311344 use case), that source() call will
# fail. This file is a drop-in replacement with the same function - either:
#   (a) edit that one source() line in rna_qc.R to point here, or
#   (b) place this file at that exact path on your system.
#
# (This is the same helper function Tcell_filter_part2.R defines
# inline for its own RNA clustering step - duplicated here only because
# rna_qc.R was intentionally left unmodified.)

Seurat_rna_umap_clustering = function(seurat_object = so, tsne.method = c("FIt-SNE", "Rtsne"), print_clusters = TRUE, reduction.name.pca = "pca_rna", reduction.name.umap = "umap_rna") {
  require(Seurat)
  seurat_object = NormalizeData(seurat_object, verbose = TRUE, normalization.method = "LogNormalize", scale.factor = 10000, assay = "RNA")
  seurat_object = FindVariableFeatures(seurat_object, selection.method = "vst", nfeatures = 2000, verbose = TRUE, assay = "RNA")
  DefaultAssay(seurat_object) = 'RNA'
  genes = rownames(seurat_object)
  Trav = genes[grepl('^Trav', genes)]; Traj = genes[grepl('^Traj', genes)]; Trac = genes[grepl('^Trac', genes)]
  Trbv = genes[grepl('^Trbv', genes)]; Trbd = genes[grepl('^Trbd', genes)]; Trbj = genes[grepl('^Trbj', genes)]
  Trd  = genes[grepl('^Trd', genes)];  Trg  = genes[grepl('^Trg', genes)]
  TCR_genes = c(Trav, Traj, Trac, Trbv, Trbd, Trbj, Trd, Trg)
  sex_specific_genes = c('Ddx3y', 'Uty', 'Xist', 'Eif2s3y', 'Kdm5d', 'Tsix')
  remove_var_genes = c(TCR_genes, sex_specific_genes)
  seurat_object@assays$RNA@var.features = seurat_object@assays$RNA@var.features[!seurat_object@assays$RNA@var.features %in% remove_var_genes]

  seurat_object = ScaleData(seurat_object, assay = "RNA")
  seurat_object = RunPCA(seurat_object, npcs = 100, ndims.print = 1:5, nfeatures.print = 5, reduction.name = reduction.name.pca, assay = "RNA", var.features = seurat_object@assays$RNA@var.features)

  x = cumsum((seurat_object@reductions[[reduction.name.pca]]@stdev**2 / sum(seurat_object@reductions[[reduction.name.pca]]@stdev**2)))
  ndims = min(which(x >= 0.8))

  message("Running UMAP")
  seurat_object = RunUMAP(seurat_object, dims = 1:ndims, reduction = reduction.name.pca, reduction.name = reduction.name.umap)

  message("RNA clusters...")
  seurat_object = FindNeighbors(seurat_object, reduction = reduction.name.pca, dims = 1:ndims, k.param = 20, verbose = TRUE)
  for (res in c(0.25, 0.5, 1, 1.5, 2, 3, 4)) {
    seurat_object = FindClusters(seurat_object, resolution = res, n.start = 10, algorithm = 1, n.iter = 10)
  }

  if (print_clusters) {
    tmp = data.frame(seurat_object@meta.data, dim1 = seurat_object@reductions[[reduction.name.umap]]@cell.embeddings[, 1], dim2 = seurat_object@reductions[[reduction.name.umap]]@cell.embeddings[, 2])
    for (i in grep("RNA_snn_res.", colnames(seurat_object@meta.data))) {
      p = ggplot(tmp) + geom_point(aes(dim1, dim2, color = tmp[, i]), size = I(1), alpha = I(1)) + theme_bw() +
        theme(axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15), legend.text = element_text(size = 20),
              axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20), legend.title = element_blank()) +
        ggtitle(sprintf("%s, %s", reduction.name.umap, colnames(seurat_object@meta.data)[i]))
      print(LabelClusters(p, id = colnames(seurat_object@meta.data)[i], size = 5))
    }
  }

  return(seurat_object)
}
