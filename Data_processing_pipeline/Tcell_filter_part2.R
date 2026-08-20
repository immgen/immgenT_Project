# Tcell_filt_part2.R
# Classify cells as T vs non-T using the module-score cutoffs in T_cutoffs.csv,
# drop non-T cells, and (if ADT data exists) cluster the surviving T cells on
# both RNA and ADT. No secondary nCount_ADT cutoff here - low-ADT-count cells
# are already dropped upstream in adt_qc.R.
#
# Usage:
#   Rscript Tcell_filt_part2.R <seurat_object.Rds> <T_cutoffs.csv> <EXP.txt> <ADT YES/NO>

args = commandArgs(TRUE)
sc_IGT_path = args[1]
cutoffs_path = args[2]
EXP_path     = args[3]
ADT          = toupper(args[4]) == "YES"

libs = c("Seurat", "ggplot2", "inflection", "grid", "gridExtra", "UCell",
         "reshape2", "dplyr", "fitdistrplus", "ggExtra", "ggrastr", "gatepoints", "rafalib", "scales")
sapply(libs, function(x) suppressMessages(suppressWarnings(library(x, character.only = TRUE, quietly = T, warn.conflicts = F))))

sc_IGT = readRDS(sc_IGT_path)
# No HTO flag is passed to this script, so detect it directly from the
# object - nCount_HTO only exists when an HTO assay was actually created.
HTO_present = "HTO" %in% Assays(sc_IGT)

EXP = read.table(EXP_path)
EXP = as.character(EXP)

cutoffs = read.csv(cutoffs_path, header = FALSE)

cutoff_T   = cutoffs[1, 2]
cutoff_B   = cutoffs[2, 2]
cutoff_ILC = cutoffs[3, 2]
cutoff_MNP = cutoffs[4, 2]

DotPlotHeatmap = function(data = tmp, title = "ADT_seurat_clusters", so_obj) {
  tmp = melt(data, variable.name = "SYMBOL")
  tmp$adt_thr = so_obj@assays$ADT@meta.features[tmp$SYMBOL, "cutoff_crl_norm"]
  tmp2 = tmp %>% group_by(cluster, SYMBOL) %>% summarize(mean = mean(value), freq = mean(value > adt_thr) * 100, n_cells = n()) %>% as.data.frame()
  tmp3 = dcast(data = tmp2, SYMBOL ~ cluster, value.var = "mean")
  rownames(tmp3) = tmp3$SYMBOL
  tmp3 = tmp3[, -1]
  tmp3 = tmp3[rowSums(tmp3[]) > 0, ]
  hc_col = hclust(dist(1 - cor(tmp3, method = "pearson")))
  hc_row = hclust(dist(1 - cor(t(tmp3), method = "pearson")))
  tmp2$cluster = factor(tmp2$cluster, levels = hc_col$labels[hc_col$order])
  tmp2$SYMBOL = factor(tmp2$SYMBOL, levels = hc_row$labels[hc_row$order])
  p = ggplot(tmp2) + geom_point(aes(x = SYMBOL, y = cluster, color = mean, size = freq, alpha = freq)) +
    scale_color_gradient(low = "blue", high = "red") + ggtitle(label = title) + theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.text.x = element_text(angle = 90, vjust = 0, hjust = 1, size = 10), axis.text.y = element_text(angle = 0, size = 15),
          axis.title.x = element_blank(), axis.title.y = element_blank())
  print(p)
}

subset_keep = sc_IGT

percentages = rbind(subset_keep$signature_1T.enriched,
                     subset_keep$signature_1mnp.enriched,
                     subset_keep$signature_1B.enriched,
                     subset_keep$signature_1ILC.enriched)
row.names(percentages) = c("T Percentages", "MNP Percentages", "B Percentages", "ILC Percentages")
percentages = as.data.frame(t(percentages))

p1_T = ggplot() + geom_point(data = percentages, aes(`T Percentages`, `B Percentages`), color = 'blue') +
  xlim(0, max(percentages[, 1])) + ylim(0, max(percentages[, 3])) + theme(text = element_text(size = 10)) +
  ggtitle("ADT T cells") + scale_y_continuous(limits = c(0, max(percentages[, 3])), breaks = seq(0, max(percentages[, 3]), by = 0.01)) +
  geom_vline(xintercept = cutoff_T, linetype = "dashed", color = "red") + geom_hline(yintercept = cutoff_B, linetype = "dashed", color = "red")

p2_T = ggplot() + geom_point(data = percentages, aes(`T Percentages`, `ILC Percentages`), color = 'blue') +
  xlim(0, max(percentages[, 1])) + ylim(0, max(percentages[, 3])) + theme(text = element_text(size = 10)) +
  ggtitle("ADT T cells") + scale_y_continuous(limits = c(0, max(percentages[, 3])), breaks = seq(0, max(percentages[, 3]), by = 0.01)) +
  geom_vline(xintercept = cutoff_T, linetype = "dashed", color = "red") + geom_hline(yintercept = cutoff_ILC, linetype = "dashed", color = "red")

p3_T = ggplot() + geom_point(data = percentages, aes(`T Percentages`, `MNP Percentages`), color = 'blue') +
  xlim(0, max(percentages[, 1])) + theme(text = element_text(size = 10)) +
  ggtitle("ADT T cells") + scale_y_continuous(limits = c(0, max(percentages[, 3])), breaks = seq(0, max(percentages[, 3]), by = 0.01)) +
  geom_vline(xintercept = cutoff_T, linetype = "dashed", color = "red") + geom_hline(yintercept = cutoff_MNP, linetype = "dashed", color = "red")

png(file = "T_labeled_RNA_CellScoring-1.png", width = 1800, height = 1050)
pushViewport(viewport(layout = grid.layout(nrow = 4, ncol = 3)))
define_region = function(row, col) viewport(layout.pos.row = row, layout.pos.col = col)
print(p1_T, vp = define_region(row = 1, col = 1))
print(p2_T, vp = define_region(row = 1, col = 2))
print(p3_T, vp = define_region(row = 1, col = 3))
dev.off()

## Classify T cells vs non-T cells from the module scores

percentages_Tcells = percentages[percentages[["T Percentages"]] > cutoff_T, ]
percentages_Tcells = percentages_Tcells[percentages_Tcells[["MNP Percentages"]] < cutoff_MNP, ]
percentages_Tcells = percentages_Tcells[percentages_Tcells[["B Percentages"]] < cutoff_B, ]
percentages_Tcells = percentages_Tcells[percentages_Tcells[["ILC Percentages"]] < cutoff_ILC, ]

percentages_nonTcells = percentages[percentages[["T Percentages"]] < cutoff_T, ]
percentages_nonTcells = percentages_nonTcells[percentages_nonTcells[["MNP Percentages"]] < cutoff_MNP, ]
percentages_nonTcells = percentages_nonTcells[percentages_nonTcells[["B Percentages"]] < cutoff_B, ]
percentages_nonTcells = percentages_nonTcells[percentages_nonTcells[["ILC Percentages"]] < cutoff_ILC, ]

percentages_nonTcells2 = percentages[percentages[["T Percentages"]] > cutoff_T, ]
percentages_nonTcells2 = percentages_nonTcells2[percentages_nonTcells2[["B Percentages"]] > cutoff_B, ]

removed_cells1_preTCR = rownames(percentages)[!rownames(percentages) %in% rownames(percentages_Tcells)]
removed_cells1 = removed_cells1_preTCR

subset_filter1_keptcells = subset(sc_IGT, cells = removed_cells1, invert = TRUE)

Tcells = subset_filter1_keptcells

nonTcells = subset(sc_IGT, cells = colnames(subset_filter1_keptcells), invert = TRUE)
nonT_cells_removed = table(nonTcells@meta.data$HTO_classification.simplified, nonTcells@meta.data$orig.ident)

## Final QC stats

if (ADT & HTO_present) {
  tmp = data.frame(Tcells@meta.data) %>%
    group_by(HTO_classification.simplified) %>%
    summarize(ncells_postqc = n(), mean_nCount_RNA_postqc = round(mean(nCount_RNA)), mean_nFeature_RNA_postqc = round(mean(nFeature_RNA)),
              mean_nCount_HTO_postqc = round(mean(nCount_HTO)), mean_nCount_ADT_postqc = round(mean(nCount_ADT))) %>% as.data.frame()
} else if (ADT) {
  tmp = data.frame(Tcells@meta.data) %>%
    group_by(HTO_classification.simplified) %>%
    summarize(ncells_postqc = n(), mean_nCount_RNA_postqc = round(mean(nCount_RNA)), mean_nFeature_RNA_postqc = round(mean(nFeature_RNA)),
              mean_nCount_ADT_postqc = round(mean(nCount_ADT))) %>% as.data.frame()
} else if (HTO_present) {
  tmp = data.frame(Tcells@meta.data) %>%
    group_by(HTO_classification.simplified) %>%
    summarize(ncells_postqc = n(), mean_nCount_RNA_postqc = round(mean(nCount_RNA)), mean_nFeature_RNA_postqc = round(mean(nFeature_RNA)),
              mean_nCount_HTO_postqc = round(mean(nCount_HTO))) %>% as.data.frame()
} else {
  tmp = data.frame(Tcells@meta.data) %>%
    group_by(HTO_classification.simplified) %>%
    summarize(ncells_postqc = n(), mean_nCount_RNA_postqc = round(mean(nCount_RNA)), mean_nFeature_RNA_postqc = round(mean(nFeature_RNA))) %>% as.data.frame()
}

tmp2 = data.frame(nonTcells@meta.data) %>% group_by(HTO_classification.simplified) %>% summarize(ncells_nonTcells = n()) %>% as.data.frame()

missing = tmp[, 1][!tmp[, 1] %in% (tmp2[, 1])]
if (length(missing) > 0) {
  missing = as.data.frame(cbind(missing, 0))
  colnames(missing) = colnames(tmp2)
  tmp2 = rbind(tmp2, missing)
}

final_table = merge(tmp2, tmp)
write.csv(final_table, 'seuratobject_singlet_postRNAfiltering_postADTfiltering_postTfiltering.csv')

if (ADT & HTO_present) {
  EXP_final_table = t(as.data.frame(c(round(mean(Tcells$nCount_RNA)), round(mean(Tcells$nFeature_RNA)), round(mean(Tcells$nCount_HTO)), round(mean(Tcells$nCount_ADT)))))
  colnames(EXP_final_table) = c('mean_nCount_RNA_postqc', 'mean_nFeature_RNA_postqc', 'mean_nCount_HTO_postqc', 'mean_nCount_ADT_postqc')
} else if (ADT) {
  EXP_final_table = t(as.data.frame(c(round(mean(Tcells$nCount_RNA)), round(mean(Tcells$nFeature_RNA)), round(mean(Tcells$nCount_ADT)))))
  colnames(EXP_final_table) = c('mean_nCount_RNA_postqc', 'mean_nFeature_RNA_postqc', 'mean_nCount_ADT_postqc')
} else if (HTO_present) {
  EXP_final_table = t(as.data.frame(c(round(mean(Tcells$nCount_RNA)), round(mean(Tcells$nFeature_RNA)), round(mean(Tcells$nCount_HTO)))))
  colnames(EXP_final_table) = c('mean_nCount_RNA_postqc', 'mean_nFeature_RNA_postqc', 'mean_nCount_HTO_postqc')
} else {
  EXP_final_table = t(as.data.frame(c(round(mean(Tcells$nCount_RNA)), round(mean(Tcells$nFeature_RNA)))))
  colnames(EXP_final_table) = c('mean_nCount_RNA_postqc', 'mean_nFeature_RNA_postqc')
}
rownames(EXP_final_table) = EXP
write.csv(EXP_final_table, 'seuratobject_EXP_singlet_postRNAfiltering_postADTfiltering_postTfiltering.csv')

## Label every cell in the full object with its inferred lineage

percentages_nonT = rbind(nonTcells$signature_1T.enriched, nonTcells$signature_1mnp.enriched, nonTcells$signature_1B.enriched, nonTcells$signature_1ILC.enriched)
row.names(percentages_nonT) = c("T Percentages", "MNP Percentages", "B Percentages", "ILC Percentages")
cell_index = apply(percentages_nonT, 2, which.max)
cell_index = as.data.frame(cell_index)
cell_index[cell_index == 2] = "MNP"
cell_index[cell_index == 3] = "B"
cell_index[cell_index == 4] = "ILC"
nonTcells$cell_type = cell_index
Idents(nonTcells) = 'cell_type'

sc_IGT$cell_type = 'n/a'
sc_IGT$Tcell = 'n/a'
sc_IGT$Tcell[colnames(sc_IGT) %in% colnames(Tcells)] = 'Y'
sc_IGT$Tcell[colnames(sc_IGT) %in% colnames(nonTcells)] = 'N'

sc_IGT$cell_type[colnames(sc_IGT) %in% colnames(Tcells)] = 'T'
if ('B' %in% levels(nonTcells))   sc_IGT$cell_type[colnames(sc_IGT) %in% WhichCells(nonTcells, idents = 'B')]   = 'B'
if ('MNP' %in% levels(nonTcells)) sc_IGT$cell_type[colnames(sc_IGT) %in% WhichCells(nonTcells, idents = 'MNP')] = 'MNP'
if ('ILC' %in% levels(nonTcells)) sc_IGT$cell_type[colnames(sc_IGT) %in% WhichCells(nonTcells, idents = 'ILC')] = 'ILC'
sc_IGT$cell_type[colnames(sc_IGT) %in% rownames(percentages_nonTcells)]  = 'NA'
sc_IGT$cell_type[colnames(sc_IGT) %in% rownames(percentages_nonTcells2)] = 'B-T like cells'

Idents(sc_IGT) = 'cell_type'
Idents(sc_IGT) = 'Tcell'

table_column = c('n_RNA_Tcells', 'n_RNA_nonTcells', 'n_final_Tcells')
table_summary = as.data.frame(t(c(dim(percentages_Tcells)[1], length(removed_cells1_preTCR), dim(Tcells)[2])))
colnames(table_summary) = table_column

pdf(file = "Tcell_cleanup_final_numbers.pdf", width = 10, height = 5, useDingbats = FALSE)
grid.table(table_summary)
dev.off()

## RNA clustering of the final T cells (identical for every configuration)

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

DefaultAssay(Tcells) = "RNA"
Final_Tcells = Seurat_rna_umap_clustering(seurat_object = Tcells, tsne.method = "FIt-SNE", print_clusters = TRUE, reduction.name.pca = "pca_rna", reduction.name.umap = "umap_rna")

if (ADT) {

  ## ADT clustering of the final T cells, plus the RNA-vs-ADT backgate plots for the removed non-T cells

  data = as.data.frame(t(as.matrix(sc_IGT@assays$ADT@counts)))
  data = log2(data + 1)
  data = data[, order(colnames(data))]

  data2 = data[rownames(data) %in% colnames(Tcells), ]
  data  = data[!rownames(data) %in% colnames(Tcells), ]

  proX = "THY1.2"
  p1 = ggplot(data2, aes(data2[, proX], data2[, "CD3"])) + geom_point(color = 'gray') +
    geom_point(data = data, aes(data[, proX], data[, "CD3"]), color = 'red') +
    xlab(proX) + ylab("CD3") + xlim(0, max(data2[, proX])) + ylim(0, max(data2[, proX])) + theme_bw() + theme(text = element_text(size = 12)) + ggtitle("RNA non-T cells")
  p2 = ggplot(data2, aes(data2[, proX], data2[, "CD2"])) + geom_point(color = 'gray') +
    geom_point(data = data, aes(data[, proX], data[, "CD2"]), color = 'red') +
    xlab(proX) + ylab("CD2") + xlim(0, max(data2[, proX])) + ylim(0, max(data2[, proX])) + theme_bw() + theme(text = element_text(size = 12)) + ggtitle("RNA non-T cells")
  p4 = ggplot(data2, aes(data2[, proX], data2[, "TCRGD"])) + geom_point(color = 'gray') +
    geom_point(data = data, aes(data[, proX], data[, "TCRGD"]), color = 'red') +
    xlab(proX) + ylab("TCRGD") + xlim(0, max(data2[, proX])) + ylim(0, max(data2[, proX])) + theme_bw() + theme(text = element_text(size = 12)) + ggtitle("RNA non-T cells")
  p5 = ggplot(data2, aes(data2[, proX], data2[, "TCRB"])) + geom_point(color = 'gray') +
    geom_point(data = data, aes(data[, proX], data[, "TCRB"]), color = 'red') +
    xlab(proX) + ylab("TCRB") + xlim(0, max(data2[, proX])) + ylim(0, max(data2[, proX])) + theme_bw() + theme(text = element_text(size = 12)) + ggtitle("RNA non-T cells")

  p1_T = ggplot(data, aes(data[, proX], data[, "CD3"])) + geom_point(color = 'gray') +
    geom_point(data = data2, aes(data2[, proX], data2[, "CD3"]), color = 'red') +
    xlab(proX) + ylab("CD3") + xlim(0, max(data2[, proX])) + ylim(0, max(data2[, proX])) + theme_bw() + theme(text = element_text(size = 12)) + ggtitle("RNA T cells")
  p2_T = ggplot(data, aes(data[, proX], data[, "CD2"])) + geom_point(color = 'gray') +
    geom_point(data = data2, aes(data2[, proX], data2[, "CD2"]), color = 'red') +
    xlab(proX) + ylab("CD2") + xlim(0, max(data2[, proX])) + ylim(0, max(data2[, proX])) + theme_bw() + theme(text = element_text(size = 12)) + ggtitle("RNA T cells")
  p4_T = ggplot(data, aes(data[, proX], data[, "TCRGD"])) + geom_point(color = 'gray') +
    geom_point(data = data2, aes(data2[, proX], data2[, "TCRGD"]), color = 'red') +
    xlab(proX) + ylab("TCRGD") + xlim(0, max(data2[, proX])) + ylim(0, max(data2[, proX])) + theme_bw() + theme(text = element_text(size = 12)) + ggtitle("RNA T cells")
  p5_T = ggplot(data, aes(data[, proX], data[, "TCRB"])) + geom_point(color = 'gray') +
    geom_point(data = data2, aes(data2[, proX], data2[, "TCRB"]), color = 'red') +
    xlab(proX) + ylab("TCRB") + xlim(0, max(data2[, proX])) + ylim(0, max(data2[, proX])) + theme_bw() + theme(text = element_text(size = 12)) + ggtitle("RNA T cells")

  png(file = "Tcell_cleanup_ADT_post_filtering.png", width = 1000, height = 1200)
  pushViewport(viewport(layout = grid.layout(nrow = 4, ncol = 2)))
  print(p1, vp = define_region(row = 1, col = 1)); print(p2, vp = define_region(row = 1, col = 2))
  print(p4, vp = define_region(row = 2, col = 1)); print(p5, vp = define_region(row = 2, col = 2))
  print(p1_T, vp = define_region(row = 3, col = 1)); print(p2_T, vp = define_region(row = 3, col = 2))
  print(p4_T, vp = define_region(row = 4, col = 1)); print(p5_T, vp = define_region(row = 4, col = 2))
  dev.off()

  reduction.name.pca = "pca_adt"
  reduction.name.umap = "umap_adt"
  Seurat_adt_umap_clustering = function(seurat_object = so, VariableFeatures = rownames(so[["ADT"]])[!grepl(rownames(so[["ADT"]]), pattern = "unmapped|Isotype|TCRV")], reduction.name.pca = "pca_adt", reduction.name.umap = "umap_adt") {
    so = seurat_object
    DefaultAssay(so) = "ADT"
    so = NormalizeData(so, assay = "ADT", normalization.method = "CLR")
    VariableFeatures(so, assay = "ADT") = VariableFeatures
    so = ScaleData(so, assay = "ADT")
    so = RunPCA(object = so, assay = "ADT", reduction.name = reduction.name.pca)
    x = cumsum((so@reductions[[reduction.name.pca]]@stdev**2 / sum(so@reductions[[reduction.name.pca]]@stdev**2)))
    ndims = min(which(x >= 0.8))
    so = FindNeighbors(so, reduction = "pca_adt", dims = 1:ndims, k.param = 20, verbose = TRUE)
    so = RunUMAP(so, dims = 1:ndims, reduction = reduction.name.pca, reduction.name = reduction.name.umap)
    so = FindClusters(so, resolution = 1, n.start = 10, algorithm = 1, n.iter = 10)
    return(so)
  }

  DefaultAssay(Final_Tcells) = "ADT"
  Final_Tcells = Seurat_adt_umap_clustering(seurat_object = Final_Tcells, VariableFeatures = rownames(Final_Tcells[["ADT"]])[!grepl(rownames(Final_Tcells[["ADT"]]), pattern = "unmapped|Isotype|TCRV")], reduction.name.pca = reduction.name.pca, reduction.name.umap = reduction.name.umap)

  tmp = data.frame(t(as.matrix(Final_Tcells@assays$ADT@data)), cluster = Final_Tcells$seurat_clusters)
  png(file = "adt_qc_postTfilt_Dotplot_ADT_ADTclusters.png", width = 1900, height = 600)
  DotPlotHeatmap(data = tmp, title = "ADT_seurat_clusters", so_obj = Final_Tcells)
  dev.off()
}

Idents(Final_Tcells) = 'sample_name'
saveRDS(Final_Tcells, "seuratobject_singlet_postRNAfiltering_postADTfiltering_postTfiltering.Rds")
saveRDS(sc_IGT, "seuratobject_singlet_postRNAfiltering_postADTfiltering.Rds")
