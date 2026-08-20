# David Zemmour
# R/4.0.1
# usage: Rscript rna_qc.R [PROJECT] [path to seuratobject_singlet.Rds]

args = commandArgs(TRUE)
path_to_so = args[1]
th_nFeature_RNA_lo =  as.numeric(args[2])
th_nFeature_RNA_hi =  as.numeric(args[3])
th_percent_mito =  as.numeric(args[4])

libs = c("Seurat", "ggplot2", "scales", "RColorBrewer", "dplyr")
sapply(libs, function(x) suppressMessages(suppressWarnings(library(x, character.only = TRUE, quietly = T, warn.conflicts  = F))))

n = 70
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
mypal = unique(unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))))
mypal = mypal[-4]

message("Read Seurat Object")
so = readRDS(path_to_so)
# No HTO flag is passed to this script, so detect it directly from the
# object - nCount_HTO only exists when an HTO assay was actually created.
HTO_present = "HTO" %in% Assays(so)

message("Run UMAP")
so = NormalizeData(so, verbose = T, normalization.method = "LogNormalize", scale.factor = 10000)
so = FindVariableFeatures(so, selection.method = "vst", nfeatures = 2000, verbose = T)
so = ScaleData(so, assay = "RNA")
so = RunPCA(so, npcs = 100, ndims.print = 1:5, nfeatures.print = 5, reduction.name = "pca_preRNAqc")
x = cumsum((so@reductions$pca_preRNAqc@stdev**2 / sum(so@reductions$pca_preRNAqc@stdev**2)))
ndims = min(which(x >= 0.8))
so = RunUMAP(so, dims = 1:ndims, reduction = "pca_preRNAqc", reduction.name = "umap_preRNAqc")

message(sprintf("Flag outliers_nGenes with < %s  or > %s genes", th_nFeature_RNA_lo, th_nFeature_RNA_hi))
so$outliers_nGenes = so$nFeature_RNA < th_nFeature_RNA_lo | so$nFeature_RNA > th_nFeature_RNA_hi
message(sprintf("Cells outliers_nGenes with < %s genes: %s of %s total cells", th_nFeature_RNA_lo, length(which(so$nFeature_RNA < th_nFeature_RNA_lo)), ncol(so)))
message(sprintf("Cells outliers_nGenes with > %s genes: %s of %s total cells", th_nFeature_RNA_hi, length(which(so$nFeature_RNA > th_nFeature_RNA_hi)), ncol(so)))

message(sprintf("Flag outliers_deadcells with percent_mito > %s", th_percent_mito))
so$percent_mito = Matrix::colSums(so@assays$RNA@counts[rownames(so) %in% rownames(so)[grep("^MT\\-|^Mt\\-|^mt\\-", rownames(so))], ]) / Matrix::colSums(so@assays$RNA@counts) *100
so$outliers_deadcells = so$percent_mito > th_percent_mito
message(sprintf("Flag outliers_deadcells with percent_mito > %s percent: %s / %s total cells", th_percent_mito, length(which(so$outliers_deadcells ==  T)), ncol(so) ))

message("Plot QC")
png("seuratobject_singlet_preRNAfiltering-1.png",width=1600, height=1600)
tmp = data.frame(so@meta.data)
p = ggplot(tmp) + geom_point(aes(nCount_RNA, nFeature_RNA, color = HTO_classification), alpha = I(0.5), size = I(1)) +
  scale_color_manual(values = mypal) +
  scale_x_continuous(trans = log_trans(base = 10), limits = c(100,max(tmp$nCount_RNA))) +
  scale_y_continuous(trans = log_trans(base = 10), limits = c(100,max(tmp$nFeature_RNA)))  +
  annotation_logticks(sides = "bl") +
  geom_hline(yintercept = th_nFeature_RNA_lo, color = "red", linetype = "dashed") +
  geom_hline(yintercept = th_nFeature_RNA_hi, color = "red", linetype = "dashed") +
  ggtitle(paste0('cutoff_nFeature_RNA_lo=',th_nFeature_RNA_lo,'/cutoff_nFeature_RNA_lo=',th_nFeature_RNA_hi)) +
  theme_bw()
#print(p)
print(p + facet_wrap(~HTO_classification, ncol = 2))
dev.off()

png("seuratobject_singlet_preRNAfiltering-2.png",width=1600, height=1600)
tmp = data.frame(so@meta.data, dim1 = so@reductions$umap_preRNAqc@cell.embeddings[,1], dim2 = so@reductions$umap_preRNAqc@cell.embeddings[,2])
p = ggplot(tmp) + geom_point(aes(dim1, dim2, color = outliers_nGenes), size = I(0.5)) +
  theme_bw() + ggtitle(label = sprintf("umap_preRNAqc, %s", "outliers_nGenes")) + scale_color_manual(values = c("grey", "red"))
print(p)
#print(p + facet_wrap(~HTO_classification, ncol = 2))
dev.off()

png("seuratobject_singlet_preRNAfiltering-3.png",width=1600, height=1600)
tmp = data.frame(so@meta.data)
p = ggplot(tmp) + geom_point(aes(nCount_RNA, percent_mito, color = HTO_classification), alpha = I(0.5), size = I(1)) +
  scale_color_manual(values = mypal) +
  scale_x_continuous(trans = log_trans(10), limits = c(100,max(tmp$nCount_RNA))) +
  scale_y_continuous(limits = c(0,100))  +
  annotation_logticks(sides = "b") +
  geom_hline(yintercept = th_percent_mito, color = "red", linetype = "dashed") +
  ggtitle(paste0('cutoff_percent_mito=',th_percent_mito)) +
  theme_bw() #+ geom_density_2d(aes(nCount_RNA, percent_mito)) + ggtitle(label = "All data")
#print(p)
print(p + facet_wrap(~HTO_classification, ncol = 2))
dev.off()

png("seuratobject_singlet_preRNAfiltering-4.png",width=1600, height=1600)
tmp = data.frame(so@meta.data, dim1 = so@reductions$umap_preRNAqc@cell.embeddings[,1], dim2 = so@reductions$umap_preRNAqc@cell.embeddings[,2])
p = ggplot(tmp) + geom_point(aes(dim1, dim2, color = outliers_nGenes), size = I(0.5)) +
  theme_bw() + ggtitle(label = sprintf("umap_preRNAqc, %s", "outliers_nGenes")) + scale_color_manual(values = c("grey", "red"))
#print(p)
#print(p + facet_wrap(~HTO_classification, ncol = 2))

p = ggplot(tmp) + geom_point(aes(dim1, dim2, color = outliers_deadcells), size = I(0.5)) + theme_bw() + ggtitle(label = sprintf("umap_preRNAqc, %s", "outliers_deadcells")) + scale_color_manual(values = c("grey", "red"))
print(p)
#print(p + facet_wrap(~HTO_classification, ncol = 2))

dev.off()

ncells_prefiltering = ncol(so)
tmp = data.frame(so@meta.data)
if (HTO_present) {
  tmp = tmp %>% group_by(HTO_classification.simplified) %>% summarize(ncells = n(),  mean_nCount_RNA = round(mean(nCount_RNA)), mean_nFeature_RNA = round(mean(nFeature_RNA)), mean_nCount_HTO = round(mean(nCount_HTO)),  ncells_outliers_nGenes = length(which(outliers_nGenes)), ncells_outliers_deadcells = length(which(outliers_deadcells))) %>% as.data.frame()
} else {
  tmp = tmp %>% group_by(HTO_classification.simplified) %>% summarize(ncells = n(),  mean_nCount_RNA = round(mean(nCount_RNA)), mean_nFeature_RNA = round(mean(nFeature_RNA)),  ncells_outliers_nGenes = length(which(outliers_nGenes)), ncells_outliers_deadcells = length(which(outliers_deadcells))) %>% as.data.frame()
}
write.table(x = tmp, file = "seuratobject_singlet_preRNAfiltering_QC_stats.txt", quote = F, sep = "\t", row.names = F, col.names = T)

saveRDS(so, file = "seuratobject_singlet_preRNAfiltering.Rds")

message("Filter out outliers")
if (length(which(!(so$outliers_deadcells|so$outliers_nGenes))) == 0) {
  message("All cells filtered out: bad quality data")
} else {
  message(sprintf("%s out of %s cells kept", length(which(!(so$outliers_deadcells|so$outliers_nGenes))), ncol(so)))
  so = so[,!(so$outliers_deadcells|so$outliers_nGenes)]
  source("/n/groups/cbdm_lab/immgen_t/old/seurat_rna_umap_clustering.R")
  pdf("seuratobject_singlet_postRNAfiltering.pdf", 10, 10, useDingbats=FALSE)
  so = Seurat_rna_umap_clustering(seurat_object = so, tsne.method = "FIt-SNE", print_clusters = T, reduction.name.pca = "pca_rna", reduction.name.umap = "umap_rna")
  dev.off()

}

tmp = data.frame(so@meta.data)
if (HTO_present) {
  tmp = tmp %>% group_by(HTO_classification.simplified) %>% summarize(ncells = n(),  mean_nCount_RNA = round(mean(nCount_RNA)), mean_nFeature_RNA = round(mean(nFeature_RNA)), mean_nCount_HTO = round(mean(nCount_HTO)),  ncells_outliers_nGenes = length(which(outliers_nGenes)), ncells_outliers_deadcells = length(which(outliers_deadcells))) %>% as.data.frame()
} else {
  tmp = tmp %>% group_by(HTO_classification.simplified) %>% summarize(ncells = n(),  mean_nCount_RNA = round(mean(nCount_RNA)), mean_nFeature_RNA = round(mean(nFeature_RNA)),  ncells_outliers_nGenes = length(which(outliers_nGenes)), ncells_outliers_deadcells = length(which(outliers_deadcells))) %>% as.data.frame()
}
write.table(x = tmp, file = "seuratobject_singlet_postRNAfiltering_QC_stats.txt", quote = F, sep = "\t", row.names = F, col.names = T)

saveRDS(so, file = "seuratobject_singlet_postRNAfiltering.Rds")

message(sprintf("Number of cells post filtering: %s of %s total cells", ncol(so),ncells_prefiltering ))
