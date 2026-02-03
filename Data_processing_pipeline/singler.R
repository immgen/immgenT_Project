# David Zemmour
# R/4.0.1
# usage: Rscript singler.R [path to $SRR_ACCN_seuratobject_postfiltering_NferenceSeuratR.Rda]

args = commandArgs(TRUE)
path_to_so = args[1] #path_to_so = "seuratobject_withHTOADT_singlet_postRNAfiltering.Rds"
reduction.name.umap = "umap_rna"

libs = c("Seurat", "SingleR", "ggplot2", "celldex") 
sapply(libs, function(x) suppressMessages(suppressWarnings(library(x, character.only = TRUE, quietly = T, warn.conflicts  = F))))

library(RColorBrewer)
n = 70
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
mypal = unique(unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))))

message("Read Seurat Object")
so = readRDS(file = path_to_so)

message("SingleR Immgen annotation reference")
ref = ImmGenData()
nsubsets = 10
subs = split(1:ncol(so), ceiling(seq_along(1:ncol(so))/(ceiling(ncol(so)/nsubsets))))
which(duplicated(unlist(subs)))
tmp_list = list()
for (i in 1:nsubsets) { 
  print(i)
  tmp_list[[i]] = SingleR(method = "single", test = so@assays$RNA@data[,subs[[i]]], ref = ref, labels = ref$label.main, genes = "de", fine.tune = T, sd.thres = 1)
}
x = do.call(rbind,tmp_list )
singler_immgen = x
all(rownames(x) == colnames(so))
so$singler_immgen = factor(x$pruned.labels)

message("SingleR annotation mousernaseq reference")
# SingleR() expects reference datasets to be normalized and log-transformed.
ref = MouseRNAseqData()
nsubsets = 10
subs = split(1:ncol(so), ceiling(seq_along(1:ncol(so))/(ceiling(ncol(so)/nsubsets))))
which(duplicated(unlist(subs)))
tmp_list = list()
for (i in 1:nsubsets) { 
  print(i)
  tmp_list[[i]] = SingleR(method = "single", test = so@assays$RNA@data[,subs[[i]]], ref = ref, labels = ref$label.main, genes = "de", fine.tune = T, sd.thres = 1)
}
x = do.call(rbind,tmp_list )
singler_mousernaseq = x
all(rownames(x) == colnames(so))
so$singler_mousernaseq = factor(x$pruned.labels)

save(singler_mousernaseq, file = "singler_results.Rda")

message("Print annotation on RNA UMAP")
pdf("seuratobject_withHTOADT_singlet_postRNAfiltering_withSingleR.pdf", 15, 15, useDingbats=FALSE) 
tmp = data.frame(so@meta.data, dim1 = so@reductions[[reduction.name.umap]]@cell.embeddings[,1], dim2 = so@reductions[[reduction.name.umap]]@cell.embeddings[,2])
tmp$singler_immgen = factor(tmp$singler_immgen)
p = ggplot(tmp) + geom_point(aes(dim1, dim2, color = singler_immgen), size = I(0.5)) + theme_bw() + ggtitle(label = sprintf("%s, %s", reduction.name.umap, "singler_immgen")) + scale_color_manual(values = mypal)
print(LabelClusters(p, id = "singler_immgen", size = 5))
print(LabelClusters(p + facet_wrap(~HTO_classification), id = "singler_immgen", split.by = "HTO_classification"))

tmp = data.frame(so@meta.data, dim1 = so@reductions[[reduction.name.umap]]@cell.embeddings[,1], dim2 = so@reductions[[reduction.name.umap]]@cell.embeddings[,2])
tmp$singler_mousernaseq = factor(tmp$singler_mousernaseq)
p = ggplot(tmp) + geom_point(aes(dim1, dim2, color = singler_mousernaseq), size = I(0.5)) + theme_bw() + ggtitle(label = sprintf("%s, %s", reduction.name.umap, "singler_mousernaseq")) + scale_color_manual(values = mypal)
print(LabelClusters(p, id = "singler_mousernaseq"))
print(LabelClusters(p + facet_wrap(~HTO_classification), id = "singler_mousernaseq", split.by = "HTO_classification"))

dev.off()

saveRDS(so, file = "seuratobject_withHTOADT_singlet_postRNAfiltering_withSingleR.Rds")