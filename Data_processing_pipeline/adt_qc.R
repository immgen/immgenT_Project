### David Zemmour
# R/4.0.1
# usage: Rscript adt_qc.R [path to seuratobject_withHTOADT_singlet_postRNAfiltering_withSingleR.Rds] [path_to_metadata_adt] [th_nCount_ADT_lo] [n_isotype_ctrl_signal_to_flag]

args = commandArgs(TRUE)
path_to_so = args[1]
path_to_metadata_adt = args[2]
th_nCount_ADT_lo = as.numeric(args[3])
n_isotype_ctrl_signal_to_flag = as.numeric(args[4])
#path_to_so = "seuratobject_withHTOADT_singlet_postRNAfiltering_withSingleR.Rds"
#path_to_metadata_adt = "adt_hash_seq_file_April2021_panel_metadata.csv"#"adt_hash_seq_file_feb20panel.csv"
#th_nCount_ADT_lo = 500
#n_isotype_ctrl_signal_to_flag = 2
reduction.name.pca = "pca_adt"
reduction.name.umap = "umap_adt"

libs = c("Seurat", "ggplot2", "reshape2", "dplyr", "fitdistrplus", "ggExtra", "ggrastr", "gridExtra", "gatepoints", "rafalib", "scales") 
sapply(libs, function(x) suppressMessages(suppressWarnings(library(x, character.only = TRUE, quietly = T, warn.conflicts  = F))))

library(RColorBrewer)
n = 70
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
mypal = unique(unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))))
mypal = mypal[-4]
ColorRamp = colorRampPalette(brewer.pal(n = 9,name = "YlOrRd"))(10)

##Functions

FindThresholdADT = function(adt_counts = so@assays$ADT@counts, clusters = so$ADT_snn_res.1, use_n_min_clusters = 2, quantile = 99, remove_top_percent = 0.005) {
  
  data = as.matrix(adt_counts)
  data = melt(data)
  data$cluster = clusters[match(data$Var2, names(clusters))]
  
  data_summarized = data %>% group_by(Var1, cluster) %>% summarize_at(.vars = c("value"), .funs = c(mean = mean)) %>% as.data.frame()
  #data_summarized2 = dcast(data_summarized, Var1 ~cluster)
  #rownames(data_summarized2) = data_summarized2[,1]
  #data_summarized2 = data_summarized2[,-1]
  #data_summarized2
  data_summarized2 = data_summarized %>% group_by(Var1) %>% top_n(use_n_min_clusters, -mean) %>% as.data.frame()
  #as.character(data_summarized2[data_summarized2$Var1 %in% adt,"cluster"])
  
  cutoff = c()
  cutoff_crl_norm = c()
  for (adt in rownames(adt_counts)) {
    min_cluster_values = as.numeric(adt_counts[adt, clusters %in% as.character(data_summarized2[data_summarized2$Var1 %in% adt,"cluster"])])
    min_cluster_values = sort(min_cluster_values)[1:floor(length(min_cluster_values)*(1-remove_top_percent))]
    if (range(min_cluster_values)[2] == 0) {
      cutoff[adt] = 0
      cutoff_crl_norm[adt] = 0
      message(sprintf("Cutoff for %s : %s reads, %s crl_norm", adt, cutoff[adt],  cutoff_crl_norm[adt]))
    } else {
      fit = fitdist(data = min_cluster_values, distr = "nbinom", method = "mle") #control=list(trace=1, REPORT=1)
      cutoff[adt] = as.numeric(x = quantile(x = fit, probs = 0.99)$quantiles[1])
      cutoff_crl_norm[adt] = cutoff[adt] / exp(mean(log(as.numeric(adt_counts[adt,])+1))) 
      message(sprintf("Cutoff for %s : %s reads, %s crl_norm", adt, cutoff[adt],  cutoff_crl_norm[adt]))
    }
  }
  
  return(data.frame(cutoff = cutoff, cutoff_crl_norm = cutoff_crl_norm))
  
}

Plot2Markers = function(so, x_adt, y_adt) {
  x = as.numeric(so@assays$ADT@data[x_adt,])
  y = as.numeric(so@assays$ADT@data[y_adt,])
  tmp = data.frame(x, y, so@meta.data)
  p = ggplot(tmp) + geom_point_rast(aes(x, y), size = 0.5, alpha = 0.5, raster.dpi = 100) +
    geom_hline(yintercept = so@assays$ADT@meta.features[x_adt,"cutoff_crl_norm"], color = "red", linetype = "dashed") + 
    geom_vline(xintercept = so@assays$ADT@meta.features[y_adt,"cutoff_crl_norm"], color = "red", linetype = "dashed") + 
    xlab(sprintf("%s (CRL)", x_adt)) +
    ylab(sprintf("%s (CRL)", y_adt)) +
    theme_classic()
  
  return(p)
}

Backgate = function(so = so, x_adt = "CD3", y_adt = "B220", x_adt_bckgplots = "CD3", size_dot_bckgplots = 1, plot_marginals = T, gate = NULL) {
require(gatepoints)
require(ggplot2)
require(ggrastr)

    x = as.numeric(so@assays$ADT@data[x_adt,])
    y = as.numeric(so@assays$ADT@data[y_adt,])
  if(is.null((gate))) {
    tmp = data.frame(x, y)
    plot(x,y, xlab = x_adt, ylab = y_adt)
    abline(v = so@assays$ADT@meta.features[x_adt,"cutoff_crl_norm"], col = "red")
    abline(h = so@assays$ADT@meta.features[y_adt,"cutoff_crl_norm"], col = "red")
    selectedPoints = fhs(tmp[,c("x", "y")], mark = TRUE)
    points(x[as.numeric(selectedPoints)], y[as.numeric(selectedPoints)], col = "red")
    gate = colnames(so@assays$ADT@data)[as.numeric(selectedPoints)]
  } else {
    message("Gated cells provided, plotting...")
  }


p_list = list()

message("Plot original gate in first plot")
tmp = data.frame(x, y, gate = colnames(so@assays$ADT@data) %in% gate)
p_list[[1]] = ggplot(tmp) + geom_point_rast(aes(x, y), size = 0.5, alpha = 0.5, raster.dpi = 100) +
    geom_hline(yintercept = so@assays$ADT@meta.features[y_adt,"cutoff_crl_norm"], color = "red", linetype = "dashed") + 
    geom_vline(xintercept = so@assays$ADT@meta.features[x_adt,"cutoff_crl_norm"], color = "red", linetype = "dashed") + 
    geom_point(data = tmp[tmp$gate == T, ], aes(x,y), color = "blue", size = size_dot_bckgplots) +
    ggtitle(label = "Original gate") +
    xlab(sprintf("%s (CRL)", x_adt)) +
    ylab(sprintf("%s (CRL)", y_adt)) +
    theme_classic()

for (adt in sort(rownames(so@assays$ADT@data))) {
  print(adt)
  x = as.numeric(so@assays$ADT@data[x_adt_bckgplots,])
  y = as.numeric(so@assays$ADT@data[adt,])
  tmp = data.frame(x, y, gate = colnames(so@assays$ADT@data) %in% gate)
  p_list[[adt]] = ggplot(tmp) + geom_point_rast(aes(x, y), size = 0.5, alpha = 0.5, raster.dpi = 100) +
    geom_hline(yintercept = so@assays$ADT@meta.features[adt,"cutoff_crl_norm"], color = "red", linetype = "dashed") + 
    geom_vline(xintercept = so@assays$ADT@meta.features[x_adt_bckgplots,"cutoff_crl_norm"], color = "red", linetype = "dashed") + 
    geom_point(data = tmp[tmp$gate == T, ], aes(x,y), color = "red", size = size_dot_bckgplots) +
    xlab(sprintf("%s (CRL)", x_adt_bckgplots)) +
    ylab(sprintf("%s (CRL)", adt)) +
    theme_classic()
}

ret = list(gate, p_list)

return(ret)
}

DotPlotHeatmap = function(data = tmp, title = "ADT_snn_res.1") {
  tmp = melt(data, variable.name = "SYMBOL")
  tmp$adt_thr = so@assays$ADT@meta.features[tmp$SYMBOL,"cutoff_crl_norm"]
  tmp2 = tmp %>% group_by(cluster, SYMBOL) %>% summarize(mean = mean(value), freq = mean(value > adt_thr)*100, n_cells = n()) %>% as.data.frame()
  ColorRamp = rev(rainbow(10, end = 4/6))
  tmp3 = dcast(data = tmp2, SYMBOL~cluster, value.var = "mean")
  rownames(tmp3) = tmp3$SYMBOL
  tmp3 = tmp3[,-1]
  tmp3 = tmp3[rowSums(tmp3[])>0,]
  hc_col = hclust(dist(1-cor(tmp3, method = "pearson")))
  hc_row = hclust(dist(1-cor(t(tmp3), method = "pearson")))

  tmp2$cluster = factor(tmp2$cluster, levels = hc_col$labels[hc_col$order])
  tmp2$SYMBOL = factor(tmp2$SYMBOL, levels = hc_row$labels[hc_row$order])

  p = ggplot(tmp2) + geom_point(aes(x = SYMBOL, y = cluster, color = mean, size = freq, alpha = freq)) + scale_color_gradient(low = "blue", high = "red") + ggtitle(label = title) + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x  = element_text(angle=90, vjust=0,  hjust=1, size=10), axis.text.y  = element_text(angle=0, size=15), axis.title.x = element_blank(), axis.title.y = element_blank())
  print(p)
}

Seurat_adt_umap_clustering = function(seurat_object = so, VariableFeatures =  rownames(so[["ADT"]])[!grepl(rownames(so[["ADT"]]), pattern = "unmapped|Isotype|TCRV")],reduction.name.pca = "pca_adt", reduction.name.umap = "umap_adt") {
  so = seurat_object
  DefaultAssay(so) = "ADT"
  message("CLR normalization")
  so = NormalizeData(so, assay = "ADT", normalization.method = "CLR") 
  message("VariableFeatures:")
  print(VariableFeatures)
  VariableFeatures(so, assay = "ADT") = VariableFeatures
  so = ScaleData(so, assay = "ADT")
  so = RunPCA(object = so, assay = "ADT", reduction.name = reduction.name.pca)
  x = cumsum((so@reductions[[reduction.name.pca]]@stdev**2 / sum(so@reductions[[reduction.name.pca]]@stdev**2)))
  ndims = min(which(x >= 0.8))
  so = FindNeighbors(so, reduction = "pca_adt", dims = 1:ndims, k.param = 20, verbose = T)
  message("UMAP ADT")
  so = RunUMAP(so, dims = 1:ndims, reduction = reduction.name.pca, reduction.name = reduction.name.umap)
  message("ADT clustering")
  so = FindClusters(so, resolution = 1, n.start = 10, algorithm = 1, n.iter = 10)
  return(so)
}

##Pipeline
message("Loading the data...")
so = readRDS(path_to_so)

tmp = sort(log10(rowSums(so@assays$ADT@counts)), decreasing = T)
tmp = tmp[!is.na(tmp) & !is.infinite(tmp)]



pdf("adt_qc_plot1_barplot_count.pdf", height = 10, width = 20, useDingbats = F)
barplot(tmp, las =2, horiz = F, cex.names = 0.5)
dev.off()

message(sprintf("Adding the ADT metadata from %s...",path_to_metadata_adt ))
meta.features = read.csv(file = path_to_metadata_adt, header = T, row.names=NULL, sep = ",")
meta.features$Protein_Symbol = gsub(pattern = "_", replacement = "-", meta.features$Protein_Symbol)
all(rownames(so@assays$ADT@meta.features) %in% meta.features$Protein_Symbol)
fd = meta.features[match(rownames(so@assays$ADT@meta.features), meta.features$Protein_Symbol),]
# fd <- na.omit(fd)
rownames(fd) = fd$Protein_Symbol
so@assays$ADT@meta.features = fd

##Prefiltering analysis
message("UMAP and clustering prefiltering...")
so = Seurat_adt_umap_clustering(seurat_object = so, VariableFeatures =  rownames(so[["ADT"]])[!grepl(rownames(so[["ADT"]]), pattern = "unmapped|Isotype|TCRV")],reduction.name.pca = reduction.name.pca, reduction.name.umap = reduction.name.umap)

message("Find ADT cutoffs prefiltering...")
adt_thrs = FindThresholdADT()
so@assays$ADT@meta.features$cutoff = adt_thrs$cutoff[match(rownames(so@assays$ADT@meta.features), rownames(adt_thrs))]
so@assays$ADT@meta.features$cutoff_crl_norm = adt_thrs$cutoff_crl_norm[match(rownames(so@assays$ADT@meta.features), rownames(adt_thrs))]
pdf("adt_qc_plot2_adt_thrs_prefiltering.pdf", width = 5, height = 5, useDingbats = F)
for (adt in sort(rownames(so@assays$ADT@counts))) {
  print(adt)
  x_adt = "TCRB"
  x = as.numeric(so@assays$ADT@data[x_adt,])
  #shist(x, unit = 0.5)
  y = as.numeric(so@assays$ADT@data[adt,])
  #shist(y, unit = 0.5, main = adt)
  #abline(v = adt_thrs[adt,"cutoff_crl_norm"], col = "red")
  tmp = data.frame(x, y)
  p = ggplot(tmp) + geom_point_rast(aes(x, y), size = 0.5, alpha = 0.5, raster.dpi = 100) +
    geom_hline(yintercept = so@assays$ADT@meta.features[adt,"cutoff_crl_norm"], color = "red", linetype = "dashed") + 
    geom_vline(xintercept = so@assays$ADT@meta.features[x_adt,"cutoff_crl_norm"], color = "red", linetype = "dashed") + 
    xlab(sprintf("%s (CRL)", x_adt)) +
    ylab(sprintf("%s (CRL)", adt)) +
    theme_classic()
  q = ggMarginal(p, type = "density")
  print(q,  newpage = TRUE) 
  
  tmp = data.frame(so@meta.data, dim1 = so@reductions[["umap_adt"]]@cell.embeddings[,1], dim2 =so@reductions[["umap_adt"]]@cell.embeddings[,2], sub_cells =  as.numeric(so@assays$ADT@data[adt,]) > so@assays$ADT@meta.features[adt,"cutoff_crl_norm"], size = so@assays$ADT@data[adt,]) 
  p = ggplot(tmp) +  geom_point_rast(aes(dim1, dim2, color = size), size = 0.5, alpha = 0.5, raster.dpi = 100) + scale_color_gradientn(colours = ColorRamp) + 
    theme(axis.text.x  = element_text(size=15), axis.text.y  = element_text(size=15), legend.text=element_text(size=20), axis.title.x = element_text(size=20) , axis.title.y = element_text(size=20), legend.title = element_blank()) + ggtitle(sprintf("%s, %s, cutoff %s","UMAP ADT", adt, so@assays$ADT@meta.features[adt,"cutoff_crl_norm"]))
  print(p)
  
  tmp = data.frame(so@meta.data, dim1 = so@reductions[["umap_rna"]]@cell.embeddings[,1], dim2 =so@reductions[["umap_rna"]]@cell.embeddings[,2], sub_cells =  as.numeric(so@assays$ADT@data[adt,]) > so@assays$ADT@meta.features[adt,"cutoff_crl_norm"], size = so@assays$ADT@data[adt,]) 
  p =  ggplot(tmp) +  geom_point_rast(aes(dim1, dim2, color = size), size = 0.5, alpha = 0.5, raster.dpi = 100) + scale_color_gradientn(colours = ColorRamp) + 
    scale_fill_gradientn(colours = ColorRamp) + theme( axis.text.x  = element_text(size=15), axis.text.y  = element_text(size=15), legend.text=element_text(size=20), axis.title.x = element_text(size=20) , axis.title.y = element_text(size=20), legend.title = element_blank()) + ggtitle(sprintf("%s, %s, cutoff %s","UMAP RNA", adt, so@assays$ADT@meta.features[adt,"cutoff_crl_norm"]))
  print(p)
}
dev.off()

png("adt_qc_plot3_UMAP_clustering_prefiltering-1.png",width=1600, height=1600)
#tmp = data.frame(so@meta.data, dim1 = so@reductions[[reduction.name.umap]]@cell.embeddings[,1], dim2 =so@reductions[[reduction.name.umap]]@cell.embeddings[,2] ) 
#p = ggplot(tmp) + geom_point(aes(dim1, dim2, color = ADT_snn_res.1), size = I(1), alpha = I(1)) + scale_color_manual(values = mypal) + theme_bw() + theme(axis.text.x  = element_text(size=15), axis.text.y  = element_text(size=15), legend.text=element_text(size=20), axis.title.x = element_text(size=20) , axis.title.y = element_text(size=20), legend.title = element_blank()) + ggtitle(sprintf("%s, %s",reduction.name.umap, "ADT_snn_res.1"))
#print(LabelClusters(p, id = "ADT_snn_res.1", size = 5))
#print(p + facet_wrap(~HTO_classification))

#p = ggplot(tmp) + geom_point(aes(dim1, dim2, color = singler_immgen), size = I(1), alpha = I(1)) + scale_color_manual(values = mypal) + theme_bw() + theme(axis.text.x  = element_text(size=15), axis.text.y  = element_text(size=15), legend.text=element_text(size=20), axis.title.x = element_text(size=20) , axis.title.y = element_text(size=20), legend.title = element_blank()) + ggtitle(sprintf("%s, %s",reduction.name.umap, "singler_immgen"))
#LabelClusters(p, id = "ADT_snn_res.1", size = 5)
#print(p + facet_wrap(~HTO_classification))

message(sprintf("Flag cells with ADT coverage less than %s", th_nCount_ADT_lo))
so$lowCountADT = so$nCount_ADT < th_nCount_ADT_lo
message(sprintf("Cells with < %s nCount_ADT: %s of %s total cells", th_nCount_ADT_lo, length(which(so$lowCountADT)), ncol(so)))

tmp = data.frame(so@meta.data, dim1 = so@reductions[[reduction.name.umap]]@cell.embeddings[,1], dim2 =so@reductions[[reduction.name.umap]]@cell.embeddings[,2] ) 
p = ggplot(tmp) + geom_point(aes(nCount_RNA, nCount_ADT, color = HTO_classification), alpha = I(0.5), size = I(1)) + 
  scale_color_manual(values = mypal) +
  scale_x_continuous(trans = log_trans(base = 10), limits = c(10,max(tmp$nCount_RNA))) + 
  scale_y_continuous(trans = log_trans(base = 10), limits = c(10,max(tmp$nCount_ADT)))  + 
  annotation_logticks(sides = "bl") + 
  geom_hline(yintercept = th_nCount_ADT_lo, color = "red", linetype = "dashed") + 
  theme_bw()
#print(p)
print(p + facet_wrap(~HTO_classification))
#print(p + facet_wrap(~singler_immgen))
dev.off()

png("adt_qc_plot3_UMAP_clustering_prefiltering-2.png",width=1600, height=1600)
p = ggplot(tmp) +  geom_point_rast(aes(dim1, dim2), size = 0.5, alpha = 0.5, raster.dpi = 100) +
  geom_point(data = tmp[tmp$lowCountADT == T,], aes(dim1, dim2), size = 5, alpha = 0.5, color = "red") + theme_bw() + theme(axis.text.x  = element_text(size=15), axis.text.y  = element_text(size=15), legend.text=element_text(size=20), axis.title.x = element_text(size=20) , axis.title.y = element_text(size=20), legend.title = element_blank()) + ggtitle(sprintf("%s, ADT counts <  %s",reduction.name.umap, th_nCount_ADT_lo))
print(p)
#print(p + facet_wrap(~HTO_classification, ncol = 2))
dev.off()

png("adt_qc_plot3_UMAP_clustering_prefiltering-3.png",width=1600, height=1600)
message(sprintf("Flag non-specific binding cells with cells that have strong signal for >= %s isotype contols ...", n_isotype_ctrl_signal_to_flag))
isotypes = rownames(so@assays$ADT@data)[grepl("Isotype", rownames(so@assays$ADT@data))]
ct = sapply(isotypes, function(iso)  so@assays$ADT@data[iso,] > so@assays$ADT@meta.features[iso,"cutoff_crl_norm"])
ct = t(ct)
autofluo = names(which(colSums(ct) >= n_isotype_ctrl_signal_to_flag))
so$autofluo = colnames(so) %in% autofluo
message(sprintf("Cells with strong signal >= %s isotype controls: %s of %s total cells", n_isotype_ctrl_signal_to_flag, length(which(so$autofluo)), ncol(so)))
tmp = data.frame(so@meta.data, dim1 = so@reductions[[reduction.name.umap]]@cell.embeddings[,1], dim2 =so@reductions[[reduction.name.umap]]@cell.embeddings[,2] ) 

p = ggplot(tmp) +  geom_point_rast(aes(dim1, dim2), size = 0.5, alpha = 0.5, raster.dpi = 100) +
  geom_point(data = tmp[tmp$autofluo == T,], aes(dim1, dim2), size = 5, alpha = 0.5, color = "red") + theme_bw() + theme(axis.text.x  = element_text(size=15), axis.text.y  = element_text(size=15), legend.text=element_text(size=20), axis.title.x = element_text(size=20) , axis.title.y = element_text(size=20), legend.title = element_blank()) + ggtitle(sprintf("%s, non-specific binders",reduction.name.umap))
print(p)
#print(p + facet_wrap(~HTO_classification))
dev.off()

bckgate = Backgate(so = so, x_adt = "Isotype.Mouse.IgG1k", y_adt = "Isotype.Mouse.IgG2ak", x_adt_bckgplots = "TCRB", size_dot_bckgplots = 1, plot_marginals = T, gate = autofluo)
gate = bckgate[[1]]
p_list = bckgate[[2]]
#q_list = lapply(p_list, function(p) { ggMarginal(p, type = "density")})
n = ceiling(sqrt(length(p_list)))
pdf("adt_qc_plot4_backgate_autofluorescence.pdf", 30, 30)
grid.arrange(grobs = p_list, widths = rep(1, n),heights = rep(1, n), nrow = n, ncol = n)
tmp = data.frame(so@meta.data, dim1 = so@reductions$umap_adt@cell.embeddings[,1], dim2 =so@reductions$umap_adt@cell.embeddings[,2], gate = colnames(so) %in% gate) 
p = ggplot(tmp) + geom_point_rast(aes(dim1, dim2), size = I(1), alpha = I(0.5), color = "black", raster.dpi = 100) +
  geom_point(data = tmp[tmp$gate == T, ], aes(dim1, dim2), color = "red", size = 2) +
  theme_bw() + theme(axis.text.x  = element_text(size=15), axis.text.y  = element_text(size=15), legend.text=element_text(size=20), axis.title.x = element_text(size=20) , axis.title.y = element_text(size=20), legend.title = element_blank()) + ggtitle("Gated cells")
print(p)
p = ggplot(tmp) + geom_point(aes(dim1, dim2, color = singler_immgen), size = I(1), alpha = I(0.5)) +
  scale_color_manual(values = mypal) +
  theme_bw() + theme(axis.text.x  = element_text(size=15), axis.text.y  = element_text(size=15), legend.text=element_text(size=50), axis.title.x = element_text(size=20) , axis.title.y = element_text(size=20), legend.title = element_blank()) + ggtitle("singler_immgen")
print(LabelClusters(p, id = "singler_immgen", size = 20))
dev.off()

ncells_prefiltering = ncol(so)
tmp = data.frame(so@meta.data)
tmp = tmp %>% group_by(HTO_classification) %>% summarize(ncells = n(),  mean_nCount_RNA = round(mean(nCount_RNA)), mean_nFeature_RNA = round(mean(nFeature_RNA)), mean_nCount_HTO = round(mean(nCount_HTO)), mean_nCount_ADT = round(mean(nCount_ADT)),  ncells_outliers_lowCountADT = length(which(lowCountADT)), ncells_outliers_autofluorescence = length(which(autofluo))) %>% as.data.frame()
write.table(x = tmp, file = "seuratobject_withHTOADT_singlet_postRNAfiltering_withSingleR_preADTfiltering.txt", quote = F, sep = "\t", row.names = F, col.names = T)

so@assays$qc_stats_4_preADTfiltering = tmp
saveRDS(so, file = "seuratobject_withHTOADT_singlet_postRNAfiltering_withSingleR_preADTfiltering.Rds")

##Postfiltering analysis

message("UMAP and clustering postfiltering...")
so = so[,!(so$autofluo | so$lowCountADT)]
so = Seurat_adt_umap_clustering(seurat_object = so, VariableFeatures =  rownames(so[["ADT"]])[!grepl(rownames(so[["ADT"]]), pattern = "unmapped|Isotype|TCRV")],reduction.name.pca = reduction.name.pca, reduction.name.umap = reduction.name.umap)

message("Find ADT cutoffs postfiltering...")
adt_thrs = FindThresholdADT()
adt_thrs$cutoff[match(rownames(so@assays$ADT@meta.features), rownames(adt_thrs))]
so@assays$ADT@meta.features$cutoff = adt_thrs$cutoff[match(rownames(so@assays$ADT@meta.features), rownames(adt_thrs))]
so@assays$ADT@meta.features$cutoff_crl_norm = adt_thrs$cutoff_crl_norm[match(rownames(so@assays$ADT@meta.features), rownames(adt_thrs))]
pdf("adt_qc_plot5_adt_thrs_postfiltering.pdf", width = 5, height = 5, useDingbats = F)
for (adt in sort(rownames(so@assays$ADT@counts))) {
  print(adt)
  x_adt = "TCRB"
  x = as.numeric(so@assays$ADT@data[x_adt,])
  #shist(x, unit = 0.5)
  y = as.numeric(so@assays$ADT@data[adt,])
  #shist(y, unit = 0.5, main = adt)
  #abline(v = adt_thrs[adt,"cutoff_crl_norm"], col = "red")
  tmp = data.frame(x, y)
  p = ggplot(tmp) + geom_point_rast(aes(x, y), size = 0.5, alpha = 0.5, raster.dpi = 100) +
    geom_hline(yintercept = so@assays$ADT@meta.features[adt,"cutoff_crl_norm"], color = "red", linetype = "dashed") + 
    geom_vline(xintercept = so@assays$ADT@meta.features[x_adt,"cutoff_crl_norm"], color = "red", linetype = "dashed") + 
    xlab(sprintf("%s (CRL)", x_adt)) +
    ylab(sprintf("%s (CRL)", adt)) +
    theme_classic()
  q = ggMarginal(p, type = "density")
  print(q,  newpage = TRUE) 
  
  tmp = data.frame(so@meta.data, dim1 = so@reductions[["umap_adt"]]@cell.embeddings[,1], dim2 =so@reductions[["umap_adt"]]@cell.embeddings[,2], sub_cells =  as.numeric(so@assays$ADT@data[adt,]) > so@assays$ADT@meta.features[adt,"cutoff_crl_norm"], size = so@assays$ADT@data[adt,]) 
  p = ggplot(tmp) +  geom_point_rast(aes(dim1, dim2, color = size), size = 0.5, alpha = 0.5, raster.dpi = 100) + scale_color_gradientn(colours = ColorRamp) + 
    scale_fill_gradientn(colours = ColorRamp) + theme(axis.text.x  = element_text(size=15), axis.text.y  = element_text(size=15), legend.text=element_text(size=20), axis.title.x = element_text(size=20) , axis.title.y = element_text(size=20), legend.title = element_blank()) + ggtitle(sprintf("%s, %s, cutoff %s","UMAP ADT", adt, so@assays$ADT@meta.features[adt,"cutoff_crl_norm"]))
  print(p)
  
  tmp = data.frame(so@meta.data, dim1 = so@reductions[["umap_rna"]]@cell.embeddings[,1], dim2 =so@reductions[["umap_rna"]]@cell.embeddings[,2], sub_cells =  as.numeric(so@assays$ADT@data[adt,]) > so@assays$ADT@meta.features[adt,"cutoff_crl_norm"], size = so@assays$ADT@data[adt,]) 
  p =  ggplot(tmp) +  geom_point_rast(aes(dim1, dim2, color = size), size = 0.5, alpha = 0.5, raster.dpi = 100) + scale_color_gradientn(colours = ColorRamp) + 
    scale_fill_gradientn(colours = ColorRamp) + theme( axis.text.x  = element_text(size=15), axis.text.y  = element_text(size=15), legend.text=element_text(size=20), axis.title.x = element_text(size=20) , axis.title.y = element_text(size=20), legend.title = element_blank()) + ggtitle(sprintf("%s, %s, cutoff %s","UMAP RNA", adt, so@assays$ADT@meta.features[adt,"cutoff_crl_norm"]))
  print(p)
}
dev.off()

pdf("adt_qc_plot6_UMAP_clustering_postfiltering.pdf", 15, 15, useDingbats=FALSE) 
tmp = data.frame(so@meta.data, dim1 = so@reductions[[reduction.name.umap]]@cell.embeddings[,1], dim2 =so@reductions[[reduction.name.umap]]@cell.embeddings[,2] ) 
p = ggplot(tmp) + geom_point(aes(dim1, dim2, color = ADT_snn_res.1), size = I(1), alpha = I(1)) + scale_color_manual(values = mypal) + theme_bw() + theme(axis.text.x  = element_text(size=15), axis.text.y  = element_text(size=15), legend.text=element_text(size=20), axis.title.x = element_text(size=20) , axis.title.y = element_text(size=20), legend.title = element_blank()) + ggtitle(sprintf("%s, %s",reduction.name.umap, "ADT_snn_res.1"))
LabelClusters(p, id = "ADT_snn_res.1", size = 5)
print(p + facet_wrap(~HTO_classification))

p = ggplot(tmp) + geom_point(aes(dim1, dim2, color = singler_immgen), size = I(1), alpha = I(1)) + scale_color_manual(values = mypal) + theme_bw() + theme(axis.text.x  = element_text(size=15), axis.text.y  = element_text(size=15), legend.text=element_text(size=20), axis.title.x = element_text(size=20) , axis.title.y = element_text(size=20), legend.title = element_blank()) + ggtitle(sprintf("%s, %s",reduction.name.umap, "singler_immgen"))
LabelClusters(p, id = "singler_immgen", size = 5)
print(p + facet_wrap(~HTO_classification))
dev.off()

ncells_postfiltering = ncol(so)
tmp = data.frame(so@meta.data)
tmp = tmp %>% group_by(HTO_classification) %>% summarize(ncells = n(),  mean_nCount_RNA = round(mean(nCount_RNA)), mean_nFeature_RNA = round(mean(nFeature_RNA)), mean_nCount_HTO = round(mean(nCount_HTO)), mean_nCount_ADT = round(mean(nCount_ADT))) %>% as.data.frame()
write.table(x = tmp, file = "seuratobject_withHTOADT_singlet_postRNAfiltering_withSingleR_postADTfiltering.txt", quote = F, sep = "\t", row.names = F, col.names = T)

so@assays$qc_stats_4_postADTfiltering = tmp
saveRDS(so, file = "seuratobject_withHTOADT_singlet_postRNAfiltering_withSingleR_postADTfiltering.Rds")

#DotPlot
tmp = data.frame(t(as.matrix(so@assays$ADT@data)), cluster = so$ADT_snn_res.1)
pdf(file = "adt_qc_plot7_Dotplot_ADT_ADTclusters.pdf", 25, 7, useDingbats=FALSE)
DotPlotHeatmap(data = tmp, title = "ADT_snn_res.1")
dev.off()

tmp = data.frame(t(as.matrix(so@assays$ADT@data)), cluster = so$RNA_snn_res.1)
pdf(file = "adt_qc_plot8_Dotplot_ADT_RNAclusters.pdf", 25, 7, useDingbats=FALSE)
DotPlotHeatmap(data = tmp, title = "RNA_snn_res.1")
dev.off()

tmp = data.frame(t(as.matrix(so@assays$ADT@data)), cluster = so$singler_immgen)
pdf(file = "adt_qc_plot9_Dotplot_ADT_SingleRImmgen.pdf", 25, 7, useDingbats=FALSE)
DotPlotHeatmap(data = tmp, title = "singler_immgen")
dev.off()
