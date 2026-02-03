args = commandArgs(TRUE)
sc_IGT = args[1] 
cutoffs = args[2]
IGT<- args[3]  # Manully set your ncount_ADT cutoff for IGT dataset
spleen<- args[4] #hashtag for spleen control to remove from t cells


libs = c("Seurat", "ggplot2","inflection","grid","gridExtra","UCell","Seurat", "ggplot2", "reshape2", "dplyr", "fitdistrplus", "ggExtra", "ggrastr", "gridExtra", "gatepoints", "rafalib", "scales")
sapply(libs, function(x) suppressMessages(suppressWarnings(library(x, character.only = TRUE, quietly = T, warn.conflicts  = F))))


# Load in seurat_object
sc_IGT = readRDS(sc_IGT)


IGT<- read.table('IGT.txt')
IGT<- as.character(IGT)

print('hi')
spleen<- read.table('spleen.txt')
rownames(spleen) <- spleen[,1]
#spleen<- table(c('Negative'))
#rownames(spleen)<-'Negative'
print('bye')

cutoffs<- read.csv('T_cutoffs.csv', header = F)

DotPlotHeatmap = function(data = tmp, title = "ADT_seurat_clusters") {
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



#________________________

cutoff_T<- cutoffs[1,2]
cutoff_B<- cutoffs[2,2]
cutoff_ILC<- cutoffs[3,2]
cutoff_MNP<-cutoffs[4,2]
ncountADT_cutoff<- cutoffs[5,2]

subset_keep<-sc_IGT
sc_IGTM<- as.data.frame(cbind(sc_IGT$nCount_RNA,sc_IGT$nCount_ADT))
colnames(sc_IGTM) <- c('nCount_RNA', 'nCount_ADT')


percentages<- rbind(subset_keep$signature_1T.enriched,
                    subset_keep$signature_1mnp.enriched,
                    subset_keep$signature_1B.enriched,
                    subset_keep$signature_1ILC.enriched)
row.names(percentages) <- c("T Percentages","MNP Percentages","B Percentages","ILC Percentages")

percentages<-as.data.frame(t(percentages))


p1_T<- ggplot()+ 
  geom_point(data=percentages ,aes(`T Percentages`, `B Percentages`),color='blue' )+ 
  xlim(0, max(percentages[,1])) + ylim(0, max(percentages[,3])) +theme(text = element_text(size=10))+
  ggtitle("ADT T cells")+ scale_y_continuous(limits = c(0, max(percentages[,3])),breaks = (seq(0,max(percentages[,3]), by = 0.01)))+
  geom_vline(xintercept=cutoff_T, linetype="dashed", color = "red")+ geom_hline(yintercept=cutoff_B, linetype="dashed", color = "red")

p2_T<-ggplot()+ 
  geom_point(data=percentages ,aes(`T Percentages`, `ILC Percentages`),color='blue' )+
  xlim(0, max(percentages[,1])) + ylim(0, max(percentages[,3])) +theme(text = element_text(size=10))+
  ggtitle("ADT T cells")+ scale_y_continuous(limits = c(0, max(percentages[,3])),breaks = (seq(0,max(percentages[,3]), by = 0.01)))+
  geom_vline(xintercept=cutoff_T, linetype="dashed", color = "red")+ geom_hline(yintercept=cutoff_ILC, linetype="dashed", color = "red")

p3_T<-ggplot()+ 
  geom_point(data=percentages ,aes(`T Percentages`, `MNP Percentages`),color='blue' )+
  xlim(0, max(percentages[,1])) +theme(text = element_text(size=10))+
  ggtitle("ADT T cells")+ scale_y_continuous(limits = c(0, max(percentages[,3])),breaks = (seq(0,max(percentages[,3]), by = 0.01)))+
  geom_vline(xintercept=cutoff_T, linetype="dashed", color = "red")+ geom_hline(yintercept=cutoff_MNP, linetype="dashed", color = "red")

png(file = paste0('T_labeled_RNA_CellScoring-1.png'),width=1800, height=1050)
# Create layout : nrow = 3, ncol = 2
pushViewport(viewport(layout = grid.layout(nrow = 4, ncol = 3)))
define_region <- function(row, col){
  viewport(layout.pos.row = row, layout.pos.col = col)
} 

print(p1_T, vp = define_region(row = 1, col = 1))
print(p2_T, vp = define_region(row = 1, col = 2))
print(p3_T, vp = define_region(row = 1, col = 3))
dev.off()

png(file = paste0('T_labeled_RNA_CellScoring-2.png'),width = 1000,height = 1300)
ggplot(sc_IGTM, aes(nCount_RNA, nCount_ADT)) + geom_point(color='gray')+
  scale_y_continuous(breaks = round(seq(min(sc_IGTM[2]),max(sc_IGTM[2]), by = 250),1))+
   geom_hline(yintercept=ncountADT_cutoff, linetype="dashed", color = "red")



dev.off()


percentages<- rbind(subset_keep$signature_1T.enriched,
                    subset_keep$signature_1mnp.enriched,
                    subset_keep$signature_1B.enriched,
                    subset_keep$signature_1ILC.enriched)
row.names(percentages) <- c("T Percentages","MNP Percentages","B Percentages","ILC Percentages")
percentages_Tcells<- percentages[,percentages[1,]>cutoff_T]
percentages_Tcells<- percentages_Tcells[,percentages_Tcells[2,]<cutoff_MNP]
percentages_Tcells<- percentages_Tcells[,percentages_Tcells[3,]<cutoff_B]
percentages_Tcells<- percentages_Tcells[,percentages_Tcells[4,]<cutoff_ILC]



percentages_nonTcells<- percentages[,percentages[1,]<cutoff_T]
percentages_nonTcells<- percentages_nonTcells[,percentages_nonTcells[2,]<cutoff_MNP]
percentages_nonTcells<- percentages_nonTcells[,percentages_nonTcells[3,]<cutoff_B]
percentages_nonTcells<- percentages_nonTcells[,percentages_nonTcells[4,]<cutoff_ILC]

percentages_nonTcells2<- percentages[,percentages[1,]>cutoff_T]
percentages_nonTcells2<- percentages_nonTcells2[,percentages_nonTcells2[3,]>cutoff_B]







removed_cells1_preTCR<- colnames(percentages)[!colnames(percentages) %in% colnames(percentages_Tcells)]
removed_cells1<- removed_cells1_preTCR

subset_filter1_keptcells<-subset(sc_IGT, cells =removed_cells1, invert=T)
removed_cells2<- which(subset_filter1_keptcells$nCount_ADT<ncountADT_cutoff)
subset_filter2_keptcells<-subset(subset_filter1_keptcells, cells =removed_cells2, invert=T)


Tcells<- subset_filter2_keptcells
print("hi")
if (rownames(spleen) == 'Negative') {
Tcells<- Tcells
} else{
Idents(sc_IGT)<- 'HTO_classification.simplified'
Idents(Tcells)<- 'HTO_classification.simplified'
spleen_sc<- subset(x= sc_IGT, idents= rownames(spleen))
Tcells<- subset(x= Tcells, idents= rownames(spleen), invert = TRUE)
saveRDS(spleen_sc, file = "spleen.Rds")
}
print('bye')


nonTcells<-subset(sc_IGT, cells =colnames(subset_filter1_keptcells), invert=T)
nonT_cells_removed <-table(nonTcells@meta.data$HTO_classification.simplified, nonTcells@meta.data$orig.ident)
nonTcells<-subset(sc_IGT, cells =colnames(subset_filter2_keptcells), invert=T)
total_cells_removed <-table(nonTcells@meta.data$HTO_classification.simplified, nonTcells@meta.data$orig.ident)
#adt_secondary_filter_cells_removed<- total_cells_removed - nonT_cells_removed
ncells_final_Tcells <-table(Tcells@meta.data$HTO_classification.simplified, Tcells@meta.data$orig.ident)
ncells_final_Tcells <-table(Tcells@meta.data$HTO_classification.simplified, Tcells@meta.data$orig.ident)
tmp = data.frame(Tcells@meta.data)
tmp = tmp %>% group_by(HTO_classification.simplified) %>% summarize(ncells_postqc = n(),  mean_nCount_RNA_postqc = round(mean(nCount_RNA)), mean_nFeature_RNA_postqc = round(mean(nFeature_RNA)), mean_nCount_HTO_postqc = round(mean(nCount_HTO)), mean_nCount_ADT_postqc = round(mean(nCount_ADT))) %>% as.data.frame()
tmp2 = data.frame(nonTcells@meta.data)
tmp2 = tmp2 %>% group_by(HTO_classification.simplified) %>% summarize(ncells_nonTcells = n()) %>% as.data.frame()
#tmp2$adt_secondary_filter_cells_removed <- as.vector(adt_secondary_filter_cells_removed)
missing<-tmp[,1][!tmp[,1] %in% (tmp2[,1])]
if (length(missing) > 0 ){
missing<-as.data.frame(cbind(missing,0))
colnames(missing)<- colnames(tmp2)
tmp2<-rbind(tmp2,missing)
}

#final_table<-cbind(tmp2,tmp[,-1])
final_table<-merge(tmp2,tmp)
write.csv(final_table,'seuratobject_withHTOADT_singlet_postRNAfiltering_postADTfiltering_postTfiltering.csv')
IGT_final_table<-t(as.data.frame(c(round(mean(Tcells$nCount_RNA)),round(mean(Tcells$nFeature_RNA)),round(mean(Tcells$nCount_HTO)),round(mean(Tcells$nCount_ADT)))))
colnames(IGT_final_table)<-c('mean_nCount_RNA_postqc','mean_nFeature_RNA_postqc','mean_nCount_HTO_postqc','mean_nCount_ADT_postqc')
rownames(IGT_final_table)<-IGT
write.csv(IGT_final_table,'seuratobject_IGT_singlet_postRNAfiltering_postADTfiltering_postTfiltering.csv')

percentages_nonT<- rbind(nonTcells$signature_1T.enriched,
                    nonTcells$signature_1mnp.enriched,
                    nonTcells$signature_1B.enriched,
                    nonTcells$signature_1ILC.enriched)
row.names(percentages_nonT) <- c("T Percentages","MNP Percentages","B Percentages","ILC Percentages")
cell_index <- apply(percentages_nonT,2,which.max)
cell_index<- as.data.frame(cell_index)
cell_index[cell_index==2]<-"MNP"
cell_index[cell_index==3]<-"B"
cell_index[cell_index==4]<-"ILC"
nonTcells$cell_type <-cell_index 
Idents(nonTcells)<- 'cell_type'

sc_IGT$cell_type<- 'n/a'
sc_IGT$Tcell<- 'n/a'

sc_IGT$Tcell[colnames(sc_IGT) %in% colnames(Tcells)]<- 'Y'
sc_IGT$Tcell[colnames(sc_IGT) %in% colnames(nonTcells)]<- 'N'


Idents(nonTcells)<- 'cell_type'
sc_IGT$cell_type[colnames(sc_IGT) %in% colnames(Tcells)]<- 'T'

if ('B' %in% levels(nonTcells)){
sc_IGT$cell_type[colnames(sc_IGT) %in% WhichCells(nonTcells,idents = 'B')]<- 'B'}

if ('MNP' %in% levels(nonTcells)){
sc_IGT$cell_type[colnames(sc_IGT) %in% WhichCells(nonTcells,idents = 'MNP')]<- 'MNP'}

if ('ILC' %in% levels(nonTcells)){
sc_IGT$cell_type[colnames(sc_IGT) %in% WhichCells(nonTcells,idents = 'ILC')]<- 'ILC'}


sc_IGT$cell_type[colnames(sc_IGT) %in% colnames(percentages_nonTcells)]<- 'NA'
sc_IGT$cell_type[colnames(sc_IGT) %in% colnames(percentages_nonTcells2)]<- 'B-T like cells'







Idents(sc_IGT)<- 'cell_type'
levels(sc_IGT)

Idents(sc_IGT)<- 'Tcell'
levels(sc_IGT)


data <- as.data.frame(t(as.matrix(sc_IGT@assays$ADT@counts)))
data <- log2(data + 1)
data <- data[, order(colnames(data))]

data2<- data[rownames(data) %in% colnames(Tcells),]
data<- data[!rownames(data) %in% colnames(Tcells),]

proX<- "THY1.2"
proY1<- "CD3"

p1<- ggplot(data2, aes(data2[, proX], data2[, proY1])) + geom_point(color='gray') +
  geom_point(data=data ,aes(data[, proX], data[, proY1]),color='red' )+ 
  xlab(proX)  +
  ylab(proY1) + xlim(0, max(data2[, proX])) + ylim(0, max(data2[, proX])) + theme_bw() +theme(text = element_text(size=12))+
  ggtitle("RNA non-T cells")


proY2<- "CD2"

p2<- ggplot(data2, aes(data2[, proX], data2[, proY2])) + geom_point(color='gray') +
  geom_point(data=data ,aes(data[, proX], data[, proY2]),color='red' )+ 
  xlab(proX)  +
  ylab(proY2) + xlim(0, max(data2[, proX])) + ylim(0, max(data2[, proX])) + theme_bw() +theme(text = element_text(size=12))+
  ggtitle("RNA non-T cells")

proY4<- "TCRGD"

p4<-ggplot(data2, aes(data2[, proX], data2[, proY4])) + geom_point(color='gray') +
  geom_point(data=data ,aes(data[, proX], data[, proY4]),color='red' )+ 
  xlab(proX)  +
  ylab(proY4) + xlim(0, max(data2[, proX])) + ylim(0, max(data2[, proX])) + theme_bw() +theme(text = element_text(size=12))+
  ggtitle("RNA non-T cells")

proY5<- "TCRB"

p5<- ggplot(data2, aes(data2[, proX], data2[, proY5])) + geom_point(color='gray') +
  geom_point(data=data ,aes(data[, proX], data[, proY5]),color='red' )+ 
  xlab(proX) +
  ylab(proY5) + xlim(0, max(data2[, proX])) + ylim(0, max(data2[, proX])) + theme_bw() +theme(text = element_text(size=12))+
  ggtitle("RNA non-T cells")




proX<- "THY1.2"
proY1<- "CD3"

p1_T<- ggplot(data, aes(data[, proX], data[, proY1])) + geom_point(color='gray') +
  geom_point(data=data2 ,aes(data2[, proX], data2[, proY1]),color='red' )+ 
  xlab(proX)  +
  ylab(proY1) + xlim(0, max(data2[, proX])) + ylim(0, max(data2[, proX])) + theme_bw() +theme(text = element_text(size=12))+
  ggtitle("RNA T cells")


proY2<- "CD2"

p2_T<- ggplot(data, aes(data[, proX], data[, proY2])) + geom_point(color='gray') +
  geom_point(data=data2 ,aes(data2[, proX], data2[, proY2]),color='red' )+ 
  xlab(proX)  +
  ylab(proY2) + xlim(0, max(data2[, proX])) + ylim(0, max(data2[, proX])) + theme_bw() +theme(text = element_text(size=12))+
  ggtitle("RNA T cells")

proY4<- "TCRGD"

p4_T<-ggplot(data, aes(data[, proX], data[, proY4])) + geom_point(color='gray') +
  geom_point(data=data2 ,aes(data2[, proX], data2[, proY4]),color='red' )+ 
  xlab(proX)  +
  ylab(proY4) + xlim(0, max(data2[, proX])) + ylim(0, max(data2[, proX])) + theme_bw() +theme(text = element_text(size=12))+
  ggtitle("RNA T cells")

proY5<- "TCRB"

p5_T<- ggplot(data, aes(data[, proX], data[, proY5])) + geom_point(color='gray') +
  geom_point(data=data2 ,aes(data2[, proX], data2[, proY5]),color='red' )+ 
  xlab(proX)+
  ylab(proY5) + xlim(0, max(data2[, proX])) + ylim(0, max(data2[, proX])) + theme_bw() +theme(text = element_text(size=12))+
  ggtitle("RNA T cells")

table_column<-c('n_RNA_Tcells','n_RNA_nonTcells','n_low_ADT_T_cells (to remove)','n_final_Tcells')
table<-as.data.frame(t(c(dim(percentages_Tcells)[2],length(removed_cells1_preTCR),length(removed_cells2),dim(Tcells)[2])))
colnames(table)<-table_column 



png(file = paste0('Tcell_cleanup_ADT_post_filtering.png'),width=1000, height=1200)
pushViewport(viewport(layout = grid.layout(nrow = 4, ncol = 2)))
define_region <- function(row, col){
  viewport(layout.pos.row = row, layout.pos.col = col)
} 
print(p1, vp = define_region(row = 1, col = 1))   # Span over two columns
print(p2, vp = define_region(row = 1, col = 2))
print(p4, vp = define_region(row = 2, col = 1))
print(p5, vp = define_region(row = 2, col = 2))
print(p1_T, vp = define_region(row = 3, col = 1))   # Span over two columns
print(p2_T, vp = define_region(row = 3, col = 2))
print(p4_T, vp = define_region(row = 4, col = 1))
print(p5_T, vp = define_region(row = 4, col = 2))

dev.off()

pdf(file = paste0('Tcell_cleanup_final_numbers.pdf'), width = 10, height = 5,useDingbats=FALSE)
print(grid.table(table),vp = define_region(row = 4, col = 1))
dev.off()

#Idents(Tcells)<- 'hash.ID'

#Tcells<- subset(x= Tcells, idents= spleen, invert = TRUE)
Seurat_rna_umap_clustering = function(seurat_object = so, tsne.method = c("FIt-SNE", "Rtsne"), print_clusters = T, reduction.name.pca = "pca_rna", reduction.name.umap = "umap_rna") {
  require(Seurat)
  seurat_object = NormalizeData(seurat_object, verbose = T, normalization.method = "LogNormalize", scale.factor = 10000,assay = "RNA")
  seurat_object = FindVariableFeatures(seurat_object, selection.method = "vst", nfeatures = 2000, verbose = T,assay = "RNA")
  DefaultAssay(seurat_object)<- 'RNA'
  genes<- rownames(seurat_object)
  Trav<- genes[grepl('^Trav',genes)]
  Traj<-genes[grepl('^Traj',genes)]
  Trac<-genes[grepl('^Trac',genes)]
  Trbv<-genes[grepl('^Trbv',genes)]
  Trbd<-genes[grepl('^Trbd',genes)]
  Trbj<-genes[grepl('^Trbj',genes)]
  Trd<-genes[grepl('^Trd',genes)]
  Trg<-genes[grepl('^Trg',genes)]
  TCR_genes<- (c(Trav,Traj,Trac,Trbv,Trbd,Trbj,Trd,Trg))
  sex_specific_genes<- c('Ddx3y','Uty','Xist','Eif2s3y','Kdm5d','Tsix')
  remove_var_genes<- c(TCR_genes,sex_specific_genes)
  seurat_object@assays$RNA@var.features<- seurat_object@assays$RNA@var.features[!seurat_object@assays$RNA@var.features %in% remove_var_genes]
  
  seurat_object = ScaleData(seurat_object, assay = "RNA")
  
  seurat_object = RunPCA(seurat_object, npcs = 100, ndims.print = 1:5, nfeatures.print = 5, reduction.name = reduction.name.pca,assay = "RNA", var.features=seurat_object@assays$RNA@var.features)
  
  x = cumsum((seurat_object@reductions[[reduction.name.pca]]@stdev**2 / sum(seurat_object@reductions[[reduction.name.pca]]@stdev**2)))
  ndims = min(which(x >= 0.8))
  ElbowPlot(object = seurat_object, ndims = 100, reduction = reduction.name.pca) + geom_vline(mapping = aes(xintercept = ndims), color = "red", linetype="dashed") + ggtitle(label = sprintf("ndims PCA for further analysis = %s", ndims))
  DimHeatmap(seurat_object, dims = c(1:3, 70:75), cells = 500, balanced = TRUE, reduction = reduction.name.pca)
  
  #message("Running tSNE...")
  #seurat_object = RunTSNE(seurat_object, dims = 1:ndims, tsne.method = tsne.method, nthreads = 8, max_iter = 2000)
  message("Running UMAP")
  seurat_object = RunUMAP(seurat_object, dims = 1:ndims, reduction = reduction.name.pca, reduction.name = reduction.name.umap)
  
  message("RNA clusters...")
  seurat_object = FindNeighbors(seurat_object, reduction = reduction.name.pca, dims = 1:ndims, k.param = 20, verbose = T)
  
  seurat_object = FindClusters(seurat_object, resolution = 0.25, n.start = 10, algorithm = 1, n.iter = 10)
  seurat_object = FindClusters(seurat_object, resolution = 0.5, n.start = 10, algorithm = 1, n.iter = 10)
  seurat_object = FindClusters(seurat_object, resolution = 1, n.start = 10, algorithm = 1, n.iter = 10)
  seurat_object = FindClusters(seurat_object, resolution = 1.5, n.start = 10, algorithm = 1, n.iter = 10)
  seurat_object = FindClusters(seurat_object, resolution = 2, n.start = 10, algorithm = 1, n.iter = 10)
  seurat_object = FindClusters(seurat_object, resolution = 3, n.start = 10, algorithm = 1, n.iter = 10)
  seurat_object = FindClusters(seurat_object, resolution = 4, n.start = 10, algorithm = 1, n.iter = 10)
  #colnames(seurat_object@meta.data) = gsub(x = colnames(seurat_object@meta.data), pattern = "^(.*?)_snn_res.", replacement = "RNACluster_Res")
  
  if (print_clusters == T) {
    message("Printing clusters...")
    tmp = data.frame(seurat_object@meta.data, dim1 = seurat_object@reductions[[reduction.name.umap]]@cell.embeddings[,1], dim2 =seurat_object@reductions[[reduction.name.umap]]@cell.embeddings[,2] )
    for( i in grep("RNA_snn_res.", colnames(seurat_object@meta.data))) {
      p = ggplot(tmp) + geom_point(aes(dim1, dim2, color = tmp[,i]), size = I(1), alpha = I(1)) + theme_bw() + theme(axis.text.x  = element_text(size=15), axis.text.y  = element_text(size=15), legend.text=element_text(size=20), axis.title.x = element_text(size=20) , axis.title.y = element_text(size=20), legend.title = element_blank()) + ggtitle(sprintf("%s, %s",reduction.name.umap, colnames(seurat_object@meta.data)[i]))
      p2 = LabelClusters(p, id =  colnames(seurat_object@meta.data)[i], size = 5)
      print(p2)
    }
  }
  
  
  return(seurat_object)
}

DefaultAssay(Tcells)<-"RNA"
Final_Tcells= Seurat_rna_umap_clustering(seurat_object = Tcells, tsne.method = "FIt-SNE", print_clusters = T, reduction.name.pca = "pca_rna", reduction.name.umap = "umap_rna")
reduction.name.pca = "pca_adt"
reduction.name.umap = "umap_adt"                                
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

DefaultAssay(Final_Tcells)<-"ADT"
Final_Tcells= Seurat_adt_umap_clustering(seurat_object = Final_Tcells, VariableFeatures =  rownames(Final_Tcells[["ADT"]])[!grepl(rownames(Final_Tcells[["ADT"]]), pattern = "unmapped|Isotype|TCRV")],reduction.name.pca = reduction.name.pca, reduction.name.umap = reduction.name.umap)


Idents(Final_Tcells)<- 'sample_name'
levels(Final_Tcells)
saveRDS(Final_Tcells,"seuratobject_withHTOADT_singlet_postRNAfiltering_withSingleR_postADTfiltering_postTfiltering.Rds")
saveRDS(sc_IGT,"seuratobject_withHTOADT_singlet_postRNAfiltering_withSingleR_postADTfiltering.Rds")

so<- Final_Tcells
cluster <- Final_Tcells$seurat_clusters

tmp = data.frame(t(as.matrix(Final_Tcells@assays$ADT@data)), cluster = Final_Tcells$seurat_clusters)
png(file = "adt_qc_postTfilt__Dotplot_ADT_ADTclusters.png",width=1900, height=600)
DotPlotHeatmap(data = tmp, title = "ADT_seurat_clusters")
dev.off()


