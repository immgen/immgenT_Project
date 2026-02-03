args = commandArgs(TRUE)
sc_IGT = args[1] 
cutoffs = args[2]
IGT<- args[3]  # Manully set your ncount_ADT cutoff for IGT dataset
spleen<- args[4] #hashtag for spleen control to remove from t cells


libs = c("Seurat", "ggplot2","inflection","grid","gridExtra","UCell","Seurat", "ggplot2", "reshape2", "dplyr", "fitdistrplus", "ggExtra", "ggrastr", "gridExtra", "gatepoints", "rafalib", "scales")
sapply(libs, function(x) suppressMessages(suppressWarnings(library(x, character.only = TRUE, quietly = T, warn.conflicts  = F))))


# Load in seurat_object
sc_IGT = readRDS(sc_IGT)


IGT<- read.table(IGT)
IGT<- as.character(IGT)

#spleen<- read.table(spleen)
#spleen<- as.character(spleen)

spleen<- read.table('spleen.txt')
rownames(spleen) <- spleen[,1]

cutoffs<- read.csv('T_cutoffs.csv', header = F)





#________________________

cutoff_T<- cutoffs[1,2]
cutoff_B<- cutoffs[2,2]
cutoff_ILC<- cutoffs[3,2]
cutoff_MNP<-cutoffs[4,2]

#sc_IGT@assays$ADT@counts@x <- 0
subset_keep<-sc_IGT

sc_IGTM<- as.data.frame(cbind(sc_IGT$nCount_RNA))
colnames(sc_IGTM) <- c('nCount_RNA')


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

pdf(file = paste0('T_labeled_RNA_CellScoring.pdf'), width = 19, height = 20)
# Create layout : nrow = 3, ncol = 2
pushViewport(viewport(layout = grid.layout(nrow = 4, ncol = 3)))
define_region <- function(row, col){
  viewport(layout.pos.row = row, layout.pos.col = col)
} 

print(p1_T, vp = define_region(row = 1, col = 1))
print(p2_T, vp = define_region(row = 1, col = 2))
print(p3_T, vp = define_region(row = 1, col = 3))


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
removed_cells2<- 0
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
tmp = tmp %>% group_by(HTO_classification.simplified) %>% summarize(ncells_postqc = n(),  mean_nCount_RNA_postqc = round(mean(nCount_RNA)), mean_nFeature_RNA_postqc = round(mean(nFeature_RNA)), mean_nCount_HTO_postqc = round(mean(nCount_HTO))) %>% as.data.frame() # , mean_nCount_ADT_postqc = round(mean(nCount_ADT))
tmp2 = data.frame(nonTcells@meta.data)
tmp2 = tmp2 %>% group_by(HTO_classification.simplified) %>% summarize(ncells_nonTcells = n()) %>% as.data.frame()

missing<-tmp[,1][!tmp[,1] %in% (tmp2[,1])]
if (length(missing) > 0 ){
missing<-as.data.frame(cbind(missing,0))
colnames(missing)<- colnames(tmp2)
tmp2<-rbind(tmp2,missing)
}


#tmp2$adt_secondary_filter_cells_removed <- as.vector(adt_secondary_filter_cells_removed)
#final_table<-cbind(tmp2,tmp[,-1])
final_table<-merge(tmp2,tmp)
write.csv(final_table,'seuratobject_withHTOADT_singlet_postRNAfiltering_postADTfiltering_postTfiltering.csv')
IGT_final_table<-t(as.data.frame(c(round(mean(Tcells$nCount_RNA)),round(mean(Tcells$nFeature_RNA)),round(mean(Tcells$nCount_HTO))))) #,round(mean(Tcells$nCount_ADT)))))
colnames(IGT_final_table)<-c('mean_nCount_RNA_postqc','mean_nFeature_RNA_postqc','mean_nCount_HTO_postqc') #,'mean_nCount_ADT_postqc')
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



Idents(Final_Tcells)<- 'sample_name'
levels(Final_Tcells)
saveRDS(Final_Tcells,"seuratobject_withHTOADT_singlet_postRNAfiltering_withSingleR_postADTfiltering_postTfiltering.Rds")
saveRDS(sc_IGT,"seuratobject_withHTOADT_singlet_postRNAfiltering_withSingleR_postADTfiltering.Rds")



