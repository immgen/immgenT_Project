# David Zemmour
# R/4.0.1
# usage: make_seurat_rna_adt_hto.R [path to project_description.csv] [path to RNA matrix.mtx.gz, genes.tsv.gz, barcodes.tsv.gz] [path to HTO+ADT matrix.mtx.gz, genes.tsv.gz, barcodes.tsv.gz] [hash or adt rows to remove from the analysis]

args = commandArgs(TRUE)
#PROJECT = args[1] # PROJECT = "project_description.csv" setwd("~/googledrive/ImmgenT/20220321_P01_mammarygland/")
path_to_rna = args[1] #path to RNA
path_to_hto_adt = args[2] #path_to_hto_adt 
path_to_sample_names = args[3] #path_to_hto_adt_remove_from_analysis = "hash_adt_rows_to_remove_from_analysis.txt"
spleen<- args[4] #hashtag for spleen control to remove from t cells

#project_description=read.csv(PROJECT,header=T,as.is=T)
#print(project_description)

print('hi')
spleen<- read.table('spleen.txt')
rownames(spleen) <- spleen[,1]
#spleen<- table(c('Negative'))
#rownames(spleen)<-'Negative'
print('bye')

libs = c("Seurat", "ggplot2", "gridExtra", "dplyr") 
sapply(libs, function(x) suppressMessages(suppressWarnings(library(x, character.only = TRUE, quietly = T, warn.conflicts  = F))))

library(RColorBrewer)
n = 70
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
mypal = unique(unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))))

message("Load RNA matrix")
umis = Read10X(data.dir = 'GEX/filtered_feature_bc_matrix', gene.column = 2)
cell_barcode = gsub(pattern = "\\-1", replacement = "", colnames(umis))
colnames(umis) = cell_barcode

message("Load HTO ADT matrix")
htos_adt = Read10X(data.dir = 'FBC/filtered_feature_bc_matrix', gene.column = 1)
cell_barcode = gsub(pattern = "\\-1", replacement = "", colnames(htos_adt))
colnames(htos_adt) = cell_barcode

message("Create seurat object")
length(intersect(colnames(htos_adt), colnames(umis)))
joint.bcs = intersect(colnames(htos_adt), colnames(umis))

message(sprintf("Joint RNA-HT/ADT barcodes %s", length(intersect(colnames(htos_adt), colnames(umis)))))
umis = umis[,joint.bcs]
htos = htos_adt[rownames(htos_adt)[grepl("Hash|HT", rownames(htos_adt))],joint.bcs]
adt = htos_adt[rownames(htos_adt)[!grepl("Hash|HT|unmapped", rownames(htos_adt))],joint.bcs]
message("These are the hashes: ")
print(htos)



hto_table<-data.frame(matrix(ncol = 1, nrow =nrow(htos) ))
for (n in 1:nrow(htos)) {
  hto_table[n,]<- sum(htos[n,])

}
rownames(hto_table)<-rownames(htos)
colnames(hto_table)<- 'hto_summed_counts'
hto_table<-t(hto_table)
write.csv(hto_table,'hto_counts_table.csv')


hashing<- read.csv(paste0(path_to_sample_names,"/samples.csv"),header=1)
hashing <- hashing[(1:28),]
hashing <- na.omit(hashing)

print(hashing[1])
rownames(hashing)<-hashing[,1]
htos <- htos[rownames(htos) %in% rownames(hashing),]

so = CreateSeuratObject(counts = umis,strip.suffix = T)
so$cell_barcode = colnames(so)
so$Cell_ID = paste("cell",1:length(colnames(so)),sep="_")

#so@assays$Experiment_details = as.character(project_description[1,])
#names(so@assays$Experiment_details) = colnames(project_description)

#htos <- htos[c("HT1","HT2","HT3","HT4","HT5","HT6","HT7","HT8", "HT9", "HT10"),]
so[["HTO"]] = CreateAssayObject(counts = htos + 1)
so[["ADT"]] = CreateAssayObject(counts = adt)

message("Print rownames and colnames so@RNA and so$ADT")
head(rownames(so@assays$RNA@counts), 50)
# head(rownames(so@assays$ADT@counts), 200)
head(colnames(so), 50)

message(sprintf("Cells with no Hash counts %s", names(which(colSums(htos) ==0))))
print(names(which(colSums(htos) ==0)))
write.table(x = names(which(colSums(htos) ==0)), file = "cells_with_no_hashtag_counts.txt",quote = F,row.names = F, col.names = F)
so = so[,colSums(htos) != 0]




message("Demultiplex cells based on HTO")
so = NormalizeData(so, assay = "HTO", normalization.method = "CLR")
so = HTODemux(so, assay = "HTO", positive.quantile = 0.999, verbose = T)

so$HTO_classification_orig=so$HTO_classification

for(i in which(so$HTO_classification=="Negative")){
  print(i)
  max_hash<-max(so@assays$HTO@counts[,i])
  max2_hash<-max(so@assays$HTO@counts[,i][-which.max(so@assays$HTO@counts[,i])])
  if(max_hash>2*max2_hash & max_hash>10){
    so$HTO_classification[i]<-names(which.max(so@assays$HTO@counts[,i]))
  }
}

png("make_seurat_rna_adt_hto-1.png",width=1400, height=1600)
RidgePlot(so, assay = "HTO", features = rownames(so[["HTO"]]), ncol = 2)
dev.off()


png("make_seurat_rna_adt_hto-2.png",width=1400, height=1600)

panel.cor <- function(x, y, digits = 2, cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  # correlation coefficient
  r <- cor(x, y)
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste("r= ", txt, sep = "")
  text(0.5, 0.6, txt)

  # p-value calculation
  p <- cor.test(x, y)$p.value
  txt2 <- format(c(p, 0.123456789), digits = digits)[1]
  txt2 <- paste("p= ", txt2, sep = "")
  if(p<0.01) txt2 <- paste("p= ", "<0.01", sep = "")
  text(0.5, 0.4, txt2)
}

pairs(t(so@assays$HTO@data), cex = 4, upper.panel = panel.cor, pch = ".")
dev.off()

png("make_seurat_rna_adt_hto-3.png",width=1400, height=1200)




panel.cor <- function(x, y, digits = 2, cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  # correlation coefficient
  r <- cor(x, y)
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste("r= ", txt, sep = "")
  text(0.5, 0.6, txt)

  # p-value calculation
  p <- cor.test(x, y)$p.value
  txt2 <- format(c(p, 0.123456789), digits = digits)[1]
  txt2 <- paste("p= ", txt2, sep = "")
  if(p<0.01) txt2 <- paste("p= ", "<0.01", sep = "")
  text(0.5, 0.4, txt2)
}


pairs(t(so@assays$HTO@data), cex = 4, upper.panel = panel.cor, pch = ".")

so = ScaleData(so, assay = "HTO", features = rownames(so@assays$HTO), verbose = FALSE)
so = RunPCA(so, assay ="HTO", features = rownames(so@assays$HTO), reduction.name = "pca_htos", approx = FALSE)
so = RunTSNE(so, assay ="HTO", reduction = "pca_htos", reduction.name = "tsne_htos", dims = 1:ncol(so@reductions$pca_htos@cell.embeddings), perplexity = 100, verbose = T, check_duplicates = FALSE, max_iter = 500)

so$HTO_classification.simplified = as.character(so$HTO_classification)
so$HTO_classification.simplified[so$HTO_classification.global == "Doublet"] = "Doublet"
tmp = data.frame(so@meta.data)
tmp$dim1 = so@reductions$tsne_htos@cell.embeddings[,1]
tmp$dim2 = so@reductions$tsne_htos@cell.embeddings[,2]
#tmp$HTO_classification.simplified = as.character(tmp$HTO_classification)
#tmp$HTO_classification.simplified[tmp$HTO_classification.global == "Doublet"] = "Doublet"

#Adding Rescured category in HTO_classification.global
so$HTO_classification.global[so$HTO_classification.global == 'Negative' & so$HTO_classification.simplified != 'Negative'] = "Rescued"

theme_my = theme(axis.text.x=element_text(size=15), axis.text.y=element_text(size=15), legend.text=element_text(size=10), axis.title.x = element_text(size=20) , axis.title.y = element_text(size=20), legend.title = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 

p = ggplot(tmp) + geom_point(aes(dim1, dim2, color = HTO_classification.simplified), alpha = I(1), size = 1) + scale_color_manual(values = mypal) + theme_bw() + theme_my + facet_wrap(~HTO_classification.global, nrow = 2)

q = HTOHeatmap(so, assay = "HTO", ncells = 5000)


grid.arrange(
  grobs = list(p,q),
  widths = c(1,1),
  heights = c(1,1,1),
  nrow = 3,
  ncol = 2
  #layout_matrix = rbind(c(1,0),
  #                      c(1,0))
  )

dev.off()

message("Saving seuratobject_withHTOADT_alldata.Rda and seuratobject_withHTOADT_singlet.Rda")

tmp = data.frame(so@meta.data)
# tmp = tmp %>% group_by(HTO_classification.simplified) %>% summarize(ncells = n(),  mean_nCount_RNA = round(mean(nCount_RNA)), mean_nFeature_RNA = round(mean(nFeature_RNA)), mean_nCount_HTO = round(mean(nCount_HTO)), mean_nCount_ADT = round(mean(nCount_ADT))) %>% as.data.frame()
tmp = tmp %>% group_by(HTO_classification.simplified) %>% summarize(ncells = n(),  mean_nCount_RNA = round(mean(nCount_RNA)), mean_nFeature_RNA = round(mean(nFeature_RNA)), mean_nCount_HTO = round(mean(nCount_HTO))) %>% as.data.frame()
write.table(x = tmp, file = "seuratobject_withHTOADT_alldata_QCstats.txt", quote = F, sep = "\t", row.names = F, col.names = T)
so@assays$qc_stats_1 = tmp

#tmp = so@assays$Experiment_details["Hashtag"]
#tmp = unlist(strsplit(x = tmp, split= ","))
#tmp = unlist(strsplit(x = tmp, split= ";"))
print('hiii')
hts = unique(so$HTO_classification.simplified)
hts = hts[!hts %in% c("Doublet", "Negative")]
so$sample_name = NA
print('bywwww')
# for (i in hts) {
#   so$sample_name[so$HTO_classification.simplified == i] = tmp[grep(pattern = paste(i,"[^0-9]", sep = ""), tmp)]
# }



#so$HTO_classification.global <- factor(so$HTO_classification.global, levels = c('Singlet', 'Doublet','Negative','Rescued'))
Idents(so)<- 'HTO_classification.global' 
table<- t((table(Idents(so))))
#colnames(table)<- c('nSinglets','nDoublets','nNegatives','nRecovered')
write.table(x = table, file = "Demultiplexing_QC_stats.txt", quote = F, sep = "\t", row.names = F, col.names = T)


saveRDS(so, file = "seuratobject_withHTOADT_alldata.Rds")

so = so[,so$HTO_classification.global %in% c("Singlet", "Rescued")]


# print("Debugging")
# print(Idents(so))

Idents(so) <- 'HTO_classification'


saveRDS(so, file = "seuratobject_withHTOADT_singlet.Rds")
tmp = data.frame(so@meta.data)
# tmp = tmp %>% group_by(HTO_classification.simplified) %>% summarize(ncells = n(),  mean_nCount_RNA = round(mean(nCount_RNA)), mean_nFeature_RNA = round(mean(nFeature_RNA)), mean_nCount_HTO = round(mean(nCount_HTO)), mean_nCount_ADT = round(mean(nCount_ADT))) %>% as.data.frame()
tmp = tmp %>% group_by(HTO_classification.simplified) %>% summarize(ncells = n(),  mean_nCount_RNA = round(mean(nCount_RNA)), mean_nFeature_RNA = round(mean(nFeature_RNA)), mean_nCount_HTO = round(mean(nCount_HTO))) %>% as.data.frame()
write.table(x = tmp, file = "seuratobject_withHTOADT_singlet_QCstats.txt", quote = F, sep = "\t", row.names = F, col.names = T)
