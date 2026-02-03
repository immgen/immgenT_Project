args = commandArgs(TRUE)
seuratobject_withHTOADT_singlet_preRNAfiltering_QC_stats = args[1]
seuratobject_withHTOADT_singlet_postRNAfiltering_withSingleR_preADTfiltering = args[2]
seuratobject_withHTOADT_singlet_postRNAfiltering_postADTfiltering_postTfiltering = args[3]
seuratobject_IGT_singlet_postRNAfiltering_postADTfiltering_postTfiltering<- args[4]  # Manully set your ncount_ADT cutoff for IGT dataset

libs = c("Seurat", "ggplot2","inflection","grid","gridExtra","UCell")
sapply(libs, function(x) suppressMessages(suppressWarnings(library(x, character.only = TRUE, quietly = T, warn.conflicts  = F))))


RNAfiltering<-read.table(seuratobject_withHTOADT_singlet_preRNAfiltering_QC_stats,header = T)
ADTfiltering<- read.table(seuratobject_withHTOADT_singlet_postRNAfiltering_withSingleR_preADTfiltering,header = T)
tmp<-cbind(RNAfiltering,ADTfiltering)
keeps<- c('HTO_classification.simplified','ncells','ncells_outliers_nGenes','ncells_outliers_deadcells','ncells_outliers_lowCountADT','ncells_outliers_autofluorescence')
tmp = tmp[keeps]
table2<-read.csv(seuratobject_withHTOADT_singlet_postRNAfiltering_postADTfiltering_postTfiltering,header = T)
table2<-table2[,-1]
#final_table<-cbind(tmp,table2)
final_table<-merge(tmp,table2)
write.csv(final_table,'seuratobject_withHTOADT_singlet_postRNAfiltering_postADTfiltering_postTfiltering_FinalTable.csv',row.names = F)
table3<-read.csv(seuratobject_IGT_singlet_postRNAfiltering_postADTfiltering_postTfiltering,header = T, row.names = 1)
final_table_IGT<-as.data.frame(t(colSums(final_table[,-c(1,9,10,11,12,13)])))
rownames(final_table_IGT)<-rownames(table3)
final_table_IGT<-cbind(final_table_IGT,table3)
write.csv(final_table_IGT,'seuratobject_IGT_singlet_postRNAfiltering_postADTfiltering_postTfiltering_FinalTable.csv',row.names = T)
