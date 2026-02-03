args = commandArgs(TRUE)
sc_object = args[1]
path_to_sample_names = args[2]
IGT = args[3]

libs = c("Seurat", "ggplot2","inflection","grid","gridExtra","UCell")
sapply(libs, function(x) suppressMessages(suppressWarnings(library(x, character.only = TRUE, quietly = T, warn.conflicts  = F))))



sc<- readRDS(sc_object)
hashing<- read.csv(paste0(path_to_sample_names,"/samples.csv"),header=0)
hashing <- hashing[(1:18),]
hashing <- na.omit(hashing)
print(hashing[1])
print(hashing[2])
for (n in 1:length(sc@meta.data$HTO_classification.simplified)){
  for (l in 1:(dim(hashing)[1])){
  if (sc@meta.data$HTO_classification.simplified[n] == as.character(hashing[l,1])) {
    sc@meta.data$sample_name[n]<- as.character(hashing[l,2])
  }
  }
}

#sc@meta.data$nFeature_ADT<-0
#sc@meta.data$nCount_ADT<-0
sc_IGT_cell_metrics<- as.data.frame(cbind(sc@meta.data$nCount_RNA,sc@meta.data$nFeature_RNA)) #,sc@meta.data$nCount_ADT,sc@meta.data$nFeature_ADT))
colnames(sc_IGT_cell_metrics) <- c('nCount_RNA', 'nFeature_RNA') #,'nCount_ADT','nFeature_ADT')
write.csv(sc_IGT_cell_metrics,"cell_metrics.csv")
  
saveRDS(sc, file = "seuratobject_withHTOADT_singlet_postRNAfiltering_withSingleR.Rds")
