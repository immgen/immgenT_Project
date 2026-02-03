args = commandArgs(TRUE)
sc_object = args[1]
sc_object2 = args[2] 
path_to_sample_names = args[3]
IGT = args[4]

libs = c("Seurat", "ggplot2","inflection","grid","gridExtra","UCell")
sapply(libs, function(x) suppressMessages(suppressWarnings(library(x, character.only = TRUE, quietly = T, warn.conflicts  = F))))



sc<- readRDS(sc_object)
hashing<- read.csv(paste0(path_to_sample_names,"/samples.csv"),header=0)
hashing <- hashing[(1:28),]
hashing <- na.omit(hashing)
print(hashing[1])
print(hashing[2])
for (n in 1:length(sc$HTO_classification.simplified)){
  for (l in 1:(dim(hashing)[1])){
  if (sc$HTO_classification.simplified[n] == as.character(hashing[l,1])) {
    sc$sample_name[n]<- as.character(hashing[l,2])
  }
  }
}

  
saveRDS(sc, file = "dataset.Rds")


sc<- readRDS(sc_object2)
hashing<- read.csv(paste0(path_to_sample_names,"/samples.csv"),header=0)
hashing <- hashing[(1:28),]
hashing <- na.omit(hashing)
print(hashing[1])
print(hashing[2])
for (n in 1:length(sc$HTO_classification.simplified)){
  for (l in 1:(dim(hashing)[1])){
  if (sc$HTO_classification.simplified[n] == as.character(hashing[l,1])) {
    sc$sample_name[n]<- as.character(hashing[l,2])
  }
  }
}


saveRDS(sc, file = "dataset_clean.Rds")
