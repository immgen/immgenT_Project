library(Seurat)
library(seqinr)

args = commandArgs(TRUE)

experiment_name <- args[1]
seurat_obj_name <- args[2]

experiment_name<- read.table(experiment_name)
experiment_name<- as.character(experiment_name)

fasta_name <- paste(experiment_name, "allCells_fullyproductiveTCR_trv_cdr3", sep = "_")

sc <- readRDS(seurat_obj_name)

tcr <- read.csv("TCR/outs/filtered_contig_annotations.csv", sep=",")

# Remove the -1 at the end of each barcode.

tcr$barcode <- gsub("-1", "", tcr$barcode)


tcr_a <- tcr[tcr$chain=="TRA",]
tcr_a$barcode<- make.names(tcr_a$barcode, unique = TRUE)
tcr_a$barcode<- gsub('\\b.1\\b', '.contig2', tcr_a$barcode, ignore.case = TRUE)
tcr_a$barcode<- gsub('\\b.2\\b', '.contig3', tcr_a$barcode, ignore.case = TRUE)
tcr_a$barcode<- gsub('\\b.3\\b', '.contig4', tcr_a$barcode, ignore.case = TRUE)

tcr_b <- tcr[tcr$chain=="TRB",]
tcr_b$barcode<- make.names(tcr_b$barcode, unique = TRUE)
tcr_b$barcode<- gsub('\\b.1\\b', '.contig2', tcr_b$barcode, ignore.case = TRUE)
tcr_b$barcode<- gsub('\\b.2\\b', '.contig3', tcr_b$barcode, ignore.case = TRUE)
tcr_b$barcode<- gsub('\\b.3\\b', '.contig4', tcr_b$barcode, ignore.case = TRUE)

# Add alpha and beta sequences to fully-productive TCRs
tcr_fa <- read.fasta("TCR/outs/all_contig.fasta")



tcr_a$full_a_seq <- toupper(unlist(getSequence(tcr_fa[names(tcr_fa) %in% tcr_a$contig_id], as.string = TRUE)))
tcr_b$full_b_seq <- toupper(unlist(getSequence(tcr_fa[names(tcr_fa) %in% tcr_b$contig_id], as.string = TRUE)))
rownames(tcr_a) <- tcr_a$barcode
rownames(tcr_b) <- tcr_b$barcode
tcr_a <- tcr_a[,33, drop=FALSE]
tcr_b <- tcr_b[,33, drop=FALSE]
rownames(tcr_a)<- paste0(experiment_name, ".", rownames(tcr_a), sep="")
colnames(tcr_a)<-c('full_a_seq_contig1')
paste(rownames(tcr_a))
paste(colnames(sc))

new_sc <- AddMetaData(object=sc, metadata=tcr_a)
rownames(tcr_b)<- paste0(experiment_name, ".", rownames(tcr_b), sep="")



tcr_a_contig2<- (tcr_a)
tcr_a_contig2[,2]<-rownames(tcr_a)
tcr_a_contig2<-tcr_a_contig2[(grepl('.contig2', rownames(tcr_a_contig2))),]
rownames_tcr_a_contig2<-tcr_a_contig2$V2
tcr_a_contig2<- as.data.frame(tcr_a_contig2[,1])
rownames(tcr_a_contig2)<-rownames_tcr_a_contig2
rownames(tcr_a_contig2) <- gsub(".contig2", "", rownames(tcr_a_contig2))
colnames(tcr_a_contig2)<-c('full_a_seq_contig2')

new_sc <- AddMetaData(object=new_sc, metadata=tcr_a_contig2)
colnames(tcr_b)<-c('full_b_seq_contig1')

new_sc <- AddMetaData(object=new_sc, metadata=tcr_b)

tcr_b_contig2<- (tcr_b)
tcr_b_contig2[,2]<-rownames(tcr_b)
tcr_b_contig2<-tcr_b_contig2[(grepl('.contig2', rownames(tcr_b_contig2))),]
rownames_tcr_b_contig2<-tcr_b_contig2$V2
tcr_b_contig2<- as.data.frame(tcr_b_contig2[,1])
rownames(tcr_b_contig2)<-rownames_tcr_b_contig2
rownames(tcr_b_contig2) <- gsub(".contig2", "", rownames(tcr_b_contig2))
colnames(tcr_b_contig2)<-c('full_b_seq_contig2')


new_sc <- AddMetaData(object=new_sc, metadata=tcr_b_contig2)


meta <- new_sc@meta.data
filt_meta <- meta[!with(meta,is.na(full_a_seq_contig1) | is.na(full_b_seq_contig1)),]
tcr_cells <- rownames(filt_meta)
sc_fullprod_tcr <- subset(new_sc, cells = tcr_cells)




 
# Add clonotypic information to fully-productive TCRs
tcr_a <- tcr[tcr$chain=="TRA",]
tcr_a$barcode<- make.names(tcr_a$barcode, unique = TRUE)
tcr_a$barcode<- gsub('\\b.1\\b', '.contig2', tcr_a$barcode, ignore.case = TRUE)


tcr_b <- tcr[tcr$chain=="TRB",]
tcr_b$barcode<- make.names(tcr_b$barcode, unique = TRUE)
tcr_b$barcode<- gsub('\\b.1\\b', '.contig2', tcr_b$barcode, ignore.case = TRUE)



tcr_a <- tcr_a[,c("barcode", "raw_clonotype_id","v_gene","j_gene","c_gene",
                'fwr1_nt', 'cdr1_nt', 'fwr2_nt', 'cdr2_nt', 'fwr3_nt', 
                'cdr3', 'cdr3_nt', 'fwr4_nt','reads','umis')]
tcr_b <- tcr_b[,c("barcode", "raw_clonotype_id","v_gene","d_gene","j_gene","c_gene",
                'fwr1_nt', 'cdr1_nt', 'fwr2_nt', 'cdr2_nt', 'fwr3_nt', 'cdr3',
                'cdr3_nt', 'fwr4_nt','reads','umis')]
names(tcr_a)[names(tcr_a) == "raw_clonotype_id"] <- "clonotype_id"
names(tcr_b)[names(tcr_b) == "raw_clonotype_id"] <- "clonotype_id"

# Reorder so barcodes are first column and set them as rownames.
rownames(tcr_a) <- tcr_a[,1]
tcr_a[,1] <- NULL
rownames(tcr_b) <- tcr_b[,1]
tcr_b[,1] <- NULL



tcr_a_contig2<- (tcr_a)
tcr_a_contig2<-tcr_a_contig2[(grepl('.contig2', rownames(tcr_a_contig2))),]
rownames(tcr_a_contig2) <- gsub(".contig2", "", rownames(tcr_a_contig2))

tcr_a_contig1<-tcr_a[!(grepl('.contig2', rownames(tcr_a))),]


colnames(tcr_a_contig2)[2:14] <- paste("alpha_ctg2", colnames(tcr_a_contig2)[2:14], sep = "_")
colnames(tcr_a_contig1)[2:14] <- paste("alpha_ctg1", colnames(tcr_a_contig1)[2:14], sep = "_")


tcr_b_contig2<- (tcr_b)
tcr_b_contig2<-tcr_b_contig2[(grepl('.contig2', rownames(tcr_b_contig2))),]
rownames(tcr_b_contig2) <- gsub(".contig2", "", rownames(tcr_b_contig2))

tcr_b_contig1<-tcr_b[!(grepl('.contig2', rownames(tcr_b))),]







# Reorder so barcodes are first column and set them as rownames.

colnames(tcr_b_contig1)[2:15] <- paste("beta_ctg1", colnames(tcr_b_contig1)[2:15], sep = "_")   
colnames(tcr_b_contig2)[2:15] <- paste("beta_ctg2", colnames(tcr_b_contig2)[2:15], sep = "_")   


rownames(tcr_a_contig1)<- paste0(experiment_name, ".", rownames(tcr_a_contig1), sep="")
rownames(tcr_a_contig2)<- paste0(experiment_name, ".", rownames(tcr_a_contig2), sep="")
rownames(tcr_b_contig1)<- paste0(experiment_name, ".", rownames(tcr_b_contig1), sep="")
rownames(tcr_b_contig2)<- paste0(experiment_name, ".", rownames(tcr_b_contig2), sep="")


# Add to the Seurat object's metadata.
clono_seurat <- AddMetaData(object=sc_fullprod_tcr, metadata=tcr_a_contig1)
clono_seurat <- AddMetaData(object=clono_seurat, metadata=tcr_a_contig2)
clono_seurat <- AddMetaData(object=clono_seurat, metadata=tcr_b_contig1)
clono_seurat <- AddMetaData(object=clono_seurat, metadata=tcr_b_contig2)

clono_seurat$v_cd3r_clono_ctg1 <- paste(clono_seurat$alpha_ctg1_v_gene, clono_seurat$alpha_ctg1_cdr3, clono_seurat$beta_ctg1_v_gene, clono_seurat$beta_ctg1_cdr3, sep = ".")
clono_seurat$v_cd3r_clono_ctg2 <- paste(clono_seurat$alpha_ctg2_v_gene, clono_seurat$alpha_ctg2_cdr3, clono_seurat$beta_ctg2_v_gene, clono_seurat$beta_ctg2_cdr3, sep = ".")

clono_seurat$IGT_cellID <- Cells(clono_seurat)
saveRDS(clono_seurat, "seuratobject_withHTOADT_singlet_fullTCRinfo_productiveonly.Rds")

# Create Alpha and Beta FASTQ files
#library(seqRFLP)
clono_meta <- clono_seurat@meta.data
clono_meta_ctg2_alpha <- clono_meta[!with(clono_meta,is.na(full_a_seq_contig2)),]
clono_meta_ctg2_beta <- clono_meta[!with(clono_meta,is.na(full_b_seq_contig2)),]

alpha_contig1 <- data.frame(paste(paste0(clono_meta$IGT_cellID,'.ctg1'), paste0("V-Gene:", clono_meta$alpha_ctg1_v_gene), paste0("CDR3:", clono_meta$alpha_ctg1_cdr3), sep = " "), clono_meta$full_a_seq_contig1)
alpha_contig2 <- data.frame(paste(paste0(clono_meta_ctg2_alpha$IGT_cellID,'.ctg2'), paste0("V-Gene:", clono_meta_ctg2_alpha$alpha_ctg2_v_gene), paste0("CDR3:", clono_meta_ctg2_alpha$alpha_ctg2_cdr3), sep = " "), clono_meta_ctg2_alpha$full_a_seq_contig2)
colnames(alpha_contig1)<- c('cell_V-Gene_CDR3_alpha','contig')
colnames(alpha_contig2)<- c('cell_V-Gene_CDR3_alpha','contig')
alpha<- rbind(alpha_contig1,alpha_contig2)
alpha<-alpha[order(alpha$`cell_V-Gene_CDR3_alpha`), ]
write.csv(alpha, paste0(fasta_name, "alpha.csv"), row.names = T, quote = F)

beta_contig1 <- data.frame(paste(paste0(clono_meta$IGT_cellID,'.ctg1'), paste0("V-Gene:", clono_meta$beta_ctg1_v_gene), paste0("CDR3:", clono_meta$beta_ctg1_cdr3), sep = " "), clono_meta$full_b_seq_contig1)
beta_contig2 <- data.frame(paste(paste0(clono_meta_ctg2_beta$IGT_cellID,'.ctg2'), paste0("V-Gene:", clono_meta_ctg2_beta$beta_ctg2_v_gene), paste0("CDR3:", clono_meta_ctg2_beta$beta_ctg2_cdr3), sep = " "), clono_meta_ctg2_beta$full_b_seq_contig2)
colnames(beta_contig1)<- c('cell_V-Gene_CDR3_beta','contig')
colnames(beta_contig2)<- c('cell_V-Gene_CDR3_beta','contig')
beta<- rbind(beta_contig1,beta_contig2)
beta<-beta[order(beta$`cell_V-Gene_CDR3_beta`), ]
write.csv(beta, paste0(fasta_name, "beta.csv"), row.names = T, quote = F)


## Write fasta via seqinr
alpha_seqs <- as.list(alpha$contig)
names(alpha_seqs) <- alpha$'cell_V-Gene_CDR3_alpha'
write.fasta(sequences = alpha_seqs, names = names(alpha_seqs), file.out = paste0(fasta_name, "alpha.fasta"))


beta_seqs <- as.list(beta$contig)
names(beta_seqs) <- beta$'cell_V-Gene_CDR3_beta'
write.fasta(sequences = beta_seqs, names = names(beta_seqs), file.out = paste0(fasta_name, "beta.fasta"))


## Custom function
#dataframe2fas <- function(df, file = "output.fasta", chain) {
#  #if (!all(c("cell_V-Gene_CDR3_beta", "contig") %in% colnames(df))) {
#  #  stop("Data frame must have columns 'cell_V-Gene_CDR3_beta' and 'contig'")
#  #}
#
#  con <- file(file, open = "w")
#  on.exit(close(con))
#
#  for (i in seq_len(nrow(df))) {
#	  if (chain == "alpha") {
#		      writeLines(paste0(">", df$"cell_V-Gene_CDR3_alpha"[i]), con)
#	  }
#	  if (chain == "beta") {
#                      writeLines(paste0(">", df$"cell_V-Gene_CDR3_beta"[i]), con)
#          }
#
#
#    writeLines(df$contig[i], con)
#  }
#}

#dataframe2fas(alpha, paste0(fasta_name, "alpha.fasta"), chain  = "alpha")
#dataframe2fas(beta, paste0(fasta_name, "beta.fasta"), chain  = "beta")
