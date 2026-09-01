# David Zemmour
# R
# usage: Rscript limma_make_tmm.R [path_to_seurat_object] [output_dir] [tmm_file_name]

options(max.print=1000)

args = commandArgs(TRUE)
if (length(args) < 3) {
    stop("Usage: Rscript limma_make_tmm.R [path_to_seurat_object] [output_dir] [tmm_file_name]")
}
path_to_seurat_object = args[1] 
output_dir = args[2] 
tmm_file_name = args[3]

# path_to_wd = "/project/zemmour/david/ImmgenT/analysis/data_integration/IGT1_96/CD4"
# path_to_seurat_object = "/project/zemmour/david/ImmgenT/analysis/data_integration/IGT1_96/CD4/igt1_96_CD4_20241113.Rds" #path to RNA: "count/outs/filtered_feature_bc_matrix/"

# Validate inputs
if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    cat("Folder created:", output_dir, "\n")
} else {
    cat("Folder already exists:", output_dir, "\n")
}

if (!file.exists(path_to_seurat_object)) {
    stop("The specified Seurat object file does not exist: ", path_to_seurat_object)
}

# setwd(output_dir)

message("loading R libraries and custom functions")
libs = c("limma", "edgeR", "Seurat","BPCells", "dplyr", "rlang","reshape2") 
sapply(libs, function(x) suppressMessages(suppressWarnings(library(x, character.only = TRUE, quietly = T, warn.conflicts  = F))))

# source("/project/jfkfloor2/zemmourlab/david/immgent/immgent_integration_git/custom_functions_david.R")

message("loading seurat object")
so = readRDS(file = path_to_seurat_object)

message("Removing TCR, mt, ribo, Gm, and Rik genes and genes that are not expressed")
tcr_genes = grepl(x = rownames(so), pattern = "Trbv|Trbd|Trbj|Trbc|Trav|Traj|Trac|Trgv|Trgd|Trgj|Trgc|Trdv|Trdj|Trdc")
gm_rik_genes = grepl(x = rownames(so), pattern = "Gm|Rik$|\\-ps$")
ribo_genes = grepl(x = rownames(so), pattern = "Rpl|Rps|Mrpl|Mrps|Rsl")
mt_genes = grepl(x = rownames(so), pattern = "^mt-")
genes_not_expressed = rowSums(so[["RNA"]]$counts) == 0
genes_to_keep = !(tcr_genes | gm_rik_genes | ribo_genes | mt_genes | genes_not_expressed)
cat("Number of genes to keep:", sum(genes_to_keep), "\n")

#Make tmm objet for limma
count = so[["RNA"]]$counts[genes_to_keep, ]
dge = DGEList(count)
dge = calcNormFactors(dge)
tmm = new("EList")
message("TMM normalization")
tmm$E = edgeR::cpm(dge, log = TRUE, prior.count = 0.1)

saveRDS(tmm, file = sprintf("%s/%s",output_dir,tmm_file_name))
message("Tmm object saved to: ", sprintf("%s/%s",output_dir,tmm_file_name))