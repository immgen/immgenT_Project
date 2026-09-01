# David Zemmour
# R
# usage: Rscript limma_fit_level2.IGTHT.R [path_to_seurat_object] [path_to_tmm_object] [output_dir] [fit_file_name]

options(max.print=1000)

# Parse arguments
args = commandArgs(TRUE)
if (length(args) < 4) {
    stop("Usage: Rscript limma_fit_level2.IGTHT.R [path_to_seurat_object] [path_to_tmm_object] [output_dir] [fit_file_name]")
}
path_to_seurat_object = args[1]
path_to_tmm_object = args[2]
output_dir = args[3] 
fit_file_name = args[4]

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

if (!file.exists(path_to_tmm_object)) {
    stop("The specified TMM object file does not exist: ", path_to_tmm_object)
}

message("loading R libraries and custom functions")
libs = c("limma", "edgeR", "Seurat","BPCells", "dplyr", "rlang","reshape2") 
sapply(libs, function(x) suppressMessages(suppressWarnings(library(x, character.only = TRUE, quietly = T, warn.conflicts  = F))))

# source("/project/jfkfloor2/zemmourlab/david/immgent/immgent_integration_git/custom_functions_david.R")

message("loading seurat object")
so = readRDS(file = path_to_seurat_object)

message("loading tmm object")
tmm = readRDS(file = path_to_tmm_object)

message("Fitting")
design = data.frame(so@meta.data[,c("annotation_level2","IGTHT")]) 
design$annotation_level2.IGTHT = paste(design$annotation_level2, design$IGTHT, sep = ".")
if (!all(rownames(design) == colnames(so))) {
    stop("Row names of 'design' do not match column names of the Seurat object. Check your input data.")
}
design = model.matrix(~ 0 + annotation_level2.IGTHT, data=design) 
colnames(design) = gsub("annotation_level2.IGTHT", "", colnames(design))
fit = lmFit(tmm, design = design)
saveRDS(fit, file = sprintf("%s/%s",output_dir,fit_file_name))
message("Fit object saved to: ", sprintf("%s/%s",output_dir,fit_file_name))