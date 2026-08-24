# filter_empty_droplets.R
# Empty-droplet removal for scRNA-seq GEX data, with OPTIONAL Feature Barcode
# (ADT and/or HTO) data.
#
# ADT and HTO are both delivered through the same 10x "Feature Barcode / Antibody Capture"
# library lane, so from THIS script's point of view they collapse into one internal
# flag:
#     FBC = (ADT == "YES") | (HTO == "YES")
#
# IMPORTANT - HTO and ADT are independent of each other here:
#   - ADT=YES, HTO=YES -> FBC matrix contains both antibody and hashtag rows
#   - ADT=YES, HTO=NO  -> FBC matrix contains only antibody rows
#   - ADT=NO,  HTO=YES -> FBC matrix contains only hashtag rows (fully
#                          filtered and used - this script does not require
#                          ADT to be present to process/keep HTO data)
#   - ADT=NO,  HTO=NO  -> no Feature Barcode step at all; RNA-only
#
# This script itself never needs to tell ADT rows apart from HTO rows - it
# just filters/writes whatever "Antibody Capture" feature-type rows exist as
# one combined FBC pool. The split into a dedicated ADT assay and/or a
# dedicated HTO assay happens downstream, in build_seurat_object.R.
#
# Usage:
#   Rscript filter_empty_droplets.R <raw_feature_bc_matrix> <filtered_feature_bc_matrix> \
#                      <output_dir> <EXP.txt> <Method> <RNA_cutoff> <FBC_cutoff> \
#                      <ADT YES/NO> <HTO YES/NO>
#
#   raw_feature_bc_matrix       combined_sample/outs/raw_feature_bc_matrix
#   filtered_feature_bc_matrix  combined_sample/outs/filtered_feature_bc_matrix
#   output_dir                  directory QC plots get written to (usually ".")
#   EXP.txt                     one-line file with this run's dataset/EXP ID
#   Method                      "Automatic" (DropletUtils knee/inflection) or "Manual" (use RNA_cutoff/FBC_cutoff)
#   RNA_cutoff                  manual UMI cutoff for GEX (only used if Method == "Manual")
#   FBC_cutoff                  manual UMI cutoff for FBC - i.e. ADT and/or HTO, whichever exist (only used if Method == "Manual" and FBC data exists)
#   ADT                         YES/NO - does this dataset have an antibody-derived tag (ADT) library
#   HTO                         YES/NO - does this dataset have a hashtag oligo (HTO) library
#
# Outputs:
#   GEX/filtered_feature_bc_matrix            always
#   FBC/filtered_feature_bc_matrix(2)          only if ADT=YES and/or HTO=YES
#   rna_elbow_plot.png / fbc_elbow_plot.png    only under Method == "Automatic" (fbc_elbow_plot.png only if FBC data exists)
#   FBC_RNA_cutoff_plot.png                    only if FBC data exists
#   RNA_cutoff_plot.png                        only if FBC data does NOT exist

args = commandArgs(TRUE)
path_to_data     = args[1]
path_to_data2    = args[2]
main_output_dir  = args[3]
EXP              = args[4]
Method           = args[5]
RNA_cutoff       = as.numeric(args[6])
FBC_cutoff       = as.numeric(args[7])
ADT              = toupper(args[8]) == "YES"
HTO              = toupper(args[9]) == "YES"
FBC              = (ADT | HTO)

libs = c("Seurat", "ggplot2", "gridExtra", "dplyr", "DropletUtils")
sapply(libs, function(x) suppressMessages(suppressWarnings(library(x, character.only = TRUE, quietly = T, warn.conflicts = F))))

library(RColorBrewer)
n = 70
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual', ]
mypal = unique(unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))))

libs = c("Seurat", "ggplot2", "reshape2", "dplyr", "fitdistrplus", "ggExtra", "scales")
sapply(libs, function(x) suppressMessages(suppressWarnings(library(x, character.only = TRUE, quietly = T, warn.conflicts = F))))

EXP <- read.table(EXP)
EXP <- as.character(EXP)

message(sprintf("ADT = %s | HTO = %s | FBC (Feature Barcode: ADT and/or HTO) data expected = %s", ADT, HTO, FBC))

sce <- read10xCounts(path_to_data)
sce$Barcode <- paste0(EXP, '.', sce$Barcode)

filtered_data <- read10xCounts(path_to_data2)
filtered_data$Barcode <- paste0(EXP, '.', filtered_data$Barcode)

filtered_data_seurat <- Read10X(path_to_data2)

if (FBC) {
  filtered_data_seurat_rna <- CreateSeuratObject(counts = filtered_data_seurat$`Gene Expression`, assay = "RNA")
  filtered_data_seurat_fbc <- CreateSeuratObject(counts = filtered_data_seurat$`Antibody Capture`, assay = "FBC")
} else {
  filtered_data_seurat_rna <- CreateSeuratObject(counts = filtered_data_seurat, assay = "RNA")
}

fbc_index <- which(sce@rowRanges@elementMetadata$Type == 'Antibody Capture')

if (FBC & length(fbc_index) > 0) {
  fbc.counts  <- counts(sce[fbc_index, ])
  rna.counts  <- counts(sce[-fbc_index, ])
  fbc.counts2 <- sce[fbc_index, ]
  rna.counts2 <- sce[-fbc_index, ]
} else {
  if (FBC) message("ADT and/or HTO was requested but no 'Antibody Capture' features were found in the raw matrix - continuing as RNA-only.")
  FBC         <- FALSE
  rna.counts  <- counts(sce)
  rna.counts2 <- sce
}

dir.create("GEX", showWarnings = FALSE)
if (FBC) dir.create("FBC", showWarnings = FALSE)

write_fbc_matrix = function(cells_of_interest_fbcM, EXP) {
  write10xCounts(
    "FBC/filtered_feature_bc_matrix",
    counts(cells_of_interest_fbcM),
    barcodes = cells_of_interest_fbcM$Barcode,
    gene.id = rownames(cells_of_interest_fbcM),
    gene.symbol = rownames(cells_of_interest_fbcM),
    overwrite = TRUE
  )
  new_barcodes <- gsub(paste0(EXP, '.'), '', cells_of_interest_fbcM$Barcode)
  cells_of_interest_fbcM$Barcode <- new_barcodes
  write10xCounts(
    "FBC/filtered_feature_bc_matrix2",
    counts(cells_of_interest_fbcM),
    barcodes = cells_of_interest_fbcM$Barcode,
    gene.id = rownames(cells_of_interest_fbcM),
    gene.symbol = rownames(cells_of_interest_fbcM),
    overwrite = TRUE
  )
}

write_gex_matrix = function(cells_of_interest_rnaM) {
  write10xCounts(
    "GEX/filtered_feature_bc_matrix",
    counts(cells_of_interest_rnaM),
    barcodes = cells_of_interest_rnaM$Barcode,
    gene.id = cells_of_interest_rnaM@rowRanges@elementMetadata$Symbol,
    gene.symbol = cells_of_interest_rnaM@rowRanges@elementMetadata$Symbol,
    overwrite = TRUE
  )
}

if (Method == 'Manual') {

  br.out_rna <- barcodeRanks(rna.counts)
  cells_of_interest_rna  <- which(br.out_rna$total > RNA_cutoff)
  cells_of_interest_rnaM <- rna.counts2[, cells_of_interest_rna]
  # Same "must also be in CellRanger's own filtered call" check as Automatic
  # mode below, so that Manual and Automatic give identical results whenever
  # their cutoffs happen to coincide.
  cells_of_interest_rnaM <- cells_of_interest_rnaM[, cells_of_interest_rnaM$Barcode %in% filtered_data$Barcode]

  if (FBC) {
    br.out_fbc <- barcodeRanks(fbc.counts)
    cells_of_interest_fbc  <- which(br.out_fbc$total > FBC_cutoff)
    cells_of_interest_fbcM <- fbc.counts2[, cells_of_interest_fbc]

    cells_of_interest_rnaM <- cells_of_interest_rnaM[, cells_of_interest_rnaM$Barcode %in% cells_of_interest_fbcM$Barcode]
    cells_of_interest_fbcM <- cells_of_interest_fbcM[, cells_of_interest_fbcM$Barcode %in% cells_of_interest_rnaM$Barcode]
    cells_of_interest_fbcM <- cells_of_interest_fbcM[, cells_of_interest_fbcM$Barcode %in% filtered_data$Barcode]

    write_fbc_matrix(cells_of_interest_fbcM, EXP)
  }

  write_gex_matrix(cells_of_interest_rnaM)

  x_cutoff <- RNA_cutoff
  y_cutoff <- if (FBC) FBC_cutoff else NA

} else {

  ## Automatic: knee/inflection from DropletUtils::barcodeRanks

  png(file.path(main_output_dir, "rna_elbow_plot.png"), width = 900, height = 900)
  br.out_rna <- barcodeRanks(rna.counts)
  plot(br.out_rna$rank, br.out_rna$total, log = "xy", xlab = "Rank", ylab = "Total")
  o <- order(br.out_rna$rank)
  lines(br.out_rna$rank[o], br.out_rna$fitted[o], col = "red")
  abline(h = metadata(br.out_rna)$knee, col = "dodgerblue", lty = 2)
  abline(h = metadata(br.out_rna)$inflection, col = "forestgreen", lty = 2)
  legend("bottomleft", lty = 2, col = c("dodgerblue", "forestgreen"), legend = c("knee", "inflection"))
  mtext(paste0('RNA cutoff=', br.out_rna@metadata$inflection), side = 3)
  dev.off()

  cells_of_interest_rna  <- which(br.out_rna$total > br.out_rna@metadata$inflection)
  cells_of_interest_rnaM <- rna.counts2[, cells_of_interest_rna]
  cells_of_interest_rnaM <- cells_of_interest_rnaM[, cells_of_interest_rnaM$Barcode %in% filtered_data$Barcode]

  if (min(filtered_data_seurat_rna$nCount_RNA) > br.out_rna@metadata$inflection) {
    x_cutoff <- min(filtered_data_seurat_rna$nCount_RNA)
  } else {
    x_cutoff <- br.out_rna@metadata$inflection
  }

  if (FBC) {
    png(file.path(main_output_dir, "fbc_elbow_plot.png"), width = 900, height = 900)
    br.out_fbc <- barcodeRanks(fbc.counts)
    plot(br.out_fbc$rank, br.out_fbc$total, log = "xy", xlab = "Rank", ylab = "Total")
    o <- order(br.out_fbc$rank)
    lines(br.out_fbc$rank[o], br.out_fbc$fitted[o], col = "red")
    abline(h = metadata(br.out_fbc)$knee, col = "dodgerblue", lty = 2)
    abline(h = metadata(br.out_fbc)$inflection, col = "forestgreen", lty = 2)
    legend("bottomleft", lty = 2, col = c("dodgerblue", "forestgreen"), legend = c("knee", "inflection"))
    mtext(paste0('FBC (ADT/HTO) cutoff=', br.out_fbc@metadata$inflection), side = 3)
    dev.off()

    cells_of_interest_fbc  <- which(br.out_fbc$total > br.out_fbc@metadata$inflection)
    cells_of_interest_fbcM <- fbc.counts2[, cells_of_interest_fbc]
    cells_of_interest_rnaM <- cells_of_interest_rnaM[, cells_of_interest_rnaM$Barcode %in% cells_of_interest_fbcM$Barcode]
    cells_of_interest_fbcM <- cells_of_interest_fbcM[, cells_of_interest_fbcM$Barcode %in% cells_of_interest_rnaM$Barcode]
    cells_of_interest_fbcM <- cells_of_interest_fbcM[, cells_of_interest_fbcM$Barcode %in% filtered_data$Barcode]

    write_fbc_matrix(cells_of_interest_fbcM, EXP)

    if (min(filtered_data_seurat_fbc$nCount_FBC) > br.out_fbc@metadata$inflection) {
      y_cutoff <- min(filtered_data_seurat_fbc$nCount_FBC)
    } else {
      y_cutoff <- br.out_fbc@metadata$inflection
    }
  } else {
    y_cutoff <- NA
  }

  write_gex_matrix(cells_of_interest_rnaM)
}

## Diagnostic scatter plot: cells vs. empty droplets, RNA UMIs vs FBC UMIs (or RNA UMIs vs rank if no FBC)

data1 <- Read10X(data.dir = path_to_data)
data1_rna <- if (is.list(data1)) data1$`Gene Expression` else data1

if (FBC) {

  data3 <- Read10X(data.dir = "FBC/filtered_feature_bc_matrix2")

  umicount_gene <- as.matrix(colSums(data1$`Gene Expression`))
  umicount_fbc  <- as.matrix(colSums(data1$`Antibody Capture`))
  data <- as.data.frame(cbind(umicount_gene, umicount_fbc))
  data[data == 0] <- NA
  data <- na.omit(data)

  add_data <- data[rownames(data) %in% colnames(data3), ]
  add_data <- cbind(add_data, 'Cells')
  colnames(add_data)[3] <- "Classification"

  subtract_data <- data[!rownames(data) %in% colnames(data3), ]
  subtract_data <- cbind(subtract_data, 'Empty Droplets')
  colnames(subtract_data)[3] <- "Classification"

  data_total <- rbind(add_data, subtract_data)

  fbc_label <- if (ADT & HTO) "ADT+HTO" else if (ADT) "ADT" else "HTO"

  png(file.path(main_output_dir, "FBC_RNA_cutoff_plot.png"), width = 1600, height = 1600)
  print(
    ggplot(data_total, aes(x = V1, y = V2, col = Classification)) + geom_point() +
      xlab("Counts (RNA)") + ylab(sprintf("Counts (%s)", fbc_label)) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black")) +
      theme(axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15),
            legend.text = element_text(size = 10), axis.title.x = element_text(size = 20),
            axis.title.y = element_text(size = 20)) +
      scale_color_manual(values = c('Red', 'Black')) +
      geom_vline(xintercept = x_cutoff, linetype = "dotted", color = "red") +
      geom_hline(yintercept = y_cutoff, linetype = "dotted", color = "red") +
      ggtitle(paste0('RNA_cutoff=', x_cutoff, '/FBC_cutoff=', y_cutoff)) +
      scale_x_continuous(trans = log_trans(base = 10), limits = c(10, max(data_total$V1))) +
      scale_y_continuous(trans = log_trans(base = 10), limits = c(10, max(data_total$V2))) +
      annotation_logticks(sides = "bl")
  )
  dev.off()

} else {

  umicount_gene <- as.matrix(colSums(data1_rna))
  data <- as.data.frame(umicount_gene)
  colnames(data) <- "V1"
  data$rank <- rank(-data$V1)

  png(file.path(main_output_dir, "RNA_cutoff_plot.png"), width = 1600, height = 900)
  print(
    ggplot(data, aes(x = rank, y = V1)) + geom_point(alpha = 0.5, size = 0.5) +
      xlab("Barcode rank") + ylab("Counts (RNA)") +
      scale_x_continuous(trans = "log10") + scale_y_continuous(trans = "log10") +
      annotation_logticks(sides = "bl") +
      geom_hline(yintercept = x_cutoff, linetype = "dotted", color = "red") +
      ggtitle(paste0('RNA_cutoff=', x_cutoff)) +
      theme_bw()
  )
  dev.off()
}

message(sprintf(
  "Done. RNA cells written to GEX/filtered_feature_bc_matrix%s",
  if (FBC) " and FBC cells written to FBC/filtered_feature_bc_matrix." else "."
))
