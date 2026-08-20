# build_seurat_object.R
# Build the Seurat object from RNA (GEX) + optional Feature Barcode (ADT/HTO)
# matrices, and demultiplex on HTO if hashtags are present.
#
# Usage:
#   Rscript build_seurat_object.R <GEX filtered_feature_bc_matrix dir> <FBC filtered_feature_bc_matrix dir> \
#                          <dir containing samples.csv> <ADT YES/NO> <HTO YES/NO>
#
#   The FBC path is ignored if ADT=NO and HTO=NO.
#   samples.csv is only read if HTO=YES (see README for its format).
#
# When HTO=NO, this script still creates HTO_classification / .global /
# .simplified columns and fills them with a constant placeholder ("no_hashing"
# / "Singlet") so every downstream script that groups or facets by
# HTO_classification keeps working without any changes.

args = commandArgs(TRUE)
path_to_rna           = args[1]
path_to_hto_adt        = args[2]
path_to_sample_names   = args[3]
ADT                    = toupper(args[4])
HTO                    = toupper(args[5])
ADT                    = (ADT == "YES")
HTO                    = (HTO == "YES")
FBC                    = (ADT | HTO)

libs = c("Seurat", "ggplot2", "gridExtra", "dplyr")
sapply(libs, function(x) suppressMessages(suppressWarnings(library(x, character.only = TRUE, quietly = T, warn.conflicts = F))))

library(RColorBrewer)
n = 70
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual', ]
mypal = unique(unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))))

message(sprintf("ADT = %s | HTO = %s | FBC data expected = %s", ADT, HTO, FBC))

message("Load RNA matrix")
umis = Read10X(data.dir = path_to_rna, gene.column = 2)
colnames(umis) = gsub(pattern = "\\-1", replacement = "", colnames(umis))

if (FBC) {

  message("Load FBC (ADT/HTO) matrix")
  htos_adt = Read10X(data.dir = path_to_hto_adt, gene.column = 1)
  colnames(htos_adt) = gsub(pattern = "\\-1", replacement = "", colnames(htos_adt))

  joint.bcs = intersect(colnames(htos_adt), colnames(umis))
  message(sprintf("Joint RNA-FBC barcodes: %s", length(joint.bcs)))
  umis = umis[, joint.bcs]

  if (HTO) {
    htos = htos_adt[rownames(htos_adt)[grepl("Hash|HT", rownames(htos_adt))], joint.bcs]
    message("These are the hashes: ")
    print(htos)
  }
  if (ADT) {
    adt = htos_adt[rownames(htos_adt)[!grepl("Hash|HT|unmapped", rownames(htos_adt))], joint.bcs]
  }
}

message("Create Seurat object")
so = CreateSeuratObject(counts = umis)
# Force the RNA assay back to the classic Seurat v3/v4 "Assay" class. Under
# Seurat v5, CreateSeuratObject() makes an "Assay5" object by default, and
# every downstream script in this pipeline uses v4-style slot access
# (so@assays$RNA@counts, @var.features, etc.), which doesn't exist on Assay5.
# CreateAssayObject() always returns a legacy Assay regardless of Seurat
# version, so this keeps the rest of the pipeline working unchanged. (ADT and
# HTO below are already created with CreateAssayObject(), so they don't need
# this - only the RNA assay, which CreateSeuratObject() builds internally.)
so[["RNA"]] = CreateAssayObject(counts = so@assays$RNA$counts)
so$cell_barcode = colnames(so)
so$Cell_ID = paste("cell", seq_len(ncol(so)), sep = "_")

if (ADT) {
  so[["ADT"]] = CreateAssayObject(counts = adt)
}

if (HTO) {

  hto_table = data.frame(hto_summed_counts = rowSums(htos))
  write.csv(t(hto_table), 'hto_counts_table.csv')

  hashing = read.csv(file.path(path_to_sample_names, "samples.csv"), header = TRUE)
  hashing = na.omit(hashing)
  rownames(hashing) = hashing[, 1]
  htos = htos[rownames(htos) %in% rownames(hashing), ]

  so[["HTO"]] = CreateAssayObject(counts = htos + 1)

  message(sprintf("Cells with no hashtag counts: %s", length(which(colSums(htos) == 0))))
  write.table(x = names(which(colSums(htos) == 0)), file = "cells_with_no_hashtag_counts.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
  so = so[, colSums(htos) != 0]

  message("Demultiplex cells based on HTO")
  so = NormalizeData(so, assay = "HTO", normalization.method = "CLR")
  so = HTODemux(so, assay = "HTO", positive.quantile = 0.999, verbose = TRUE)
  so$HTO_classification_orig = so$HTO_classification

  for (i in which(so$HTO_classification == "Negative")) {
    max_hash  = max(so@assays$HTO@counts[, i])
    max2_hash = max(so@assays$HTO@counts[, i][-which.max(so@assays$HTO@counts[, i])])
    if (max_hash > 2 * max2_hash & max_hash > 10) {
      so$HTO_classification[i] = names(which.max(so@assays$HTO@counts[, i]))
    }
  }

  png("make_seurat-1_HTO_ridgeplot.png", width = 1400, height = 1600)
  print(RidgePlot(so, assay = "HTO", features = rownames(so[["HTO"]]), ncol = 2))
  dev.off()

  panel.cor = function(x, y, digits = 2, cex.cor, ...) {
    usr = par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r = cor(x, y)
    txt = format(c(r, 0.123456789), digits = digits)[1]
    txt = paste("r= ", txt, sep = "")
    text(0.5, 0.6, txt)
    p = cor.test(x, y)$p.value
    txt2 = format(c(p, 0.123456789), digits = digits)[1]
    txt2 = paste("p= ", txt2, sep = "")
    if (p < 0.01) txt2 = paste("p= ", "<0.01", sep = "")
    text(0.5, 0.4, txt2)
  }

  png("make_seurat-2_HTO_pairs.png", width = 1400, height = 1200)
  pairs(t(so@assays$HTO@data), cex = 4, upper.panel = panel.cor, pch = ".")
  dev.off()

  so = ScaleData(so, assay = "HTO", features = rownames(so@assays$HTO), verbose = FALSE)
  so = RunPCA(so, assay = "HTO", features = rownames(so@assays$HTO), reduction.name = "pca_htos", approx = FALSE)
  so = RunTSNE(so, assay = "HTO", reduction = "pca_htos", reduction.name = "tsne_htos",
               dims = 1:ncol(so@reductions$pca_htos@cell.embeddings), perplexity = 100,
               verbose = TRUE, check_duplicates = FALSE, max_iter = 500)

  so$HTO_classification.simplified = as.character(so$HTO_classification)
  so$HTO_classification.simplified[so$HTO_classification.global == "Doublet"] = "Doublet"
  so$HTO_classification.global[so$HTO_classification.global == 'Negative' & so$HTO_classification.simplified != 'Negative'] = "Rescued"

  tmp = data.frame(so@meta.data, dim1 = so@reductions$tsne_htos@cell.embeddings[, 1], dim2 = so@reductions$tsne_htos@cell.embeddings[, 2])
  theme_my = theme(axis.text.x = element_text(size = 15), axis.text.y = element_text(size = 15), legend.text = element_text(size = 10),
                    axis.title.x = element_text(size = 20), axis.title.y = element_text(size = 20), legend.title = element_blank(),
                    panel.grid.major = element_blank(), panel.grid.minor = element_blank())

  p = ggplot(tmp) + geom_point(aes(dim1, dim2, color = HTO_classification.simplified), alpha = I(1), size = 1) +
    scale_color_manual(values = mypal) + theme_bw() + theme_my + facet_wrap(~HTO_classification.global, nrow = 2)
  q = HTOHeatmap(so, assay = "HTO", ncells = 5000)

  png("make_seurat-3_HTO_tsne_heatmap.png", width = 1400, height = 1200)
  grid.arrange(grobs = list(p, q), widths = c(1, 1), heights = c(1, 1, 1), nrow = 3, ncol = 2)
  dev.off()

} else {

  message("HTO = NO: skipping hashtag demultiplexing. Labeling all cells with a placeholder group ('no_hashing') so downstream scripts that group/facet by HTO_classification continue to work unchanged.")
  so$HTO_classification            = "no_hashing"
  so$HTO_classification.global     = "Singlet"
  so$HTO_classification.simplified = "no_hashing"
  so$HTO_classification_orig       = "no_hashing"
}

message("Saving seuratobject_alldata.Rds and seuratobject_singlet.Rds")

qc_summary = function(so) {
  tmp = data.frame(so@meta.data)
  if (HTO) {
    tmp %>% group_by(HTO_classification.simplified) %>%
      summarize(ncells = n(), mean_nCount_RNA = round(mean(nCount_RNA)), mean_nFeature_RNA = round(mean(nFeature_RNA)),
                mean_nCount_HTO = round(mean(nCount_HTO))) %>% as.data.frame()
  } else {
    tmp %>% group_by(HTO_classification.simplified) %>%
      summarize(ncells = n(), mean_nCount_RNA = round(mean(nCount_RNA)), mean_nFeature_RNA = round(mean(nFeature_RNA))) %>% as.data.frame()
  }
}

write.table(x = qc_summary(so), file = "seuratobject_alldata_QCstats.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

so$sample_name = NA

Idents(so) = 'HTO_classification.global'
demux_table = t(table(Idents(so)))
write.table(x = demux_table, file = "Demultiplexing_QC_stats.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

saveRDS(so, file = "seuratobject_alldata.Rds")

so = so[, so$HTO_classification.global %in% c("Singlet", "Rescued")]
Idents(so) = 'HTO_classification'

saveRDS(so, file = "seuratobject_singlet.Rds")

write.table(x = qc_summary(so), file = "seuratobject_singlet_QCstats.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
