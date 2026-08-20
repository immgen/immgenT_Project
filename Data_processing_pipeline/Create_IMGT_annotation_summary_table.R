# Create_IMGT_annotation_summary_table.R
# Joins IMGT High-V Quest's TRA/TRB annotation output back onto the
# TCR-annotated Seurat object's metadata, producing a per-cell clonotype
# summary table.
#
# Run this AFTER:
#   1. tcr_annotate_all_chains.R has written the alpha/beta FASTAs (see README).
#   2. You've submitted those FASTAs to IMGT High-V Quest
#      (https://www.imgt.org/HighV-QUEST/) as two separate jobs (alpha,
#      beta). This requires an IMGT account - create/log into one before
#      submitting. IMGT emails you when each job finishes (typically ~6
#      hours to a day).
#   3. You've downloaded each job's result archive and unzipped it with
#      `tar xvf`, e.g.:
#        mkdir -p TCR_annotations/TRA TCR_annotations/TRB
#        tar xvf <alpha_job_results>.txz -C TCR_annotations/TRA
#        tar xvf <beta_job_results>.txz -C TCR_annotations/TRB
#      Each folder now has IMGT's numbered output files - this script only
#      reads `6_Junction.txt` from each.
#
# Usage:
#   Rscript Create_IMGT_annotation_summary_table.R <TRA_dir> <TRB_dir> \
#     <path_to_tcr_object> <path_to_final_object> <samples.csv> <EXP.txt> \
#     [output_dir]
#
#   TRA_dir               directory with IMGT's unzipped TRA output (has 6_Junction.txt)
#   TRB_dir               directory with IMGT's unzipped TRB output (has 6_Junction.txt)
#   path_to_tcr_object    seuratobject_singlet_fullTCRinfo_productiveANDnonproductive_pairedANDunpaired.Rds
#                         (written by tcr_annotate_all_chains.R)
#   path_to_final_object  seuratobject_singlet_postRNAfiltering_postADTfiltering_postTfiltering.Rds
#                         (written by Tcell_filter_part2.R - supplies the
#                         Enhanced table's cluster/UMAP columns)
#   samples.csv           hashtag -> sample_name map, same file used
#                         throughout the pipeline
#   EXP.txt               one-line file with this run's dataset/EXP ID
#   output_dir            where to write the output tables (default ".")
#
# Outputs (all prefixed with this run's EXP ID):
#   <EXP>_TRA.txt / <EXP>_TRB.txt           - reduced IMGT junction columns, as-is
#   <EXP>_summary_table.tsv / _xl.csv       - one row per cell: hashtag, sample_name,
#                                             alpha/beta V/J/CDR3 calls (both contigs
#                                             where present), and inferred clonotype IDs
#   <EXP>_Enhanced_summary_table_xl.csv     - the above + RNA cluster and UMAP coordinates

library(dplyr)
library(stringr)
library(Seurat)

args = commandArgs(TRUE)
path_to_TRA_dir      = args[1]
path_to_TRB_dir      = args[2]
path_to_tcr_object   = args[3]
path_to_final_object = args[4]
path_to_samples      = args[5]
path_to_EXP          = args[6]
output_dir           = if (length(args) >= 7) args[7] else "."

EXP <- read.table(path_to_EXP)
EXP <- as.character(EXP)

# The barcode/contig info IMGT's Sequence.ID carries forward is exactly what
# tcr_annotate_all_chains.R wrote as each FASTA header's cell ID:
#   <EXP>.<16bp 10x barcode>.<ctg1/ctg2>
# The prefix length depends on how long EXP is for this run, so it's
# computed here rather than hardcoded - change EXP and this still works.
prefix_len   = nchar(EXP) + 1      # +1 for the "." separator after EXP
barcode_len  = 16                  # standard 10x cell barcode length
contig_start = prefix_len + barcode_len + 2   # +2 skips the "." before ctgN
contig_end   = contig_start + 3               # "ctg1"/"ctg2" is 4 characters

message("Read IMGT Junction files")
TRA <- read.delim(file.path(path_to_TRA_dir, "6_Junction.txt"), header = TRUE)
TRB <- read.delim(file.path(path_to_TRB_dir, "6_Junction.txt"), header = TRUE)

TRA_reduced <- dplyr::select(TRA, c("Sequence.number", "Sequence.ID", "V.DOMAIN.Functionality", "V.GENE.and.allele", "J.GENE.and.allele", "JUNCTION..AA.", "JUNCTION", "N.REGION"))
TRB_reduced <- dplyr::select(TRB, c("Sequence.number", "Sequence.ID", "V.DOMAIN.Functionality", "V.GENE.and.allele", "J.GENE.and.allele", "JUNCTION..AA.", "JUNCTION", "N1.REGION", "N2.REGION"))

colnames(TRA_reduced) <- c("id", "barcode", "alpha.functionality", "alpha.vgene", "alpha.jgene", "alpha.junction", "alpha.junction.nt", "alpha.vj.nregion")
colnames(TRB_reduced) <- c("id", "barcode", "beta.functionality", "beta.vgene", "beta.jgene", "beta.junction", "beta.junction.nt", "beta.vd.nregion", "beta.dj.nregion")

TRA_reduced <- TRA_reduced[,-1]
TRB_reduced <- TRB_reduced[,-1]

write.table(TRA_reduced, file.path(output_dir, paste0(EXP, "_TRA.txt")), sep = "\t", row.names = F, col.names = T, quote = T, na = "")
write.table(TRB_reduced, file.path(output_dir, paste0(EXP, "_TRB.txt")), sep = "\t", row.names = F, col.names = T, quote = T, na = "")

message("Read TCR-annotated Seurat object and build hashtag -> sample_name map")
sc <- readRDS(path_to_tcr_object)
sc_metadata <- sc@meta.data %>%
  mutate(cellID = EXP_cellID, hashtag = HTO_classification.simplified) %>%
  dplyr::select(cellID, hashtag)

# Same hashtag -> sample_name lookup as add_sample_names.R - works whether
# HTO=YES (real hashtags) or HTO=NO (samples.csv has one row mapping the
# placeholder "no_hashing" to your sample name).
hashing <- read.csv(path_to_samples, header = TRUE)
hashing <- na.omit(hashing)
sc_metadata$sample_name <- hashing[match(sc_metadata$hashtag, hashing[, 1]), 2]

message("Format TRA/TRB calls and split by contig")
TRA_formatted <- TRA_reduced %>%
  mutate(contig = substr(barcode, contig_start, contig_end)) %>%
  mutate(
    barcode = substr(barcode, prefix_len + 1, prefix_len + barcode_len),
    alpha.vgene = str_remove_all(unlist(rapply(str_extract_all(alpha.vgene, "(TRA[0-9A-Z-/]+)\\*"), function(x) paste(x, collapse = "|"), how = "replace")), "\\*"),
    alpha.jgene = str_remove_all(unlist(rapply(str_extract_all(alpha.jgene, "(TRA[0-9A-Z-/]+)\\*"), function(x) paste(x, collapse = "|"), how = "replace")), "\\*"),
    alpha.junction.nt = toupper(alpha.junction.nt),
    alpha.vj.nregion = toupper(alpha.vj.nregion)
  )
TRA_contig1 <- TRA_formatted %>% filter(contig == 'ctg1') %>% dplyr::select(-contig)
TRA_contig2 <- TRA_formatted %>% filter(contig == 'ctg2') %>% dplyr::select(-contig)
TRA_out <- TRA_contig1 %>% left_join(TRA_contig2, by = "barcode")

TRB_formatted <- TRB_reduced %>%
  mutate(contig = substr(barcode, contig_start, contig_end)) %>%
  mutate(
    barcode = substr(barcode, prefix_len + 1, prefix_len + barcode_len),
    beta.vgene = str_remove_all(unlist(rapply(str_extract_all(beta.vgene, "(TRB[0-9A-Z-/]+)\\*"), function(x) paste(x, collapse = "|"), how = "replace")), "\\*"),
    beta.jgene = str_remove_all(unlist(rapply(str_extract_all(beta.jgene, "(TRB[0-9A-Z-/]+)\\*"), function(x) paste(x, collapse = "|"), how = "replace")), "\\*"),
    beta.junction.nt = toupper(beta.junction.nt),
    beta.vd.nregion = toupper(beta.vd.nregion),
    beta.dj.nregion = toupper(beta.dj.nregion)
  )
TRB_contig1 <- TRB_formatted %>% filter(contig == 'ctg1') %>% dplyr::select(-contig)
TRB_contig2 <- TRB_formatted %>% filter(contig == 'ctg2') %>% dplyr::select(-contig)
TRB_out <- TRB_contig1 %>% left_join(TRB_contig2, by = "barcode")

tcr_out <- TRA_out %>% full_join(TRB_out, by = "barcode")
colnames(tcr_out) <- c("cellID","alpha.functionality","alpha.vgene","alpha.jgene","alpha.junction","alpha.junction.nt","alpha.vj.nregion","alpha2.functionality","alpha2.vgene","alpha2.jgene","alpha2.junction","alpha2.junction.nt","alpha2.vj.nregion","beta.functionality","beta.vgene","beta.jgene","beta.junction","beta.junction.nt","beta.vd.nregion","beta.dj.nregion","beta2.functionality","beta2.vgene","beta2.jgene","beta2.junction","beta2.junction.nt","beta2.vd.nregion","beta2.dj.nregion")

sc_metadata$cellID <- rownames(sc_metadata)
sc_metadata$cellID <- substr(sc_metadata$cellID, prefix_len + 1, prefix_len + barcode_len)

message("Join IMGT calls onto the Seurat metadata and infer clonotype IDs")
summary_table <- sc_metadata %>% left_join(tcr_out, by = "cellID") %>%
  mutate(clonotype_alpha_beta = case_when((!is.na(alpha.vgene))&(!is.na(beta.vgene)) ~ paste(alpha.vgene, alpha.jgene, alpha.junction, beta.vgene, beta.jgene, beta.junction, sep = "."))) %>%
  mutate(clonotype_alpha_beta2 = case_when((!is.na(alpha.vgene))&(!is.na(beta2.vgene)) ~ paste(alpha.vgene, alpha.jgene, alpha.junction, beta2.vgene, beta2.jgene, beta2.junction, sep = "."))) %>%
  mutate(clonotype_alpha2_beta = case_when((!is.na(alpha2.vgene))&(!is.na(beta.vgene)) ~ paste(alpha2.vgene, alpha2.jgene, alpha2.junction, beta.vgene, beta.jgene, beta.junction, sep = "."))) %>%
  mutate(clonotype_alpha2_beta2 = case_when((!is.na(alpha2.vgene))&(!is.na(beta2.vgene)) ~ paste(alpha2.vgene, alpha2.jgene, alpha2.junction, beta2.vgene, beta2.jgene, beta2.junction, sep = ".")))

summary_table$cellID <- paste0(EXP, ".", summary_table$cellID)

write.table(summary_table, file.path(output_dir, paste0(EXP, "_summary_table.tsv")), sep = "\t", row.names = F, col.names = T, quote = T, na = "")
write.csv(summary_table, file.path(output_dir, paste0(EXP, "_summary_table_xl.csv")), row.names = F, col.names = T, quote = T, na = "")

message("Build Enhanced summary table (adds RNA cluster + UMAP coordinates)")
so <- readRDS(path_to_final_object)
so_metadata <- so@meta.data[,c("cell_barcode", "RNA_snn_res.0.5")]
colnames(so_metadata) <- c("cellID", "clusterID")
so_metadata$UMAP_x <- so@reductions$umap_rna@cell.embeddings[,1]
so_metadata$UMAP_y <- so@reductions$umap_rna@cell.embeddings[,2]

enhanced_summary_table <- summary_table[, c(1:3)] %>%
  merge(so_metadata, by = "cellID") %>%
  merge(summary_table[, -c(2:3)], by = "cellID")

write.csv(enhanced_summary_table, file.path(output_dir, paste0(EXP, "_Enhanced_summary_table_xl.csv")), row.names = F, col.names = T, quote = T, na = "")

message("Done.")
