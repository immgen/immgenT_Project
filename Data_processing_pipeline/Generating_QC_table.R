# Generating_QC_table.R
# Stitch together the per-stage QC stats files into the two final summary
# tables (per-hashtag and per-IGT/run).
#
# This replaces both Generating_QC_table.R (RNA+ADT) and Generating_QC_table_RNA.R
# (RNA only), driven by ADT=YES/NO.
#
# Usage:
#   Rscript Generating_QC_table.R <preRNAfiltering QC txt> <preADTfiltering QC txt> \
#                                  <postTfiltering csv> <IGT postTfiltering csv> <ADT YES/NO>
#
#   <preADTfiltering QC txt> is only read if ADT == YES. Pass any placeholder
#   (e.g. "NA") for it if ADT == NO.

args = commandArgs(TRUE)
preRNAfiltering    = args[1]
preADTfiltering    = args[2]
postTfiltering     = args[3]
IGT_postTfiltering = args[4]
ADT                = toupper(args[5]) == "YES"

libs = c("Seurat", "ggplot2", "inflection", "grid", "gridExtra", "UCell")
sapply(libs, function(x) suppressMessages(suppressWarnings(library(x, character.only = TRUE, quietly = T, warn.conflicts = F))))

RNAfiltering = read.table(preRNAfiltering, header = TRUE)

if (ADT) {
  ADTfiltering = read.table(preADTfiltering, header = TRUE)
  tmp = cbind(RNAfiltering, ADTfiltering)
  keeps = c('HTO_classification.simplified', 'ncells', 'ncells_outliers_nGenes', 'ncells_outliers_deadcells',
            'ncells_outliers_lowCountADT', 'ncells_outliers_autofluorescence')
  tmp = tmp[keeps]
  table2 = read.csv(postTfiltering, header = TRUE)
  table2 = table2[, -1]
  final_table = merge(tmp, table2)
} else {
  keeps = c('HTO_classification.simplified', 'ncells', 'ncells_outliers_nGenes', 'ncells_outliers_deadcells')
  tmp = RNAfiltering[keeps]
  table2 = read.csv(postTfiltering, header = TRUE)
  final_table = cbind(tmp, table2)
}

write.csv(final_table, 'seuratobject_singlet_postRNAfiltering_postADTfiltering_postTfiltering_FinalTable.csv', row.names = FALSE)

table3 = read.csv(IGT_postTfiltering, header = TRUE, row.names = 1)

## Sum every numeric QC column across hashtags for the one-row-per-run table;
## drop the ID/count columns that shouldn't be summed.
non_summable = c('HTO_classification.simplified', names(final_table)[grepl('^ncells', names(final_table))])
sum_cols = setdiff(names(final_table), non_summable)

final_table_IGT = as.data.frame(t(colSums(final_table[, sum_cols, drop = FALSE])))
rownames(final_table_IGT) = rownames(table3)
final_table_IGT = cbind(final_table_IGT, table3)

write.csv(final_table_IGT, 'seuratobject_IGT_singlet_postRNAfiltering_postADTfiltering_postTfiltering_FinalTable.csv', row.names = TRUE)
