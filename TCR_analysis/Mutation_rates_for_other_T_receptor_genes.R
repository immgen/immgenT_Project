library(readr)
library(stringr)
library(dplyr)
library(tidyr)

IGT_number = "IGT10" # IGT10, IGT11, IGT12

gene_coding_regions_all <- read_delim("./gene_coding_regions.txt", 
                                      delim = "\t", escape_double = FALSE, 
                                      trim_ws = TRUE)
# gene_coding_regions.txt format:
# gene	chr	gene_start	gene_end	exon_start	exon_end
# Traj17	14	54201775	54201837	54201775	54201837

################################################################################
# Function to extract regions, handling overlapping regions
extract_regions <- function(data, regions) {
  # Get unique positions that fall within any region
  positions <- integer()
  
  # Loop through each region
  for(i in 1:nrow(regions)) {
    # Get all positions in current region
    positions <- union(positions, 
                       regions$exon_start[i]:regions$exon_end[i])
  }
  
  # Sort positions
  positions <- sort(unique(positions))
  
  # Extract rows that match any of these positions
  result <- data[data$position %in% positions, ]
  
  # Sort by position
  result <- result[order(result$position), ]
  
  return(result)
}

# Function to convert quality to numerical score
mean_phred <- function(phred_string) {
  if(phred_string == '*' | phred_string == '') {
    return(0)
  } else {
    phred_scores <- utf8ToInt(phred_string) - 33
    return(mean(phred_scores))
  }
}

extract_mutations <- function(df) {
  df %>%
    mutate(
      reads = strsplit(read.bases.modified, ""),
      qualities = strsplit(read.base.qualities, ""),
      UMIs = strsplit(UMIs, ",")
    ) %>%
    filter(
      lengths(reads) == lengths(qualities), 
      lengths(reads) == lengths(UMIs),
      total.count > 5
    ) %>%
    unnest(cols = c(reads, qualities, UMIs)) %>%
    # unnest(cols = c(reads, qualities)) %>%
    filter(reads != reference.base, !grepl("[<>]", reads)) %>%  # Remove '<' and '>'
    rename(mutation = reads, quality = qualities, UMIs = UMIs) %>%
    # rename(mutation = reads, quality = qualities) %>%
    mutate(quality = sapply(quality, mean_phred))  # Convert quality to numerical score
}

extract_non_mutations <- function(df) {
  df %>%
    mutate(
      reads = strsplit(read.bases.modified, ""),
      qualities = strsplit(read.base.qualities, ""),
      UMIs = strsplit(UMIs, ",")
    ) %>%
    filter(
      lengths(reads) == lengths(qualities), 
      lengths(UMIs) == lengths(UMIs),
      total.count > 5
    ) %>%
    unnest(cols = c(reads, qualities, UMIs)) %>%
    # unnest(cols = c(reads, qualities)) %>%
    filter(reads == reference.base, !grepl("[<>]", reads)) %>%  # Remove '<' and '>'
    rename(mutation = reads, quality = qualities, UMIs = UMIs) %>%
    # rename(mutation = reads, quality = qualities) %>%
    mutate(quality = sapply(quality, mean_phred))  # Convert quality to numerical score
}

# Function to extract unique UMIs and their counts
get_unique_umis_counts <- function(umis) {
  if (umis == "") return(0)
  
  umi_values <- unlist(strsplit(umis, ","))
  unique_umis <- unique(umi_values)

  return(length(unique_umis))
}

# Function to extract UMIs based on position
extract_UMIs <- function(read_bases, ref_base, mut_base, umi_str) {
  bases <- unlist(strsplit(read_bases, ""))
  umi_list <- unlist(strsplit(umi_str, ","))
  
  ref_indices <- which(bases == ref_base)
  mut_indices <- which(bases == mut_base)
  
  ref_UMIs <- paste(umi_list[ref_indices], collapse = ",")
  mut_UMIs <- paste(umi_list[mut_indices], collapse = ",")
  
  return(c(ref_UMIs, mut_UMIs))
}

################################################################################
################################################################################

other_genes = c("Cd4", "Cd8a", "Cd28", "Cd3e")
mutation_rates = c()

for(genename in other_genes){
  gene_coding_regions <- gene_coding_regions_all %>% filter(gene == genename) %>% select(exon_start, exon_end)
  
  samtools_outputs <- read_delim(paste0("./",IGT_number,"_samtools_outputs_minBQ_0_",genename,".txt"), delim = "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
  colnames(samtools_outputs) <- c("chr", "position","reference.base", "total.mapped.reads", "read.bases", "read.base.qualities","UMIs")
  
  samtools_outputs <- samtools_outputs %>%
    mutate(read.bases.modified = str_remove_all(read.bases, "-[0-9]+[ACGTNacgtn]+|\\+[0-9]+[ACGTNacgtn*#]+|\\*")) %>% 
    mutate(read.bases.modified = toupper(str_replace_all(read.bases.modified, ",", reference.base))) %>% 
    mutate(read.bases.modified = toupper(str_replace_all(read.bases.modified, "\\.", reference.base)))
  
  samtools_outputs <- samtools_outputs %>% rowwise() %>% 
    mutate(total.count = str_count(read.bases.modified, paste0("[^><]"))) %>% ungroup()
  
  # Extract regions
  samtools_outputs_filtered <- extract_regions(samtools_outputs, gene_coding_regions)
  
  mutations_df <- extract_mutations(samtools_outputs_filtered)
  mutations_df_filtered <- mutations_df %>% filter(quality > 30) %>% select(chr, position, reference.base, mutation, quality, UMIs) %>% add_count(chr, position, reference.base, mutation, quality, UMIs, name = "mutation_read_count") %>% distinct()
  mutations_df_filtered_UMIs_over1 <- mutations_df_filtered %>% group_by(position) %>% filter(n() > 1) %>% ungroup()
  
  mutation_rates = c(mutation_rates, sum(mutations_df_filtered_UMIs_over1$mutation_read_count)/nrow(samtools_outputs_filtered)/(sum(samtools_outputs_filtered$total.count)/90)*100)
}

out_tbl = data.frame(genesymbol = other_genes, mutation_rate = mutation_rates)
write.table(out_tbl, file = paste0(IGT_number, "_mutation_rates_corrected_UMIs_over_1.txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
