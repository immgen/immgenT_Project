#### V gene ####
library(readr)
library(stringr)
library(dplyr)
library(tidyr)
library(readxl)

Mutations_Lpart_pureJax <- read_excel("./Mutations_Lpart_pureJax.xlsx")
Mutations_longStrings_pureJax <- read_csv("./Mutations_longStrings_pureJax.csv", col_types = cols_only(`Sequence ID` = col_guess()))
filter_out_contigs <- unique(c(Mutations_Lpart_pureJax$`Sequence ID`, Mutations_longStrings_pureJax$`Sequence ID`))

file_names <- list.files(path = "./samtools_outputs_Vgene/", full.names = FALSE)

IGTs <- sapply(file_names, function(file_name) {
  parts <- str_split(file_name, "_")[[1]]
  return(parts[1])
})

unique_IGTs <- unique(IGTs)

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

# Function to filter qualities and corresponding UMIs
filter_qualities_umis <- function(qualities, umis, qual_threshold) {
  if (qualities == "") return("")
  
  qual_values <- as.numeric(unlist(strsplit(qualities, ",")))
  umi_values <- unlist(strsplit(umis, ","))
  
  if (length(qual_values) != length(umi_values)) return("")
  
  filtered_umis <- umi_values[qual_values > qual_threshold]
  paste(filtered_umis, collapse = ",")
}

# Function to extract unique UMIs and their counts
get_unique_umis_counts <- function(umis) {
  if (umis == "") return(c("", ""))
  
  umi_values <- unlist(strsplit(umis, ","))
  unique_umis <- unique(umi_values)
  umi_counts <- sapply(unique_umis, function(u) sum(umi_values == u))
  
  unique_umis_str <- paste(unique_umis, collapse = ",")
  umi_counts_str <- paste(umi_counts, collapse = ",")
  
  return(c(unique_umis_str, umi_counts_str))
}

# Function to sum UMI counts
sum_umi_counts <- function(counts) {
  if (counts == "") return(0)
  sum(as.numeric(unlist(strsplit(counts, ","))))
}

phred_conversion <- function(phred_string) {
  if(phred_string == '*' | phred_string == '') {
    return("")
  } else {
    phred_scores <- utf8ToInt(phred_string) - 33
    return(paste(phred_scores, collapse = ","))
  }
}

mean_quality <- function(qualities) {
  if (qualities == "") return(0)
  
  qual_values <- as.numeric(unlist(strsplit(qualities, ",")))
  
  return(round(mean(qual_values)))
}

# A helper function to safely count UMIs, handling empty strings
count_umis <- function(umi_string) {
  if (is.na(umi_string) || nchar(umi_string) == 0) {
    return(0)
  } else {
    return(length(strsplit(umi_string, ",")[[1]]))
  }
}



for(i in 1:length(unique_IGTs)){
  IGT_number = unique_IGTs[i]
  # IGT_number = "IGT2"
  mutations_details <- read_delim(paste0("./samtools_outputs_Vgene/",IGT_number,"_mutations_tmp.tsv"), delim = "\t", escape_double = FALSE, trim_ws = TRUE)
  mutations_details <- mutations_details %>% mutate(
    # Take the part before the first comma
    first_part = str_split(IMGT.annotation, ",", n = 2, simplify = TRUE)[,1],
    # Get everything after the last ">"
    after_arrow = str_extract(first_part, "[^>]+$"),
    # Take the last character, convert to uppercase
    mutation.base = str_to_upper(str_sub(after_arrow, -1))
  ) %>% select(-c(first_part, after_arrow))
  samtools_outputs <- read_delim(paste0("./samtools_outputs_Vgene/",IGT_number,"_samtools_outputs_minBQ_0.txt"), delim = "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
  colnames(samtools_outputs) <- c("contig.id", "v.mut.contig.pos","mutation.base.tmp", "total.mapped.reads", "read.bases", "read.base.qualities","UMIs")
  
  
  mutations_details_out <- mutations_details %>% 
    left_join(samtools_outputs, by = c("contig.id", "v.mut.contig.pos")) %>% 
    select(-c(v.insertions, v.deletions, v.ins.contig)) %>% 
    mutate(read.bases.modified = str_remove_all(read.bases, "-[0-9]+[ACGTNacgtn]+|\\+[0-9]+[ACGTNacgtn*#]+|\\*")) %>% 
    mutate(read.bases.modified = toupper(str_replace_all(read.bases.modified, ",", mutation.base.tmp)))
  
  mutations_details_out <- mutations_details_out %>% rowwise() %>% 
    mutate(IMGT.reference.base.qualities = paste0(str_sub(read.base.qualities, str_locate_all(read.bases.modified, IMGT.reference.base)[[1]][,1], str_locate_all(read.bases.modified, IMGT.reference.base)[[1]][,1]), collapse = "")) %>% 
    mutate(mutation.base.qualities = paste0(str_sub(read.base.qualities, str_locate_all(read.bases.modified, mutation.base)[[1]][,1], str_locate_all(read.bases.modified, mutation.base)[[1]][,1]), collapse = "")) %>%
    mutate(IMGT.reference.UMIs = extract_UMIs(read.bases.modified, IMGT.reference.base, mutation.base, UMIs)[1],) %>% 
    mutate(mutation.base.UMIs = extract_UMIs(read.bases.modified, IMGT.reference.base, mutation.base, UMIs)[2])
  
  mutations_details_out <- mutations_details_out %>% mutate(IMGT.reference.base.qualities.modified = phred_conversion(IMGT.reference.base.qualities), mutation.base.qualities.modified = phred_conversion(mutation.base.qualities), functionality = str_replace_all(functionality,"\\s\\(see\\scomment\\)",""), vgene = str_replace_all(vgene,"\\s\\(see\\scomment\\)",""), jgene = str_replace_all(jgene,"\\s\\(see\\scomment\\)",""))
  
  mutations_details_out <- mutations_details_out %>% 
    mutate(
      IMGT.reference.base.qualities.modified.filtered = lapply(strsplit(IMGT.reference.base.qualities.modified, ","), 
                                                               function(x) paste(x[as.numeric(x) > 30], collapse = ",")),
      mutation.base.qualities.modified.filtered = lapply(strsplit(mutation.base.qualities.modified, ","), 
                                                         function(x) paste(x[as.numeric(x) > 30], collapse = ",")),
      IMGT.reference.UMIs.filtered = mapply(filter_qualities_umis, IMGT.reference.base.qualities.modified, IMGT.reference.UMIs, 20),
      mutation.base.UMIs.filtered = mapply(filter_qualities_umis, mutation.base.qualities.modified, mutation.base.UMIs, 30),
      IMGT.reference.unique.UMIs = mapply(function(x) get_unique_umis_counts(x)[1], IMGT.reference.UMIs.filtered),
      IMGT.reference.UMI.read.counts = mapply(function(x) get_unique_umis_counts(x)[2], IMGT.reference.UMIs.filtered),
      mutation.base.unique.UMIs = mapply(function(x) get_unique_umis_counts(x)[1], mutation.base.UMIs.filtered),
      mutation.base.UMI.read.counts = mapply(function(x) get_unique_umis_counts(x)[2], mutation.base.UMIs.filtered),
      IMGT.reference.base.count = mapply(sum_umi_counts, IMGT.reference.UMI.read.counts),
      mutation.base.count = mapply(sum_umi_counts, mutation.base.UMI.read.counts),
      total.mapped.reads.filtered = IMGT.reference.base.count + mutation.base.count,
      IMGT.referece.mean.quality = mean_quality(IMGT.reference.base.qualities.modified.filtered),
      mutation.base.mean.quality = mean_quality(mutation.base.qualities.modified.filtered)
    ) %>%
    select(Sequence.ID, TR, functionality, vgene, jgene, v.mutations, v.mut.gene, v.mut.contig, contig.id, contig.length, v.mut.gene.pos, v.mut.contig.pos, IMGT.annotation, IMGT.reference.base, mutation.base, total.mapped.reads, total.mapped.reads.filtered, IMGT.reference.base.count, mutation.base.count, IMGT.referece.mean.quality, mutation.base.mean.quality, IMGT.reference.unique.UMIs, IMGT.reference.UMI.read.counts, mutation.base.unique.UMIs, mutation.base.UMI.read.counts)
  

  mutations_details_out <- mutations_details_out %>%
    rowwise() %>%
    mutate(
      # Use the helper function for robust counting
      IMGT.reference.UMI.counts = count_umis(IMGT.reference.unique.UMIs),
      mutation.base.UMI.counts = count_umis(mutation.base.unique.UMIs),
      
      # Identify the shared UMIs (this part works correctly even with empty strings)
      shared.UMIs = list(intersect(
        strsplit(IMGT.reference.unique.UMIs, ",")[[1]],
        strsplit(mutation.base.unique.UMIs, ",")[[1]]
      )),
      
      # Count the number of shared UMIs
      shared.UMI.counts = length(shared.UMIs),
      
      # Create a formatted string for shared UMI read counts
      shared.UMI.read.counts = {
        # Proceed only if there are shared UMIs to process
        if (length(shared.UMIs) > 0) {
          ref_umis <- strsplit(IMGT.reference.unique.UMIs, ",")[[1]]
          ref_counts <- strsplit(IMGT.reference.UMI.read.counts, ",")[[1]]
          mut_umis <- strsplit(mutation.base.unique.UMIs, ",")[[1]]
          mut_counts <- strsplit(mutation.base.UMI.read.counts, ",")[[1]]
          
          ref_map <- setNames(ref_counts, ref_umis)
          mut_map <- setNames(mut_counts, mut_umis)
          
          sapply(shared.UMIs, function(umi) {
            paste0(
              umi,
              "(",
              ref_map[umi],
              ":",
              mut_map[umi],
              ")"
            )
          }) %>% paste(collapse = ";")
        } else {
          # Return an empty string if there are no shared UMIs
          ""
        }
      }
    ) %>%
    select(-shared.UMIs)
  

  #productive
  mutations_details_out_productive <- mutations_details_out %>% filter(functionality == "productive") %>% filter(!(Sequence.ID %in% filter_out_contigs))
  
  mutations_details_out_filtered_contig <- mutations_details_out_productive %>% filter(v.mut.contig.pos > 10)
  
  mutations_details_out_candidate_contig <- mutations_details_out_filtered_contig %>% filter((mutation.base.count > 5) & !str_detect(vgene, str_c(c("TRAV9N-4","TRAV9-1", "TRAV9D-1", "TRAV13D3"), collapse = "|")))
  write_csv(mutations_details_out_candidate_contig, paste0("./",IGT_number,"_mutations_details_candidates_filtered_by_contig_pos.csv"), na = "")

}

