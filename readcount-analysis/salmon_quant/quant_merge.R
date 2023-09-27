################################################################################
#
# Description: This script combines salmon gene quantification output files in a folder.
# Note, only genes intersection across all per-sample files will be reported in the merged matrix.
#
################################################################################

# Load the necessary libraries
library(dplyr)
library(tidyr)
library(purrr)
library(here)

# List all the quant.gene.sf files in the directory
file_list <- list.files(path = here("UMCCR/research/projects/Nuoyi_RNAseq/quant"), pattern = "quant.genes.sf", full.names = TRUE)

# Initialize an empty data frame to store the combined data
combined_data <- data.frame(Gene = character())

# Function to read and process a single file
read_and_process_file <- function(file_path) {
  file_data <- read.delim(file_path, header = TRUE, stringsAsFactors = FALSE)
  file_name <- tools::file_path_sans_ext(basename(file_path))
  gene_data <- file_data %>%
    dplyr::select(Name, NumReads) %>%
    dplyr::filter(!grepl('PAR_Y', Name)) %>%
    dplyr::mutate(Name = sub("\\..*", "", Name)) %>%
    dplyr::rename(Gene = Name) %>%
    dplyr::rename(!!file_path := NumReads)
  return(gene_data)
}

# Use purrr::map to read and process all files and combine them into one data frame
# Use reduce to perform a full join on the "Name" column 
combined_data <- purrr::map(file_list, read_and_process_file) %>%
  purrr::reduce(full_join, by = "Gene")

# Rename the columns to remove the full path and file extension
colnames(combined_data) <- gsub(".*/(.*).quant.genes.sf", "\\1", colnames(combined_data))


# Write the combined data frame to a CSV file
write.table(combined_data, file = "/Users/kanwals/UMCCR/research/projects/Nuoyi_RNAseq/quant/FBXO.counts.matrix.txt", sep = "\t", quote = FALSE, row.names = FALSE)
