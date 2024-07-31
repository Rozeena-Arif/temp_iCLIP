
library(tidyverse)

# Read command line arguments
args <- commandArgs(trailingOnly = TRUE)
SINV_annotation <- args[1]
output_file <- args[2]

# Read the annotation file, modify it, and write the output
read_tsv(SINV_annotation, skip = 3, col_names = FALSE) %>%
  mutate(X9 = str_replace(X9, "original_biotype", "gene_biotype")) %>%
  mutate(X9 = case_when(
    str_detect(X9, "UTR") & !str_detect(X9, "gene_biotype") ~ paste0(X9, " gene_biotype \"rna\";"),
    TRUE ~ X9
  )) %>%
  mutate(X9 = case_when(
    str_detect(X9, "gene_biotype") ~ X9,
    TRUE ~ paste0(X9, " gene_biotype \"protein_coding\";")
  )) %>%
  mutate(X9 = str_remove_all(X9, "\"Genbank:.*\" ")) %>%
  write_tsv(output_file, col_names = FALSE, escape = "none")

