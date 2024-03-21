# 00_making_df

library(here)
library(tidyverse)

data_dir <- here("bulkRNAseq" ,"data")

list.files(data_dir)

counts_df <- read_tsv(file.path(data_dir, "transcript_counts_all_samples.tsv"))
tpm_df <- read_tsv(file.path(data_dir, "transcript_tpm_all_samples.tsv"))

counts_df %>%
  colnames()

new_counts_df <- counts_df %>%
  dplyr::select(target_id, length, str_subset(colnames(.), "CX3CR6Neg|Slamf6"))
new_tpm_df <- tpm_df %>%
  dplyr::select(target_id, length, str_subset(colnames(.), "CX3CR6Neg|Slamf6"))

write_tsv(new_counts_df, file.path(data_dir, "transcript_counts.tsv"))
write_tsv(new_tpm_df, file.path(data_dir, "transcript_tpm.tsv"))
