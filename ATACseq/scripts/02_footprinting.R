# 02_footprinting.R

# Purpose: To identify which motifs are binding within the enhancer region
# utilizing the output of TOBIAS BINDetect (see footprinting.sh)

library(here)
library(dplyr)
library(readr)
library(purrr)
library(stringr)
library(tibble)
library(tidyr)

footprinting_dir <- here("ATACseq", "data", "footprinting")
bindetect_dir <- file.path(footprinting_dir, "BINDetect_output")

motif_dirs <- list.files(bindetect_dir) %>%
  str_subset("MA")

bindetect_df <- motif_dirs %>%
  set_names() %>%
  imap_dfr( ~ file.path(bindetect_dir, ..1, sprintf("%s_overview.txt", ..1)) %>%
              read_tsv(),
            .id = "motif")

# chr1_sequence_94073947-94075951
motifs_binding_in_enh <- bindetect_df %>%
  filter(TFBS_chr == "chr1") %>%
  filter(TFBS_start > 94073947) %>%
  filter(TFBS_end < 94075951)

# # If desired, you can also export all motifs
# saveRDS(
#   bindetect_df,
#   here(
#     "ATACseq",
#     "data",
#     "processed_data_tables",
#     "02_footprinting.TFBS.rds"
#   )
# )

# Export motifs in enhancer
write_tsv(
  motifs_binding_in_enh,
  here(
    "ATACseq",
    "data",
    "processed_data_tables",
    "02_footprinting.TFBS_in_enh.tsv"
  )
)
