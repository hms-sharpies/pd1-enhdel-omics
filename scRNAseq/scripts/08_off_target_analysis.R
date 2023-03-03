# 08_off_target_analysis

library(here)
library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(stringr)
library(tibble)
library(readr)
library(Seurat)
library(SeuratObject)
library(ggrepel)

theme_set(theme_classic())

# purpose ------

# To answer the question, "What are the transcriptional regulators of
# the exhaustion associated enhancer and how are they expressed in the
# different subpopulations?"

# read in data -----

so <- readRDS(here("scRNAseq", "data", "processed_data_objects", "00_seurat_workflow.so.rds"))
exh_enh_test_fimo <- read_tsv(here("scRNAseq", "data", "exh_enh_test_fimo.tsv")) %>%
  filter(!is.na(sequence_name)) %>%
  mutate(motif_gene = motif_alt_id %>% str_extract("\\.[:alnum:]*$")) %>%
  mutate(motif_gene = str_sub(motif_gene, start = 2)) %>%
  dplyr::relocate(motif_gene, .before = motif_id)

genes_df <- read_tsv(here("scRNAseq", "data", "genes.gtf"), skip = 5, col_names = FALSE) %>%
  filter(X3 == "gene") %>%
  separate(X9, into = c("gene_id", "gene_version", "gene_name", "gene_source", "gene_biotype"),
           sep = "; ") %>%
  mutate(gene_name = str_sub(gene_name, start = 12, end = -2)) %>%
  dplyr::relocate(gene_name, .before = X1)
