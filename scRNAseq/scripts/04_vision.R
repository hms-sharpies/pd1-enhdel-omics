# 04_vision

library(here)
library(dplyr)
library(ggplot2)
library(purrr)
library(stringr)
library(tibble)
library(readr)
library(Seurat)
library(SeuratObject)
library(VISION)

# read in data -----

so <- readRDS(here("scRNAseq", "data", "processed_data_objects", "00_seurat_workflow.so.rds"))

LCMVexh_signature_path <- here("scRNAseq", "data", "gene_signatures", "LCMVGeneSets.gmt")
vo <- Vision(so, signatures = LCMVexh_signature_path)
vo <- analyze(vo)

out_path <- here("scRNAseq", "data", "processed_data_objects", "04_vision.vo.rds")
saveRDS(vo, out_path)

# extract signature scores -----

sscores_df <- vo@SigScores %>%
  as_tibble(rownames = "cell")

out_path <- here("scRNAseq", "data", "processed_data_tables", "04_vision.SigScores.tsv")
readr::write_tsv(sscores_df, out_path)
