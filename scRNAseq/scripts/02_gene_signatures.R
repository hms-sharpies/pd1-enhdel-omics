# 02_gene_signatures

library(here)
library(dplyr)
library(ggplot2)
library(purrr)
library(stringr)
library(tibble)
library(readr)

library(msigdbr)

hallmarks_gs <- msigdbr(species = "Mus musculus", category = "H")
c2_gs <- msigdbr(species = "Mus musculus", category = "C2")
c7_gs <- msigdbr(species = "Mus musculus", category = "C7")
gs_df <- bind_rows(
  hallmarks_gs,
  c2_gs,
  c7_gs
)
write_tsv(gs_df, here("scRNAseq", "data", "gene_signatures",
                      "02_gene_signatures.msigdbr_H-C2-C7.tsv"))

