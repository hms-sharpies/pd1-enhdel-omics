# 03_pseudobulk

library(here)
library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(stringr)
library(tibble)
library(readr)
library(SeuratObject)

library(SingleCellExperiment)
library(scater)
library(Matrix.utils)
library(DESeq2)
library(Rsc)

detach("package:SeuratObject", unload = T)

set.seed(42)

# read in data -----

so <- readRDS(here("scRNAseq", "data", "processed_data_objects", "00_seurat_workflow.so.rds"))
metadata_df <- so@meta.data %>%
  as_tibble(rownames = "cell")
smallest_cluster_cell_count <- metadata_df %>%
  filter(Clusters != "dividing") %>%
  group_by(Clusters, genotype) %>%
  summarise(count = n()) %>%
  .$count %>% min()
cells <- metadata_df %>%
  filter(Clusters != "dividing") %>%
  group_by(Clusters, genotype) %>%
  slice_sample(n = smallest_cluster_cell_count) %>%
  .$cell
so_down <- subset(so, cells = cells)

# make pseudobulk ------

pb_list <- make_pseudobulk(so_down, "Clusters", "orig.ident", "genotype")
deseq_list <- map(pb_list, ~pb_to_deseq(..1, sample_var = "orig.ident", vars = "genotype"))
rlog_list <- map(deseq_list, rlog, blind = TRUE)

out_path <- here("scRNAseq", "data", "processed_data_objects", "05_pseudobulk.deseq_list.rds")
saveRDS(deseq_list, out_path)

counts <- purrr::map_dfr(
  deseq_list,
  ~counts(..1, normalized = FALSE) %>%
    as_tibble(rownames = "gene"),
  .id = "cluster"
) %>%
  pivot_longer(cols = 3:ncol(.),
               names_to = "sample",
               values_to = "norm_count") %>%
  pivot_wider(id_cols = gene,
              names_from = c(sample, cluster),
              values_from = norm_count,
              values_fill = 0)
write_tsv(
  counts,
  here("scRNAseq", "data", "processed_data_tables",
       "05_pseudobulk.counts.tsv")
)

norm_counts <- purrr::map_dfr(
    deseq_list,
    ~counts(..1, normalized = TRUE) %>%
      as_tibble(rownames = "gene"),
    .id = "cluster"
  ) %>%
  pivot_longer(cols = 3:ncol(.),
               names_to = "sample",
               values_to = "norm_count") %>%
  pivot_wider(id_cols = gene,
              names_from = c(sample, cluster),
              values_from = norm_count,
              values_fill = 0)
write_tsv(
  norm_counts,
  here("scRNAseq", "data", "processed_data_tables",
       "05_pseudobulk.norm_counts.tsv")
)

# pca loadings -----

get_pca_genes <- function(deseq_obj, ntop = 500) {
  rv <- deseq_obj %>%
    assay() %>%
    rowVars()
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  pca <- prcomp(t(assay(deseq_obj)[select,]))
  loadings <- as_tibble(pca$rotation, rownames = "gene")
  return(loadings)
}
pca_loadings <- map(rlog_list, ~get_pca_genes(..1, ntop = dim(..1)[1])) %>%
  imap(~mutate(..1, cluster_id = ..2)) %>%
  purrr::reduce(bind_rows)
out_path <- here("scRNAseq", "data", "processed_data_tables", "05_pseudobulk.pca_loadings.tsv")
readr::write_tsv(pca_loadings, out_path)
