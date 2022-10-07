# 01_dge

library(here)
library(Seurat)
library(dplyr)
library(ggplot2)
library(purrr)
library(stringr)
library(tibble)
library(readr)
library(SeuratObject)

# read in data -----

so <- readRDS(here("scRNAseq", "data", "processed_data_objects", "00_seurat_workflow.so.rds"))

# cluster comparisons -----

dge_cluster_df <- FindAllMarkers(so %>% set_idents("Clusters")) %>%
  as_tibble()
out_path <- here("scRNAseq", "data", "processed_data_tables", "01_dge.cluster_one-v-all.tsv")
readr::write_tsv(dge_cluster_df, out_path)

# genotype comparisons within cluster -----

findmarkers_genotype_one.v.one_in_cluster <- function(so, g1, g2, cluster) {
  message(paste0("Clusters: ", cluster, "; genotypes: ", g1, " ", g2))
  Idents(so) <- "Clusters"
  fm_df <- FindMarkers(
      so,
      assay = "RNA",
      ident.1 = g1, ident.2 = g2, group.by = "genotype",
      subset.ident = cluster,
      # set the following parameters to get the complete list of ranked genes
      logfc.threshold = -Inf,
      min.pct = -Inf,
      min.diff.pct = -Inf,
      return.thresh = 1
    ) %>%
    tibble::rownames_to_column("gene") %>%
    mutate(contrast = paste0(g1, "_", g2)) %>%
    mutate(Clusters = cluster)
  return(fm_df)
}
genotypes <- so@meta.data$genotype %>% unique()
genotype_combn2 <- combn(genotypes, 2, simplify = FALSE)
clusters <- so@meta.data$Clusters %>% unique() %>% as.character()
dge_df <- cross(list(genotype_combn2, clusters)) %>%
  map_dfr(~findmarkers_genotype_one.v.one_in_cluster(so, ..1[[1]][1], ..1[[1]][2], ..1[[2]]))

out_path <- here("scRNAseq", "data", "processed_data_tables", "01_dge.genotype_one-v-one_per-cluster.tsv")
readr::write_tsv(dge_df, out_path)

