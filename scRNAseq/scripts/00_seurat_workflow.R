# 00_seurat_workflow

library(here)
library(Seurat)
library(dplyr)
library(ggplot2)
library(purrr)
library(stringr)
library(SeuratObject)

set.seed(8675309)

# read in data -----

counts_mtx_dir <- here("scRNAseq", "data", "counts_matrices")
so_list <- counts_mtx_dir %>%
  list.files() %>% rlang::set_names() %>%
  imap(~{
    file.path(counts_mtx_dir, ..1) %>%
      Seurat::Read10X() %>%
      Seurat::CreateSeuratObject(project = ..2, min.cells = 3, min.features = 200)
  })

so_raw <- merge(so_list[[1]], y = so_list[2:length(so_list)],
                add.cell.ids = names(so_list))

# add metadata -----

so_raw[["percent.mt"]] <- PercentageFeatureSet(so_raw, pattern = "^mt-")
so_raw[["percent.ribo"]] <- PercentageFeatureSet(so_raw, pattern = "^Rpl|^Rps")

genotype_metadata <- so_raw@meta.data %>%
  mutate(genotype = str_extract(orig.ident, ".")) %>%
  select(genotype)
so_raw <- AddMetaData(so_raw, metadata = genotype_metadata)

# filter by QC features ------

so <- subset(so_raw, subset = nFeature_RNA > 1000 & nFeature_RNA < 3500 & percent.mt < 6)

# filter and rerun clustering workflow -----

so <- so %>%
  NormalizeData() %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 3000) %>%
  ScaleData(features = rownames(.)) %>%
  RunPCA()

so <- so %>%
  JackStraw(num.replicate = 100) %>%
  ScoreJackStraw(dims = 1:20)

so <- so %>%
  FindNeighbors(dims = 1:15) %>%
  RunUMAP(dims = 1:15)

resolutions <- c(0.15, 0.175, 0.18, 0.8)
for (res in resolutions) {
  so <- FindClusters(so, resolution = res)
}

# add cluster metadata -----

cluster_ids_key <- c(
  "4" = "progenitor exhausted",
  "1" = "effector-like",
  "2" = "transitory exhausted",
  "0" = "terminally exhausted",
  "3" = "dividing"
)

cluster_num_key <- c(
  "4" = "progenitor",
  "1" = "effector-like",
  "2" = "transitory",
  "0" = "terminal",
  "3" = "dividing"
)

cluster_ids_metadata <- so@meta.data %>%
  mutate(cluster_ids = factor(cluster_ids_key[as.character(RNA_snn_res.0.175)],
                              levels = cluster_ids_key)) %>%
  mutate(Clusters = factor(cluster_num_key[as.character(RNA_snn_res.0.175)],
                           levels = cluster_num_key)) %>%
  dplyr::select(cluster_ids, Clusters)

so <- AddMetaData(so, cluster_ids_metadata)

saveRDS(so, here("scRNAseq", "data", "processed_data_objects", "00_seurat_workflow.so.rds"))

# in case needed, run sc clustering workflow on full dataset -----

so_raw <- so_raw %>%
  NormalizeData() %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 3000) %>%
  ScaleData(features = rownames(.)) %>%
  RunPCA()

so_raw <- so_raw %>%
  JackStraw(num.replicate = 100) %>%
  ScoreJackStraw(dims = 1:40)

so_raw <- so_raw %>%
  FindNeighbors(dims = 1:15) %>%
  RunUMAP(dims = 1:15)

resolutions <- c(0.1, 0.2, 0.8)
for (res in resolutions) {
  so_raw <- FindClusters(so_raw, resolution = res)
}

saveRDS(so_raw, here("scRNAseq", "data", "processed_data_objects", "00_seurat_workflow.so_raw.rds"))
