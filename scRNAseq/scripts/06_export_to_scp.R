# 06_export_to_scp

library(here)
library(Seurat)
library(dplyr)
library(ggplot2)
library(purrr)
library(stringr)
library(tibble)
library(readr)
library(SeuratObject)

scp_dir <- here("scRNAseq", "data", "single_cell_portal")
if (!dir.exists(scp_dir)) dir.create(scp_dir)

# read in data -----

so <-
  readRDS(here(
    "scRNAseq",
    "data",
    "processed_data_objects",
    "00_seurat_workflow.so.rds"
  ))

# export counts and data matrices ------

write_tsv(
  so@assays$RNA@counts %>% as_tibble(rownames = "GENE"),
  file.path(scp_dir, "scp_counts.txt.gz")
)
write_tsv(
  so@assays$RNA@data %>% as_tibble(rownames = "GENE"),
  here(
    "scRNAseq",
    "data",
    "single_cell_portal",
    "scp_lognorm_counts.txt.gz"
  )
)

# export metadata -------

scp_export_metadata_file(
  so,
  output_dir = scp_dir,
  biosample_id = orig.ident,
  donor_id = str_extract(biosample_id, "[:digit:]"),
  species = "NCBITaxon_10090",
  species__ontology_label = "Mus musculus",
  disease = "MONDO_0001449",
  disease__ontology_label = "lymphocytic choriomeningitis",
  organ = "UBERON_0002106",
  organ__ontology_label = "spleen",
  library_preparation_protocol = "EFO_0009899",
  library_preparation_protocol__ontology_label = "10x 3' v2",
  str_subset_regex_group = "genotype",
  str_subset_regex_numeric = "percent|nFeature|nCount"
)

# export clustering -----

scp_export_clustering_file(
  so,
  "Clusters",
  "_RNA",
  "umap",
  output_dir = scp_dir,
  filename = "scp_clustering"
)
