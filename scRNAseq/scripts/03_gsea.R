# 03_gsea

library(here)
library(dplyr)
library(ggplot2)
library(purrr)
library(stringr)
library(tibble)
library(readr)

library(liger)
source(here("helper.R"))

# read in data -----

gs_df <- read_tsv(here("scRNAseq", "data", "gene_signatures",
                       "02_gene_signatures.msigdbr_H-C2-C7.tsv"))
gs_df_key <- gs_df %>%
  select(gs_name, gs_cat)
gs_list <- gs_df %>%
  named_group_split(gs_name) %>%
  map(~..1$gene_symbol)

dge_path <- here("scRNAseq", "data", "processed_data_tables", "01_dge.genotype_one-v-one_per-cluster.tsv")
dge_df <- readr::read_tsv(dge_path)

# run GSEA using liger functions -----

make_rank_list <- function(.data,
                           gene = gene,
                           rank_by = log2FoldChange) {
  rank_by <- enquo(rank_by)
  gene <- enquo(gene)
  data_fc <- .data %>%
    arrange(desc(!!rank_by))
  data_fc %>%
    .[[rlang::as_name(rank_by)]] %>%
    set_names(data_fc[[rlang::as_name(gene)]])
}

rank_lists_one.v.all_ordered.log2FC <- dge_df %>%
  named_group_split(Clusters, contrast) %>%
  map(~make_rank_list(..1, rank_by = avg_logFC, gene = gene))
gsea_res_one.v.all_ordered.log2FC <- map_dfr(
  rank_lists_one.v.all_ordered.log2FC,
  ~liger::iterative.bulk.gsea(
    values = ..1,
    set.list = gs_list,
    power = 0
  ) %>% as_tibble(rownames = "gs_name"),
  .id = "cluster_contrast"
)
gsea_res_one.v.all_ordered.log2FC_tbl <- gsea_res_one.v.all_ordered.log2FC %>%
  separate(cluster_contrast, into = c("Clusters", "contrast"), sep = " / ")

write_tsv(
  gsea_res_one.v.all_ordered.log2FC_tbl,
  here("scRNAseq", "data", "processed_data_tables", "03_gsea.genotype_one-v-one_per-cluster_log2FC.tsv")
)

# rerun with more permutations for just oxphos and eff exh ------

# gsea_specific_df <- rank_lists_one.v.all_ordered.log2FC %>%
#   map_dfr(
#     ~ liger::bulk.gsea(
#       ..1,
#       set.list = gs_list[c("HALLMARK_OXIDATIVE_PHOSPHORYLATION",
#                            "GSE9650_EFFECTOR_VS_EXHAUSTED_CD8_TCELL_UP")],
#       n.rand = 1e6,
#       power = 0
#     ) %>% as_tibble(rownames = "gs_name"),
#     .id = "Clusters_contrast"
#   )
#
# gsea_specific_df <- gsea_specific_df %>%
#   separate(Clusters_contrast, into = c("Clusters", "contrast"), sep = " / ")
#
# write_tsv(
#   gsea_specific_df,
#   here("scRNAseq", "data", "processed_data_tables", "03_gsea.genotype_one-v-one_per-cluster_log2FC_specific.tsv")
# )

# also with Treg and Apoptosis --------

# gsea_specific2_df <- rank_lists_one.v.all_ordered.log2FC %>%
#   map_dfr(
#     ~ liger::bulk.gsea(
#       ..1,
#       set.list = gs_list[c("HAMAI_APOPTOSIS_VIA_TRAIL_UP",
#                            "GSE14350_TREG_VS_TEFF_UP")],
#       n.rand = 1e6,
#       power = 0
#     ) %>% as_tibble(rownames = "gs_name"),
#     .id = "Clusters_contrast"
#   )
#
# gsea_specific2_df <- gsea_specific2_df %>%
#   separate(Clusters_contrast, into = c("Clusters", "contrast"), sep = " / ")
#
# write_tsv(
#   gsea_specific2_df,
#   here("scRNAseq", "data", "processed_data_tables", "03_gsea.genotype_one-v-one_per-cluster_log2FC_specific2.tsv")
# )
