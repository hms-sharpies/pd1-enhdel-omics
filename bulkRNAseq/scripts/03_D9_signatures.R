# 05_D9_signatures.R

library(here)
library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(stringr)
library(tibble)
library(readr)
library(forcats)
library(ggrepel)
library(Matrix.utils)

theme_set(theme_classic())
source(here("helper.R"))

# read in -------

gene_count_D9_df_filtered <- read_tsv(here("bulkRNAseq", "processed_data_tables", "01_preprocessing.gene_count_D9_df_filtered.txt"))%>%
  mutate_at(2:ncol(.), as.integer)
sample_metadata_df <- read_tsv(here("bulkRNAseq", "processed_data_tables", "01_preprocessing.sample_metadata_D9_df.txt"))

gs_df <- read_tsv(here("scRNAseq", "data", "gene_signatures",
                       "02_gene_signatures.msigdbr_H-C2-C7.tsv"))

gs_subset_bp_df <- gs_df %>%
  filter(str_detect(gs_name, "HALLMARK|REACTOME"))
length(unique(gs_subset_bp_df$gs_name))
# 1665
gs_subset_bp_list <- gs_subset_bp_df %>%
  named_group_split(gs_name) %>%
  map(~..1$gene_symbol)

gs_subset_tcell_df <- gs_df %>%
  filter(str_detect(gs_name, "_T_CELL|TCELL|CD8|CD4|EXH")) %>%
  filter(str_detect(gs_name, "HALLMARK|REACTOME", negate = T)) %>%
  filter(str_detect(gs_name, "NEUTROPHIL|MONO|ERYTHROBLAST|DC|NKT|LIN|BCELL|LSK|CD40|CD44|MAST", negate = T))
length(unique(gs_subset_tcell_df$gs_name))
# 1446
gs_subset_tcell_list <- gs_subset_tcell_df %>%
  named_group_split(gs_name) %>%
  map(~..1$gene_symbol)

dge_list <- list()
dge_list[["progenitor"]] <- read_tsv(
  here("bulkRNAseq", "processed_data_tables", "02_D9_analysis.dge_progexh_genotype.Enhdel.v.WT.txt")
) %>% mutate(cluster_id = "progexh")
dge_list[["terminal"]] <- read_tsv(
  here("bulkRNAseq", "processed_data_tables", "02_D9_analysis.dge_termexh_genotype.Enhdel.v.WT.txt")
) %>% mutate(cluster_id = "termexh")
dge_df <- dge_list %>%
  purrr::reduce(bind_rows) %>%
  mutate(signed_neglog10_p = ifelse(log2FoldChange > 0, -log10(pvalue), log10(pvalue))) %>%
  drop_na(log2FoldChange, pvalue)

# run GSEA using liger functions -----

dge_df %>%
  group_by(cluster_id) %>%
  arrange(desc(signed_neglog10_p)) %>%
  mutate(rank = row_number()) %>%
  ggplot() +
  aes(rank, signed_neglog10_p) +
  geom_line(aes(color = cluster_id)) +
  facet_grid(~cluster_id) +
  geom_hline(yintercept = c(-log10(0.05), log10(0.05)),
             linetype = "dotted")

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

rank_lists <- dge_df %>%
  named_group_split(cluster_id) %>%
  map(~make_rank_list(..1, rank_by = signed_neglog10_p, gene = gene))

gsea_res_bp <- map_dfr(
  rank_lists,
  ~liger::iterative.bulk.gsea(
    values = ..1,
    set.list = gs_subset_bp_list,
    power = 0
  ) %>% as_tibble(rownames = "gs_name"),
  .id = "cluster_id"
)
write_tsv(
  gsea_res_bp,
  here("bulkRNAseq", "processed_data_tables", "03_D9_signatures.gsea_bp_D9_genotype.Enhdel.v.WT.txt")
)
gsea_res_bp_wide <- gsea_res_bp %>%
  pivot_wider(id_cols = gs_name,
              names_from = cluster_id,
              values_from = c(p.val, q.val, sscore, edge))
write_tsv(
  gsea_res_bp_wide,
  here("bulkRNAseq", "processed_data_tables", "03_D9_signatures.gsea_bp_wide_D9_genotype.Enhdel.v.WT.txt")
)

gsea_res_tcell <- map_dfr(
  rank_lists,
  ~liger::iterative.bulk.gsea(
    values = ..1,
    set.list = gs_subset_tcell_list,
    power = 0
  ) %>% as_tibble(rownames = "gs_name"),
  .id = "cluster_id"
)
write_tsv(
  gsea_res_tcell,
  here("bulkRNAseq", "processed_data_tables", "03_D9_signatures.gsea_tcell_D9_genotype.Enhdel.v.WT.txt")
)
gsea_res_tcell_wide <- gsea_res_tcell %>%
  pivot_wider(id_cols = gs_name,
              names_from = cluster_id,
              values_from = c(p.val, q.val, sscore, edge))
write_tsv(
  gsea_res_tcell_wide,
  here("bulkRNAseq", "processed_data_tables", "03_D9_signatures.gsea_tcell_wide_D9_genotype.Enhdel.v.WT.txt")
)

