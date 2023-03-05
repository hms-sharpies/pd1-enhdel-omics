# 08_off_target_analysis

library(here)
library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(stringr)
library(forcats) # for fct_reorder
library(tibble)
library(readr)
library(Seurat)
library(SeuratObject)
library(ggrepel)

source(here("helper.R"))
theme_set(theme_classic())

# purpose ------

# To answer the question, "How were off-target effects of CRISPR editing
# tested for and excluded?"

# read in data -----

so <- readRDS(here("scRNAseq", "data", "processed_data_objects", "00_seurat_workflow.so.rds"))

genes_df <- read_tsv(here("scRNAseq", "data", "genes.gtf"), skip = 5, col_names = FALSE) %>%
  filter(X3 == "gene") %>%
  separate(X9, into = c("gene_id", "gene_version", "gene_name", "gene_source", "gene_biotype"),
           sep = "; ") %>%
  mutate(gene_name = str_sub(gene_name, start = 12, end = -2)) %>%
  dplyr::relocate(gene_name, .before = X1)

off_target_genes <-
  read_tsv(here("scRNAseq", "data", "off_target_CFD_0.3.txt"),
           col_names = FALSE)$X1

off_target_genes_present <- off_target_genes[off_target_genes %in% rownames(so@assays$RNA@counts)]

off_target_genes %>% length()
# there are 70 genes
off_target_genes_present %>% length()
# 37 of the 70 are present in the data
37 / 70
# about 53%

gene_expr_tidy <- get_metadata_from_so(so, genes = off_target_genes_present) %>%
  pivot_longer(cols = off_target_genes_present,
               names_to = "gene")

zscore_mean_df <- gene_expr_tidy %>%
  group_by(orig.ident, gene) %>%
  summarise(mean_value = mean(value)) %>%
  group_by(gene) %>%
  mutate(zscore_mean_value = as.vector(scale(mean_value))) %>%
  ungroup()
zscore_mean_df %>%
  ggplot() +
  aes(orig.ident, fct_reorder(gene, mean_value), fill = mean_value) +
  geom_tile() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_viridis()

# permutation testing -----

sample_n_genes <- function(so, n, seed = 42) {
  set.seed(seed)
  so@assays$RNA@counts %>%
    rownames() %>%
    sample(size = n)
}

get_lm_est <- function(so, seed = 42) {
  sampled_genes <- sample_n_genes(so, n = 37, seed = seed)
  cor_tbl <- get_metadata_from_so(so, genes = sampled_genes) %>%
    pivot_longer(cols = sampled_genes,
                 names_to = "gene") %>%
    group_by(orig.ident, gene) %>%
    summarise(mean_value = mean(value)) %>%
    group_by(gene) %>%
    mutate(zscore_mean_value = as.vector(scale(mean_value))) %>%
    ungroup() %>%
    mutate(genotype = str_extract(orig.ident, ".")) %>%
    group_by(genotype, gene) %>%
    summarise(avg_mean_value = mean(mean_value)) %>%
    ungroup() %>%
    pivot_wider(id_cols = gene,
                names_from = genotype,
                values_from = avg_mean_value) %>%
    lm(W ~ D - 1, data = .) %>%
    broom::tidy()
}

lm_est_out <- map_dfr(1:100, ~get_lm_est(so, seed = ..1))
lm_est_out2 <- map_dfr(101:200, ~get_lm_est(so, seed = ..1))

zscore_mean_df %>%
  mutate(genotype = str_extract(orig.ident, ".")) %>%
  group_by(genotype, gene) %>%
  summarise(avg_mean_value = mean(mean_value)) %>%
  ungroup() %>%
  pivot_wider(id_cols = gene,
              names_from = genotype,
              values_from = avg_mean_value) %>%
  ggplot() +
  aes(D, W) +
  geom_point() +
  stat_regline_equation() +
  stat_smooth(method = "lm") +
  labs(x = "Expression in EnhDel", y = "Expression in WT",
       title = "Expression of predicted off-target genes")

zscore_mean_df %>%
  mutate(genotype = str_extract(orig.ident, ".")) %>%
  group_by(genotype, gene) %>%
  summarise(avg_mean_value = mean(mean_value)) %>%
  ungroup() %>%
  pivot_wider(id_cols = gene,
              names_from = genotype,
              values_from = avg_mean_value) %>%
  lm(W ~ D - 1, data = .) %>%
  broom::tidy()

lm_est_out_TOTAL <-
  lm_est_out %>%
  bind_rows(lm_est_out2)
lm_est_out_TOTAL %>%
  filter(term == "D") %>%
  ggplot() +
  aes(estimate) +
  geom_histogram(bins = 15) +
  annotate("segment",
           y = 60, yend = 60,
           x = 1.03, xend = 1.03 + 0.1*0.3,
           arrow = arrow()) +
  annotate("text",
           y = 60,
           x = 1.03 + 0.1*0.3,
           label = "30% of simulations",
           hjust = 0) +
  annotate("segment",
           y = 60, yend = 60,
           x = 1.03, xend = 1.03 - 0.1*0.7,
           arrow = arrow()) +
  annotate("text",
           y = 60,
           x = 1.03 - 0.1*0.7,
           label = "70% of simulations",
           hjust = 1) +
  geom_vline(xintercept = 1.03, color = "red") +
  labs(x = "Estimate of EnhDel Effect", y = "Frequency",
       subtitle = "Red line denotes estimate of predicted off-target genes") +
  ggtitle("Distribution from random samplings of 37 genes (200 trials)")

sum(lm_est_out_TOTAL$estimate < 1.03) / 200
sum(lm_est_out_TOTAL$estimate > 1.03) / 200
