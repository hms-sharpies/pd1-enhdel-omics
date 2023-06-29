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
source(here("aesthetics.R"))
theme_set(theme_classic())

# purpose ------

# To answer the question, "How were off-target effects of CRISPR editing
# tested for and excluded?"

# read in data -----

so <- readRDS(here("scRNAseq", "data", "processed_data_objects", "00_seurat_workflow.so.rds"))

# off target Tier 1 ------

off_target_genes <- c("Cmtr2", "Snx12", "Lrrc23", "Perp", "2900026A02Rik", "Zc3h13", "Vgll4")
off_target_genes_present <- off_target_genes[off_target_genes %in% rownames(so@assays$RNA@counts)]

length(off_target_genes) == length(off_target_genes_present)
# all off_target genes are present

get_metadata_from_so(so, genes = off_target_genes) %>%
  pivot_longer(cols = off_target_genes,
               names_to = "gene",
               values_to = "expr") %>%
  mutate(cluster_ids = factor(str_wrap(cluster_ids, width = 10),
                              levels = unique(str_wrap(
                                names(palette_cluster_ids), width = 10
                              )))) %>%
  ggplot() +
  aes(genotype, expr) +
  geom_violin(aes(fill = genotype),
              scale = "width",
              draw_quantiles = 0.5) +
  facet_grid(gene ~ cluster_ids) +
  stat_compare_means(
    comparisons = combn(c("D", "W", "K"), 2, simplify = FALSE),
    label = "p.signif",
    size = TEXT_SIZE * GGPLOT_TEXT_SCALE_FACTOR
  ) +
  scale_fill_manual(values = palette_genotype) +
  scale_y_continuous(expand = c(0, 0), limits = ~ c(0, ..1[2] * 1.2)) +
  remove_x_spine() +
  theme(
    strip.background = element_blank(),
    legend.position = "none",
    strip.text.y = element_text(
      angle = 0,
      hjust = 0,
      vjust = 0.5
    )
  )

# off target CFD 0.3 ------

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

# chromosome 1 ------

# get genes from chromosome 1
chr1_genes_df <- genes_df %>%
  filter(X1 == "1")

all_genes_in_so <- rownames(so@assays$RNA@counts)
length(all_genes_in_so)
# there are 14,462 genes in the Seurat object

chr1_genes_in_so <- all_genes_in_so[all_genes_in_so %in% chr1_genes_df$gene_name]
length(chr1_genes_in_so)
# there are 859 genes on Chr1 found

chr1_off_target_genes <- off_target_genes[off_target_genes %in% chr1_genes_in_so]
length(chr1_off_target_genes)
# of those, 5 are off target genes

dge_df <- read_tsv(
  here(
    "scRNAseq",
    "data",
    "processed_data_tables",
    "01_dge.genotype_one-v-one_per-cluster.tsv"
  )
) %>% filter(contrast == "D_W")

dge_df %>%
  filter(gene %in% chr1_genes_in_so) %>%
  group_by(Clusters) %>%
  mutate(signed_loglog_p = ifelse(avg_logFC > 0,
                                    log10(-log10(p_val)),
                                    -log10(-log10(p_val)))) %>%
  arrange(signed_loglog_p) %>%
  mutate(rank = row_number()) %>%
  ggplot() +
  aes(rank, signed_loglog_p) +
  # geom_point(size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  # geom_hline(yintercept = c(log10(-log10(0.05)),-log10(-log10(0.05))),
  #            color = "red") +
  # geom_line() +
  geom_point(size = 0.1, alpha = 0.5) +
  geom_point(data = ~..1 %>% filter(gene %in% chr1_off_target_genes),
             size = 2) +
  geom_text_repel(aes(label = gene),
                  data = ~..1 %>% filter(gene %in% chr1_off_target_genes)) +
  facet_wrap(~Clusters, ncol = 5, scale = "free_y") +
  remove_x_spine() +
  theme(strip.background = element_blank(),
        axis.text.x = element_blank()) +
  labs(title = "Genes differentially expressed between EnhDel vs WT on Chr1",
       subtitle = "Labelled genes are predicted off-target",
       y = "(WT) <-- Directional log(-log(p-value)) --> (EnhDel)")

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
