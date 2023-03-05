# 07_motif_analysis

library(here)
library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(stringr)
library(tibble)
library(readr)
library(Seurat)
library(SeuratObject)
library(ggrepel)

source(here("helper.R"))
theme_set(theme_classic())

# purpose ------

# To answer the question, "What are the transcriptional regulators of
# the exhaustion associated enhancer and how are they expressed in the
# different subpopulations?"

# read in data -----

so <- readRDS(here("scRNAseq", "data", "processed_data_objects", "00_seurat_workflow.so.rds"))
exh_enh_test_fimo <- read_tsv(here("scRNAseq", "data", "exh_enh_test_fimo.tsv")) %>%
  filter(!is.na(sequence_name)) %>%
  mutate(motif_gene = motif_alt_id %>% str_extract("\\.[:alnum:]*$")) %>%
  mutate(motif_gene = str_sub(motif_gene, start = 2)) %>%
  dplyr::relocate(motif_gene, .before = motif_id)

genes_df <- read_tsv(here("scRNAseq", "data", "genes.gtf"), skip = 5, col_names = FALSE) %>%
  filter(X3 == "gene") %>%
  separate(X9, into = c("gene_id", "gene_version", "gene_name", "gene_source", "gene_biotype"),
           sep = "; ") %>%
  mutate(gene_name = str_sub(gene_name, start = 12, end = -2)) %>%
  dplyr::relocate(gene_name, .before = X1)

# some sanity checks  ------

genes_df$gene_name %>% unique() %>% length()
# there are 31017 unique genes from the mouse transcriptome

rownames(so@assays$RNA@counts) %>% length()
# there are 14462 genes found in the Seurat object

genes_in_gtf <- rownames(so@assays$RNA@counts) %in% (unique(genes_df$gene_name))
sum(genes_in_gtf)
# for some reason, only 13923 genes from the so are in the gtf file??? Seems sus

rownames(so@assays$RNA@counts)[!(genes_in_gtf)]
# Not sure which genes these are...
# Seems like a pretty random assortment of genes to me. Perhaps had to do with
# that initial filtering step where I restricted to only things labelled "gene"
# Maybe I'm using the wrong gtf file?
# It's probably good enough for now so I'm going to ignore it



# genes present -----

intersect(unique(genes_df$gene_name), exh_enh_test_fimo$motif_gene) %>% length()
# There are 39 genes that the motif analysis found that are in the mouse genome

intersect(rownames(so@assays$RNA), exh_enh_test_fimo$motif_gene) %>% length()
# There are 22 genes that the motif analysis found in the Seurat object

# genes identified by motif analysis -------

exh_enh_test_fimo_sig <- exh_enh_test_fimo %>%
  group_by(motif_gene) %>%
  summarise(at_least_1_sig = max(`q-value`) < 0.05)

exh_enh_test_fimo %>%
  filter(`q-value` < 0.05) %>%
  filter(motif_gene %in% unique(genes_df$gene_name)) %>%
  .$motif_gene %>% unique() %>% length()
# 16 significant of the 39 identified.

exh_enh_test_fimo %>%
  filter(`q-value` < 0.05) %>%
  filter(motif_gene %in% rownames(so@assays$RNA)) %>%
  .$motif_gene %>% unique() %>% length()
# There are 8 significant genes with gene expression.


# gene expression -------

motif_genes <- intersect(rownames(so@assays$RNA), exh_enh_test_fimo$motif_gene)

get_metadata_from_so(so, genes = motif_genes) %>%
  filter(genotype == "W") %>%
  pivot_longer(cols = motif_genes,
               names_to = "gene") %>%
  plot_violin(group = cluster_ids, var = value) +
  facet_grid(gene~.) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        strip.background = element_blank(),
        strip.text = element_text(angle = 0))


mean_expr_df <- get_metadata_from_so(so, genes = motif_genes) %>%
  filter(genotype == "W") %>%
  pivot_longer(cols = motif_genes,
               names_to = "gene") %>%
  group_by(cluster_ids, gene) %>%
  summarise(mean_expr = mean(value)) %>%
  group_by(gene) %>%
  mutate(zscore_mean_expr = as.vector(scale(mean_expr)),
         scale_max_mean_expr = mean_expr / max(mean_expr))

mean_expr_df %>%
  ggplot() +
  aes(cluster_ids, gene, fill = zscore_mean_expr) +
  geom_tile() +
  geom_text(aes(label = format(mean_expr, digits = 2))) +
  scale_fill_viridis_c() +
  remove_spines() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

mean_expr_df %>%
  ggplot() +
  aes(cluster_ids, gene, fill = scale_max_mean_expr) +
  geom_tile() +
  geom_text(aes(label = format(mean_expr, digits = 2))) +
  scale_fill_viridis_c() +
  remove_spines() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


mean_expr_df %>%
  left_join(exh_enh_test_fimo_sig, by = c("gene" = "motif_gene")) %>%
  filter(at_least_1_sig) %>%
  ggplot() +
  aes(cluster_ids, gene, fill = zscore_mean_expr) +
  geom_tile() +
  geom_text(aes(label = format(mean_expr, digits = 2))) +
  scale_fill_viridis_c() +
  remove_spines() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_text(hjust = 0.5)) +
  labs(x = "Cluster", y = "gene", fill = "Z-score\nmean expr.")

mean_expr_df %>%
  ungroup() %>%
  # filter(cluster_ids != "dividing") %>%
  left_join(exh_enh_test_fimo_sig, by = c("gene" = "motif_gene")) %>%
  ggplot() +
  aes(cluster_ids, zscore_mean_expr) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_line(aes(group = gene, color = at_least_1_sig), alpha = 0.3) +
  geom_text_repel(
    aes(x = 1, y = zscore_mean_expr, label = gene),
    hjust = 1,
    nudge_x = -1,
    color = "black",
    direction = "y",
    data = ~ ..1 %>%
      filter(cluster_ids == "progenitor exhausted")
  ) +
  geom_point(aes(size = mean_expr, color = at_least_1_sig)) +
  scale_color_manual(values = c("black", "red")) +
  remove_x_spine() +
  ggtitle("Relative Expression of EnhExh-associated TFs in Exh. Subsets") +
  labs(x = "Clusters", y = "Relative Mean Expr.", size = "Mean Expr.",
       subtitle = "In WT CD8+ T cells only")


mean_expr_df %>%
  ungroup() %>%
  # filter(cluster_ids != "dividing") %>%
  left_join(exh_enh_test_fimo_sig, by = c("gene" = "motif_gene")) %>%
  filter(at_least_1_sig) %>%
  ggplot() +
  aes(cluster_ids, zscore_mean_expr) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_line(aes(group = gene), alpha = 0.3) +
  geom_text_repel(
    aes(x = 1, y = zscore_mean_expr, label = gene),
    hjust = 1,
    nudge_x = -1,
    color = "black",
    direction = "y",
    data = ~ ..1 %>%
      filter(cluster_ids == "progenitor exhausted")
  ) +
  geom_point(aes(size = mean_expr)) +
  scale_color_manual(values = c("black", "red")) +
  remove_x_spine() +
  ggtitle("Relative Expression of EnhExh-associated TFs in Exh. Subsets") +
  labs(x = "Clusters", y = "Relative Mean Expr.", size = "Mean Expr.",
       subtitle = "In WT CD8+ T cells only")


# TFs from Sen et al. Science 2015 CRISPR tiling screen ------

# relevant TFs are SOX3, TBX21, RAR
# associated genes should be Sox3, Tbx21, Rara, Rarg

get_metadata_from_so(so, genes = c("Tbx21", "Rara", "Rarg")) %>%
  filter(genotype == "W") %>%
  pivot_longer(cols = c("Tbx21", "Rara", "Rarg"),
               names_to = "gene") %>%
  plot_violin(group = cluster_ids, var = value) +
  facet_grid(gene~.) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

