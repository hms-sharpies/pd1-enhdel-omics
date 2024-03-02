# 01_preprocessing

library(here)
library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)
library(stringr)
library(tibble)
library(readr)
library(ggrepel)
library(Matrix.utils)
library(DESeq2)
library(Rsc)

set.seed(42)
theme_set(theme_classic() + theme(
  strip.background = element_blank()
))

raw_count_df <- read_tsv(here("bulkRNAseq", "data", "transcript_counts_all_samples.tsv"))
nrow(raw_count_df)

sample_cols <- raw_count_df %>%
  colnames() %>%
  .[c(-1, -2)]

# convert ensembl names to gene symbols ------

gtf_df <- read_tsv(here("bulkRNAseq", "data", "Mus_musculus.GRCm39.111.gtf.gz"),
                   skip = 5, col_names = FALSE)
gene_id_to_gene_name_df <- gtf_df %>%
  mutate(transcript_id = str_extract(X9, "ENSMUST[:digit:]*")) %>%
  mutate(gene_name = str_extract(X9, "gene_name \"[:alnum:]*") %>%
           str_remove("gene_name \"")) %>%
  dplyr::select(transcript_id, gene_name) %>%
  drop_na(transcript_id) %>%
  unique()

count_df <- raw_count_df %>%
  mutate(transcript_id = str_remove(target_id, "\\..*$")) %>%
  dplyr::relocate(transcript_id, .after = target_id) %>%
  left_join(gene_id_to_gene_name_df, by = "transcript_id") %>%
  dplyr::relocate(gene_name, .after = transcript_id)

# NA transcripts? -----

# there are 150 transcripts with no associated gene
count_df %>%
  filter(is.na(gene_name))

# make gene_count_df ------

# we'll drop NA genes

gene_count_df <- count_df %>%
  drop_na(gene_name) %>% # drop na genes
  group_by(gene_name) %>%
  summarise_at(sample_cols, sum)

# filter out genes present in only one sample --------

dim(gene_count_df)
# 34065 37

gene_count_tidy <- gene_count_df %>%
  pivot_longer(cols = all_of(sample_cols),
               names_to = "sample",
               values_to = "count")

rare_genes <- gene_count_tidy %>%
  group_by(gene_name) %>%
  summarise(nneg_n_samples = sum(count > 0)) %>%
  filter(nneg_n_samples <= 1) %>%
  .$gene_name
length(rare_genes)
# 9904

gene_count_df_filtered <- gene_count_df %>%
  filter(!(gene_name %in% rare_genes))
dim(gene_count_df_filtered)
#  24161    37

gene_count_tidy_filtered <- gene_count_df_filtered %>%
  pivot_longer(cols = all_of(sample_cols),
               names_to = "sample",
               values_to = "count")

# PCA ------

pca_out <- gene_count_df_filtered %>%
  column_to_rownames("gene_name") %>%
  t() %>%
  prcomp(scale = T)

pca_df <- pca_out$x %>%
  as_tibble(rownames = "sample") %>%
  separate(sample,
           into = c("day", "mice", "genotype", "exh_subset", "sample"),
           sep = "_",
           remove = FALSE)
pca_df %>%
  ggplot() +
  aes(PC1, PC2) +
  geom_point(aes(color = paste(exh_subset, day), shape = genotype),
             size = 3) +
  scale_color_manual(values = c("red4", "red1",
                                "mediumpurple4", "mediumpurple1",
                                "aquamarine4","aquamarine1")) +
  geom_text_repel(aes(label = sample))

# from the PCA, it makes sense to drop S25, which is an outlier.
# not sure what's up with S24, so we'll keep it for now.

gene_count_df_filtered2 <- gene_count_df_filtered %>%
  dplyr::select(-str_subset(colnames(.), "S25"))

# investigate PC2 loadings -------

# pca_out$rotation %>%
#   as_tibble(rownames = "gene") %>%
#   arrange(desc(PC2)) %>%
#   head(20) %>%
#   View()

saveRDS(pca_out, here("bulkRNAseq", "results", "raw_pca_out.rds"))

# Number of genes expressed? -------

gene_count_tidy_filtered2 <- gene_count_df_filtered2 %>%
  pivot_longer(cols = 2:ncol(.),
               names_to = "sample",
               values_to = "count")

gene_count_tidy_filtered2 %>%
  ggplot() +
  aes(sample, log10(count + 1)) +
  geom_violin(scale = "width", fill = "grey") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

# filter gene count for D9 and D30 separately -----

gene_count_D30_df <- gene_count_df_filtered2 %>%
  dplyr::select(gene_name, str_subset(colnames(.), "D30"))

gene_count_D30_tidy <- gene_count_D30_df %>%
  pivot_longer(cols = 2:ncol(.),
               names_to = "sample",
               values_to = "count")
rare_genes_D30 <- gene_count_D30_tidy %>%
  group_by(gene_name) %>%
  summarise(nneg_n_samples = sum(count > 0)) %>%
  filter(nneg_n_samples <= 1) %>%
  .$gene_name
length(rare_genes_D30)
# 5178

gene_count_D30_df_filtered <- gene_count_D30_df %>%
  filter(!(gene_name %in% rare_genes_D30)) %>%
  mutate_at(2:ncol(.), as.integer)
nrow(gene_count_D30_df_filtered)
# 18983


gene_count_D9_df <- gene_count_df_filtered2 %>%
  dplyr::select(gene_name, str_subset(colnames(.), "D9"))

gene_count_D9_tidy <- gene_count_D9_df %>%
  pivot_longer(cols = 2:ncol(.),
               names_to = "sample",
               values_to = "count")
rare_genes_D9 <- gene_count_D9_tidy %>%
  group_by(gene_name) %>%
  summarise(nneg_n_samples = sum(count > 0)) %>%
  filter(nneg_n_samples <= 1) %>%
  .$gene_name
length(rare_genes_D9)
# 906

gene_count_D9_df_filtered <- gene_count_D9_df %>%
  filter(!(gene_name %in% rare_genes_D9)) %>%
  mutate_at(2:ncol(.), as.integer)
nrow(gene_count_D9_df_filtered)
# 23255

# write gene count tables -----

saveRDS(gene_count_D30_df_filtered, here("bulkRNAseq", "results", "gene_count_D30_df_filtered.rds"))
saveRDS(gene_count_D9_df_filtered, here("bulkRNAseq", "results", "gene_count_D9_df_filtered.rds"))

sample_metadata_D30_df <-
  tibble(sample = colnames(gene_count_D30_df_filtered)[2:ncol(gene_count_D30_df_filtered)]) %>%
  separate(sample,
           into = c("day", "mice", "genotype", "exh_subset", "sample_num"),
           sep = "_",
           remove = FALSE)
write_tsv(sample_metadata_D30_df, here("bulkRNAseq", "results", "sample_metadata_D30_df.txt"))

sample_metadata_D9_df <-
  tibble(sample = colnames(gene_count_D9_df_filtered)[2:ncol(gene_count_D9_df_filtered)]) %>%
  separate(sample,
           into = c("day", "mice", "genotype", "exh_subset", "sample_num"),
           sep = "_",
           remove = FALSE)
write_tsv(sample_metadata_D9_df, here("bulkRNAseq", "results", "sample_metadata_D9_df.txt"))

# Also remove sample 24, re-filter, re-save D30 -------

gene_count_df_filtered3 <- gene_count_df_filtered2 %>%
  dplyr::select(-str_subset(colnames(.), "S24"))

gene_count_D30_df <- gene_count_df_filtered3 %>%
  dplyr::select(gene_name, str_subset(colnames(.), "D30"))

gene_count_D30_tidy <- gene_count_D30_df %>%
  pivot_longer(cols = 2:ncol(.),
               names_to = "sample",
               values_to = "count")
rare_genes_D30 <- gene_count_D30_tidy %>%
  group_by(gene_name) %>%
  summarise(nneg_n_samples = sum(count > 0)) %>%
  filter(nneg_n_samples <= 1) %>%
  .$gene_name
length(rare_genes_D30)
# 5329

gene_count_D30_df_filtered <- gene_count_D30_df %>%
  filter(!(gene_name %in% rare_genes_D30)) %>%
  mutate_at(2:ncol(.), as.integer)
nrow(gene_count_D30_df_filtered)
# 18832

saveRDS(gene_count_D30_df_filtered, here("bulkRNAseq", "results", "gene_count_D30_df_filtered_removeS24.rds"))

sample_metadata_D30_df <-
  tibble(sample = colnames(gene_count_D30_df_filtered)[2:ncol(gene_count_D30_df_filtered)]) %>%
  separate(sample,
           into = c("day", "mice", "genotype", "exh_subset", "sample_num"),
           sep = "_",
           remove = FALSE)
write_tsv(sample_metadata_D30_df, here("bulkRNAseq", "results", "sample_metadata_D30_df_removeS24.txt"))


# Redo PCA and save loadings ---------

gene_count_df_filtered3_tidy <- gene_count_df_filtered3 %>%
  pivot_longer(cols = 2:ncol(.),
               names_to = "sample",
               values_to = "count")
rare_genes <- gene_count_df_filtered3_tidy %>%
  group_by(gene_name) %>%
  summarise(nneg_n_samples = sum(count > 0)) %>%
  filter(nneg_n_samples <= 1) %>%
  .$gene_name
length(rare_genes)
# 84

gene_count_df_filtered4 <- gene_count_df_filtered3 %>%
  filter(!(gene_name %in% rare_genes)) %>%
  mutate_at(2:ncol(.), as.integer)
nrow(gene_count_df_filtered4)
# 24077
#
# pca_out <- gene_count_df_filtered4 %>%
#   column_to_rownames("gene_name") %>%
#   t() %>%
#   prcomp(scale = T)
#
# pca_df <- pca_out$x %>%
#   as_tibble(rownames = "sample") %>%
#   separate(sample,
#            into = c("day", "mice", "genotype", "exh_subset", "sample"),
#            sep = "_",
#            remove = FALSE)
# pca_df %>%
#   ggplot() +
#   aes(PC1, PC2) +
#   geom_point(aes(color = paste(exh_subset, day), shape = genotype),
#              size = 3) +
#   scale_color_manual(values = c("red4", "red1",
#                                 "mediumpurple4", "mediumpurple1",
#                                 "aquamarine4","aquamarine1")) +
#   geom_text_repel(aes(label = sample))
