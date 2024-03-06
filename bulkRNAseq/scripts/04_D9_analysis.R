# 04_D9_analysis.R

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
library(scattermore)

theme_set(theme_classic())

# read in -------

gene_count_D9_df_filtered <- readRDS(here("bulkRNAseq", "results", "gene_count_D9_df_filtered.rds"))
sample_metadata_df <- read_tsv(here("bulkRNAseq", "results", "sample_metadata_D9_df.txt"))

# WT, Prog vs Term ---------

counts_mtx <- gene_count_D9_df_filtered %>%
  dplyr::select(gene_name, str_subset(colnames(.), "CD4512")) %>%
  column_to_rownames("gene_name") %>%
  as.matrix()
de <- DESeqDataSetFromMatrix(counts_mtx,
                             colData = sample_metadata_df %>%
                               filter(str_detect(sample, "CD4512")) %>%
                               column_to_rownames("sample"),
                             design = formula("~ exh_subset"))
de <- DESeq2::DESeq(de)
dge_exh.prog.vs.term <-
  results(de, contrast = c("exh_subset",
                           "Slamf6",
                           "CX3CR6Neg")) %>%
  as_tibble(rownames = "gene")

# dge_exh.prog.vs.term %>%
#   filter(padj < 0.05) %>%
#   arrange(log2FoldChange) %>%
#   View()

dge_exh.prog.vs.term %>%
  ggplot() +
  aes(log2FoldChange, -log10(padj)) +
  geom_scattermore() +
  geom_point(size = 3,
             data = ~..1 %>% filter(gene %in% c("Slamf6", "Havcr2"))) +
  geom_text_repel(aes(label = gene),
                  data = ~..1 %>% filter(gene %in% c("Havcr2", "Slamf6")))
# okay great, slamf6 and havcr2 appear in the correct directions here

write_tsv(dge_exh.prog.vs.term,
          here("bulkRNAseq", "results", "dge_D9_exh.prog.vs.term.txt"))


# EnhDel, Prog vs Term ---------

counts_mtx <- gene_count_D9_df_filtered %>%
  dplyr::select(gene_name, str_subset(colnames(.), "CD451_")) %>%
  column_to_rownames("gene_name") %>%
  as.matrix()
de <- DESeqDataSetFromMatrix(counts_mtx,
                             colData = sample_metadata_df %>%
                               filter(str_detect(sample, "CD451_")) %>%
                               column_to_rownames("sample"),
                             design = formula("~ exh_subset"))
de <- DESeq2::DESeq(de)
dge_exh.prog.vs.term <-
  results(de, contrast = c("exh_subset",
                           "Slamf6",
                           "CX3CR6Neg")) %>%
  as_tibble(rownames = "gene")

# dge_exh.prog.vs.term %>%
#   filter(padj < 0.05) %>%
#   arrange(log2FoldChange) %>%
#   View()

dge_exh.prog.vs.term %>%
  ggplot() +
  aes(log2FoldChange, -log10(padj)) +
  geom_scattermore() +
  geom_point(size = 3,
             data = ~..1 %>% filter(gene %in% c("Slamf6", "Havcr2"))) +
  geom_text_repel(aes(label = gene),
                  data = ~..1 %>% filter(gene %in% c("Havcr2", "Slamf6")))

write_tsv(dge_exh.prog.vs.term,
          here("bulkRNAseq", "results", "dge_D9_enhdel_exh.prog.vs.term.txt"))

# Within subset, enhdel vs WT -------

counts_mtx_list <- list()
counts_mtx_list[["progenitor"]] <- gene_count_D9_df_filtered %>%
  dplyr::select(gene_name, str_subset(colnames(.), "Slamf6")) %>%
  column_to_rownames("gene_name") %>%
  as.matrix()
counts_mtx_list[["transitory"]] <- gene_count_D9_df_filtered %>%
  dplyr::select(gene_name, str_subset(colnames(.), "CX3CR6Pos")) %>%
  column_to_rownames("gene_name") %>%
  as.matrix()
counts_mtx_list[["terminal"]] <- gene_count_D9_df_filtered %>%
  dplyr::select(gene_name, str_subset(colnames(.), "CX3CR6Neg")) %>%
  column_to_rownames("gene_name") %>%
  as.matrix()

de_list <- counts_mtx_list %>%
  map(~{
    DESeqDataSetFromMatrix(..1,
                           colData = sample_metadata_df %>%
                             filter(sample %in% colnames(..1)) %>%
                             column_to_rownames("sample"),
                           design = formula("~ genotype"))
  })
de_list <- map(de_list, DESeq2::DESeq)
results_list <- map(de_list, ~{
  results(..1, name = "genotype_CD4512_vs_CD451") %>%
    as_tibble(rownames = "gene")
})

results_list$progenitor %>%
  # filter(gene %in% c("Pdcd1", "Ifit1", "Cxcl10", "Plac8"))
  ggplot() +
  aes(log2FoldChange, -log10(padj)) +
  geom_scattermore() +
  geom_point(size = 2,
             data = ~..1 %>% filter(gene %in% c("Pdcd1", "Ifit1", "Cxcl10",
                                                "Plac8",
                                                "Itgb7"))) +
  geom_text_repel(aes(label = gene),
                  data = ~..1 %>% filter(gene %in% c("Pdcd1", "Ifit1", "Cxcl10",
                                                     "Plac8", "Itgb7"))) +
  labs(caption = "Pdcd1 not significant")


results_list$transitory %>%
  # filter(gene %in% c("Pdcd1", "Plac8", "Cxcl10"))
  ggplot() +
  aes(log2FoldChange, -log10(padj)) +
  geom_scattermore() +
  geom_point(size = 2,
             data = ~..1 %>% filter(gene %in% c("Pdcd1", "Plac8", "Cxcl10"))) +
  geom_text_repel(aes(label = gene),
                  data = ~..1 %>% filter(gene %in% c("Pdcd1", "Plac8", "Cxcl10")))


results_list$terminal %>%
  # filter(gene %in% c("Pdcd1", "Gzmb", "Plac8", "Klrd1", "Jun"))
  ggplot() +
  aes(log2FoldChange, -log10(padj)) +
  geom_scattermore() +
  geom_point(size = 2,
             data = ~..1 %>% filter(gene %in% c("Pdcd1", "Gzmb", "Plac8",
                                                "Klrd1", "Jun"))) +
  geom_text_repel(aes(label = gene),
                  data = ~..1 %>% filter(gene %in% c("Pdcd1", "Gzmb", "Plac8",
                                                     "Klrd1", "Jun")))

write_tsv(results_list$progenitor,
          here("bulkRNAseq", "results", "dge_D9_progexh_genotype.Enhdel.v.WT.txt"))
write_tsv(results_list$transitory,
          here("bulkRNAseq", "results", "dge_D9_transexh_genotype.Enhdel.v.WT.txt"))
write_tsv(results_list$terminal,
          here("bulkRNAseq", "results", "dge_D9_termexh_genotype.Enhdel.v.WT.txt"))

# results_list$terminal%>%
#   filter(padj < 0.05) %>%
#   View()

# read_tsv(here("bulkRNAseq", "results", "dge_D9_progexh_genotype.Enhdel.v.WT.txt")) %>%
#   filter(padj < 0.05) %>%
  # filter(gene %>% str_detect("Actb")) %>%
  # View()
