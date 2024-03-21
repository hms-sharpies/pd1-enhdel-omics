# 02_D9_analysis.R

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

# read in -------

gene_count_D9_df_filtered <- read_tsv(here("bulkRNAseq", "processed_data_tables", "01_preprocessing.gene_count_D9_df_filtered.txt")) %>%
  mutate_at(2:ncol(.), as.integer)
sample_metadata_df <- read_tsv(here("bulkRNAseq", "processed_data_tables", "01_preprocessing.sample_metadata_D9_df.txt"))

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

norm_counts <- as.data.frame(counts(de, normalized = TRUE))
write_tsv(norm_counts, here("bulkRNAseq", "processed_data_tables", "02_D9_analysis.deseq_norm_counts_exh.prog.vs.term.txt"))

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

write_tsv(dge_exh.prog.vs.term, here("bulkRNAseq", "processed_data_tables", "02_D9_analysis.dge_exh.prog.vs.term.txt"))

# Within subset, enhdel vs WT -------

counts_mtx_list <- list()
counts_mtx_list[["progenitor"]] <- gene_count_D9_df_filtered %>%
  dplyr::select(gene_name, str_subset(colnames(.), "Slamf6")) %>%
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

norm_counts_df <- de_list %>%
  map(counts, normalized = TRUE) %>%
  map(as.data.frame) %>%
  purrr::reduce(bind_cols)
write_tsv(norm_counts_df, here("bulkRNAseq", "processed_data_tables", "02_D9_analysis.deseq_norm_counts_all_clusters_genotype.Enhdel.v.WT.txt"))

write_tsv(results_list$progenitor,
          here("bulkRNAseq", "processed_data_tables", "02_D9_analysis.dge_progexh_genotype.Enhdel.v.WT.txt"))
write_tsv(results_list$terminal,
          here("bulkRNAseq", "processed_data_tables", "02_D9_analysis.dge_termexh_genotype.Enhdel.v.WT.txt"))

