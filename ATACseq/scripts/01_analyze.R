# 01_analyze

library(here)
library(dplyr)
library(readr)
library(purrr)
library(stringr)
library(tibble)
library(tidyr)

library(sva)
library(DESeq2)

# create count_tbl -----

counts_dir <- here("ATACseq", "data", "counts")
counts_files <- list.files(counts_dir) %>%
  file.path(counts_dir, .) %>%
  set_names(str_extract(., "[:alpha:]*[:digit:]*_[:alpha:]*\\.") %>% str_sub(end = -2))
count_list <- imap(counts_files, ~{
  count_df <- read_tsv(..1, col_names = FALSE)
  colnames(count_df) <- c("Chr", "Start", "End", ..2)
  return(count_df)
})
count_tbl <- purrr::reduce(count_list, left_join, by = c("Chr", "Start", "End"))
count_df <- count_tbl %>%
  mutate(peak_id = paste0(Chr, ":", Start, "-", End)) %>%
  select(-Chr, -Start, -End) %>%
  column_to_rownames("peak_id")
count_mtx <- as.matrix(count_df)

# create condition_df ------

condition_df <- data.frame(sample = colnames(count_mtx)) %>%
  mutate(cell_type = str_extract(sample, "_[:alpha:]*") %>% str_sub(start = 2),
         genotype = str_extract(sample, "^[:alpha:]"),
         replicate = str_extract(sample, "\\d\\d|\\d"),
         group = paste0(genotype, "_", cell_type),
         batch = case_when(
           replicate == "1" | replicate == "2" ~ "batch 1",
           str_detect(replicate, "3|4") ~ "batch 2"
         ))

# run ComBatSeq -----

batch <- condition_df$batch
adjusted <- ComBat_seq(count_mtx, batch=batch, group = NULL)

# filter peaks with less than 5 reads across all conditions -----

count_mtx_adj <- as.matrix(adjusted[(apply(adjusted, 1, max) > 5),])

# deseq -----

dedsATAC <- DESeqDataSetFromMatrix(countData = count_mtx_adj,
                                   colData = condition_df,
                                   design = ~ group)
deATAC <- DESeq(dedsATAC, modelMatrixType = "standard")

dfSizeFactors <- as.data.frame(sizeFactors(deATAC))
write.table(dfSizeFactors, here("ATACseq", "data", "processed_data_tables",
                                "01_analyze.SizeFactors.txt"))

dfNormalizedCounts <- counts(deATAC, normalized = TRUE)
write.table(dfNormalizedCounts, here("ATACseq", "data", "processed_data_tables",
                                     "01_analyze.NormCountsATAC.txt"))

# results ------

results_W_S.v.CN <- results(deATAC, contrast = c("group", "W_S", "W_CN")) %>%
  as.data.frame() %>%
  mutate(p_BH = p.adjust(pvalue, method = "BH")) %>%
  as_tibble(rownames = "peak_id")
write_tsv(results_W_S.v.CN, here("ATACseq", "data", "processed_data_tables",
                                 "01_analyze.deseq_results_W_S.v.CN.csv"))


results_list <- cross2(list("S", "CN"), combn(c("D", "K", "W"), m = 2, simplify = FALSE)) %>%
  {set_names(., nm = map(., ~paste0(..1[[2]][1], ".v.", ..1[[2]][2], "_", ..1[[1]])))} %>%
  map(function(cross) {
    cluster <- cross[[1]]
    genotype1 <- cross[[2]][1]; genotype2 <- cross[[2]][2]
    group1 <- paste0(genotype1, "_", cluster)
    group2 <- paste0(genotype2, "_", cluster)
    results(deATAC, contrast = c("group", group1, group2))
  })

results_tidy <- results_list %>%
  map(as.data.frame) %>%
  map(mutate, p_BH = p.adjust(pvalue, method = "BH")) %>%
  map_dfr(as_tibble, rownames = "peak_id", .id = "comparison")
write_tsv(results_tidy, here("ATACseq", "data", "processed_data_tables",
                             "01_analyze.deseq_results_tidy.tsv"))

results_megatable <- results_tidy %>%
  pivot_wider(id_cols = peak_id,
              names_from = comparison,
              values_from = c(log2FoldChange, pvalue, padj, p_BH),
              names_sort = FALSE,
              names_vary = "slowest")
write_tsv(results_megatable, here("ATACseq", "data", "processed_data_tables",
                                  "01_analyze.deseq_results_megatable.tsv"))
