# 00_qc

library(ChIPseeker)
library(here)
library(stringr)
library(ggplot2)
library(purrr)
library(dplyr)
library(readr)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
library(clusterProfiler)

bedfiles <- list.files(here("ATACseq", "data", "peaks")) %>%
  str_subset("bed") %>%
  str_subset("summits") %>%
  here("ATACseq", "data", "peaks", .)
bedfile.labels <- gsub("_0.001.summits.bed", "", basename(bedfiles))

peaks <- map(bedfiles %>% set_names(bedfile.labels), readPeakFile)
promoter <- getPromoters(TxDb=txdb, upstream=1000, downstream=1000)
tagMatrixList <- map(peaks, ~getTagMatrix(..1, windows = promoter))
tss_plots <- tagMatrixList %>%
  imap(~{
    plotAvgProf(..1, xlim=c(-1000, 1000),
                xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency") +
      ggtitle(..2)
  })

tss_df <- map_dfr(tss_plots, ~..1$data, .id = "sample") %>%
  arrange(sample)
write_tsv(tss_df, here("ATACseq", "data", "processed_data_tables", "00_qc.tss.tsv"))
