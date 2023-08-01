library(biomaRt)
hmart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
mmart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")

so <- readRDS(here("scRNAseq", "data", "processed_data_objects", "00_seurat_workflow.so.rds"))
mgenes <- rownames(so)

mouse_human_genes <- getLDS(attributes = c("mgi_symbol"), mart = mmart, values = mgenes,
                            attributesL = c("hgnc_symbol"), martL = hmart,
                            verbose = TRUE)

saveRDS(mouse_human_genes, here("scRNAseq", "data", "mouse_human_genes.rds"))
