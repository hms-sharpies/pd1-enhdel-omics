# Driver script to generate all processed data
# WARNING (***): noted scripts will take a long time to run on PC (potentially >1 hr)
# See README.md for instructions on where to download the outputs of this script instead

library(here)

scRNA_script_dir <- here("scRNAseq", "scripts")
source(file.path(scRNA_script_dir, "00_seurat_workflow.R"))
source(file.path(scRNA_script_dir, "01_dge.R")) # ***WARNING
source(file.path(scRNA_script_dir, "02_gene_signatures.R"))
source(file.path(scRNA_script_dir, "03_gsea.R")) # ***WARNING, this took about 15-ish hours to run, I haven't timed it though
source(file.path(scRNA_script_dir, "04_vision.R")) # ***WARNING
source(file.path(scRNA_script_dir, "05_pseudobulk.R"))

ATAC_script_dir <- here("ATACseq", "scripts")
source(file.path(ATAC_script_dir, "00_qc.R"))
source(file.path(ATAC_script_dir, "01_analyze.R"))
