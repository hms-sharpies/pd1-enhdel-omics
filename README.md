
<!-- README.md is generated from README.Rmd. Please edit that file -->

# pd1-enhdel-omics

<!-- badges: start -->

<!-- badges: end -->

This repository provides code to reproduce results and figures from
**Deletion of an exhaustion-specific PD-1 enhancer improves T cell
function and viral control**.

The associated GEO repository is GSE212507. Additionally, the single
cell counts matrices and other associated data tables have been uploaded
to the Single Cell Portal under the [PD-1 enhancer
deletion](https://singlecell.broadinstitute.org/single_cell/study/SCP1772/pd-1-enhancer-deletion).

## Contents

  - `scRNAseq`: directory containing the data and scripts for
    reproducing scRNA-seq analyses and figures
  - `ATACseq`: directory containing the data and scripts for reproducing
    ATAC-seq analyses and figures
  - `figures`: directory containing the scripts to reproduce the figures
  - `renv`: directory containing all information for tracking R package
    versioning and downloads. See [Introduction to
    renv](https://rstudio.github.io/renv/articles/renv.html) for details
  - `aesthetics.R`: file containing aesthetics and themes for the
    project
  - `helper.R`: file containing libraries and helper functions
  - `driver.R`: optional script to regenerate all intermediate results.
    May take a long time (\> 24 hours) to execute.

## Instructions

To regenerate figures from raw data, take the following steps:

1.  Clone the GitHub.
2.  Restore the renv environment from the renv.lock file, using
    `renv::restore()`. The [VISION
    package](https://github.com/YosefLab/VISION), a signature scoring
    method developed by the Yosef lab, may have to be manually installed
    from their GitHub repository. The package can be installed using
    `devtools` by running the code below:

<!-- end list -->

    library(devtools)
    devtools::install_github("YosefLab/VISION")

3.  Download the single cell counts matrices from the Single Cell Portal
    and place them in the `scRNAseq/data/counts_matrices`.
4.  Regenerate all intermediate data tables and objects by running the
    script `driver.R`. Warning: It may take a long time, \>24 hours
    depending on compute power, to execute.
5.  Use Rmarkdown to knit the `.Rmd` files in the `figures/` directory
    to output the figures in either HTML or PDF format.

Alternatively, to bypass the long compute time in step 3, follow the
instructions below to generate/acquire processed data files, then knit
the `.Rmd` files in the `figures/` directory.

1.  Run the script at `scRNAseq/scripts/02_gene_signatures` to generate
    the
    `scRNAseq/data/gene_signatures/02_gene_signatures.msigdbr_H-C2-C7.tsv`
2.  Download the processed data files from the Study Files section
    [Single Cell
    Portal](https://singlecell.broadinstitute.org/single_cell/study/SCP1772/pd-1-enhancer-deletion#study-download)
    and place them at the following paths.

<!-- end list -->

  - scRNAseq/data/counts\_matrices.zip -\> decompress
  - scRNAseq/data/processed\_data\_tables.zip -\> decompress
  - scRNAseq/data/processed\_data\_objects/00\_seurat\_workflow.so.rds
  - scRNAseq/data/processed\_data\_objects/05\_pseudobulk.deseq\_list.rds

## Contact

Please contact a corresponding author for any questions, comments, or
concerns regarding the paper in general.

## Compute environment and version

``` r
sessionInfo()
#> R version 4.0.3 (2020-10-10)
#> Platform: x86_64-apple-darwin17.0 (64-bit)
#> Running under: macOS Big Sur 10.16
#> 
#> Matrix products: default
#> BLAS:   /Library/Frameworks/R.framework/Versions/4.0/Resources/lib/libRblas.dylib
#> LAPACK: /Library/Frameworks/R.framework/Versions/4.0/Resources/lib/libRlapack.dylib
#> 
#> locale:
#> [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices datasets  utils     methods   base     
#> 
#> loaded via a namespace (and not attached):
#>  [1] digest_0.6.29       magrittr_2.0.3      evaluate_0.15      
#>  [4] rlang_1.0.6         stringi_1.7.6       cli_3.4.1          
#>  [7] renv_0.15.4         rstudioapi_0.13     rmarkdown_2.14     
#> [10] tools_4.0.3         stringr_1.4.0       xfun_0.30          
#> [13] yaml_2.3.5          fastmap_1.1.0       compiler_4.0.3     
#> [16] BiocManager_1.30.17 htmltools_0.5.2     knitr_1.39
```

### A note on reproducibility

Anyone running the code should use the same R and R package versions as
defined by `renv`. However, I cannot guarantee that one will get the
exact same results in the single cell clustering unless one also ensures
that they are using the same processor and operating system. See [this
discussion](https://github.com/scverse/scanpy/issues/2014) about
problems in reproducibility of single cell clustering analyses
persisting in Scanpy.
