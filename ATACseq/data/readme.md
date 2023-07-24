# ATACseq data directory

## Contents

 * `counts`: counts matrices as output by ATACseq alignment workflow (described in Methods of paper)
 * `peaks`: peaks called for each sample
 * `processed_data_tables`: data tables used for figures
 
## Subdirectory Contents

### `peaks`

 * `peaks/all_merged_peaks_for_quantitation.merged.bed`: bedfile containing all merged peaks
 * `peaks/merged_peaks_annotated_header.txt`: text file with header definition of `all_merged_peakers_for_quantitation.merged.bed`

### `footprinting`

 * `footprinting/mm10.blacklist.bed`: bedfile containing updated mm10 blacklist regions
 * `footprinting/genome.fa`: fasta file for mm10 genome
 * `footprinting/JASPAR2022_CORE_non-redundant_pfms.meme`: all vertebrate motifs, downloaded from the [JASPAR database website](https://jaspar.genereg.net/downloads/)
