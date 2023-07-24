# footprinting.sh

# Purpose: to regenerate processed data tables from footprinting analysis
# Set up:
  # Install tobias3.7_env with YML file through conda
  # Install samtools (https://www.htslib.org/)

##############
# MERGE BAMS #
##############

samtools merge KW_CN.bam K1_CN.bam K2_CN.bam K3_CN.bam K4_CN.bam W1_CN.bam W2_CN.bam
samtools merge KW_S.bam K1_S.bam K2_S.bam K3_S.bam K4_S.bam W1_S.bam W2_S.bam

############################
# FOOTPRINTING WITH TOBIAS #
############################

conda activate tobias3.7_env

BAMDIR=ATACseq/data/bams
PEAKDIR=ATACseq/data/peaks
OUTDIR=ATACseq/data/footprinting

# Term exh samples....

SAMPLE=KW_CN

TOBIAS ATACorrect \
--bam ${BAMDIR}/${SAMPLE}.bam \
--genome ${OUTDIR}/genome.fa \
--peaks ${PEAKDIR}/all_merged_peaks_for_quantitation.merged.bed \
--blacklist ${OUTDIR}/mm10.blacklist.bed \
--outdir ${OUTDIR}/ATACorrect_${SAMPLE} \
--cores 2

TOBIAS FootprintScores \
--signal ${OUTDIR}/ATACorrect_${SAMPLE}/${SAMPLE}_corrected.bw \
--regions ${PEAKDIR}/all_merged_peaks_for_quantitation.merged.bed \
--output ${OUTDIR}/${SAMPLE}_footprints.bw \
--cores 2

TOBIAS ScoreBed \
--bed ${PEAKDIR}/all_merged_peaks_for_quantitation.merged.bed \
--bigwigs ${OUTDIR}/ATACorrect_${SAMPLE}/${SAMPLE}_corrected.bw \
--output ${OUTDIR}/merged_peaks_${SAMPLE}_scored.bed

# Prog exh samples ...

SAMPLE=KW_S

TOBIAS ATACorrect \
--bam ${BAMDIR}/${SAMPLE}.bam \
--genome ${OUTDIR}/genome.fa \
--peaks ${PEAKDIR}/all_merged_peaks_for_quantitation.merged.bed \
--blacklist ${OUTDIR}/mm10.blacklist.bed \
--outdir ${OUTDIR}/ATACorrect_${SAMPLE} \
--cores 2

TOBIAS FootprintScores \
--signal ${OUTDIR}/ATACorrect_${SAMPLE}/${SAMPLE}_corrected.bw \
--regions ${PEAKDIR}/all_merged_peaks_for_quantitation.merged.bed \
--output ${OUTDIR}/${SAMPLE}_footprints.bw \
--cores 2

TOBIAS ScoreBed \
--bed ${PEAKDIR}/all_merged_peaks_for_quantitation.merged.bed \
--bigwigs ${OUTDIR}/ATACorrect_${SAMPLE}/${SAMPLE}_corrected.bw \
--output ${OUTDIR}/merged_peaks_${SAMPLE}_scored.bed

# Predicting binding from footprinting profile...

TOBIAS BINDetect \
--motifs ${OUTDIR}/JASPAR2022_CORE_non-redundant_pfms.meme \
--signals ${OUTDIR}/KW_CN_footprints.bw ${OUTDIR}/KW_S_footprints.bw \
--genome ${OUTDIR}/genome.fa \
--peaks ${PEAKDIR}/all_merged_peaks_for_quantitation.merged.bed \
--peak_header ${PEAKDIR}/merged_peaks_annotated_header.txt \
--outdir ${OUTDIR}/BINDetect_output \
--cond_names KW_CN KW_S \
--cores 2

##############################
# MAKE PROCESSED DATA TABLES #
##############################

Rscript ATACseq/scripts/02_footprinting.R
