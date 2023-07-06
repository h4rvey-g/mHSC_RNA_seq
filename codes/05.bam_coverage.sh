#!/bin/bash
# This script is running inside the container, njaved/deeptools, so run the following command to get into the container:
# docker run -it -v $(pwd):/mHSC_RNA_seq njaved/deeptools /bin/bash
mkdir -p /mHSC_RNA_seq/results/05.bam_coverage
# index all bam files first
for bamfile in /mHSC_RNA_seq/data/04.STARalign/*/*.bam; do
    if [ ! -f "${bamfile}.bai" ]; then
        samtools index "${bamfile}"
    fi
done
# use multiBamSummary to get the coverage of all bam files in /mHSC_RNA_seq/data/04.STARalign
# multiBamSummary bins \
#     --bamfiles /mHSC_RNA_seq/data/04.STARalign/*/*.bam \
#     --labels $(ls /mHSC_RNA_seq/data/04.STARalign/*/*.bam | xargs -n 1 basename | sed 's/Aligned\.sortedByCoord\.out\.bam$//') \
#     --numberOfProcessors 30 \
#     --outRawCounts /mHSC_RNA_seq/results/05.bam_coverage/STARalign_clean_counts.tab \
#     --outFileName /mHSC_RNA_seq/results/05.bam_coverage/STARalign_clean.npz \
#     --binSize 10000
# use plotCorrelation to plot the correlation between samples
plotCorrelation \
    --corData /mHSC_RNA_seq/results/05.bam_coverage/STARalign_clean.npz \
    --corMethod pearson \
    --whatToPlot heatmap \
    --labels $(ls /mHSC_RNA_seq/data/04.STARalign/*/*.bam | xargs -n 1 basename | sed 's/Aligned\.sortedByCoord\.out\.bam$//') \
    --plotFile /mHSC_RNA_seq/results/05.bam_coverage/STARalign_clean_heatmap.pdf \
    --outFileCorMatrix /mHSC_RNA_seq/results/05.bam_coverage/STARalign_clean_correlation.tab \
    --removeOutliers \
    --plotNumbers \
    --colorMap RdYlBu \
    --plotTitle "STARalign_clean Spearman Correlation" \
    --plotHeight 10 \
    --plotWidth 10
# plot scatter plot
plotCorrelation \
    --corData /mHSC_RNA_seq/results/05.bam_coverage/STARalign_clean.npz \
    --corMethod pearson \
    --whatToPlot scatterplot \
    --plotFile /mHSC_RNA_seq/results/05.bam_coverage/STARalign_clean_scatterplot.pdf \
    --labels $(ls /mHSC_RNA_seq/data/04.STARalign/*/*.bam | xargs -n 1 basename | sed 's/Aligned\.sortedByCoord\.out\.bam$//') \
    --outFileCorMatrix /mHSC_RNA_seq/results/05.bam_coverage/STARalign_clean_correlation.tab \
    --removeOutliers \
    --plotTitle "STARalign_clean Spearman Correlation"
# use plotPCA to plot the PCA of samples
plotPCA \
    --corData /mHSC_RNA_seq/results/05.bam_coverage/STARalign_clean.npz \
    --plotFile /mHSC_RNA_seq/results/05.bam_coverage/STARalign_clean_PCA.pdf \
    --labels $(ls /mHSC_RNA_seq/data/04.STARalign/*/*.bam | xargs -n 1 basename | sed 's/Aligned\.sortedByCoord\.out\.bam$//') \
    --plotTitle "STARalign_clean PCA" \
    --plotHeight 10 \
    --plotWidth 10
