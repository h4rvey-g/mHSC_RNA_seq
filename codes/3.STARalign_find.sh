#!/bin/bash
# Align the reads to the genome using the QCed fq files in 2.trim_galore

# input read len = 150

getAlign() {
    index_dir="./STARindex_150"
    indir=0.cleandata
    outdir=3.STARalign_find # use the basename of current .sh file to create a outdir
    # make dir at current working directionary
    #	mkdir -p $outdir #STAR can make folder by itself
    line=($1)
    sampleID="${line[0]}"
    # group=${line[2]}
    # echo "Processing: "$sampleID::$group
    STAR \
        --runMode alignReads \
        --twopassMode Basic \
        --genomeDir ${index_dir} \
        --runThreadN 24 \
        --readFilesIn $indir/${sampleID}_1.clean.fq.gz $indir/${sampleID}_2.clean.fq.gz \
        --readFilesCommand zcat \
        --outFileNamePrefix $outdir/$sampleID/$sampleID \
        --outSAMtype BAM SortedByCoordinate \
        --outBAMsortingThreadN 5 \
        --outFilterType BySJout \
        --quantMode GeneCounts \
        --outReadsUnmapped Fastx
}
export -f getAlign
# date
parallel --bar getAlign :::: samples.txt.tmp
# date
# /data0/apps/anaconda3/bin/multiqc -i 3.STARalign_tmp -o 3.STARalign_tmp 3.STARalign_tmp
