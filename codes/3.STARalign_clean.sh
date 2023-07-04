#!/bin/bash
# Align the reads to the genome using the QCed fq files in 2.trim_galore

# input read len = 135

getAlign() {
	index_dir="./STARindex"
	indir=2.trim_galore
	outdir=3.STARalign_clean # use the basename of current .sh file to create a outdir
	# make dir at current working directionary
	#	mkdir -p $outdir #STAR can make folder by itself
	line=($1)
	sampleID="${line[0]}-${line[2]}"
	# group=${line[2]}
	# echo "Processing: "$sampleID::$group
	STAR \
		--runMode alignReads \
		--twopassMode Basic \
		--genomeDir ${index_dir} \
		--runThreadN 10 \
		--readFilesIn $indir/${sampleID}_val_1.fq.gz $indir/${sampleID}_val_2.fq.gz \
		--readFilesCommand zcat \
		--outFileNamePrefix $outdir/$sampleID/$sampleID \
		--outSAMtype BAM SortedByCoordinate \
		--outBAMsortingThreadN 5 \
		--outFilterType BySJout \
		--quantMode GeneCounts
}
export -f getAlign
# date
parallel --bar getAlign :::: samples.txt
# date
/data0/apps/anaconda3/bin/multiqc -i 3.STARalign_clean -o 3.STARalign_clean 3.STARalign_clean
