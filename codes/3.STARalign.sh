#!/bin/bash
# Align the reads to the genome using the QCed fq files in 2.trim_galore

# input read len = 150

getAlign() {
	index_dir="./STARindex_150"
	indir=2.sortmerna/results
	outdir=3.STARalign_norrna # use the basename of current .sh file to create a outdir
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
		--runThreadN 2 \
		--readFilesIn $indir/${sampleID}_fwd.fq.gz $indir/${sampleID}_rev.fq.gz \
		--readFilesCommand zcat \
		--outFileNamePrefix $outdir/$sampleID/$sampleID \
		--outSAMtype BAM SortedByCoordinate \
		--outBAMsortingThreadN 5 \
		--outFilterType BySJout \
		--quantMode GeneCounts 
		# --outFilterScoreMinOverLread 0 \
		# --outFilterMatchNminOverLread 0
}
export -f getAlign
# date
parallel --progress --keep-order --line-buffer getAlign ::: Lin1-WT Lin2-WT Lin3-WT Lin3-KO LSK2
# date
/data0/apps/anaconda3/bin/multiqc -i 3.STARalign_norrna -o 3.STARalign_norrna 3.STARalign_norrna
