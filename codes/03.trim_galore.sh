#bin/bash
# based on the report at 1.cleandata_QC/1.cleandata_QC_multiqc_report.html, need to hard-cut the first 9bp according to the failed "Per Base Sequence Content"

# input read len = 150
# cut 5' len = 15
# length cutoff 150-15=135

getTrim() {
	indir=./data/02.CleanData
	outdir_trim=./data/03.TrimData
	outdir_QC=./results/03.trim_galore_QC
	# make dir at current working directionary
	mkdir -p $outdir_trim
	mkdir -p $outdir_QC
	line1=($1)
	sampleID=${line1[0]}
	groupName=${line1[2]}
	echo "Processing: $sampleID, $groupName"
	trim_galore \
		--paired \
		--gzip \
		--clip_R1 10 \
		--clip_R2 10 \
		--stringency 10 \
		--length 140 \
		--cores 4 \
		--fastqc --fastqc_args "-o $outdir_QC" \
		-o $outdir_trim \
		--basename "$sampleID-$groupName" \
		$indir/${sampleID}/${sampleID}_1.clean.fq.gz $indir/${sampleID}/${sampleID}_2.clean.fq.gz
}
export -f getTrim
parallel --progress --keep-order --line-buffer getTrim :::: ./data/samples.txt
/data0/apps/anaconda3/bin/multiqc -i 03.trim_galore_QC -o ./results/03.trim_galore_QC ./results/03.trim_galore_QC
