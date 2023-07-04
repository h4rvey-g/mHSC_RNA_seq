#bin/bash
# based on the report at 1.cleandata_QC/1.cleandata_QC_multiqc_report.html, need to hard-cut the first 9bp according to the failed "Per Base Sequence Content"

# input read len = 150
# cut 5' len = 15
# length cutoff 150-15=135

getTrim() {
	indir=0.cleandata
	outdir_fq=2.trim_galore # use the basename of current .sh file to create a outdir
	outdir_rep=2.trim_galore_QC
	# make dir at current working directionary
	mkdir -p $outdir_fq
	mkdir -p $outdir_rep
	line1=($1)
	sampleID=${line1[0]}
	groupName=${line1[2]}
	echo "Processing: $sampleID, $groupName"
	trim_galore \
		--paired \
		--gzip \
		--clip_R1 15 \
		--clip_R2 15 \
		--stringency 10 \
		--length 135 \
		--cores 4 \
		--fastqc --fastqc_args "-o $outdir_rep" \
		-o $outdir_fq \
		--basename "$sampleID-$groupName" \
		$indir/${sampleID}_1.clean.fq.gz $indir/${sampleID}_2.clean.fq.gz
}
export -f getTrim
parallel getTrim :::: samples.txt
/data0/apps/anaconda3/bin/multiqc -i 2.trim_galore_QC -o 2.trim_galore_QC 2.trim_galore_QC
