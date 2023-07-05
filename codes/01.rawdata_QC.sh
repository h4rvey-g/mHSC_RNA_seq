#bin/bash
# get the fastqc results for all fq.gz files and integrated report by multiqc
indir=0.rawdata
outdir=$(basename $0 .sh) # use the basename of current .sh file to create a outdir
mkdir -p $outdir # make dir at current working directionary

for fq in $indir/*.fq.gz
do
	fastqc -o $outdir $fq &
done
# use multiqc to combine the reports from the folder $outdir, name the integrated report file with outdir name, put the output report file into $outdir
multiqc -i $outdir -o $outdir $outdir
