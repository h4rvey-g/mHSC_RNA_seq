#bin/bash
# get the fastqc results for all fq.gz files and integrated report by multiqc
indir=0.cleandata
outdir=$(basename $0 .sh) # use the basename of current .sh file to create a outdir
mkdir -p $outdir          # make dir at current working directionary
export PATH="/home/guozhonghao/.conda/envs/fastqc/bin:$PATH"
# for fq in $(find $indir -type f -name "*.fq.gz"); do
# 	fastqc -o $outdir $fq &
# done
# convert the for loop to a function and use parallel to run the function
getFastqc() {
	fq=$1
	/home/guozhonghao/.conda/envs/fastqc/bin/fastqc -o ./1.cleandata_QC $fq
}
export -f getFastqc
parallel getFastqc :::: <(find $indir -type f -name "*.fq.gz")
# use multiqc to combine the reports from the folder $outdir, name the integrated report file with outdir name, put the output report file into $outdir
multiqc -i $outdir -o $outdir $outdir
