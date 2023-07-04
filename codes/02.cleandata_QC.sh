#bin/bash
source ~/.envrc
echo $(pwd)
# get the fastqc results for all fq.gz files and integrated report by multiqc
indir=./data/02.CleanData
outdir=./results/02.CleanData_QC
mkdir -p $outdir          # make dir at current working directionary
# export PATH="/home/guozhonghao/.conda/envs/fastqc/bin:$PATH"
# for fq in $(find $indir -type f -name "*.fq.gz"); do
# 	fastqc -o $outdir $fq &
# done
# convert the for loop to a function and use parallel to run the function
getFastqc() {
	fq=$1
outdir=./results/02.CleanData_QC
	# mkdir -p $outdir/$(dirname $fq)
	fastqc -o $outdir $fq
}
export -f getFastqc
parallel --progress --keep-order --line-buffer getFastqc :::: <(find $indir -type f -name "*.fq.gz")
# use multiqc to combine the reports from the folder $outdir, name the integrated report file with outdir name, put the output report file into $outdir
/data0/apps/anaconda3/bin/multiqc -i $outdir -o $outdir $outdir
