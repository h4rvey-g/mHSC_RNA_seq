#!/bin/bash

workdir=./Analysis_XYG
bams_indir=$workdir/3.1.bam_clean

#1: prepare the sample group, rep1 in one group and rep2 in the other.
sampleMeta=./Analysis_XYG/samples.txt
sample_grpA=rmats_bams_goupA.txt
sample_grpB=rmats_bams_goupB.txt
sample_grpC=rmats_bams_goupC.txt
sample_grpD=rmats_bams_goupD.txt
tmp_file=$workdir/tmp.txt
#output the group 1
echo -n >$tmp_file #reset the tmp file
ls $bams_indir/*A.bam | while read -r line; do
  echo -n "$line," >>$tmp_file
done
sed '$ s/,$//g' $tmp_file >$sample_grpA #remove the last ","

echo -n >$tmp_file #reset the tmp file
ls $bams_indir/*B.bam | while read -r line; do
  echo -n "$line," >>$tmp_file
done
sed '$ s/,$//g' $tmp_file >$sample_grpB #remove the last ","

echo -n >$tmp_file #reset the tmp file
ls $bams_indir/*C.bam | while read -r line; do
  echo -n "$line," >>$tmp_file
done
sed '$ s/,$//g' $tmp_file >$sample_grpC #remove the last ","

echo -n >$tmp_file #reset the tmp file
ls $bams_indir/*D.bam | while read -r line; do
  echo -n "$line," >>$tmp_file
done
sed '$ s/,$//g' $tmp_file >$sample_grpD #remove the last ","

sample_grpA=$(ls $bams_indir/*A.bam | tr '\n' ',') sample_grpA=${sample_grpA%?} echo $sample_grpA > rmats_bams_goupA.txt

sample_grpB=$(ls $bams_indir/*B.bam | tr '\n' ',') sample_grpB=${sample_grpB%?} echo $sample_grpB > rmats_bams_goupB.txt

sample_grpC=$(ls $bams_indir/*C.bam | tr '\n' ',') sample_grpC=${sample_grpC%?} echo $sample_grpC > rmats_bams_goupC.txt

sample_grpD=$(ls $bams_indir/*D.bam | tr '\n' ',') sample_grpD=${sample_grpD%?} echo $sample_grpD > rmats_bams_goupD.txt
#

#2: do AS detection using rMATs turbo
rmats_outdir=$workdir/4.rMATs
genome_gtf="/data0/work/XuLab/references/mouse/GRCm39_GENCODE_M31/gencode.vM31.primary_assembly.annotation.gtf"

rmats \
  --b1 $sample_grpA \
  --b2 $sample_grpB \
  --gtf $genome_gtf \
  --od $rmats_outdir \
  --tmp $rmats_outdir/tmp \
  --readLength 150 \
  -t paired \
  --nthread 30 \
  --statoff \
  --novelSS

rmats \
  --b1 $sample_grpC \
  --b2 $sample_grpD \
  --gtf $genome_gtf \
  --od $rmats_outdir \
  --tmp $rmats_outdir/tmp \
  --readLength 150 \
  -t paired \
  --nthread 30 \
  --statoff \
  --novelSS