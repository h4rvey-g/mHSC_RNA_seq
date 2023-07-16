#!/bin/bash
workdir=./data
bams_indir=$workdir/04.STARalign
rmats_outdir=$workdir/06.rMATs
#1: prepare the sample group, rep1 in one group and rep2 in the other.
sampleMeta=./data/samples.txt
sample_groups=(A B C)
for group in "${sample_groups[@]}"; do
  sample_grp="$rmats_outdir/rmats_bams_group${group}.txt"
  echo -n >"$sample_grp"
  find "$bams_indir" -name "*-""$group""Aligned.sortedByCoord.out.bam" | while read -r line; do
    echo -n "$line," >>"$sample_grp"
  done
  sed -i '$ s/,$//g' "$sample_grp" #remove the last ","
done

# #2: do as detection using rmats turbo
genome_gtf="./data/00.reference/genome/gencode.vM32.chr_patch_hapl_scaff.annotation.gtf"
alias rmats='/data0/apps/anaconda3/bin/python /data0/apps/anaconda3/bin/rmats.py'
# use rmats to compare each group in a, b, c of $sample_groups, use loop
sample_groups=(A B C)
for ((i = 0; i < ${#sample_groups[@]}; i++)); do
  for ((j = i + 1; j < ${#sample_groups[@]}; j++)); do
    b1="$rmats_outdir/rmats_bams_group${sample_groups[i]}.txt"
    b2="$rmats_outdir/rmats_bams_group${sample_groups[j]}.txt"
    /data0/apps/anaconda3/bin/python /data0/apps/anaconda3/bin/rmats.py \
      --b1 $b2 \
      --b2 $b1 \
      --gtf $genome_gtf \
      --od "$rmats_outdir/${sample_groups[i]}_vs_${sample_groups[j]}" \
      --tmp "$rmats_outdir/tmp/${sample_groups[i]}_vs_${sample_groups[j]}" \
      --readLength 140 \
      -t paired \
      --nthread 30 \
      --statoff \
      --novelSS
  done
done
# remove all tmp in subfolders
find "$rmats_outdir" -name "tmp" | xargs rm -rf
# compress all files and folders in $rmat_outdir in one zip
rm -f "./results/06.rMATs.zip"
zip -r "./results/06.rMATs.zip" "$rmats_outdir"
