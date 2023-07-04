#!/bin/bash
# eval "$(mamba shell hook -s bash)"
# mamba activate sortmerna
get_sortmerna() {
    sortmerna -ref ./Mus_musculus.GRCm39.ncrna.fa -reads $1 -reads $2 \
        -workdir "./2.sortmerna/$(basename $1 _1.clean.fq.gz)" -fastx -paired_out -other "./2.sortmerna/$(basename $1 _1.clean.fq.gz)/$(basename $1 _1.clean.fq.gz)" -out2 -threads 30
}
export -f get_sortmerna
parallel --link --progress --keep-order --line-buffer \
    get_sortmerna :::: <(find 0.cleandata -type f -name "*_1.clean.fq.gz") :::: <(find 0.cleandata -type f -name "*_2.clean.fq.gz")
/data0/apps/anaconda3/bin/multiqc -i 2.sortmerna -o 2.sortmerna 2.sortmerna
# move all _fwd.fq.gz and _rev.fq.gz to 2.sortmerna/results
mkdir -p 2.sortmerna/results
find 2.sortmerna -type f -name "*_fwd.fq.gz" -exec mv {} 2.sortmerna/results \;
find 2.sortmerna -type f -name "*_rev.fq.gz" -exec mv {} 2.sortmerna/results \;
kraken --db kraken_db --threads 30 --gzip-compressed --paired --output 4.kraken/kraken.out \
--report 4.kraken/kraken.report --use-names --confidence 0.1 --report-zero-counts --fastq-input 2.sortmerna/results/*_fwd.fq.gz 2.sortmerna/results/*_rev.fq.gz
