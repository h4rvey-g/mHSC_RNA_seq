#!/bin/bash
nextflow run file:/data0/work/guozhonghao/HCC/tools/rnaseq.git \
    --input data/01.RawData/samplesheet.csv \
    --outdir data/101.RNA_workflow \
    --gtf data/reference_data/alias/mm10/gencode_gtf/default/mm10.gtf.gz \
    --fasta data/reference_data/alias/mm10/fasta/default/mm10.fa \
    --featurecounts_group_type gene_type \
    --save_reference \
    -profile docker \
    --gencode \
    -resume \
    -c ./nf.config # --salmon_index data/reference_data/alias/hg38/salmon_sa_index/default/ \

nextflow run file:/data0/work/guozhonghao/mHSC_RNA_seq/tools/rnasplice.git \
    -r dev \
    --source fastq \
    --input data/101.RNA_workflow/star_salmon/samplesheet_splice.csv \
    --outdir data/102.splice_workflow \
    --gtf data/reference_data/alias/mm10/gencode_gtf/default/mm10.gtf \
    --fasta data/reference_data/alias/mm10/fasta/default/mm10.fa \
    --contrasts data/101.RNA_workflow/star_salmon/contrast_splice.csv \
    --miso_genes 'ENSMUSG00000000544.14' \
    --save_reference \
    -profile docker \
    --gencode \
    -resume \
    -c ./nf_splice.config
