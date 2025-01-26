# group is a list of A_vs_B, A_vs_C, B_vs_C
group=(A_vs_B A_vs_C B_vs_C)
AS_type=(A3SS A5SS MXE RI SE)
gene_list=(ENSMUSG00000031785.17
    ENSMUSG00000020160.19
    ENSMUSG00000016494.10
    ENSMUSG00000032698.16
    ENSMUSG00000025492.7
    ENSMUSG00000032596.15
    ENSMUSG00000039410.17
    ENSMUSG00000031785.17
    ENSMUSG00000032596.15
    ENSMUSG00000039410.17
    ENSMUSG00000038227.16
    ENSMUSG00000046879.8
    ENSMUSG00000042817.16
)

get_sashimiplot() {
    group=$1
    AS_type=$2
    mkdir -p "data/06.rMATs/sashimiplot/${group}/${AS_type}"
    path_up="data/06.rMATs/${group}/summary/${AS_type}_up.tsv"
    path_down="data/06.rMATs/${group}/summary/${AS_type}_down.tsv"
    # filter GeneID column to get the rows of the gene_list
    gene_list=(ENSMUSG00000031785.17
        ENSMUSG00000020160.19
        ENSMUSG00000016494.10
        ENSMUSG00000032698.16
        ENSMUSG00000025492.7
        ENSMUSG00000032596.15
        ENSMUSG00000039410.17
        ENSMUSG00000031785.17
        ENSMUSG00000032596.15
        ENSMUSG00000039410.17
        ENSMUSG00000038227.16
        ENSMUSG00000046879.8
        ENSMUSG00000042817.16
    )
    echo "" >"data/06.rMATs/sashimiplot/${group}/${AS_type}.MATS.JCEC.txt"
    head -n 1 $path_up >"data/06.rMATs/sashimiplot/${group}/${AS_type}.MATS.JCEC.txt"
    for gene in ${gene_list[@]}; do
        grep $gene $path_up >>"data/06.rMATs/sashimiplot/${group}/${AS_type}.MATS.JCEC.txt"
        grep $gene $path_down >>"data/06.rMATs/sashimiplot/${group}/${AS_type}.MATS.JCEC.txt"
    done
    AS_file="data/06.rMATs/sashimiplot/${group}/${AS_type}.MATS.JCEC.txt"
    # extract the first and last letter from group
    group1=$(echo $group | cut -c1)
    group2=$(echo $group | cut -c6)
    group1_path="data/06.rMATs/rmats_bams_group${group1}.txt"
    group2_path="data/06.rMATs/rmats_bams_group${group2}.txt"
    rmats2sashimiplot --b1 $group2_path --b2 $group1_path \
        -e $AS_file --event-type $AS_type \
        --l1 "Group $group2" --l2 "Group $group1" --exon_s 1 --intron_s 5 \
        -o "data/06.rMATs/sashimiplot/${group}/${AS_type}/"
}

export -f get_sashimiplot
parallel -j 30 get_sashimiplot ::: ${group[@]} ::: ${AS_type[@]}
# zip the sashimiplot results, but don't add dirs which names has "index", output zip is results/Lin_AS_GO_KEGG.zip
find data/06.rMATs/sashimiplot -type f -not -path "*index*" | zip results/Lin_AS_sashimiplot.zip -@
