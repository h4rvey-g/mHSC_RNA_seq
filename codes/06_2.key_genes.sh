group=(A_vs_B A_vs_C B_vs_C)
AS_type=(A3SS A5SS MXE RI SE)
# extract lines from all data/06.rMATs/${group}/summary/${AS_type}_*.tsv, which "GeneID" matches the gene_list, only keep column "GeneID", "geneSymbol", "PValue", "FDR", "IncLevelDifference", save all results to results/06_2.key_genes/Lin_AS_result.tsv
echo "" >results/06_2.key_genes/Lin_AS_result.tsv
get_AS_result() {
    gene_list=(
        ENSMUSG00000006389.14
        ENSMUSG00000004043.15
        ENSMUSG00000006386.16
        ENSMUSG00000015053.15
        ENSMUSG00000052534.16
        ENSMUSG00000026815.15
        ENSMUSG00000038418.8
        ENSMUSG00000028717.13
        ENSMUSG00000038227.16
        ENSMUSG00000020160.19
        ENSMUSG00000040732.21
        ENSMUSG00000069769.14
        ENSMUSG00000069769.14
        ENSMUSG00000002111.9
        ENSMUSG00000031162.15
        ENSMUSG00000002111.9
        ENSMUSG00000034957.11
        ENSMUSG00000056501.4
        ENSMUSG00000002111.9
        ENSMUSG00000034957.11
        ENSMUSG00000056501.4
        ENSMUSG00000031162.15
        ENSMUSG00000026815.15
        ENSMUSG00000022952.18
        ENSMUSG00000022952.18
        ENSMUSG00000002111.9
        ENSMUSG00000056501.4
    )
    group=$1
    AS_type=$2
    path_up="data/06.rMATs/${group}/summary/${AS_type}_up.tsv"
    path_down="data/06.rMATs/${group}/summary/${AS_type}_down.tsv"
    tmp_file="results/06_2.key_genes/${group}_${AS_type}_tmp.tsv"
    echo "" >$tmp_file
    head -n 1 $path_up >$tmp_file
    for gene in ${gene_list[@]}; do
        # extract the lines which "GeneID" matches the gene_list
        grep $gene $path_up >>$tmp_file
        grep $gene $path_down >>$tmp_file
    done
    # find column number of "GeneID", "geneSymbol", "PValue", "FDR", "IncLevelDifference" in the header
    col_nums=$(head -n 1 $tmp_file | tr '\t' '\n' | cat -n | grep -E "GeneID|geneSymbol|PValue|FDR|IncLevelDifference" | cut -f1 | tr '\n' ',')
    # remove the space and the last comma in col_nums
    col_nums=$(echo $col_nums | sed 's/ //g' | sed 's/,$//')
    # echo $col_nums
    # extract columns "GeneID", "geneSymbol", "PValue", "FDR", "IncLevelDifference"
    cut -f$col_nums $tmp_file >$tmp_file.tmp
    # add column "group" and "AS_type" to the result, the column names are
    sed -i "s/^/$group\t$AS_type\t/g" $tmp_file.tmp
    # save the result to results/06_2.key_genes/Lin_AS_result.tsv
    cat $tmp_file.tmp >>results/06_2.key_genes/Lin_AS_result.tsv
}

export -f get_AS_result
parallel -j 30 get_AS_result ::: ${group[@]} ::: ${AS_type[@]}
sort results/06_2.key_genes/Lin_AS_result.tsv | uniq >results/06_2.key_genes/unique_Lin_AS_result.tsv
mv results/06_2.key_genes/unique_Lin_AS_result.tsv results/06_2.key_genes/Lin_AS_result.tsv
# delete the tmp files
rm results/06_2.key_genes/*tmp*
# # File path
# file="results/06_2.key_genes/Lin_AS_result.tsv"

# # Find line numbers with "GeneID"
# line_numbers=$(grep -n "GeneID" "$file" | awk -F: '{print $1}')
# # Convert line numbers to an array
# line_numbers_array=($(echo ${line_numbers}))
# # Loop through line numbers and delete all except the first occurrence
# for ((i = 1; i < ${#line_numbers_array[@]}; i++)); do
#     sed -i "${line_numbers_array[i]}d" "$file"
# done
get_sashimiplot() {
    group=$1
    AS_type=$2
    mkdir -p "results/06_2.key_genes/${group}/${AS_type}"
    path_up="data/06.rMATs/${group}/summary/${AS_type}_up.tsv"
    path_down="data/06.rMATs/${group}/summary/${AS_type}_down.tsv"
    # filter GeneID column to get the rows of the gene_list
    gene_list=(
        ENSMUSG00000006389.14
        ENSMUSG00000004043.15
        ENSMUSG00000006386.16
        ENSMUSG00000015053.15
        ENSMUSG00000052534.16
        ENSMUSG00000026815.15
        ENSMUSG00000038418.8
        ENSMUSG00000028717.13
        ENSMUSG00000038227.16
        ENSMUSG00000020160.19
        ENSMUSG00000040732.21
        ENSMUSG00000069769.14
        ENSMUSG00000069769.14
        ENSMUSG00000002111.9
        ENSMUSG00000031162.15
        ENSMUSG00000002111.9
        ENSMUSG00000034957.11
        ENSMUSG00000056501.4
        ENSMUSG00000002111.9
        ENSMUSG00000034957.11
        ENSMUSG00000056501.4
        ENSMUSG00000031162.15
        ENSMUSG00000026815.15
        ENSMUSG00000022952.18
        ENSMUSG00000022952.18
        ENSMUSG00000002111.9
        ENSMUSG00000056501.4
    )
    # make gene_list unique
    gene_list=($(echo "${gene_list[@]}" | tr ' ' '\n' | sort -u | tr '\n' ' '))
    echo "" >"results/06_2.key_genes/${group}/${AS_type}.MATS.JCEC.txt"
    head -n 1 $path_up >"results/06_2.key_genes/${group}/${AS_type}.MATS.JCEC.txt"
    for gene in ${gene_list[@]}; do
        grep $gene $path_up >>"results/06_2.key_genes/${group}/${AS_type}.MATS.JCEC.txt"
        grep $gene $path_down >>"results/06_2.key_genes/${group}/${AS_type}.MATS.JCEC.txt"
    done
    AS_file="results/06_2.key_genes/${group}/${AS_type}.MATS.JCEC.txt"
    # extract the first and last letter from group
    group1=$(echo $group | cut -c1)
    group2=$(echo $group | cut -c6)
    group1_path="data/06.rMATs/rmats_bams_group${group1}.txt"
    group2_path="data/06.rMATs/rmats_bams_group${group2}.txt"
    /home/guozhonghao/mambaforge/envs/rmats2sashimiplot/bin/rmats2sashimiplot --b1 $group2_path --b2 $group1_path \
        -e $AS_file --event-type $AS_type \
        --l1 "Group $group2" --l2 "Group $group1" --exon_s 1 --intron_s 5 \
        -o "results/06_2.key_genes/${group}/${AS_type}/"
}

export -f get_sashimiplot
parallel -j 30 get_sashimiplot ::: ${group[@]} ::: ${AS_type[@]}
# zip the sashimiplot results, but don't add dirs which names has "index", output zip is results/Lin_AS_GO_KEGG.zip
find results/06_2.key_genes -type f -not -path "*index*" | zip results/Lin_key_genes_AS.zip -@
