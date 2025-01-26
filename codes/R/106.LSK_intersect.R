Lin_LSK_intersect <- function(group1, group2, intersected_data) {
    # replace group1 and group2 with the actual group names
    # if WT, replace with "WT6"
    # if KO4, replace with "KO10"
    # if KO5, replace with "KO12"
    # intersected_data <- intersected_data_KO4_WT
    group1_LSK <- group1 %>%
        str_replace("WT", "WT6") %>%
        str_replace("KO4", "KO10") %>%
        str_replace("KO5", "KO12")
    group2_LSK <- group2 %>%
        str_replace("WT", "WT6") %>%
        str_replace("KO4", "KO10") %>%
        str_replace("KO5", "KO12")
    intersected_data_Lin <- intersected_data
    # load file from https://ghproxy.net/https://raw.githubusercontent.com/h4rvey-g/mHSC_LSK_RNA_seq/master/results/103.intersect/group1_vs_group2_intersection_p0.05_fc1.0.tsv
    intersected_data_LSK <- read_tsv(
        paste0("https://ghproxy.net/https://raw.githubusercontent.com/h4rvey-g/mHSC_LSK_RNA_seq/master/results/103.intersect/", group1_LSK, "_vs_", group2_LSK, "_intersection_p0.05_fc1.0.tsv")
    )
    # Get direction-matched intersections
    intersected_Lin_LSK <- intersect(
        intersected_data_Lin$gene_symbol[sign(intersected_data_Lin[[paste0(group1, "_vs_", group2, "_deg_log2fold")]]) == 1],
        intersected_data_LSK$gene_symbol[sign(intersected_data_LSK[[paste0(group1_LSK, "_vs_", group2_LSK, "_deg_log2fold")]]) == 1]
    ) %>%
        c(intersect(
            intersected_data_Lin$gene_symbol[sign(intersected_data_Lin[[paste0(group1, "_vs_", group2, "_deg_log2fold")]]) == -1],
            intersected_data_LSK$gene_symbol[sign(intersected_data_LSK[[paste0(group1_LSK, "_vs_", group2_LSK, "_deg_log2fold")]]) == -1]
        ))
    # create a new data frame with the intersected genes, and the log2fold values from both datasets
    intersected_Lin_LSK_data <- intersected_data_Lin %>%
        filter(gene_symbol %in% intersected_Lin_LSK) %>%
        dplyr::select(gene_symbol, geneID, paste0(group1, "_vs_", group2, "_deg_log2fold")) %>%
        distinct() %>%
        rename_at(vars(-gene_symbol), ~ paste0("Lin_", .)) %>%
        left_join(
            intersected_data_LSK %>%
                filter(gene_symbol %in% intersected_Lin_LSK) %>%
                dplyr::select(gene_symbol, paste0(group1_LSK, "_vs_", group2_LSK, "_deg_log2fold")) %>%
                distinct() %>%
                rename_at(vars(-gene_symbol), ~ paste0("LSK_", .)),
            by = "gene_symbol"
        ) %>%
        # arrange by mean of abs log2fold values from both datasets
        arrange(
            desc(
                rowMeans(abs(dplyr::select(., contains("deg_log2fold"))))
            )
        )
    dir.create("results/106.Lin_LSK_intersect", showWarnings = FALSE)
    write_tsv(
        intersected_Lin_LSK_data,
        paste0("results/106.Lin_LSK_intersect/", group1, "_vs_", group2, "_Lin_LSK_intersect.tsv")
    )
    intersected_Lin_LSK_data
}
