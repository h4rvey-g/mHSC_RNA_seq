intersect_data <- function(noiseq_res, exon_data, deg_data, p_cutoff = 0.05, log2fc_cutoff = 1, group1, group2) {
    # dtu_data <- dtu_data_KO4_WT
    # noiseq_res <- noiseq_results_KO4_WT
    # exon_data <- exon_data_KO4_WT
    # deg_data <- deg_data_KO4_WT
    # 过滤数据
    dtu_filtered <- noiseq_res %>%
        # col contains("dtu_gene") < p_cutoff
        filter_at(vars(contains("gene_padj")), all_vars(. < p_cutoff)) %>%
        filter_at(vars(contains("log2fold")), all_vars(abs(.) > log2fc_cutoff))
    exon_filtered <- exon_data %>%
        # col contains("exon_gene") < p_cutoff
        filter_at(vars(contains("gene_padj")), all_vars(. < p_cutoff))
    # filter_at(vars(contains("log2fold")), all_vars(abs(.) > log2fc_cutoff))
    deg_filtered <- deg_data %>%
        # col contains("deg_pval") < p_cutoff
        filter_at(vars(contains("pval")), all_vars(. < p_cutoff)) %>%
        filter_at(vars(contains("log2fold")), all_vars(abs(.) > log2fc_cutoff))

    # 提取基因ID
    dtu_genes <- unique(dtu_filtered$geneID)
    exon_genes <- unique(exon_filtered$geneID)
    deg_genes <- unique(deg_filtered$geneID)

    # 创建基因列表
    gene_list <- list(
        DTU = dtu_genes,
        DEX = exon_genes,
        DEG = deg_genes
    )

    # 绘制Venn图
    library(ggVennDiagram)
    p <- ggVennDiagram(
        gene_list,
        label = "both",
        label_alpha = 0.5
    ) +
        scale_fill_distiller(palette = "RdYlBu") +
        labs(title = "Intersection of DTU, DEX and DEG") +
        theme(plot.background = element_rect(fill = "white"))

    ggsave(sprintf("results/103.intersect/%s_vs_%s_intersection_p%.2f_fc%.1f.png", group1, group2, p_cutoff, log2fc_cutoff), p)
    # inner join of all three filtered data
    intersected_data <- dtu_filtered %>%
        dplyr::select(geneID, contains("log2fold"), contains("gene_padj"), contains("mean"), transcript_id) %>%
        inner_join(
            exon_filtered %>%
                dplyr::select(-contains("exon_exon_log2fold"), -contains("exon_exon_padj"), -featureID) %>%
                distinct(),
            by = "geneID"
        ) %>%
        inner_join(deg_filtered, by = "geneID")
    # use genekitr::transId to convert geneID to gene name
    conversion <- genekitr::transId(intersected_data$geneID %>% unique(), "symbol", "mouse")
    conversion <- conversion %>%
        rename(gene_symbol = symbol)
    intersected_data <- intersected_data %>%
        inner_join(conversion, by = c("geneID" = "input_id")) %>%
        # move gene_symbol to the first column
        dplyr::select(gene_symbol, everything())
    write_tsv(
        intersected_data,
        sprintf("results/103.intersect/%s_vs_%s_intersection_p%.2f_fc%.1f.tsv", group1, group2, p_cutoff, log2fc_cutoff)
    )
    intersected_data
}
