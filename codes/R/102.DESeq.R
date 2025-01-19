run_DEG <- function(se_data, group1, group2) {
    se_genes <- se_data$se_genes
    se_genes_de <- se_genes %>%
        identify_abundant(factor_of_interest = condition) %>%
        test_differential_abundance(
            ~ 0 + condition,
            method = "limma_voom",
            contrasts = c(paste0("condition", group1, " - condition", group2)),
            action = "get"
        )
    se_genes_de <- se_genes_de %>%
        # rename transcript ID to geneID, remove version number
        mutate(geneID = transcript %>% str_replace("\\.\\d+$", "")) %>%
        dplyr::select(-transcript) %>%
        # rename containing logFC to group1_vs_group2_deg_log2fold
        rename_with(~ paste0(group1, "_vs_", group2, "_deg_log2fold"), contains("logFC")) %>%
        # rename containing pval to group1_vs_group2_deg_padj
        rename_with(~ paste0(group1, "_vs_", group2, "_deg_pval"), contains("adj.P")) %>%
        dplyr::select(geneID, contains("deg"))
    se_genes_de
}
