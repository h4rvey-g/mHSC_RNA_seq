update_hash <- function() {
  # library(jsonlite)

  # Read index info.json
  index_info <- fromJSON("data/102.splice_workflow/genome/index/salmon/info.json")

  # Hash values from index
  hash_values <- list(
    index_seq_hash = index_info$SeqHash,
    index_name_hash = index_info$NameHash,
    index_seq_hash512 = index_info$SeqHash512,
    index_name_hash512 = index_info$NameHash512,
    index_decoy_seq_hash = index_info$DecoySeqHash,
    index_decoy_name_hash = index_info$DecoyNameHash
  )

  # Get all sample directories
  sample_dirs <- list.dirs(
    "data/102.splice_workflow/star_salmon/salmon",
    full.names = TRUE,
    recursive = FALSE
  )

  # Update each sample's meta_info.json
  for (dir in sample_dirs) {
    meta_info_path <- file.path(dir, "aux_info", "meta_info.json")
    if (file.exists(meta_info_path)) {
      # Read existing meta_info.json
      meta_info <- fromJSON(meta_info_path)

      # Update hash values
      meta_info[names(hash_values)] <- hash_values

      # Write updated meta_info.json
      write(
        toJSON(meta_info, auto_unbox = TRUE, pretty = TRUE),
        meta_info_path
      )
    }
  }
}
load_data <- function() {
  # library(tximeta)
  # library(tidyverse)
  # library(SummarizedExperiment)

  # Create sample metadata dataframe
  samples <- data.frame(
    names = c(
      "Lin-KO4_1", "Lin-KO4_2", "Lin-KO4_3",
      "Lin-KO5_1", "Lin-KO5_2", "Lin-KO5_3",
      "Lin-WT_1", "Lin-WT_2", "Lin-WT_3"
    ),
    condition = c(rep("KO4", 3), rep("KO5", 3), rep("WT", 3)),
    stringsAsFactors = FALSE
  )

  # Add file paths
  samples$files <- file.path(
    "data/102.splice_workflow/star_salmon/salmon",
    samples$names, "quant.sf"
  )
  samples$gene_files <- file.path(
    "data/102.splice_workflow/star_salmon/salmon",
    samples$names, "quant.genes.sf"
  )

  # Verify all files exist
  stopifnot(all(file.exists(samples$files)))
  stopifnot(all(file.exists(samples$gene_files)))

  # Define paths
  indexDir <- "data/102.splice_workflow/genome/index/salmon"
  fastaPath <- "data/reference_data/alias/mm10/fasta/default/mm10.fa"
  gtfPath <- "data/reference_data/alias/mm10/gencode_gtf/default/mm10.gtf"

  # Create linkedTxome
  makeLinkedTxome(
    indexDir = indexDir,
    source = "LocalGENCODE",
    organism = "Mus musculus",
    release = "M22",
    genome = "mm10",
    fasta = fastaPath,
    gtf = gtfPath,
    write = TRUE,
    jsonFile = file.path(dirname(indexDir), "salmon_mm10_index.json")
  )
  # Import transcript-level data
  se_tx <- tximeta(samples)

  # Import gene-level data
  se_genes <- summarizeToGene(se_tx, assignRanges = "abundant", countsFromAbundance = "lengthScaledTPM")
  se_tx_2 <- SummarizedExperiment(
    assays = assays(se_tx),
    colData = colData(se_tx),
    rowData = rowData(se_tx)
  ) %>%
    # remove version number
    mutate(
      tx_name = .feature %>% str_replace("\\.\\d+$", ""),
      gene_name = unlist(gene_id) %>% str_replace("\\.\\d+$", "")
    )
  se_genes_2 <- SummarizedExperiment(
    assays = assays(se_genes),
    colData = colData(se_genes),
    rowData = rowData(se_genes)
  ) %>%
    # remove version number
    mutate(
      gene_name = .feature %>% str_replace("\\.\\d+$", "")
    )
  list(
    se_tx = se_tx_2,
    se_genes = se_genes_2
  )
}

load_dtu <- function(group1, group2) {
  dtu <- read_tsv(
    file.path(
      "data/102.splice_workflow/star_salmon/dexseq_dtu/results/stager",
      paste0("getAdjustedPValues.", group1, "-", group2, ".tsv")
    )
  )
  # add comp_name as prefix to column names, except cols geneID and txID
  colnames(dtu) <- c("geneID", "txID", paste0(group1, "_vs_", group2, "_dtu_", colnames(dtu)[-c(1, 2)]))
  # add padj as suffix to column names containing "dtu_gene" and "dtu_transcript"
  colnames(dtu) <- colnames(dtu) %>%
    str_replace_all("dtu_gene", "dtu_gene_padj") %>%
    str_replace_all("dtu_transcript", "dtu_transcript_padj")
  dex <- read_tsv(
    paste0("data/102.splice_workflow/star_salmon/dexseq_dtu/results/dexseq/DEXSeqResults.", group1, "-", group2, ".tsv")
  ) %>%
    dplyr::select(featureID, starts_with("log2fold")) %>%
    # rename col starting with log2fold to group1_vs_group2_dex_log2fold, and reverse the value
    rename_with(~ paste0(group1, "_vs_", group2, "_dex_log2fold"), starts_with("log2fold")) %>%
    mutate_at(vars(contains("log2fold")), ~ . * -1) %>%
    mutate(feature_name = featureID %>%
      # remove version number
      str_replace("\\.\\d+$", "")) %>%
    dplyr::select(-featureID)
  dtu_results <- dtu %>%
    left_join(dex, by = c("txID" = "feature_name"))
}

load_exon <- function(group1, group2) {
  # Load DEXSeq results
  dexseq <- read_csv(
    file.path(
      "data/102.splice_workflow/star_salmon/dexseq_exon/results",
      paste0("DEXSeqResults.", group1, "_vs_", group2, ".csv")
    )
  ) %>%
    dplyr::select(groupID, featureID, contains("log2fold"), contains("padj")) %>%
    # rename contains log2fold to group1_vs_group2_exon_log2fold
    rename_with(~ paste0(group1, "_vs_", group2, "_exon_exon_log2fold"), contains("log2fold")) %>%
    mutate_at(vars(contains("log2fold")), ~ . * -1) %>%
    # rename contains padj to group1_vs_group2_exon_padj
    rename_with(~ paste0(group1, "_vs_", group2, "_exon_exon_padj"), contains("padj")) %>%
    # Remove version number in groupID, and rename to geneID
    mutate(geneID = groupID %>% str_replace("\\.\\d+$", "")) %>%
    dplyr::select(-groupID)

  # Load per-gene q-values
  gene_qval <- read_csv(
    file.path(
      "data/102.splice_workflow/star_salmon/dexseq_exon/results",
      paste0("perGeneQValue.", group1, "_vs_", group2, ".csv")
    )
  ) %>%
    # Add prefix to column names except groupID
    rename_with(
      ~ paste0(group1, "_vs_", group2, "_exon_gene_", .),
      -groupID
    ) %>%
    # Remove version number in groupID, and rename to geneID
    mutate(geneID = groupID %>% str_replace("\\.\\d+$", "")) %>%
    dplyr::select(-groupID)

  # Merge results
  exon_results <- dexseq %>%
    left_join(gene_qval, by = "geneID")

  return(exon_results)
}
