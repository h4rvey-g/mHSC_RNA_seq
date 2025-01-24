run_isoformswitch <- function(group1, group2,
                              salmon_dir,
                              gtf_file,
                              transcript_fasta,
                              alpha = 0.05,
                              dIFcutoff = 0.1,
                              intersected_data) {
    # library(IsoformSwitchAnalyzeR)
    # Create directory if it doesn't exist
    dir.create("results/104.isoformswitch", recursive = TRUE, showWarnings = FALSE)

    # Remove all files in the results/104.isoformswitch/group1_vs_group2 directory
    unlink(sprintf("results/104.isoformswitch/%s_vs_%s", group1, group2), recursive = TRUE, force = TRUE)

    isoformSwitchAnalysisPart1 <- function(
        # Arguments
        switchAnalyzeRlist,
        pathToGTF = NULL,
        pathToOutput = NULL,
        alpha = NULL,
        dIFcutoff = NULL) {
        isConditional <- switchAnalyzeRlist$sourceId != "preDefinedSwitches"

        # preFilter
        if (isConditional) {
            switchAnalyzeRlist <-
                preFilter(
                    switchAnalyzeRlist = switchAnalyzeRlist,
                    removeSingleIsoformGenes = TRUE,
                    quiet = TRUE,
                    alpha = alpha,
                    dIFcutoff = dIFcutoff
                )
        }

        # Test isoform switches
        if (isConditional) {
            if (any(switchAnalyzeRlist$conditions$nrReplicates > 5)) {
                switchAnalyzeRlist <-
                    isoformSwitchTestSatuRn(
                        switchAnalyzeRlist,
                        reduceToSwitchingGenes = TRUE,
                        quiet = TRUE,
                        alpha = alpha,
                        dIFcutoff = dIFcutoff
                    )
            } else {
                switchAnalyzeRlist <-
                    isoformSwitchTestDEXSeq(
                        switchAnalyzeRlist,
                        reduceToSwitchingGenes = FALSE,
                        quiet = TRUE,
                        alpha = alpha,
                        dIFcutoff = dIFcutoff
                    )
            }

            if (nrow(switchAnalyzeRlist$isoformSwitchAnalysis) == 0) {
                stop("No isoform switches were identified with the current cutoffs.")
            }
        }


        # Predict ORF

        if (is.null(switchAnalyzeRlist$orfAnalysis)) {
            # Add known annoation

            switchAnalyzeRlist <- addORFfromGTF(
                switchAnalyzeRlist = switchAnalyzeRlist,
                pathToGTF = pathToGTF,
                quiet = TRUE
            )

            # Predict novel once (if any are missing)

            if (any(switchAnalyzeRlist$orfAnalysis$orf_origin == "not_annotated_yet")) {
                switchAnalyzeRlist <- analyzeNovelIsoformORF(
                    switchAnalyzeRlist = switchAnalyzeRlist,
                    analysisAllIsoformsWithoutORF = TRUE,
                )
            }
        }

        # Extract and write sequences

        switchAnalyzeRlist <- extractSequence(
            switchAnalyzeRlist = switchAnalyzeRlist,
            extractNTseq = TRUE,
            extractAAseq = TRUE,
            addToSwitchAnalyzeRlist = TRUE,
            writeToFile = FALSE,
            pathToOutput = pathToOutput,
            quiet = TRUE,
            alpha = alpha,
            dIFcutoff = dIFcutoff
        )

        return(switchAnalyzeRlist)
    }

    isoformSwitchAnalysisPart2 <- function(
        # Core arguments
        switchAnalyzeRlist,
        # Analysis and output arguments
        pathToOutput = NULL,
        # Other arguments
        alpha = NULL,
        dIFcutoff = NULL) {
        # Predict intron retentions

        switchAnalyzeRlist <-
            analyzeAlternativeSplicing(
                switchAnalyzeRlist = switchAnalyzeRlist,
                quiet = TRUE,
                alpha = alpha,
                dIFcutoff = dIFcutoff
            )

        # Predict functional consequences

        switchAnalyzeRlist <-
            analyzeSwitchConsequences(
                switchAnalyzeRlist = switchAnalyzeRlist,
                consequencesToAnalyze = c("intron_retention", "ORF_seq_similarity", "NMD_status"),
                quiet = TRUE,
                alpha = alpha,
                dIFcutoff = dIFcutoff
            )

        # Make isoform switch plots

        switchPlotTopSwitches(
            switchAnalyzeRlist = switchAnalyzeRlist,
            n = Inf,
            pathToOutput = pathToOutput,
            filterForConsequences = TRUE,
            splitFunctionalConsequences = FALSE,
            quiet = TRUE,
            alpha = alpha,
            dIFcutoff = dIFcutoff
        )

        # Make overall consequences

        pdf(
            file = paste(
                pathToOutput,
                "common_switch_consequences.pdf",
                sep = ""
            ),
            width = 10,
            height = 7
        )
        print(
            extractConsequenceSummary(
                switchAnalyzeRlist = switchAnalyzeRlist,
                plot = TRUE,
                returnResult = FALSE,
                alpha = alpha,
                dIFcutoff = dIFcutoff
            )
        )
        dev.off()

        return(switchAnalyzeRlist)
    }

    isoformSwitchAnalysisCombined <- function(
        # Core arguments
        switchAnalyzeRlist,
        # Annotation arguments
        pathToGTF = NULL,
        # Analysis and output arguments
        pathToOutput = NULL,
        # Other arguments
        alpha = NULL,
        dIFcutoff = NULL) {
        # Run part 1

        switchAnalyzeRlist <-
            isoformSwitchAnalysisPart1(
                switchAnalyzeRlist = switchAnalyzeRlist,
                pathToOutput = pathToOutput,
                pathToGTF = pathToGTF,
                alpha = alpha,
                dIFcutoff = dIFcutoff
            )

        # Run part 2 without annoation

        switchAnalyzeRlist <-
            isoformSwitchAnalysisPart2(
                switchAnalyzeRlist = switchAnalyzeRlist,
                pathToOutput = pathToOutput,
                alpha = alpha,
                dIFcutoff = dIFcutoff
            )

        return(switchAnalyzeRlist)
    }
    # 导入Salmon定量结果
    salmonQuant <- importIsoformExpression(
        parentDir = salmon_dir,
        showProgress = FALSE,
        quiet = TRUE
    )
    counts <- salmonQuant$counts %>%
        column_to_rownames(var = "isoform_id")
    abundance <- salmonQuant$abundance %>%
        column_to_rownames(var = "isoform_id")

    # 构建设计矩阵
    design <- tibble(
        sampleID = colnames(salmonQuant$counts),
        condition = factor(ifelse(
            str_detect(sampleID, group1), group1,
            ifelse(str_detect(sampleID, group2), group2, NA)
        ))
    ) %>%
        filter(!is.na(condition))

    # 构建对比矩阵
    contrasts <- data.frame(
        condition_1 = group1,
        condition_2 = group2
    )

    # 创建switchAnalyzeRlist对象
    switchList <- importRdata(
        isoformCountMatrix = counts[, design$sampleID],
        isoformRepExpression = abundance[, design$sampleID],
        designMatrix = design,
        isoformExonAnnoation = gtf_file,
        isoformNtFasta = transcript_fasta,
        comparisonsToMake = contrasts,
        showProgress = FALSE,
        quiet = TRUE,
        ignoreAfterPeriod = TRUE # 忽略版本号后的部分
    )

    conversion <- genekitr::transId(switchList$isoformFeatures$gene_id %>% unique(), "symbol", "mouse", keepNA = TRUE)
    genes_to_keep <- conversion %>%
        filter(symbol %in% intersected_data$gene_symbol) %>%
        pull(input_id)

    switchList_filtered <- subsetSwitchAnalyzeRlist(
        switchList,
        subset = switchList$isoformFeatures$gene_id %in% genes_to_keep
    )
    # 运行分析
    tryCatch(
        {
            # 执行完整分析流程
            switchList_filtered_analysed <- isoformSwitchAnalysisCombined(
                switchList_filtered,
                pathToGTF = gtf_file,
                pathToOutput = sprintf("results/104.isoformswitch/%s_vs_%s", group1, group2),
                alpha = alpha,
                dIFcutoff = dIFcutoff
            )

            # 保存结果
            write_csv(
                extractSwitchSummary(switchList_filtered_analysed, filterForConsequences = FALSE),
                sprintf("results/104.isoformswitch/%s_vs_%s_summary.csv", group1, group2)
            )

            write_csv(
                switchList_filtered_analysed$isoformFeatures,
                sprintf("results/104.isoformswitch/%s_vs_%s_isoformfeatures.csv", group1, group2)
            )

            saveRDS(
                switchList_filtered_analysed,
                sprintf("results/104.isoformswitch/%s_vs_%s_switchlist.rds", group1, group2)
            )
        },
        error = function(e) {
            message(sprintf("Error in isoform switch analysis: %s", conditionMessage(e)))
            return(NULL)
        }
    )

    return(switchList_filtered)
}
