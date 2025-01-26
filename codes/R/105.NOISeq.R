run_noiseq_analysis <- function(se_data, group1, group2) {
    # Extract transcript-level counts
    se_tx <- se_data$se_tx
    counts <- assays(se_tx)$counts

    # Create factors dataframe
    samples <- data.frame(
        names = colnames(counts),
        condition = colData(se_tx)$condition,
        stringsAsFactors = FALSE
    )

    # Filter for relevant conditions
    keep_samples <- samples$condition %in% c(group1, group2)
    samples <- samples[keep_samples, ]
    counts <- counts[, keep_samples]

    # Get transcript lengths
    tx_lengths <- se_tx %>% pull(length)
    names(tx_lengths) <- rownames(counts)

    # Pre-filter with more nuanced criteria
    cpm <- edgeR::cpm(counts)

    # Filter out transcripts with zero or near-zero expression
    keep_nonzero <- rowSums(counts) > 0

    # For non-zero transcripts, keep those with at least 0.01 CPM in at least 2 samples
    keep <- rowSums(cpm >= 0.01) >= 2

    # Apply filters
    counts <- counts[keep, ]
    tx_lengths <- tx_lengths[keep]
    # Print filtering statistics
    cat(sprintf("Filtering summary:\n"))
    cat(sprintf("- Total transcripts: %d\n", length(keep)))
    cat(sprintf("- Transcripts with zero counts: %d\n", sum(!keep_nonzero)))
    cat(sprintf(
        "- Transcripts kept: %d (%.1f%%)\n",
        sum(keep), 100 * sum(keep) / length(keep)
    ))

    # Print filtering statistics
    cat(sprintf(
        "Keeping %d out of %d transcripts (%.1f%%)\n",
        sum(keep), length(keep),
        100 * sum(keep) / length(keep)
    ))

    # Create NOISeq data object
    noiseq_data <- NOISeq::readData(
        data = counts,
        length = tx_lengths,
        factors = data.frame(
            condition = samples$condition
        )
    )

    # Run NOISeqBIO analysis with CPM filtering method
    noiseq_results <- NOISeq::noiseqbio(
        noiseq_data,
        k = 0.5, # Parameter for replacing zero counts
        norm = "rpkm", # Use RPKM normalization
        factor = "condition",
        conditions = c(group1, group2),
        lc = 1, # Length correction factor
        r = 20, # Number of samples to simulate for null distribution
        adj = 1.5, # Adjustment factor for kernel density estimation
        plot = FALSE, # Don't generate diagnostic plots
        filter = 1, # Use CPM filtering method
        cpm = 1, # Minimum CPM threshold
        cv.cutoff = 100 # Maximum coefficient of variation allowed
    )

    # Extract results and add more information
    results_df <- noiseq_results@results[[1]] %>%
        as.data.frame() %>%
        tibble::rownames_to_column("transcript_id") %>%
        dplyr::mutate(
            group1 = group1,
            group2 = group2,
            comparison = paste0(group1, "_vs_", group2),
            significant = prob >= 0.95 # Using stricter threshold for biological replicates
        )
    # add geneID, based on data/102.splice_workflow/genome/mm10.tx2gene.tsv
    tx2gene <- read_tsv("data/102.splice_workflow/genome/mm10.tx2gene.tsv", col_names = c("transcript_id", "geneID", "gene_symbol"))
    results_df <- results_df %>%
        inner_join(tx2gene, by = "transcript_id") %>%
        dplyr::select(-gene_symbol) %>%
        as_tibble()
    # rename transcript_id to geneID, log2FC to log2fold
    # dplyr::rename(
    #     geneID = transcript_id,
    #     log2fold = log2FC
    # )
    # dplyr::rename(
    #     !!paste0(group1, "_mean") := "group1_mean",
    #     !!paste0(group2, "_mean") := "group2_mean"
    # )
    # Convert probability to p-value (since NOISeq gives prob = 1-FDR)
    prob2pvalue <- function(prob) {
        # Handle NA values
        p <- rep(NA_real_, length(prob))
        non_na <- !is.na(prob)
        # Convert prob to p-value: prob is essentially 1-FDR
        p[non_na] <- 1 - prob[non_na]
        return(p)
    }

    perGeneQValue <- function(results_df, method = "exact", nperm = 24) {
        # Get p-values from probabilities
        pvals <- prob2pvalue(results_df$prob)

        # Get only testable transcripts
        wTest <- which(!is.na(pvals))
        if (length(wTest) == 0) {
            stop("No valid p-values found")
        }

        # Use only those transcripts that were testable
        pvals <- pvals[wTest]
        geneID <- factor(results_df$geneID[wTest])
        geneSplit <- split(seq_along(geneID), geneID)

        # Summarize p-values of transcripts for one gene: take the minimum
        pGene <- sapply(geneSplit, function(i) min(pvals[i]))
        if (!all(is.finite(pGene))) {
            stop("Invalid p-values found after gene-level summarization")
        }

        # Determine the thetas to be used
        theta <- unique(sort(pGene))

        # Compute q-values based on method choice
        if (method == "exact") {
            q <- perGeneQValueExact(pGene, theta, geneSplit)
        } else if (method == "simulation") {
            q <- perGeneQValueBySimulation(pGene, theta, geneSplit, nperm)
        } else {
            stop("Invalid method specified. Use 'exact' or 'simulation'")
        }

        # Return a named vector of q-values per gene
        res <- rep(NA_real_, length(pGene))
        res <- q[match(pGene, theta)]
        res <- pmin(1, res)
        names(res) <- names(geneSplit)

        if (any(is.na(res))) {
            warning("NA values found in final q-values")
        }

        return(res)
    }

    perGeneQValueExact <- function(pGene, theta, geneSplit) {
        stopifnot(length(pGene) == length(geneSplit))

        # Compute the numerator
        numTranscripts <- lengths(geneSplit)
        tab <- tabulate(numTranscripts)
        notZero <- (tab > 0)

        numerator <- mapply(
            function(m, n) m * (1 - (1 - theta)^n),
            m = tab[notZero],
            n = which(notZero)
        )
        numerator <- rowSums(matrix(numerator, ncol = sum(notZero)))

        # Compute the denominator
        bins <- cut(pGene, breaks = c(-Inf, theta), right = TRUE, include.lowest = TRUE)
        counts <- tabulate(bins, nbins = nlevels(bins))
        denom <- cumsum(counts)

        stopifnot(denom[length(denom)] == length(pGene))

        return(numerator / denom)
    }

    perGeneQValueBySimulation <- function(pGene, theta, geneSplit, nperm = 24) {
        nr <- sum(lengths(geneSplit))

        # Generate random p-values and compute minimum per gene
        pRand <- apply(
            matrix(runif(nr * nperm), nrow = nr, ncol = nperm), 2,
            function(p) sapply(geneSplit, function(i) min(p[i]))
        )

        stopifnot(nrow(pRand) == length(pGene), ncol(pRand) == nperm)

        # Compute histograms
        hTest <- hist(pGene, breaks = c(theta, +Inf), plot = FALSE)
        hRand <- hist(pRand, breaks = c(theta, +Inf), plot = FALSE)

        stopifnot(
            sum(hTest$counts) == length(pGene),
            sum(hRand$counts) == length(pRand)
        )

        # Compute q-values
        numPos <- cumsum(hTest$counts)
        numFalsePos <- cumsum(hRand$counts) / nperm

        return(numFalsePos / numPos)
    }

    # Example usage:
    qvals <- perGeneQValue(results_df, method = "exact") %>%
        as.data.frame() %>%
        # colname is group1_vs_group2_noisq_qval
        rename(!!paste0(group1, "_vs_", group2, "_dtu_gene_padj") := ".") %>%
        rownames_to_column("geneID")
    results_df <- results_df %>%
        left_join(qvals, by = "geneID") %>%
        # Remove version number
        mutate(geneID = sub("\\..*", "", geneID)) %>%
        # rename log2FC to group1_vs_group2_dtu_log2fold
        rename(!!paste0(group1, "_vs_", group2, "_dtu_log2fold") := log2FC) %>%
        # add dtu_ prefix to column paste0(group1, "_mean") and paste0(group2, "_mean")
        rename(!!paste0("dtu_", group1, "_mean") := paste0(group1, "_mean")) %>%
        rename(!!paste0("dtu_", group2, "_mean") := paste0(group2, "_mean")) %>%
        # move geneID, transcript_id to the front
        dplyr::select(geneID, transcript_id, everything())
    return(results_df)
}
