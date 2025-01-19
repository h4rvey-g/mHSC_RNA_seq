library(targets)
library(crew)
library(tidyverse)
library(tarchetypes)
tar_source(c(
    "codes/R/101.load_data.R",
    "codes/R/102.DESeq.R",
    "codes/R/103.intersect.R",
    "codes/R/104.isoformswitch.R"
))
tar_option_set(
    tidy_eval = FALSE,
    packages <- c(
        "tidyverse", "patchwork", "tidyplots", "SummarizedExperiment", "tximeta", "tidySummarizedExperiment",
        "tarchetypes", "tidybulk", "genekitr", "ggVennDiagram", "IsoformSwitchAnalyzeR"
    ),
    controller = crew_controller_local(workers = 20, seconds_timeout = 6000),
    format = "qs",
    storage = "worker", retrieval = "worker",
    seed = 42
)
tar_config_set(
    script = "codes/_targets.R",
    store = "data/_targets"
)
options(max.print = 12, spe = "mouse")
group_comparison <- tibble(
    group1 = c("KO4", "KO5", "KO5"),
    group2 = c("WT", "WT", "KO4")
)
mapped <- tar_map(
    values = group_comparison,
    names = c("group1", "group2"), # Use both columns for naming
    # Target for loading DTU data
    tar_target(dtu_data, load_dtu(group1, group2)),
    tar_target(exon_data, load_exon(group1, group2)),
    tar_target(deg_data, run_DEG(se_data, group1, group2)),
    tar_target(intersected_data, intersect_data(dtu_data, exon_data, deg_data,
        p_cutoff = 0.05,
        log2fc_cutoff = 1, group1, group2
    )),
    tar_target(isoform_switch, run_isoformswitch(group1, group2, salmon_dir, gtf_file, transcript_fasta,
        alpha = 0.05, dIFcutoff = 0.1, intersected_data
    ))
)
# Combine DTU results into a named list
combined_dtu <- tar_combine(
    all_dtu_results,
    mapped[["dtu_data"]],
    command = list(!!!.x)
)
combined_exon <- tar_combine(
    all_exon_results,
    mapped[["exon_data"]],
    command = list(!!!.x)
)
combined_deg <- tar_combine(
    all_deg_results,
    mapped[["deg_data"]],
    command = list(!!!.x)
)
combined_isoform <- tar_combine(
    all_isoform_results,
    mapped[["isoform_switch"]],
    command = list(!!!.x)
)

# Return the complete pipeline
list(
    tar_target(salmon_dir, "data/102.splice_workflow/star_salmon/salmon", format = "file"),
    tar_target(gtf_file, "data/reference_data/alias/mm10/gencode_gtf/default/mm10.gtf", format = "file"),
    tar_target(transcript_fasta, "data/102.splice_workflow/genome/genome.transcripts.fa", format = "file"),
    tar_target(se_data, load_data()),
    mapped,
    combined_dtu,
    combined_exon,
    combined_deg,
    combined_isoform
)
