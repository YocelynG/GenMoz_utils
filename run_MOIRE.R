#!/usr/bin/env Rscript

# Run MOIRE analyses on all allele_*.txt files in a directory.
#
# Input files may be tab- or comma-separated, with either header format:
#   sampleID  locus  pseudo_cigar  reads  norm.reads.locus
#   "SampleID","Locus","PseudoCIGAR","Reads","WSAF"
#
# Usage:
#   Rscript run_moire.R [input_directory] [output_directory] [threads]
#
# Example:
#   Rscript run_moire.R \
#     /home/ygutierrez/Documents/CSP/HFS_data/MOIRE_Dcifer_HFS22 \
#     /home/ygutierrez/Documents/CSP/HFS_data/MOIRE_Dcifer_HFS22 \
#     20

args <- commandArgs(trailingOnly = TRUE)

default_dir <- "/home/gmatambisso/Documents/Moire_ANC2022/allele_tables_by_district"

input_dir <- if (length(args) >= 1L) args[[1L]] else default_dir
output_dir <- if (length(args) >= 2L) args[[2L]] else input_dir
n_threads <- if (length(args) >= 3L) suppressWarnings(as.integer(args[[3L]])) else 20L

if (!dir.exists(input_dir)) {
  stop("Input directory does not exist: ", input_dir, call. = FALSE)
}

if (is.na(n_threads) || n_threads < 1L) {
  stop("The number of threads must be a positive integer.", call. = FALSE)
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

required_packages <- c("dplyr", "moire")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0L) {
  stop(
    "Missing required R package(s): ",
    paste(missing_packages, collapse = ", "),
    "\nInstall them before running this script.",
    call. = FALSE
  )
}

suppressPackageStartupMessages(library(dplyr))

pops <- list.files(
  input_dir,
  pattern = "^allele_.*\\.txt$",
  full.names = TRUE
)

if (length(pops) == 0L) {
  stop(
    "No files matching '^allele_.*\\.txt$' were found in: ",
    input_dir,
    call. = FALSE
  )
}

burnin <- 1e4
num_samples <- 1e4
pt_chains <- seq(1, 0.5, length.out = n_threads)
failed_files <- character(0)

message("Found ", length(pops), " input file(s).")
message("Output directory: ", normalizePath(output_dir, mustWork = FALSE))
message("Parallel-tempering chains/threads: ", n_threads)

for (pop in pops) {
  message("\n[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] Processing: ", basename(pop))

  success <- tryCatch({
    # Detect the field separator (tab or comma) from the header line.
    header_line <- readLines(pop, n = 1L, warn = FALSE)
    field_sep <- if (grepl("\t", header_line, fixed = TRUE)) "\t" else ","

    moire_input <- read.table(
      pop,
      header = TRUE,
      sep = field_sep,
      quote = "\"",
      comment.char = "",
      stringsAsFactors = FALSE,
      check.names = FALSE
    )

    # Supported header formats, mapped to the internal column names:
    #   sampleID, locus, pseudo_cigar, reads, norm.reads.locus
    #   "SampleID","Locus","PseudoCIGAR","Reads","WSAF"
    column_aliases <- list(
      SampleID = c("SampleID", "sampleID"),
      Locus = c("Locus", "locus"),
      PseudoCIGAR = c("PseudoCIGAR", "pseudo_cigar"),
      Reads = c("Reads", "reads"),
      norm.Reads.Locus = c("WSAF", "norm.reads.locus", "norm.Reads.Locus")
    )

    matched_columns <- vapply(
      column_aliases,
      function(aliases) {
        hit <- intersect(aliases, colnames(moire_input))
        if (length(hit) > 0L) hit[[1L]] else NA_character_
      },
      character(1)
    )

    if (anyNA(matched_columns)) {
      missing_columns <- vapply(
        column_aliases[is.na(matched_columns)],
        paste,
        character(1),
        collapse = " or "
      )
      stop(
        "Missing required column(s): ",
        paste(missing_columns, collapse = "; "),
        "\n  Columns found: ",
        paste(colnames(moire_input), collapse = ", ")
      )
    }

    moire_input <- moire_input %>%
      select(all_of(matched_columns)) %>%
      mutate(
        SampleID = as.character(SampleID),
        Locus = as.character(Locus),
        PseudoCIGAR = as.character(PseudoCIGAR),
        Reads = as.numeric(Reads),
        norm.Reads.Locus = as.numeric(norm.Reads.Locus)
      )

    # Remove masking from the pseudo-CIGAR strings.
    moire_input$PseudoCIGAR <- gsub(
      "\\d+\\+[^N]*N",
      "",
      moire_input$PseudoCIGAR
    )

    # Empty pseudo-CIGAR strings represent the reference allele.
    moire_input$PseudoCIGAR <- ifelse(
      moire_input$PseudoCIGAR == "" | is.na(moire_input$PseudoCIGAR),
      ".",
      moire_input$PseudoCIGAR
    )

    # Aggregate alleles that become identical after masking is removed.
    moire_input_agg <- moire_input %>%
      group_by(SampleID, Locus, PseudoCIGAR) %>%
      summarise(
        reads = sum(Reads, na.rm = TRUE),
        norm.Reads.Locus = sum(norm.Reads.Locus, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(Allele = paste0(Locus, "__", PseudoCIGAR))

    colnames(moire_input_agg) <- tolower(colnames(moire_input_agg))
    colnames(moire_input_agg)[colnames(moire_input_agg) == "sampleid"] <- "sample_id"
    moire_input_agg <- as.data.frame(moire_input_agg)

    message(
      "  Samples: ", dplyr::n_distinct(moire_input_agg$sample_id),
      "; loci: ", dplyr::n_distinct(moire_input_agg$locus),
      "; allele records: ", nrow(moire_input_agg)
    )

    dat_filter <- moire::load_long_form_data(moire_input_agg)

    mcmc_results <- moire::run_mcmc(
      dat_filter,
      is_missing = dat_filter$is_missing,
      verbose = TRUE,
      burnin = burnin,
      samples_per_chain = num_samples,
      pt_chains = pt_chains,
      pt_num_threads = n_threads,
      thin = 10
    )

    coi_results <- moire::summarize_coi(mcmc_results)
    ecoi_results <- moire::summarize_effective_coi(mcmc_results)

    # Combine standard COI and effective COI summaries by sample.
    all_coi_results <- merge(
      coi_results,
      ecoi_results,
      by = "sample_id",
      all = TRUE,
      suffixes = c("_coi", "_effective_coi")
    )

    all_coi_results$sample_id <- as.character(all_coi_results$sample_id)
    all_coi_results$sample_id <- ifelse(
      !is.na(all_coi_results$sample_id) &
        !grepl("\\.", all_coi_results$sample_id),
      paste0(all_coi_results$sample_id, ".0"),
      all_coi_results$sample_id
    )

    af_results <- moire::summarize_allele_freqs(mcmc_results)
    he_results <- moire::summarize_he(mcmc_results)

    coi_results_filename <- sub(
      "\\.txt$",
      "_COI_Results.txt",
      basename(pop)
    )
    af_results_filename <- sub(
      "\\.txt$",
      "_AF_Results.txt",
      basename(pop)
    )
    he_results_filename <- sub(
      "\\.txt$",
      "_He_Results.txt",
      basename(pop)
    )

    write.table(
      all_coi_results,
      file.path(output_dir, coi_results_filename),
      row.names = FALSE,
      quote = FALSE,
      sep = "\t"
    )
    write.table(
      af_results,
      file.path(output_dir, af_results_filename),
      row.names = FALSE,
      quote = FALSE,
      sep = "\t"
    )
    write.table(
      he_results,
      file.path(output_dir, he_results_filename),
      row.names = FALSE,
      quote = FALSE,
      sep = "\t"
    )

    message("  Completed: ", basename(pop))
    TRUE
  }, error = function(e) {
    message("  ERROR: ", conditionMessage(e))
    FALSE
  })

  if (!success) {
    failed_files <- c(failed_files, basename(pop))
  }
}

if (length(failed_files) > 0L) {
  message(
    "\nCompleted with errors in ", length(failed_files), " file(s):\n  ",
    paste(failed_files, collapse = "\n  ")
  )
  quit(status = 1L)
}

message("\nAll MOIRE analyses completed successfully.")
