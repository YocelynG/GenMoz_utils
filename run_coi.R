#!/usr/bin/env Rscript

# ============================================================
# Read-based COI from allele tables
# ============================================================
#
# For each sample and locus, alleles are ranked by reads and kept
# until they cover <threshold> of the locus reads (the allele that
# crosses the threshold is kept). The COI of a sample is the
# largest number of kept alleles at any locus.
#
# Accepted allele table formats:
#   .txt / .tsv / .csv, comma- or tab-separated (detected from
#   the header line). Filenames must start with "allele"
#   (case-insensitive) when a directory is given.
#
# Required allele table columns (same aliases as run_moire.R):
#   SampleID / sampleID / sample_id
#   Locus / locus
#   PseudoCIGAR / pseudo_cigar
#   Reads / reads
#
#
# Usage:
#
# Rscript run_coi.R --input <allele_input> [--output <output_dir>]
#                   [--threshold <value>] [--count_reference]
#
# Options:
#
#   -i, --input          Allele table file, or a directory of allele*
#                        files. Required.
#
#   -o, --output         Output directory (default: the input directory).
#
#   -t, --threshold      Cumulative read proportion covered by the kept
#                        alleles (default 0.98).
#
#   --count_reference    Count PseudoCIGAR "." as the reference allele.
#                        Without this flag, "." is treated as a missing
#                        allele and removed before counting.
#
#   -h, --help           Show the help message.
#
#
# Outputs generated for each allele table:
#
#   *_COI_Table.txt
#     SampleID, COI, avg_alleles_per_locus (tab-separated,
#     one row per sample)
#
#
# Examples:
#
#   Rscript run_coi.R --input allele_tables_by_district --output coi_output
#   Rscript run_coi.R -i allele_tables_by_district -o coi_output -t 0.95 --count_reference


# ============================================================
# Command-line arguments
# ============================================================

if (!requireNamespace("optparse", quietly = TRUE)) {
  stop(
    "Missing required R package: optparse",
    "\nInstall it before running this script.",
    call. = FALSE
  )
}

option_list <- list(
  optparse::make_option(
    c("-i", "--input"),
    help = "Allele table file, or a directory of allele* files. Required."
  ),
  optparse::make_option(
    c("-o", "--output"),
    help = "Output directory [default: the input directory]"
  ),
  optparse::make_option(
    c("-t", "--threshold"),
    type = "double",
    default = 0.98,
    help = "Cumulative read proportion covered by the kept alleles [default: %default]"
  ),
  optparse::make_option(
    "--count_reference",
    action = "store_true",
    default = FALSE,
    help = "Count PseudoCIGAR \".\" as the reference allele instead of removing it"
  )
)

opt <- optparse::parse_args(
  optparse::OptionParser(
    usage = "Rscript run_coi.R --input <allele_input> [options]",
    option_list = option_list
  )
)

if (is.null(opt$input)) {
  stop(
    "--input is required. Run 'Rscript run_coi.R --help' for usage.",
    call. = FALSE
  )
}

allele_input <- opt$input

output_dir <- if (!is.null(opt$output)) {
  opt$output
} else if (dir.exists(allele_input)) {
  allele_input
} else {
  dirname(allele_input)
}

threshold <- opt$threshold

count_reference <- opt$count_reference


# ============================================================
# Basic checks
# ============================================================

if (!file.exists(allele_input)) {
  stop("Allele input does not exist: ", allele_input, call. = FALSE)
}

if (is.na(threshold) || threshold <= 0 || threshold > 1) {
  stop("The threshold must be a number in (0, 1].", call. = FALSE)
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)


# ============================================================
# Required packages
# ============================================================

required_packages <- c(
  "dplyr"
)

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

suppressPackageStartupMessages(
  library(dplyr)
)


# ============================================================
# Find input files
# ============================================================

pops <- if (dir.exists(allele_input)) {
  allele_files <- list.files(
    allele_input,
    pattern = "^allele.*\\.(txt|tsv|csv)$",
    full.names = TRUE,
    ignore.case = TRUE
  )

  # Result files from run_moire.R, run_coi.R and run_dcifer.R also
  # start with "allele", so skip them when they share the folder.
  allele_files[!grepl("(_Results|_Sample_Summary|_COI_Table)\\.txt$", allele_files)]
} else {
  allele_input
}

if (length(pops) == 0L) {
  stop(
    "No files matching allele*.txt, allele*.tsv, or allele*.csv ",
    "(case-insensitive) were found in: ",
    allele_input,
    call. = FALSE
  )
}

failed_files <- character(0)

message("Found ", length(pops), " allele table(s).")
message("Output directory: ", normalizePath(output_dir, mustWork = FALSE))
message("Read threshold: ", threshold)
message(
  "PseudoCIGAR \".\": ",
  if (count_reference) "counted as reference allele" else "removed as missing"
)


# ============================================================
# Helper function: normalize column names
# ============================================================

# Same rule as run_moire.R: lower case, drop everything that is
# not a letter or digit (SampleID / sample_id -> sampleid).

normalize_column_name <- function(x) {
  gsub("[^a-z0-9]", "", tolower(x))
}


# ============================================================
# Helper function: read a comma- or tab-separated table
# ============================================================

# The separator is taken from the header line rather than the
# extension, since some .txt allele tables are comma-separated.

read_delim_auto <- function(path) {
  header <- readLines(path, n = 1L, warn = FALSE)

  separator <- if (
    lengths(regmatches(header, gregexpr("\t", header))) >=
      lengths(regmatches(header, gregexpr(",", header)))
  ) {
    "\t"
  } else {
    ","
  }

  read.table(
    path,
    header = TRUE,
    sep = separator,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    quote = "\"",
    comment.char = "",
    fill = TRUE,
    colClasses = "character"
  )
}


# ============================================================
# Helper function: find a column by alias
# ============================================================

find_column <- function(df, aliases, label, path) {
  matches <- which(normalize_column_name(colnames(df)) %in% aliases)

  if (length(matches) == 0L) {
    stop(
      "Missing required column '", label, "' in ", basename(path),
      "\n\nColumns present in file:\n  ",
      paste(colnames(df), collapse = ", "),
      call. = FALSE
    )
  }

  colnames(df)[matches[1L]]
}


# ============================================================
# Helper function: missing alleles
# ============================================================

is_missing_allele <- function(x) {
  missing <- is.na(x) | x == ""

  if (!count_reference) {
    missing <- missing | x == "."
  }

  missing
}


# ============================================================
# Process each allele table
# ============================================================

for (pop in pops) {

  message(
    "\n[",
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    "] Processing: ",
    basename(pop)
  )

  success <- tryCatch({

    base_filename <- tools::file_path_sans_ext(basename(pop))


    # ========================================================
    # Allele table
    # ========================================================

    allele_raw <- read_delim_auto(pop)

    sample_col <- find_column(allele_raw, "sampleid", "SampleID", pop)
    locus_col <- find_column(allele_raw, "locus", "Locus", pop)
    cigar_col <- find_column(allele_raw, "pseudocigar", "PseudoCIGAR", pop)
    reads_col <- find_column(allele_raw, "reads", "Reads", pop)

    message(
      "  Detected column mapping:\n",
      "    SampleID <- ", sample_col, "\n",
      "    Locus <- ", locus_col, "\n",
      "    PseudoCIGAR <- ", cigar_col, "\n",
      "    Reads <- ", reads_col
    )

    df_cluster <- tibble::tibble(
      sample = allele_raw[[sample_col]],
      Locus = allele_raw[[locus_col]],
      Allele = allele_raw[[cigar_col]],
      Reads = suppressWarnings(as.numeric(allele_raw[[reads_col]]))
    ) %>%
      filter(!is.na(sample), sample != "")

    all_samples <- unique(df_cluster$sample)


    # ========================================================
    # Alleles kept per locus
    # ========================================================

    df_kept_per_locus <- df_cluster %>%
      filter(!is_missing_allele(Allele)) %>%
      group_by(sample, Locus, Allele) %>%
      summarise(Reads = sum(Reads, na.rm = TRUE), .groups = "drop") %>%
      group_by(sample, Locus) %>%
      arrange(desc(Reads), Allele, .by_group = TRUE) %>%
      mutate(
        total_reads = sum(Reads, na.rm = TRUE),
        prop = ifelse(total_reads > 0, Reads / total_reads, 0),
        cumprop = cumsum(prop),
        cumprop_before = lag(cumprop, default = 0),
        keep = cumprop_before < threshold
      ) %>%
      filter(keep) %>%
      summarise(n_alleles = n_distinct(Allele), .groups = "drop")


    # ========================================================
    # COI table
    # ========================================================

    df_coi_table <- df_kept_per_locus %>%
      group_by(sample) %>%
      summarise(
        COI = max(n_alleles),
        avg_alleles_per_locus = mean(n_alleles),
        .groups = "drop"
      )

    no_coi <- setdiff(all_samples, df_coi_table$sample)

    if (length(no_coi) > 0L) {
      message(
        "  WARNING: ", length(no_coi),
        " sample(s) have no alleles left and get no COI: ",
        paste(no_coi, collapse = ", ")
      )
    }

    message(
      "  Samples: ", nrow(df_coi_table),
      "; loci: ", n_distinct(df_kept_per_locus$Locus),
      "; COI range: ", paste(range(df_coi_table$COI), collapse = "-"),
      "; polyclonal (COI > 1): ", sum(df_coi_table$COI > 1)
    )


    # ========================================================
    # Save COI table
    # ========================================================

    coi_table_filename <- paste0(base_filename, "_COI_Table.txt")

    write.table(
      df_coi_table %>% select(SampleID = sample, COI, avg_alleles_per_locus),
      file = file.path(output_dir, coi_table_filename),
      row.names = FALSE,
      quote = FALSE,
      sep = "\t"
    )

    message("  Output file:\n    ", coi_table_filename)

    message("  Completed: ", basename(pop))


    TRUE


  }, error = function(e) {

    message("  ERROR: ", conditionMessage(e))

    FALSE
  })


  # ==========================================================
  # Record failed files
  # ==========================================================

  if (!success) {
    failed_files <- c(failed_files, basename(pop))
  }
}


# ============================================================
# Final status
# ============================================================

if (length(failed_files) > 0L) {

  message(
    "\nCompleted with errors in ",
    length(failed_files),
    " file(s):\n  ",
    paste(failed_files, collapse = "\n  ")
  )

  quit(status = 1L)
}


message("\nAll COI tables completed successfully.")
