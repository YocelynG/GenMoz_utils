#!/usr/bin/env Rscript

# ============================================================
# Run Dcifer between-host relatedness (IBD) on allele tables
# ============================================================
#
# Wraps PGEcore::dcifer_ibd_wrapper(). Each allele table is
# converted to the PGEcore input format and, optionally, paired
# with the MOIRE allele frequencies and COI produced by
# run_moire.R.
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
#
# Other columns (Reads, WSAF, ...) are ignored: Dcifer only uses
# allele presence/absence.
#
#
# Usage:
#
# Rscript run_dcifer.R --input <allele_input> --output <output_dir> [options]
#
# Options:
#
#   -i, --input         Allele table file, or a directory of allele*
#                       files. Required.
#
#   -o, --output        Output directory. Required.
#
#   -a, --allele_freq   MOIRE *_AF_Results.txt file, or a directory
#                       containing them. Omitted: Dcifer estimates
#                       allele frequencies from the allele table.
#
#   -c, --coi           COI table, or a directory containing them:
#                       MOIRE *_COI_Results.txt (sample_id, post_coi_med)
#                       or run_coi.R *_COI_Table.txt (SampleID, COI).
#                       Omitted: naive COI, the 2nd-highest allele
#                       count across loci (as dcifer::getCOI).
#
#   -m, --metadata      Metadata file with SampleID and <pop_col>.
#                       When given, allele frequencies are estimated
#                       per population and both within- and
#                       between-population pairs are run. Cannot be
#                       combined with --allele_freq (PGEcore
#                       restriction).
#
#   -p, --pop_col       Metadata column defining populations
#                       (default "district"). Only used with --metadata.
#
#   -t, --threads       Parallel workers (default 4).
#
#   -h, --help          Show the help message.
#
# When --allele_freq / --coi is a directory, the file for each allele
# table is found by name, as written by run_moire.R and run_coi.R
#
#
# Outputs generated for each allele table:
#
#   *_Dcifer_IBD_Results.txt
#     specimen_name_a, specimen_name_b, btwn_host_rel,
#     p_value, CI_lower, CI_upper
#
#   dcifer_inputs/*_allele_table.tsv, *_coi.tsv (and
#   *_allele_freq.tsv, *_metadata.tsv when used): the exact
#   files passed to PGEcore, kept for reproducibility.
#
#
# Examples:
#
# Per district, with MOIRE frequencies and COI:
#   Rscript run_dcifer.R --input allele_tables_by_district \
#     --allele_freq Moire_output --coi Moire_output \
#     --output dcifer_output --threads 4
#

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
    help = "Output directory. Required."
  ),
  optparse::make_option(
    c("-a", "--allele_freq"),
    help = paste0(
      "MOIRE *_AF_Results.txt file, or a directory containing them ",
      "[default: estimated by Dcifer]"
    )
  ),
  optparse::make_option(
    c("-c", "--coi"),
    help = paste0(
      "MOIRE *_COI_Results.txt or run_coi.R *_COI_Table.txt file, or a ",
      "directory containing them [default: naive COI from allele counts]"
    )
  ),
  optparse::make_option(
    c("-m", "--metadata"),
    help = paste0(
      "Metadata file with SampleID and --pop_col; runs within- and ",
      "between-population pairs. Cannot be combined with --allele_freq"
    )
  ),
  optparse::make_option(
    c("-p", "--pop_col"),
    default = "district",
    help = "Metadata column defining populations [default: %default]"
  ),
  optparse::make_option(
    c("-t", "--threads"),
    type = "integer",
    default = 4L,
    help = "Parallel workers [default: %default]"
  )
)

opt <- optparse::parse_args(
  optparse::OptionParser(
    usage = "Rscript run_dcifer.R --input <allele_input> --output <output_dir> [options]",
    option_list = option_list
  )
)

if (is.null(opt$input) || is.null(opt$output)) {
  stop(
    "--input and --output are required. ",
    "Run 'Rscript run_dcifer.R --help' for usage.",
    call. = FALSE
  )
}

# Options left out are NULL; NA marks "not provided" below
optional_arg <- function(x) {
  if (is.null(x)) NA_character_ else x
}

allele_input <- opt$input
moire_af_input <- optional_arg(opt$allele_freq)
moire_coi_input <- optional_arg(opt$coi)
metadata_path <- optional_arg(opt$metadata)
output_dir <- opt$output
n_threads <- opt$threads
pop_col <- opt$pop_col


# ============================================================
# Dcifer parameters
# ============================================================

# MOIRE COI summary column used as the Dcifer COI (rounded to an
# integer). Dcifer raises it per locus to the number of observed
# alleles when it is lower.
coi_column <- "post_coi_med"

# MOIRE allele frequency summary column; renormalized per locus.
af_column <- "post_allele_freqs_mean"

rnull <- 0
alpha <- 0.05


# ============================================================
# Basic checks
# ============================================================

if (!file.exists(allele_input)) {
  stop("Allele input does not exist: ", allele_input, call. = FALSE)
}

for (optional_input in c(moire_af_input, moire_coi_input, metadata_path)) {
  if (!is.na(optional_input) && !file.exists(optional_input)) {
    stop("Input does not exist: ", optional_input, call. = FALSE)
  }
}

if (is.na(n_threads) || n_threads < 1L) {
  stop("The number of threads must be a positive integer.", call. = FALSE)
}

if (!is.na(moire_af_input) && !is.na(metadata_path)) {
  stop(
    "--allele_freq and --metadata cannot be used together: with ",
    "metadata, PGEcore estimates allele frequencies per population. ",
    "Leave out one of them.",
    call. = FALSE
  )
}

input_prep_dir <- file.path(output_dir, "dcifer_inputs")

dir.create(input_prep_dir, recursive = TRUE, showWarnings = FALSE)


# ============================================================
# Required packages
# ============================================================

required_packages <- c(
  "dplyr",
  "PGEcore",
  "dcifer",
  "foreach",
  "iterators",
  "doParallel",
  "parallelly"
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

for (moire_input in c(moire_af_input, moire_coi_input)) {
  if (!is.na(moire_input) && !dir.exists(moire_input) && length(pops) > 1L) {
    stop(
      "A single MOIRE file was given for ", length(pops), " allele tables: ",
      moire_input, "\nPass the MOIRE output directory instead.",
      call. = FALSE
    )
  }
}

failed_files <- character(0)

message("Found ", length(pops), " allele table(s).")
message("Output directory: ", normalizePath(output_dir, mustWork = FALSE))
message("Parallel workers: ", n_threads)
message(
  "Allele frequencies: ",
  if (!is.na(moire_af_input)) {
    paste0("MOIRE (", af_column, ")")
  } else if (!is.na(metadata_path)) {
    paste0("estimated by Dcifer per '", pop_col, "'")
  } else {
    "estimated by Dcifer"
  }
)
message(
  "COI: ",
  if (!is.na(moire_coi_input)) {
    paste0("from ", moire_coi_input)
  } else {
    "naive, from allele counts"
  }
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
# Helper function: locate the AF / COI file for an allele table
# ============================================================

# Tries each suffix in turn and returns the first file that exists.

resolve_moire_file <- function(moire_input, base_filename, suffixes) {
  if (is.na(moire_input)) {
    return(NA_character_)
  }

  if (!dir.exists(moire_input)) {
    return(moire_input)
  }

  paths <- file.path(moire_input, paste0(base_filename, suffixes))
  found <- paths[file.exists(paths)]

  if (length(found) == 0L) {
    stop(
      "File not found, tried: ", paste(paths, collapse = ", "),
      call. = FALSE
    )
  }

  found[1L]
}


# ============================================================
# Helper function: write tab-delimited Dcifer input
# ============================================================

write_dcifer_input <- function(data, path) {
  write.table(
    data,
    file = path,
    row.names = FALSE,
    quote = FALSE,
    sep = "\t"
  )

  path
}


# ============================================================
# Helper function: stop any parallel cluster left by a failure
# ============================================================

stop_default_cluster <- function() {
  cl <- parallel::getDefaultCluster()

  if (!is.null(cl)) {
    try(parallel::stopCluster(cl), silent = TRUE)
    parallel::setDefaultCluster(NULL)
  }
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

    message(
      "  Detected column mapping:\n",
      "    SampleID <- ", sample_col, "\n",
      "    Locus <- ", locus_col, "\n",
      "    PseudoCIGAR <- ", cigar_col
    )

    # Masking is removed exactly as in run_moire.R, so alleles that
    # differ only in masked regions are collapsed and allele names
    # match the MOIRE allele frequency table ("." = reference).
    allele_table <- tibble::tibble(
      specimen_name = allele_raw[[sample_col]],
      target_name = allele_raw[[locus_col]],
      seq = gsub("\\d+\\+[^N]*N", "", allele_raw[[cigar_col]])
    ) %>%
      filter(
        !is.na(specimen_name), specimen_name != "",
        !is.na(target_name), target_name != ""
      ) %>%
      mutate(seq = ifelse(is.na(seq) | seq == "", ".", seq)) %>%
      distinct()


    # ========================================================
    # Metadata (population mode)
    # ========================================================

    metadata_file <- NA_character_

    if (!is.na(metadata_path)) {
      metadata_raw <- read_delim_auto(metadata_path)

      meta_sample_col <- find_column(
        metadata_raw, "sampleid", "SampleID", metadata_path
      )

      if (!pop_col %in% colnames(metadata_raw)) {
        stop(
          "Population column '", pop_col, "' not found in ",
          basename(metadata_path),
          call. = FALSE
        )
      }

      metadata <- tibble::tibble(
        specimen_name = metadata_raw[[meta_sample_col]],
        population = metadata_raw[[pop_col]]
      ) %>%
        filter(!is.na(population), population != "") %>%
        distinct(specimen_name, .keep_all = TRUE)

      # Samples without a population (e.g. lab controls) cannot be
      # assigned allele frequencies, so they are dropped.
      no_metadata <- setdiff(allele_table$specimen_name, metadata$specimen_name)

      if (length(no_metadata) > 0L) {
        message(
          "  Dropping ", length(no_metadata),
          " sample(s) without a '", pop_col, "' value in the metadata: ",
          paste(no_metadata, collapse = ", ")
        )

        allele_table <- allele_table %>%
          filter(!specimen_name %in% no_metadata)
      }

      metadata_file <- write_dcifer_input(
        metadata %>%
          filter(specimen_name %in% allele_table$specimen_name) %>%
          rename(!!pop_col := population),
        file.path(input_prep_dir, paste0(base_filename, "_metadata.tsv"))
      )
    }

    samples <- unique(allele_table$specimen_name)


    # ========================================================
    # COI from MOIRE or run_coi.R
    # ========================================================

    coi_path <- resolve_moire_file(
      moire_coi_input, base_filename, c("_COI_Results.txt", "_COI_Table.txt")
    )

    if (!is.na(coi_path)) {
      coi_raw <- read_delim_auto(coi_path)

      if (all(c("sample_id", coi_column) %in% colnames(coi_raw))) {
        coi_ids <- coi_raw$sample_id
        coi_values <- coi_raw[[coi_column]]
        coi_source <- coi_column
      } else if (all(c("SampleID", "COI") %in% colnames(coi_raw))) {
        coi_ids <- coi_raw$SampleID
        coi_values <- coi_raw$COI
        coi_source <- "COI"
      } else {
        stop(
          "Expected columns sample_id and ", coi_column, " (MOIRE) or ",
          "SampleID and COI (run_coi.R) in ", basename(coi_path),
          call. = FALSE
        )
      }

      # run_moire.R appends ".0" to sample IDs without a "."; undo it
      # for IDs that only match the allele table without the suffix.
      coi <- tibble::tibble(
        specimen_name = ifelse(
          coi_ids %in% samples,
          coi_ids,
          sub("\\.0$", "", coi_ids)
        ),
        coi = pmax(1L, as.integer(round(as.numeric(coi_values))))
      ) %>%
        filter(specimen_name %in% samples)

      missing_coi <- setdiff(samples, coi$specimen_name)

      if (length(missing_coi) > 0L) {
        stop(
          length(missing_coi), " sample(s) in the allele table have no COI in ",
          basename(coi_path), ": ",
          paste(head(missing_coi, 10L), collapse = ", "),
          if (length(missing_coi) > 10L) ", ...",
          call. = FALSE
        )
      }

      message("  COI: ", basename(coi_path), " (", coi_source, ")")

    } else {

      # Naive COI computed here instead of by dcifer::getCOI(), which
      # names its result "<sample>.<locus>" and so no longer matches
      # the sample IDs. Same rule: 2nd-highest allele count per sample.
      coi <- allele_table %>%
        count(specimen_name, target_name, name = "n_alleles") %>%
        group_by(specimen_name) %>%
        summarise(
          coi = sort(n_alleles, decreasing = TRUE)[min(2L, n())],
          .groups = "drop"
        )

      message("  COI: naive (2nd-highest allele count across loci)")
    }

    coi_file <- write_dcifer_input(
      coi,
      file.path(input_prep_dir, paste0(base_filename, "_coi.tsv"))
    )


    # ========================================================
    # Allele frequencies from MOIRE
    # ========================================================

    af_file <- NA_character_
    af_path <- resolve_moire_file(
      moire_af_input, base_filename, "_AF_Results.txt"
    )

    if (!is.na(af_path)) {
      af_raw <- read_delim_auto(af_path)

      if (!all(c("locus", "allele", af_column) %in% colnames(af_raw))) {
        stop(
          "Expected columns locus, allele and ", af_column, " in ",
          basename(af_path),
          call. = FALSE
        )
      }

      # MOIRE alleles are "<locus>__<masked PseudoCIGAR>"
      allele_freq <- tibble::tibble(
        target_name = af_raw$locus,
        seq = substring(af_raw$allele, nchar(af_raw$locus) + 3L),
        freq = as.numeric(af_raw[[af_column]])
      ) %>%
        group_by(target_name) %>%
        mutate(freq = freq / sum(freq)) %>%
        ungroup()

      # MOIRE leaves out loci with a single allele in the data. They
      # carry no relatedness information, so they are dropped here
      # too; a missing polymorphic locus means a mismatched file.
      alleles_per_locus <- allele_table %>%
        distinct(target_name, seq) %>%
        count(target_name, name = "n_alleles") %>%
        filter(!target_name %in% allele_freq$target_name)

      if (any(alleles_per_locus$n_alleles > 1L)) {
        stop(
          "Polymorphic loci missing from ", basename(af_path), ": ",
          paste(
            alleles_per_locus$target_name[alleles_per_locus$n_alleles > 1L],
            collapse = ", "
          ),
          call. = FALSE
        )
      }

      if (nrow(alleles_per_locus) > 0L) {
        message(
          "  Dropping ", nrow(alleles_per_locus),
          " monomorphic locus/loci not in the MOIRE frequencies"
        )

        allele_table <- allele_table %>%
          filter(!target_name %in% alleles_per_locus$target_name)
      }

      table_alleles <- paste(allele_table$target_name, allele_table$seq, sep = ":")
      af_alleles <- paste(allele_freq$target_name, allele_freq$seq, sep = ":")
      not_in_af <- setdiff(table_alleles, af_alleles)
      not_in_table <- setdiff(af_alleles, table_alleles)

      if (length(not_in_af) > 0L || length(not_in_table) > 0L) {
        stop(
          "Alleles in ", basename(pop), " and ", basename(af_path),
          " do not match (", length(not_in_af), " only in the allele table, ",
          length(not_in_table), " only in the MOIRE frequencies). ",
          "Check that the MOIRE output belongs to this allele table.",
          call. = FALSE
        )
      }

      af_file <- write_dcifer_input(
        allele_freq,
        file.path(input_prep_dir, paste0(base_filename, "_allele_freq.tsv"))
      )

      message("  Allele frequencies: ", basename(af_path))
    }

    message(
      "  Samples: ", length(samples),
      "; loci: ", n_distinct(allele_table$target_name),
      "; alleles: ", nrow(distinct(allele_table, target_name, seq))
    )

    allele_file <- write_dcifer_input(
      allele_table,
      file.path(input_prep_dir, paste0(base_filename, "_allele_table.tsv"))
    )


    # ========================================================
    # Run Dcifer
    # ========================================================

    results_filename <- paste0(base_filename, "_Dcifer_IBD_Results.txt")

    optional_file <- function(x) {
      if (is.na(x)) NULL else x
    }

    dcifer_results <- PGEcore::dcifer_ibd_wrapper(
      allele_table = allele_file,
      relatedness_output = file.path(output_dir, results_filename),
      coi_table = coi_file,
      allele_freq_table = optional_file(af_file),
      specimen_metadata = optional_file(metadata_file),
      pop_name_col = if (is.na(metadata_file)) NULL else pop_col,
      rnull = rnull,
      alpha = alpha,
      threads = n_threads
    )

    message(
      "  Pairs: ", nrow(dcifer_results),
      "; significant at alpha = ", alpha, ": ",
      sum(dcifer_results$p_value < alpha, na.rm = TRUE)
    )

    message("  Output file:\n    ", results_filename)

    message("  Completed: ", basename(pop))


    TRUE


  }, error = function(e) {

    stop_default_cluster()

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


message("\nAll Dcifer analyses completed successfully.")
