# sequence_utils.R
# FASTA parsing and sequence manipulation for TRUENORTH

library(Biostrings)
library(dplyr)

# Data paths relative to package root
DATA_DIR <- "data/fastas"

ORGANISM_FILES <- list(

human = "human/hg38-mito-and-nuclear-tRNAs.fa",
  yeast = "yeast/sacCer-mito-and-nuclear-tRNAs.fa",
  ecoli = "ecoli/ecoliK12MG1655-tRNAs.fa"
)

#' Parse a FASTA file
#'
#' @param file_path Path to FASTA file
#' @return DNAStringSet with sequences named by headers
#' @export
parse_fasta <- function(file_path) {
  if (!file.exists(file_path)) {
    stop("FASTA file not found: ", file_path)
  }

  sequences <- readDNAStringSet(file_path)

  # Convert any U to T (some files use RNA notation)
  sequences <- DNAStringSet(gsub("U", "T", as.character(sequences)))

  sequences
}

#' Parse tRNA header into components
#'
#' Handles different formats:
#' - Human/Yeast: nuc-tRNA-Ala-AGC-1-1, mito-tRNA-Phe-GAA-1-1
#' - E. coli: tRNA-Ala-GGC-1-1
#' - Yeast mito (simple): mito-tRNA-Pro-UGG
#'
#' @param header FASTA header string (without >)
#' @param organism One of "human", "yeast", "ecoli"
#' @return Named list with: compartment, amino_acid, anticodon, gene_family, copy_number
#' @export
parse_trna_header <- function(header, organism) {
  # Initialize result
  result <- list(
    compartment = NA_character_,
    amino_acid = NA_character_,
    anticodon = NA_character_,
    gene_family = NA_integer_,
    copy_number = NA_integer_
  )

  if (organism == "ecoli") {
    # E. coli format: tRNA-Ala-GGC-1-1
    # Also handles variants like tRNA-Ile2-CAT-1-1 (Ile2 = Isoleucine variant 2)
    pattern <- "^tRNA-([A-Za-z]+[0-9]*)-([A-Z]{3})-([0-9]+)-([0-9]+)$"
    match <- regmatches(header, regexec(pattern, header))[[1]]

    if (length(match) == 5) {
      result$compartment <- NA_character_  # E. coli has no compartments
      result$amino_acid <- match[2]
      result$anticodon <- match[3]
      result$gene_family <- as.integer(match[4])
      result$copy_number <- as.integer(match[5])
    }
  } else {
    # Human/Yeast format: (nuc|mito)-tRNA-Ala-AGC-1-1
    # Also handle yeast mito without gene family: mito-tRNA-Pro-UGG

    # Try full format first
    pattern_full <- "^(nuc|mito)-tRNA-([A-Za-z]+)-([A-Z]{3})-([0-9]+)-([0-9]+)$"
    match <- regmatches(header, regexec(pattern_full, header))[[1]]

    if (length(match) == 6) {
      result$compartment <- ifelse(match[2] == "nuc", "nuclear", "mitochondrial")
      result$amino_acid <- match[3]
      result$anticodon <- match[4]
      result$gene_family <- as.integer(match[5])
      result$copy_number <- as.integer(match[6])
    } else {
      # Try simple mito format: mito-tRNA-Pro-UGG
      pattern_simple <- "^(nuc|mito)-tRNA-([A-Za-z]+)-([A-Z]{3})$"
      match <- regmatches(header, regexec(pattern_simple, header))[[1]]

      if (length(match) == 4) {
        result$compartment <- ifelse(match[2] == "nuc", "nuclear", "mitochondrial")
        result$amino_acid <- match[3]
        result$anticodon <- match[4]
        result$gene_family <- 1L  # Default for simple format
        result$copy_number <- 1L
      }
    }
  }

  result
}

#' Load tRNA sequences and metadata for an organism
#'
#' @param organism One of "human", "yeast", "ecoli"
#' @param data_dir Base directory for data files (default: "data/fastas")
#' @return Data frame with tRNA metadata and sequences
#' @export
load_organism_data <- function(organism, data_dir = NULL) {
  if (!organism %in% names(ORGANISM_FILES)) {
    stop("Unknown organism: ", organism,
         ". Must be one of: ", paste(names(ORGANISM_FILES), collapse = ", "))
  }

  # Find data directory
  if (is.null(data_dir)) {
    # Try to find relative to working directory or package
    possible_paths <- c(
      DATA_DIR,
      file.path("..", DATA_DIR),
      file.path(system.file(package = "truenorth"), DATA_DIR)
    )
    data_dir <- Find(dir.exists, possible_paths)
    if (is.null(data_dir)) {
      stop("Could not find data directory. Tried: ", paste(possible_paths, collapse = ", "))
    }
  }

  file_path <- file.path(data_dir, ORGANISM_FILES[[organism]])
  sequences <- parse_fasta(file_path)

  # Parse headers and build data frame
  headers <- names(sequences)
  parsed <- lapply(headers, parse_trna_header, organism = organism)

  df <- data.frame(
    id = headers,
    organism = organism,
    compartment = sapply(parsed, `[[`, "compartment"),
    amino_acid = sapply(parsed, `[[`, "amino_acid"),
    anticodon = sapply(parsed, `[[`, "anticodon"),
    gene_family = sapply(parsed, `[[`, "gene_family"),
    copy_number = sapply(parsed, `[[`, "copy_number"),
    sequence = as.character(sequences),
    length = width(sequences),
    stringsAsFactors = FALSE
  )

  # Add composite keys
  df$family_id <- paste(df$amino_acid, df$anticodon, df$gene_family, sep = "-")
  df$isoacceptor_id <- paste(df$amino_acid, df$anticodon, sep = "-")

  df
}

#' Generate reverse complement of a DNA sequence
#'
#' The probe sequence is the reverse complement of the target region
#' (probe hybridizes to target in antiparallel orientation)
#'
#' @param sequence Character string or DNAString
#' @return Character string of reverse complement
#' @export
reverse_complement <- function(sequence) {
  if (is.character(sequence)) {
    sequence <- DNAString(sequence)
  }
  as.character(reverseComplement(sequence))
}

#' Filter tRNA data frame by various criteria
#'
#' @param trna_df Data frame from load_organism_data
#' @param compartment "nuclear", "mitochondrial", or NULL for all
#' @param amino_acid Three-letter amino acid code, or NULL for all
#' @param anticodon Three-letter anticodon, or NULL for all
#' @param family_id Family identifier (AA-Anticodon-GeneFamily), or NULL for all
#' @return Filtered data frame
#' @export
filter_trnas <- function(trna_df,
                         compartment = NULL,
                         amino_acid = NULL,
                         anticodon = NULL,
                         family_id = NULL) {
  result <- trna_df

  if (!is.null(compartment)) {
    result <- result[result$compartment == compartment | is.na(result$compartment), ]
  }

  if (!is.null(amino_acid)) {
    result <- result[result$amino_acid == amino_acid, ]
  }

  if (!is.null(anticodon)) {
    result <- result[result$anticodon == anticodon, ]
  }

  if (!is.null(family_id)) {
    result <- result[result$family_id == family_id, ]
  }

  result
}

#' Get unique amino acids in dataset
#'
#' @param trna_df Data frame from load_organism_data
#' @return Character vector of amino acid codes
#' @export
get_amino_acids <- function(trna_df) {
  sort(unique(trna_df$amino_acid))
}

#' Get unique anticodons for an amino acid
#'
#' @param trna_df Data frame from load_organism_data
#' @param amino_acid Three-letter amino acid code, or NULL for all
#' @return Character vector of anticodons
#' @export
get_anticodons <- function(trna_df, amino_acid = NULL) {
  if (!is.null(amino_acid)) {
    trna_df <- trna_df[trna_df$amino_acid == amino_acid, ]
  }
  sort(unique(trna_df$anticodon))
}

#' Get gene families for an isoacceptor
#'
#' @param trna_df Data frame from load_organism_data
#' @param amino_acid Three-letter amino acid code
#' @param anticodon Three-letter anticodon
#' @return Integer vector of gene family numbers
#' @export
get_gene_families <- function(trna_df, amino_acid, anticodon) {
  filtered <- trna_df[trna_df$amino_acid == amino_acid &
                       trna_df$anticodon == anticodon, ]
  sort(unique(filtered$gene_family))
}

#' Extract a subsequence from a tRNA
#'
#' @param sequence Character string of full tRNA sequence
#' @param start Start position (1-based)
#' @param end End position (1-based, inclusive)
#' @return Character string of subsequence
#' @export
extract_region <- function(sequence, start, end) {
  if (start < 1 || end > nchar(sequence) || start > end) {
    stop("Invalid region: start=", start, ", end=", end,
         ", sequence length=", nchar(sequence))
  }
  substr(sequence, start, end)
}
