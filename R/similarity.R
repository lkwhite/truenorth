# similarity.R
# Pre-computed similarity analysis for specificity-aware probe design

library(Biostrings)
library(dplyr)

# Check for pwalign package (required for pairwiseAlignment in BioC >= 3.19)
if (!requireNamespace("pwalign", quietly = TRUE)) {
  message("Note: pwalign package not found. Install with: BiocManager::install('pwalign')")
  message("Similarity computations will use a simpler alignment method.")
  USE_PWALIGN <- FALSE
} else {
  USE_PWALIGN <- TRUE
}

# Cache directory for pre-computed similarity data
SIMILARITY_DIR <- "data/similarity"

#' Compute pairwise alignment between two sequences
#'
#' @param seq1 Character string or DNAString
#' @param seq2 Character string or DNAString
#' @return List with: score, percent_identity, alignment object
#' @keywords internal
align_pair <- function(seq1, seq2) {
  if (is.character(seq1)) seq1 <- DNAString(seq1)
  if (is.character(seq2)) seq2 <- DNAString(seq2)

  # Try pairwiseAlignment if pwalign is available
  if (USE_PWALIGN) {
    tryCatch({
      aln <- pairwiseAlignment(
        pattern = seq1,
        subject = seq2,
        type = "global"
      )

      n_match <- nmatch(aln)
      aln_length <- nchar(aln)
      pct_identity <- (n_match / aln_length) * 100

      return(list(
        score = score(aln),
        percent_identity = pct_identity,
        n_matches = n_match,
        n_mismatches = nmismatch(aln),
        alignment_length = aln_length,
        alignment = aln
      ))
    }, error = function(e) {
      # Fall through to simple method
    })
  }

  # Simple fallback: compare sequences position by position
  # This works well for tRNAs which are similar length
  seq1_str <- as.character(seq1)
  seq2_str <- as.character(seq2)
  len1 <- nchar(seq1_str)
  len2 <- nchar(seq2_str)

  # Use shorter length for comparison
  compare_len <- min(len1, len2)
  bases1 <- strsplit(substr(seq1_str, 1, compare_len), "")[[1]]
  bases2 <- strsplit(substr(seq2_str, 1, compare_len), "")[[1]]

  n_match <- sum(bases1 == bases2)
  n_mismatch <- compare_len - n_match
  # Add penalty for length difference
  length_diff <- abs(len1 - len2)

  aln_length <- max(len1, len2)
  pct_identity <- (n_match / aln_length) * 100

  list(
    score = n_match - n_mismatch - length_diff,
    percent_identity = pct_identity,
    n_matches = n_match,
    n_mismatches = n_mismatch + length_diff,
    alignment_length = aln_length,
    alignment = NULL
  )
}

#' Compute pairwise identity matrix for all tRNAs
#'
#' This is O(n^2) so may take a minute for large datasets (human ~454 tRNAs)
#'
#' @param trna_df Data frame from load_organism_data
#' @param verbose Print progress messages
#' @return Named matrix of percent identities (symmetric)
#' @export
compute_identity_matrix <- function(trna_df, verbose = TRUE) {
  n <- nrow(trna_df)
  ids <- trna_df$id

  # Initialize matrix
  identity_matrix <- matrix(100, nrow = n, ncol = n,
                            dimnames = list(ids, ids))

  # Compute upper triangle (matrix is symmetric)
  total_pairs <- n * (n - 1) / 2
  pair_count <- 0

  if (verbose) {
    message("Computing pairwise alignments for ", n, " tRNAs (", total_pairs, " pairs)...")
  }

  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      aln <- align_pair(trna_df$sequence[i], trna_df$sequence[j])
      identity_matrix[i, j] <- aln$percent_identity
      identity_matrix[j, i] <- aln$percent_identity

      pair_count <- pair_count + 1
      if (verbose && pair_count %% 10000 == 0) {
        message("  ", pair_count, "/", total_pairs, " pairs computed...")
      }
    }
  }

  if (verbose) message("Done.")
  identity_matrix
}

#' Find tRNAs similar to a target above a threshold
#'
#' @param target_id ID of target tRNA
#' @param identity_matrix Matrix from compute_identity_matrix
#' @param threshold Minimum percent identity to include (default 70)
#' @return Data frame with id, percent_identity, sorted by similarity
#' @export
find_similar_trnas <- function(target_id, identity_matrix, threshold = 70) {
  if (!target_id %in% rownames(identity_matrix)) {
    stop("Target ID not found in identity matrix: ", target_id)
  }

  similarities <- identity_matrix[target_id, ]

  # Exclude self
  similarities <- similarities[names(similarities) != target_id]

  # Filter by threshold and sort
  above_threshold <- similarities[similarities >= threshold]
  above_threshold <- sort(above_threshold, decreasing = TRUE)

  data.frame(
    id = names(above_threshold),
    percent_identity = as.numeric(above_threshold),
    stringsAsFactors = FALSE
  )
}

#' Find positions where target diverges from similar sequences
#'
#' For each position in the target, count how many of the similar sequences
#' have a different base at that position.
#'
#' @param target_id ID of target tRNA
#' @param similar_ids Vector of IDs of similar tRNAs to compare against
#' @param trna_df Data frame from load_organism_data
#' @return Data frame with position, divergence_count, divergence_fraction
#' @export
find_divergent_positions <- function(target_id, similar_ids, trna_df) {
  target_row <- trna_df[trna_df$id == target_id, ]
  if (nrow(target_row) == 0) {
    stop("Target ID not found: ", target_id)
  }

  target_seq <- target_row$sequence
  target_len <- nchar(target_seq)
  target_bases <- strsplit(target_seq, "")[[1]]

  # Get similar sequences
  similar_df <- trna_df[trna_df$id %in% similar_ids, ]

  if (nrow(similar_df) == 0) {
    # No similar sequences - all positions are equally "divergent"
    return(data.frame(
      position = 1:target_len,
      divergence_count = 0,
      divergence_fraction = 0,
      target_base = target_bases,
      stringsAsFactors = FALSE
    ))
  }

  # For each position, count mismatches across similar sequences
  divergence_counts <- sapply(1:target_len, function(pos) {
    target_base <- target_bases[pos]

    # Count how many similar sequences have a different base at this position
    # Handle different length sequences by checking bounds
    mismatches <- sum(sapply(similar_df$sequence, function(seq) {
      if (pos > nchar(seq)) return(1)  # Position doesn't exist = mismatch
      substr(seq, pos, pos) != target_base
    }))

    mismatches
  })

  data.frame(
    position = 1:target_len,
    divergence_count = divergence_counts,
    divergence_fraction = divergence_counts / nrow(similar_df),
    target_base = target_bases,
    stringsAsFactors = FALSE
  )
}

#' Recommend probe regions based on divergence from similar sequences
#'
#' Finds contiguous regions with high divergence that could make good probe targets.
#'
#' @param target_id ID of target tRNA
#' @param trna_df Data frame from load_organism_data
#' @param similarity_data Pre-computed similarity data from compute_similarity_data
#' @param min_length Minimum region length (default 20)
#' @param max_length Maximum region length (default 25)
#' @param similarity_threshold Include tRNAs with >= this % identity as potential off-targets
#' @return Data frame of recommended regions with start, end, avg_divergence, score
#' @export
recommend_probe_regions <- function(target_id, trna_df, similarity_data,
                                    min_length = 20, max_length = 25,
                                    similarity_threshold = 70) {
  # Get similar tRNAs
  similar_df <- find_similar_trnas(target_id, similarity_data$identity_matrix,
                                   threshold = similarity_threshold)
  similar_ids <- similar_df$id

  # Get divergent positions
  div_pos <- find_divergent_positions(target_id, similar_ids, trna_df)

  target_len <- nrow(div_pos)

  # Score each possible window
  regions <- list()

  for (len in min_length:max_length) {
    for (start in 1:(target_len - len + 1)) {
      end <- start + len - 1
      window_div <- div_pos$divergence_fraction[start:end]

      # Score: higher is better (more divergent from off-targets)
      # We want regions where most positions differ from similar tRNAs
      avg_divergence <- mean(window_div)
      min_divergence <- min(window_div)
      max_divergence <- max(window_div)

      # Composite score: favor regions that are consistently divergent
      score <- avg_divergence * 0.7 + min_divergence * 0.3

      regions[[length(regions) + 1]] <- data.frame(
        start = start,
        end = end,
        length = len,
        avg_divergence = avg_divergence,
        min_divergence = min_divergence,
        max_divergence = max_divergence,
        score = score,
        stringsAsFactors = FALSE
      )
    }
  }

  result <- bind_rows(regions)

  # Sort by score (highest first)
  result <- result[order(-result$score), ]

  # Add rank
  result$rank <- 1:nrow(result)

  result
}

#' Compute full similarity data for an organism
#'
#' @param trna_df Data frame from load_organism_data
#' @param verbose Print progress messages
#' @return List with identity_matrix and metadata
#' @export
compute_similarity_data <- function(trna_df, verbose = TRUE) {
  organism <- unique(trna_df$organism)
  if (length(organism) != 1) {
    stop("trna_df should contain data for exactly one organism")
  }

  if (verbose) message("Computing similarity data for ", organism, "...")

  identity_matrix <- compute_identity_matrix(trna_df, verbose = verbose)

  list(
    organism = organism,
    n_trnas = nrow(trna_df),
    identity_matrix = identity_matrix,
    computed_at = Sys.time()
  )
}

#' Save similarity data to cache
#'
#' @param similarity_data Output from compute_similarity_data
#' @param data_dir Base data directory (default: "data")
#' @export
save_similarity_data <- function(similarity_data, data_dir = "data") {
  cache_dir <- file.path(data_dir, "similarity")
  if (!dir.exists(cache_dir)) {
    dir.create(cache_dir, recursive = TRUE)
  }

  organism <- similarity_data$organism
  file_path <- file.path(cache_dir, paste0(organism, "_similarity.rds"))

  saveRDS(similarity_data, file_path)
  message("Saved similarity data to: ", file_path)

  invisible(file_path)
}

#' Load similarity data from cache
#'
#' @param organism One of "human", "yeast", "ecoli"
#' @param data_dir Base data directory (default: "data")
#' @return Similarity data list, or NULL if not found
#' @export
load_similarity_data <- function(organism, data_dir = "data") {
  cache_dir <- file.path(data_dir, "similarity")
  file_path <- file.path(cache_dir, paste0(organism, "_similarity.rds"))

  if (!file.exists(file_path)) {
    return(NULL)
  }

  readRDS(file_path)
}

#' Get or compute similarity data for an organism
#'
#' Loads from cache if available, otherwise computes and caches.
#'
#' @param organism One of "human", "yeast", "ecoli"
#' @param trna_df Data frame from load_organism_data (only needed if computing)
#' @param data_dir Base data directory
#' @param force_recompute If TRUE, recompute even if cache exists
#' @param verbose Print progress messages
#' @return Similarity data list
#' @export
get_similarity_data <- function(organism, trna_df = NULL, data_dir = "data",
                                force_recompute = FALSE, verbose = TRUE) {
  # Try loading from cache
  if (!force_recompute) {
    cached <- load_similarity_data(organism, data_dir)
    if (!is.null(cached)) {
      if (verbose) message("Loaded cached similarity data for ", organism)
      return(cached)
    }
  }

  # Need to compute
  if (is.null(trna_df)) {
    stop("No cached data found and trna_df not provided. ",
         "Either provide trna_df or run compute_similarity_data first.")
  }

  similarity_data <- compute_similarity_data(trna_df, verbose = verbose)
  save_similarity_data(similarity_data, data_dir)

  similarity_data
}
