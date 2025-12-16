# validation.R
# Specificity checking and off-target analysis

library(Biostrings)
library(dplyr)

#' Count mismatches between probe and a target sequence at a position
#'
#' @param probe_seq Probe sequence (5' to 3')
#' @param target_seq Target tRNA sequence
#' @param position Start position in target where probe would bind (1-based)
#' @return Integer count of mismatches, or NA if position invalid
#' @export
count_mismatches_at_position <- function(probe_seq, target_seq, position) {
  probe_len <- nchar(probe_seq)
  target_len <- nchar(target_seq)

  # Check bounds
  if (position < 1 || position + probe_len - 1 > target_len) {
    return(NA_integer_)
  }

  # Get the target region
  target_region <- substr(target_seq, position, position + probe_len - 1)

  # Probe binds as reverse complement
  probe_rc <- reverse_complement(probe_seq)

  # Count mismatches
  probe_bases <- strsplit(toupper(probe_rc), "")[[1]]
  target_bases <- strsplit(toupper(target_region), "")[[1]]

  sum(probe_bases != target_bases)
}

#' Find best binding position for a probe on a target sequence
#'
#' Scans all positions and finds the one with fewest mismatches.
#'
#' @param probe_seq Probe sequence
#' @param target_seq Target tRNA sequence
#' @return List with: position, mismatches, binds (TRUE if mismatches <= threshold)
#' @export
find_best_binding <- function(probe_seq, target_seq, max_mismatches = 5) {
  probe_len <- nchar(probe_seq)
  target_len <- nchar(target_seq)

  if (probe_len > target_len) {
    return(list(position = NA, mismatches = NA, binds = FALSE))
  }

  # Check all possible positions
  best_pos <- NA
  best_mismatches <- Inf

  for (pos in 1:(target_len - probe_len + 1)) {
    mm <- count_mismatches_at_position(probe_seq, target_seq, pos)
    if (!is.na(mm) && mm < best_mismatches) {
      best_mismatches <- mm
      best_pos <- pos
    }
  }

  list(
    position = best_pos,
    mismatches = best_mismatches,
    binds = best_mismatches <= max_mismatches
  )
}

#' Check probe specificity against all tRNAs
#'
#' @param probe_seq Probe sequence (what you order, 5' to 3')
#' @param target_id The intended target tRNA ID
#' @param trna_df Data frame from load_organism_data
#' @param similarity_data Pre-computed similarity (optional, speeds up search)
#' @param max_mismatches Maximum mismatches to report (default 5)
#' @return Data frame of all hits with mismatch counts
#' @export
check_probe_specificity <- function(probe_seq, target_id, trna_df,
                                    similarity_data = NULL,
                                    max_mismatches = 5) {
  # Get target info for categorization
  target_row <- trna_df[trna_df$id == target_id, ]
  if (nrow(target_row) == 0) {
    stop("Target ID not found: ", target_id)
  }

  target_family <- target_row$family_id
  target_isoacceptor <- target_row$isoacceptor_id
  target_aa <- target_row$amino_acid

  # Determine which tRNAs to check
  # If similarity data available, prioritize similar ones
  if (!is.null(similarity_data)) {
    # Get similar tRNAs (potential off-targets)
    similar <- find_similar_trnas(target_id, similarity_data$identity_matrix, threshold = 50)
    check_ids <- c(target_id, similar$id)
    # Also include any we might have missed
    check_ids <- unique(c(check_ids, trna_df$id))
  } else {
    check_ids <- trna_df$id
  }

  # Check each tRNA
  results <- lapply(check_ids, function(id) {
    row <- trna_df[trna_df$id == id, ]
    if (nrow(row) == 0) return(NULL)

    binding <- find_best_binding(probe_seq, row$sequence, max_mismatches)

    if (is.na(binding$mismatches)) return(NULL)

    # Categorize relationship to target
    if (id == target_id) {
      category <- "target"
    } else if (row$family_id == target_family) {
      category <- "same_family"
    } else if (row$isoacceptor_id == target_isoacceptor) {
      category <- "same_isoacceptor"
    } else if (row$amino_acid == target_aa) {
      category <- "same_amino_acid"
    } else {
      category <- "other"
    }

    data.frame(
      id = id,
      amino_acid = row$amino_acid,
      anticodon = row$anticodon,
      gene_family = row$gene_family,
      copy_number = row$copy_number,
      compartment = row$compartment,
      best_position = binding$position,
      mismatches = binding$mismatches,
      category = category,
      stringsAsFactors = FALSE
    )
  })

  result_df <- bind_rows(results)

  # Filter to those within mismatch threshold and sort
  result_df <- result_df[result_df$mismatches <= max_mismatches, ]
  result_df <- result_df[order(result_df$mismatches, result_df$category), ]

  # Add rank
  result_df$rank <- 1:nrow(result_df)

  result_df
}

#' Categorize off-target hits by relationship to target
#'
#' @param hits_df Data frame from check_probe_specificity
#' @return Summary data frame with counts by category
#' @export
categorize_off_targets <- function(hits_df) {
  if (nrow(hits_df) == 0) {
    return(data.frame(
      category = character(),
      count = integer(),
      min_mismatches = integer(),
      max_mismatches = integer(),
      stringsAsFactors = FALSE
    ))
  }

  hits_df %>%
    group_by(category) %>%
    summarize(
      count = n(),
      min_mismatches = min(mismatches),
      max_mismatches = max(mismatches),
      .groups = "drop"
    ) %>%
    arrange(min_mismatches)
}

#' Generate visual alignment report for a probe vs off-targets
#'
#' Shows the probe sequence aligned to each hit with mismatches highlighted.
#'
#' @param probe_seq Probe sequence
#' @param hits_df Data frame from check_probe_specificity
#' @param trna_df Data frame from load_organism_data
#' @param max_display Maximum number of alignments to show (default 10)
#' @return Character vector of formatted alignment strings
#' @export
generate_alignment_report <- function(probe_seq, hits_df, trna_df, max_display = 10) {
  if (nrow(hits_df) == 0) {
    return("No hits to display.")
  }

  # Limit display
  display_df <- head(hits_df, max_display)

  probe_rc <- reverse_complement(probe_seq)
  probe_len <- nchar(probe_seq)

  lines <- c(
    paste0("Probe: 5'-", probe_seq, "-3'"),
    paste0("       (binds as: 3'-", probe_rc, "-5')"),
    "",
    paste0("Off-target analysis (", nrow(hits_df), " hits, showing top ", nrow(display_df), "):"),
    paste(rep("-", 60), collapse = ""),
    ""
  )

  for (i in 1:nrow(display_df)) {
    hit <- display_df[i, ]
    target_row <- trna_df[trna_df$id == hit$id, ]
    target_seq <- target_row$sequence

    # Get binding region
    pos <- hit$best_position
    if (is.na(pos)) next

    target_region <- substr(target_seq, pos, pos + probe_len - 1)

    # Build match string
    probe_bases <- strsplit(probe_rc, "")[[1]]
    target_bases <- strsplit(target_region, "")[[1]]

    match_str <- sapply(1:length(probe_bases), function(j) {
      if (j > length(target_bases)) return(" ")
      if (probe_bases[j] == target_bases[j]) "|" else "X"
    })
    match_str <- paste(match_str, collapse = "")

    # Format output
    lines <- c(lines,
               paste0(hit$id, " (", hit$category, ", ", hit$mismatches, " mismatches):"),
               paste0("  Target: 5'-", target_region, "-3'  [pos ", pos, "-", pos + probe_len - 1, "]"),
               paste0("          ", match_str),
               paste0("  Probe:  3'-", probe_rc, "-5'"),
               ""
    )
  }

  if (nrow(hits_df) > max_display) {
    lines <- c(lines,
               paste0("... and ", nrow(hits_df) - max_display, " more hits not shown."))
  }

  paste(lines, collapse = "\n")
}

#' Calculate overall specificity score for a probe
#'
#' Higher score = more specific (fewer/weaker off-target hits).
#'
#' @param hits_df Data frame from check_probe_specificity
#' @param target_id The intended target tRNA ID
#' @return Numeric score (0-1, higher is better)
#' @export
calculate_specificity_score <- function(hits_df, target_id) {
  # Filter out the target itself
  off_targets <- hits_df[hits_df$id != target_id, ]

  if (nrow(off_targets) == 0) {
    return(1.0)  # Perfect specificity
  }

  # Penalize based on off-target characteristics
  # Exact matches (0 mismatches) are worst
  # More mismatches = less penalty

  penalty <- sum(sapply(1:nrow(off_targets), function(i) {
    mm <- off_targets$mismatches[i]
    category <- off_targets$category[i]

    # Base penalty by mismatch count
    base_penalty <- 1 / (mm + 1)  # 0 mm = 1.0, 1 mm = 0.5, 2 mm = 0.33, etc.

    # Modify by category (same family less concerning)
    if (category == "same_family") {
      base_penalty <- base_penalty * 0.5  # Family members often intentional
    } else if (category == "same_isoacceptor") {
      base_penalty <- base_penalty * 0.8
    }

    base_penalty
  }))

  # Convert penalty to score (0-1)
  score <- 1 / (1 + penalty)
  score
}

#' Full validation of a probe: specificity check with summary
#'
#' @param probe_seq Probe sequence
#' @param target_id The intended target tRNA ID
#' @param trna_df Data frame from load_organism_data
#' @param similarity_data Pre-computed similarity (optional)
#' @param max_mismatches Maximum mismatches to consider (default 5)
#' @return List with: hits_df, summary, specificity_score, alignment_report
#' @export
validate_probe <- function(probe_seq, target_id, trna_df,
                           similarity_data = NULL, max_mismatches = 5) {
  # Check specificity
  hits <- check_probe_specificity(
    probe_seq, target_id, trna_df,
    similarity_data = similarity_data,
    max_mismatches = max_mismatches
  )

  # Categorize
  summary <- categorize_off_targets(hits)

  # Score
  score <- calculate_specificity_score(hits, target_id)

  # Generate report
  report <- generate_alignment_report(probe_seq, hits, trna_df)

  list(
    hits = hits,
    summary = summary,
    specificity_score = score,
    alignment_report = report,
    probe_sequence = probe_seq,
    target_id = target_id
  )
}
