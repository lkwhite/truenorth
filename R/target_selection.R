# target_selection.R
# Hierarchical target selection with divergence analysis

library(dplyr)

#' Build hierarchical tree of tRNAs for an organism
#'
#' Returns a nested structure showing:
#' Amino acid -> Isoacceptor (anticodon) -> Gene family -> Individual copies
#'
#' @param trna_df Data frame from load_organism_data
#' @param compartment Filter by "nuclear", "mitochondrial", or NULL for all
#' @return Nested list representing the hierarchy
#' @export
build_trna_hierarchy <- function(trna_df, compartment = NULL) {
  if (!is.null(compartment)) {
    trna_df <- trna_df[trna_df$compartment == compartment | is.na(trna_df$compartment), ]
  }

  hierarchy <- list()

  for (aa in sort(unique(trna_df$amino_acid))) {
    aa_df <- trna_df[trna_df$amino_acid == aa, ]
    hierarchy[[aa]] <- list(
      count = nrow(aa_df),
      isoacceptors = list()
    )

    for (ac in sort(unique(aa_df$anticodon))) {
      ac_df <- aa_df[aa_df$anticodon == ac, ]
      isoacceptor_id <- paste(aa, ac, sep = "-")
      hierarchy[[aa]]$isoacceptors[[ac]] <- list(
        count = nrow(ac_df),
        isoacceptor_id = isoacceptor_id,
        families = list()
      )

      for (fam in sort(unique(ac_df$gene_family))) {
        fam_df <- ac_df[ac_df$gene_family == fam, ]
        family_id <- paste(aa, ac, fam, sep = "-")
        hierarchy[[aa]]$isoacceptors[[ac]]$families[[as.character(fam)]] <- list(
          count = nrow(fam_df),
          family_id = family_id,
          copies = fam_df$id
        )
      }
    }
  }

  class(hierarchy) <- c("trna_hierarchy", "list")
  attr(hierarchy, "organism") <- unique(trna_df$organism)
  attr(hierarchy, "compartment") <- compartment
  attr(hierarchy, "total_trnas") <- nrow(trna_df)

  hierarchy
}

#' Print summary of tRNA hierarchy
#'
#' @param x A trna_hierarchy object
#' @param ... Additional arguments (ignored)
#' @export
print.trna_hierarchy <- function(x, ...) {
  cat("tRNA Hierarchy for", attr(x, "organism"), "\n")
  if (!is.null(attr(x, "compartment"))) {
    cat("Compartment:", attr(x, "compartment"), "\n")
  }
  cat("Total tRNAs:", attr(x, "total_trnas"), "\n")
  cat("Amino acids:", length(x), "\n\n")


  for (aa in names(x)) {
    cat(aa, "(", x[[aa]]$count, "tRNAs):\n", sep = "")
    for (ac in names(x[[aa]]$isoacceptors)) {
      iso <- x[[aa]]$isoacceptors[[ac]]
      cat("  ", aa, "-", ac, " (", iso$count, "): ", sep = "")
      fam_counts <- sapply(iso$families, function(f) f$count)
      cat(paste(names(fam_counts), "=", fam_counts, collapse = ", "), "\n")
    }
  }
}

#' Select tRNAs by hierarchical criteria
#'
#' @param trna_df Data frame from load_organism_data
#' @param amino_acids Vector of amino acid codes to include, or NULL for all
#' @param anticodons Vector of anticodons to include, or NULL for all
#' @param family_ids Vector of family IDs (AA-Anticodon-Family), or NULL for all
#' @param trna_ids Vector of specific tRNA IDs, or NULL for all
#' @param compartment "nuclear", "mitochondrial", or NULL for all
#' @param exclude_amino_acids Amino acids to exclude
#' @param exclude_anticodons Anticodons to exclude
#' @param exclude_family_ids Family IDs to exclude
#' @param exclude_trna_ids Specific tRNA IDs to exclude
#' @return Filtered data frame
#' @export
select_trnas <- function(trna_df,
                         amino_acids = NULL,
                         anticodons = NULL,
                         family_ids = NULL,
                         trna_ids = NULL,
                         compartment = NULL,
                         exclude_amino_acids = NULL,
                         exclude_anticodons = NULL,
                         exclude_family_ids = NULL,
                         exclude_trna_ids = NULL) {

  result <- trna_df

  # Apply compartment filter

if (!is.null(compartment)) {
    result <- result[result$compartment == compartment | is.na(result$compartment), ]
  }

  # Apply inclusion filters (if specified, limit to these)
  if (!is.null(amino_acids)) {
    result <- result[result$amino_acid %in% amino_acids, ]
  }

  if (!is.null(anticodons)) {
    result <- result[result$anticodon %in% anticodons, ]
  }

  if (!is.null(family_ids)) {
    result <- result[result$family_id %in% family_ids, ]
  }

  if (!is.null(trna_ids)) {
    result <- result[result$id %in% trna_ids, ]
  }

  # Apply exclusion filters
  if (!is.null(exclude_amino_acids)) {
    result <- result[!result$amino_acid %in% exclude_amino_acids, ]
  }

  if (!is.null(exclude_anticodons)) {
    result <- result[!result$anticodon %in% exclude_anticodons, ]
  }

  if (!is.null(exclude_family_ids)) {
    result <- result[!result$family_id %in% exclude_family_ids, ]
  }

  if (!is.null(exclude_trna_ids)) {
    result <- result[!result$id %in% exclude_trna_ids, ]
  }

  result
}

#' Create a target selection specifying desired and avoid groups
#'
#' @param trna_df Data frame from load_organism_data
#' @param desired_ids Vector of tRNA IDs that probe SHOULD hit
#' @param avoid_ids Vector of tRNA IDs that probe should NOT hit
#' @return List with desired_df, avoid_df, and metadata
#' @export
create_target_selection <- function(trna_df, desired_ids, avoid_ids = NULL) {
  desired_df <- trna_df[trna_df$id %in% desired_ids, ]

  if (nrow(desired_df) == 0) {
    stop("No desired targets found")
  }

  if (!is.null(avoid_ids) && length(avoid_ids) > 0) {
    avoid_df <- trna_df[trna_df$id %in% avoid_ids, ]
  } else {
    # Default: avoid everything else
    avoid_df <- trna_df[!trna_df$id %in% desired_ids, ]
  }

  selection <- list(
    desired = desired_df,
    avoid = avoid_df,
    desired_ids = desired_ids,
    avoid_ids = avoid_df$id,
    n_desired = nrow(desired_df),
    n_avoid = nrow(avoid_df)
  )

  class(selection) <- c("target_selection", "list")
  selection
}

#' Analyze sequence conservation within a group of tRNAs
#'
#' Shows how similar the sequences are to each other - important for
#' understanding if one probe can hit all desired targets.
#'
#' @param trna_df Data frame of tRNAs to analyze (subset of full data)
#' @param similarity_data Pre-computed similarity data
#' @return List with conservation statistics
#' @export
analyze_group_conservation <- function(trna_df, similarity_data = NULL) {
  n <- nrow(trna_df)

  if (n == 0) {
    return(list(
      n_sequences = 0,
      message = "No sequences to analyze"
    ))
  }

  if (n == 1) {
    return(list(
      n_sequences = 1,
      identical = TRUE,
      min_identity = 100,
      max_identity = 100,
      mean_identity = 100,
      message = "Single sequence - no comparison needed"
    ))
  }

  ids <- trna_df$id

  # Get pairwise identities
  if (!is.null(similarity_data) && all(ids %in% rownames(similarity_data$identity_matrix))) {
    # Use pre-computed matrix
    submatrix <- similarity_data$identity_matrix[ids, ids]
    # Get upper triangle (excluding diagonal)
    pairwise_ids <- submatrix[upper.tri(submatrix)]
  } else {
    # Compute on the fly
    pairwise_ids <- c()
    for (i in 1:(n-1)) {
      for (j in (i+1):n) {
        aln <- align_pair(trna_df$sequence[i], trna_df$sequence[j])
        pairwise_ids <- c(pairwise_ids, aln$percent_identity)
      }
    }
  }

  # Check for identical sequences
  sequences <- trna_df$sequence
  n_unique <- length(unique(sequences))

  list(
    n_sequences = n,
    n_unique_sequences = n_unique,
    all_identical = n_unique == 1,
    min_identity = round(min(pairwise_ids), 1),
    max_identity = round(max(pairwise_ids), 1),
    mean_identity = round(mean(pairwise_ids), 1),
    sd_identity = round(sd(pairwise_ids), 2),
    pairwise_identities = pairwise_ids
  )
}

#' Analyze divergence between desired and avoid groups
#'
#' Shows how different the avoid group is from desired - important for
#' understanding if specificity is achievable.
#'
#' @param selection Target selection from create_target_selection
#' @param similarity_data Pre-computed similarity data
#' @return List with divergence statistics
#' @export
analyze_group_divergence <- function(selection, similarity_data = NULL) {
  desired_df <- selection$desired
  avoid_df <- selection$avoid

  if (nrow(avoid_df) == 0) {
    return(list(
      n_desired = nrow(desired_df),
      n_avoid = 0,
      message = "No avoid targets - any probe region will work"
    ))
  }

  desired_ids <- desired_df$id
  avoid_ids <- avoid_df$id

  # Get cross-group identities
  if (!is.null(similarity_data)) {
    mat <- similarity_data$identity_matrix
    cross_ids <- c()
    for (d_id in desired_ids) {
      if (d_id %in% rownames(mat)) {
        for (a_id in avoid_ids) {
          if (a_id %in% colnames(mat)) {
            cross_ids <- c(cross_ids, mat[d_id, a_id])
          }
        }
      }
    }
  } else {
    # Compute on the fly
    cross_ids <- c()
    for (i in 1:nrow(desired_df)) {
      for (j in 1:nrow(avoid_df)) {
        aln <- align_pair(desired_df$sequence[i], avoid_df$sequence[j])
        cross_ids <- c(cross_ids, aln$percent_identity)
      }
    }
  }

  # Find the most similar avoid targets (hardest to discriminate)
  if (!is.null(similarity_data)) {
    # For each desired, find most similar avoid
    closest_avoid <- lapply(desired_ids, function(d_id) {
      if (!d_id %in% rownames(similarity_data$identity_matrix)) return(NULL)
      avoid_similarities <- similarity_data$identity_matrix[d_id, avoid_ids]
      avoid_similarities <- avoid_similarities[!is.na(avoid_similarities)]
      if (length(avoid_similarities) == 0) return(NULL)
      max_idx <- which.max(avoid_similarities)
      data.frame(
        desired_id = d_id,
        closest_avoid_id = names(avoid_similarities)[max_idx],
        identity = avoid_similarities[max_idx],
        stringsAsFactors = FALSE
      )
    })
    closest_avoid_df <- bind_rows(closest_avoid)
  } else {
    closest_avoid_df <- NULL
  }

  list(
    n_desired = nrow(desired_df),
    n_avoid = nrow(avoid_df),
    min_cross_identity = round(min(cross_ids), 1),
    max_cross_identity = round(max(cross_ids), 1),
    mean_cross_identity = round(mean(cross_ids), 1),
    closest_avoid_targets = closest_avoid_df,
    # Interpretation
    specificity_challenge = ifelse(max(cross_ids) > 90,
                                   "HIGH - some avoid targets very similar",
                                   ifelse(max(cross_ids) > 80,
                                          "MODERATE - some avoid targets fairly similar",
                                          "LOW - avoid targets are divergent"))
  )
}

#' Analyze position-by-position conservation/divergence for probe design
#'
#' For each position in the reference sequence, shows:
#' - Conservation within desired group (higher = easier to hit all)
#' - Divergence from avoid group (higher = easier to avoid)
#'
#' @param selection Target selection from create_target_selection
#' @param reference_id ID of reference sequence (must be in desired group)
#' @return Data frame with per-position statistics
#' @export
analyze_region_divergence <- function(selection, reference_id = NULL) {
  desired_df <- selection$desired
  avoid_df <- selection$avoid

  # Use first desired as reference if not specified
  if (is.null(reference_id)) {
    reference_id <- desired_df$id[1]
  }

  ref_row <- desired_df[desired_df$id == reference_id, ]
  if (nrow(ref_row) == 0) {
    stop("Reference ID must be in desired group: ", reference_id)
  }

  ref_seq <- ref_row$sequence
  ref_len <- nchar(ref_seq)
  ref_bases <- strsplit(ref_seq, "")[[1]]

  # For each position, calculate:
  # 1. Conservation in desired group (% matching reference)
  # 2. Divergence from avoid group (% NOT matching reference)

  other_desired <- desired_df[desired_df$id != reference_id, ]
  results <- lapply(1:ref_len, function(pos) {
    ref_base <- ref_bases[pos]

    # Conservation in desired (including reference = 100% for that position)
    if (nrow(other_desired) == 0) {
      desired_conservation <- 100
    } else {
      matches <- sapply(other_desired$sequence, function(seq) {
        if (pos > nchar(seq)) return(FALSE)
        substr(seq, pos, pos) == ref_base
      })
      desired_conservation <- (sum(matches) + 1) / (length(matches) + 1) * 100
    }

    # Divergence from avoid
    if (nrow(avoid_df) == 0) {
      avoid_divergence <- 100  # No avoid targets = fully divergent
    } else {
      mismatches <- sapply(avoid_df$sequence, function(seq) {
        if (pos > nchar(seq)) return(TRUE)  # Different length = mismatch
        substr(seq, pos, pos) != ref_base
      })
      avoid_divergence <- sum(mismatches) / length(mismatches) * 100
    }

    data.frame(
      position = pos,
      reference_base = ref_base,
      desired_conservation = round(desired_conservation, 1),
      avoid_divergence = round(avoid_divergence, 1),
      # Combined score: want high conservation AND high divergence
      selectivity_score = round((desired_conservation + avoid_divergence) / 2, 1),
      stringsAsFactors = FALSE
    )
  })

  result_df <- bind_rows(results)
  attr(result_df, "reference_id") <- reference_id
  attr(result_df, "n_desired") <- nrow(desired_df)
  attr(result_df, "n_avoid") <- nrow(avoid_df)

  result_df
}

#' Find optimal probe regions for a target selection
#'
#' Identifies regions that are:
#' - Conserved across desired targets (probe will hit all)
#' - Divergent from avoid targets (probe won't cross-hybridize)
#'
#' @param selection Target selection from create_target_selection
#' @param min_length Minimum probe length (default 20)
#' @param max_length Maximum probe length (default 25)
#' @param min_conservation Minimum conservation in desired group (default 90%)
#' @param min_divergence Minimum divergence from avoid group (default 50%)
#' @return Data frame of recommended regions
#' @export
find_selective_regions <- function(selection,
                                   min_length = 20,
                                   max_length = 25,
                                   min_conservation = 90,
                                   min_divergence = 50) {

  # Get position-by-position analysis
  pos_analysis <- analyze_region_divergence(selection)
  ref_len <- nrow(pos_analysis)

  regions <- list()

  for (len in min_length:max_length) {
    for (start in 1:(ref_len - len + 1)) {
      end <- start + len - 1

      window <- pos_analysis[start:end, ]

      # Calculate window statistics
      mean_conservation <- mean(window$desired_conservation)
      min_conservation_in_window <- min(window$desired_conservation)
      mean_divergence <- mean(window$avoid_divergence)
      min_divergence_in_window <- min(window$avoid_divergence)
      mean_selectivity <- mean(window$selectivity_score)

      # Check if region meets criteria
      if (min_conservation_in_window >= min_conservation &&
          min_divergence_in_window >= min_divergence) {
        quality <- "EXCELLENT"
      } else if (mean_conservation >= min_conservation &&
                 mean_divergence >= min_divergence) {
        quality <- "GOOD"
      } else if (mean_conservation >= 80 && mean_divergence >= 30) {
        quality <- "FAIR"
      } else {
        quality <- "POOR"
      }

      regions[[length(regions) + 1]] <- data.frame(
        start = start,
        end = end,
        length = len,
        mean_conservation = round(mean_conservation, 1),
        min_conservation = round(min_conservation_in_window, 1),
        mean_divergence = round(mean_divergence, 1),
        min_divergence = round(min_divergence_in_window, 1),
        selectivity_score = round(mean_selectivity, 1),
        quality = quality,
        stringsAsFactors = FALSE
      )
    }
  }

  result <- bind_rows(regions)

  # Sort by selectivity score
  result <- result[order(-result$selectivity_score), ]
  result$rank <- 1:nrow(result)

  # Add reference info
  attr(result, "reference_id") <- attr(pos_analysis, "reference_id")
  attr(result, "n_desired") <- attr(pos_analysis, "n_desired")
  attr(result, "n_avoid") <- attr(pos_analysis, "n_avoid")

  result
}

#' Print summary of target selection
#'
#' @param x A target_selection object
#' @param ... Additional arguments (ignored)
#' @export
print.target_selection <- function(x, ...) {
  cat("Target Selection\n")
  cat("================\n\n")

  cat("DESIRED targets (", x$n_desired, "):\n", sep = "")
  if (x$n_desired <= 10) {
    for (id in x$desired_ids) {
      row <- x$desired[x$desired$id == id, ]
      cat("  ", id, "\n", sep = "")
    }
  } else {
    for (id in head(x$desired_ids, 5)) {
      cat("  ", id, "\n", sep = "")
    }
    cat("  ... and", x$n_desired - 5, "more\n")
  }

  cat("\nAVOID targets (", x$n_avoid, "):\n", sep = "")
  if (x$n_avoid <= 5) {
    for (id in x$avoid_ids) {
      cat("  ", id, "\n", sep = "")
    }
  } else {
    cat("  (", x$n_avoid, " tRNAs to avoid)\n", sep = "")
  }
}
