# probe_design.R
# Probe generation and property calculations

library(Biostrings)

# Nearest-neighbor thermodynamic parameters (SantaLucia 1998)
# dH in kcal/mol, dS in cal/(mol·K)
NN_PARAMS <- list(
  # DNA/DNA hybridization parameters
  AA = list(dH = -7.9, dS = -22.2),
  TT = list(dH = -7.9, dS = -22.2),  # Same as AA (complement)
  AT = list(dH = -7.2, dS = -20.4),
  TA = list(dH = -7.2, dS = -21.3),
  CA = list(dH = -8.5, dS = -22.7),
  TG = list(dH = -8.5, dS = -22.7),  # Same as CA (complement)
  GT = list(dH = -8.4, dS = -22.4),
  AC = list(dH = -8.4, dS = -22.4),  # Same as GT (complement)
  CT = list(dH = -7.8, dS = -21.0),
  AG = list(dH = -7.8, dS = -21.0),  # Same as CT (complement)
  GA = list(dH = -8.2, dS = -22.2),
  TC = list(dH = -8.2, dS = -22.2),  # Same as GA (complement)
  CG = list(dH = -10.6, dS = -27.2),
  GC = list(dH = -9.8, dS = -24.4),
  GG = list(dH = -8.0, dS = -19.9),
  CC = list(dH = -8.0, dS = -19.9),  # Same as GG (complement)

  # Initiation parameters
  init_gc = list(dH = 0.1, dS = -2.8),   # G/C terminal
  init_at = list(dH = 2.3, dS = 4.1),    # A/T terminal

  # Symmetry correction (for self-complementary)
  symmetry = list(dH = 0, dS = -1.4)
)

# Gas constant in cal/(mol·K)
R_GAS <- 1.987

#' Known modification-prone regions in tRNAs
#'
#' Returns a list of regions with modification risk levels.
#' Positions are approximate and based on canonical tRNA numbering.
#'
#' @param seq_length Length of tRNA sequence (for scaling)
#' @return Data frame with start, end, risk (1-3), and description
#' @export
get_modification_zones <- function(seq_length = 76) {
  # Scale positions based on typical 76nt tRNA
  scale <- seq_length / 76

  zones <- data.frame(
    start = round(c(
      16,   # D-loop start
      34,   # Anticodon
      54    # TΨC loop
    ) * scale),
    end = round(c(
      20,   # D-loop end
      37,   # Anticodon + position 37
      56    # TΨC loop end
    ) * scale),
    risk = c(2, 3, 2),  # 1=low, 2=moderate, 3=high
    description = c(
      "D-loop (dihydrouridine)",
      "Anticodon loop (heavily modified)",
      "TΨC loop (pseudouridine, methylation)"
    ),
    stringsAsFactors = FALSE
  )

  zones
}

#' Calculate modification penalty for a probe region
#'
#' Higher values = more modification risk = worse for hybridization
#'
#' @param start Probe start position
#' @param end Probe end position
#' @param seq_length tRNA sequence length
#' @return Numeric penalty score (0 = no modifications, higher = worse)
#' @export
calculate_modification_penalty <- function(start, end, seq_length = 76) {
  zones <- get_modification_zones(seq_length)

  penalty <- 0
  for (i in 1:nrow(zones)) {
    zone <- zones[i, ]
    # Check for overlap
    overlap_start <- max(start, zone$start)
    overlap_end <- min(end, zone$end)

    if (overlap_start <= overlap_end) {
      # Has overlap - penalty proportional to overlap and risk
      overlap_length <- overlap_end - overlap_start + 1
      probe_length <- end - start + 1
      overlap_fraction <- overlap_length / probe_length

      # Penalty scales with risk level and overlap fraction
      # Risk 3 (anticodon): up to 30 points
      # Risk 2 (D-loop, TΨC): up to 15 points
      penalty <- penalty + (zone$risk * 10 * overlap_fraction)
    }
  }

  penalty
}

#' Calculate GC content of a sequence
#'
#' @param sequence DNA sequence string
#' @return Numeric GC percentage (0-100)
#' @export
calculate_gc_content <- function(sequence) {
  sequence <- toupper(sequence)
  bases <- strsplit(sequence, "")[[1]]

  gc_count <- sum(bases %in% c("G", "C"))
  total <- length(bases)

  if (total == 0) return(0)
  (gc_count / total) * 100
}

#' Calculate melting temperature using basic formula
#'
#' Uses the Wallace rule for short oligos and a modified formula for longer ones.
#'
#' @param sequence DNA sequence string
#' @return Numeric Tm in Celsius
#' @export
calculate_tm_basic <- function(sequence) {
  sequence <- toupper(sequence)
  len <- nchar(sequence)

  if (len == 0) return(NA_real_)

  # Count bases
  bases <- strsplit(sequence, "")[[1]]
  nA <- sum(bases == "A")
  nT <- sum(bases == "T")
  nG <- sum(bases == "G")
  nC <- sum(bases == "C")

  if (len < 14) {
    # Wallace rule for short oligos
    tm <- 2 * (nA + nT) + 4 * (nG + nC)
  } else {
    # Modified formula for longer oligos
    gc_content <- (nG + nC) / len
    tm <- 64.9 + 41 * (gc_content - 0.164)
  }

  tm
}

#' Calculate melting temperature using nearest-neighbor method
#'
#' More accurate for 15-30 nt probes. Uses SantaLucia (1998) unified parameters.
#'
#' @param sequence DNA sequence string
#' @param probe_conc Probe concentration in M (default 250 nM)
#' @param na_conc Sodium concentration in M (default 50 mM)
#' @return Numeric Tm in Celsius
#' @export
calculate_tm_nn <- function(sequence, probe_conc = 250e-9, na_conc = 50e-3) {
  sequence <- toupper(sequence)
  len <- nchar(sequence)

  if (len < 2) return(NA_real_)

  # Sum thermodynamic parameters for each dinucleotide
  dH_total <- 0
  dS_total <- 0

  for (i in 1:(len - 1)) {
    dinuc <- substr(sequence, i, i + 1)

    if (dinuc %in% names(NN_PARAMS)) {
      params <- NN_PARAMS[[dinuc]]
      dH_total <- dH_total + params$dH
      dS_total <- dS_total + params$dS
    } else {
      # Unknown dinucleotide (e.g., contains N) - use average
      dH_total <- dH_total + (-8.0)
      dS_total <- dS_total + (-22.0)
    }
  }

  # Add initiation parameters based on terminal bases
  first_base <- substr(sequence, 1, 1)
  last_base <- substr(sequence, len, len)

  if (first_base %in% c("G", "C")) {
    dH_total <- dH_total + NN_PARAMS$init_gc$dH
    dS_total <- dS_total + NN_PARAMS$init_gc$dS
  } else {
    dH_total <- dH_total + NN_PARAMS$init_at$dH
    dS_total <- dS_total + NN_PARAMS$init_at$dS
  }

  if (last_base %in% c("G", "C")) {
    dH_total <- dH_total + NN_PARAMS$init_gc$dH
    dS_total <- dS_total + NN_PARAMS$init_gc$dS
  } else {
    dH_total <- dH_total + NN_PARAMS$init_at$dH
    dS_total <- dS_total + NN_PARAMS$init_at$dS
  }

  # Calculate Tm
  # Tm = dH / (dS + R * ln(Ct/4)) - 273.15
  # For non-self-complementary: Ct = probe concentration
  dH_cal <- dH_total * 1000  # Convert kcal to cal
  tm <- dH_cal / (dS_total + R_GAS * log(probe_conc / 4)) - 273.15

  # Salt correction (Owczarzy et al. 2004 approximation)
  # Tm_corrected = Tm + 16.6 * log10([Na+])
  tm_corrected <- tm + 16.6 * log10(na_conc)

  tm_corrected
}

#' Generate candidate probe sequences from a target
#'
#' Creates all possible probes of specified lengths from the target sequence.
#' The returned probe_sequence is the reverse complement (what you order).
#'
#' @param sequence Target tRNA sequence
#' @param min_length Minimum probe length (default 20)
#' @param max_length Maximum probe length (default 25)
#' @param region Region constraint: "5prime", "3prime", "middle", or c(start, end)
#' @return Data frame with probe information
#' @export
generate_candidate_probes <- function(sequence, min_length = 20, max_length = 25,
                                      region = NULL) {
  seq_len <- nchar(sequence)

  # Determine region bounds
  if (is.null(region)) {
    start_bound <- 1
    end_bound <- seq_len
  } else if (is.character(region)) {
    if (region == "5prime") {
      start_bound <- 1
      end_bound <- min(35, seq_len)  # First ~35 nt
    } else if (region == "3prime") {
      start_bound <- max(1, seq_len - 34)  # Last ~35 nt
      end_bound <- seq_len
    } else if (region == "middle") {
      start_bound <- max(1, floor(seq_len * 0.25))
      end_bound <- min(seq_len, ceiling(seq_len * 0.75))
    } else {
      stop("Unknown region: ", region, ". Use '5prime', '3prime', 'middle', or c(start, end)")
    }
  } else if (is.numeric(region) && length(region) == 2) {
    start_bound <- max(1, region[1])
    end_bound <- min(seq_len, region[2])
  } else {
    stop("region must be NULL, a string, or a numeric vector of length 2")
  }

  probes <- list()

  for (len in min_length:max_length) {
    # Probe must fit within bounds
    for (start in start_bound:(end_bound - len + 1)) {
      if (start < 1) next

      end <- start + len - 1
      if (end > seq_len) next

      target_region <- substr(sequence, start, end)
      probe_seq <- reverse_complement(target_region)

      gc <- calculate_gc_content(probe_seq)
      tm_basic <- calculate_tm_basic(probe_seq)
      tm_nn <- calculate_tm_nn(probe_seq)

      probes[[length(probes) + 1]] <- data.frame(
        start = start,
        end = end,
        length = len,
        target_region = target_region,
        probe_sequence = probe_seq,
        gc_content = gc,
        tm_basic = tm_basic,
        tm_nn = tm_nn,
        stringsAsFactors = FALSE
      )
    }
  }

  if (length(probes) == 0) {
    # Return empty data frame with correct columns
    return(data.frame(
      start = integer(),
      end = integer(),
      length = integer(),
      target_region = character(),
      probe_sequence = character(),
      gc_content = numeric(),
      tm_basic = numeric(),
      tm_nn = numeric(),
      stringsAsFactors = FALSE
    ))
  }

  bind_rows(probes)
}

#' Score a probe based on various criteria
#'
#' @param probe_row Row from probe data frame (or vector with named elements)
#' @param target_gc_min Minimum acceptable GC% (default 40)
#' @param target_gc_max Maximum acceptable GC% (default 60)
#' @param target_tm Target Tm for scoring (default NULL = no Tm penalty)
#' @param specificity_score Score from validation (higher = more specific)
#' @return Numeric score (higher is better)
#' @export
score_probe <- function(probe_row, target_gc_min = 40, target_gc_max = 60,
                        target_tm = NULL, specificity_score = NULL) {
  score <- 100  # Base score

  # GC content penalty
  gc <- probe_row$gc_content
  if (gc < target_gc_min) {
    score <- score - (target_gc_min - gc) * 2
  } else if (gc > target_gc_max) {
    score <- score - (gc - target_gc_max) * 2
  }

  # Tm penalty (if target specified)
  if (!is.null(target_tm)) {
    tm_diff <- abs(probe_row$tm_nn - target_tm)
    score <- score - tm_diff * 0.5
  }

  # Specificity bonus (if provided)
  if (!is.null(specificity_score)) {
    score <- score + specificity_score * 50
  }

  # Length preference (slightly favor middle of range)
  len <- probe_row$length
  if (len == 22 || len == 23) {
    score <- score + 2  # Slight bonus for "sweet spot" length
  }

  max(0, score)  # Don't go negative
}

#' Design probes for specified tRNA target(s)
#'
#' Main entry point for probe design. Uses similarity data to recommend regions.
#'
#' @param target_ids Vector of tRNA IDs (for single mode) or family_id (for family mode)
#' @param trna_df Data frame from load_organism_data
#' @param similarity_data Pre-computed similarity data (optional but recommended)
#' @param target_mode "single" or "family"
#' @param min_length Minimum probe length (default 20)
#' @param max_length Maximum probe length (default 25)
#' @param region Region constraint (default NULL = use recommended regions)
#' @param top_n Return top N probes (default 20)
#' @return Data frame of candidate probes with scores
#' @export
design_probes <- function(target_ids, trna_df, similarity_data = NULL,
                          target_mode = "single",
                          min_length = 20, max_length = 25,
                          region = NULL, top_n = 20) {

  # Get target sequence(s)
  if (target_mode == "family") {
    # target_ids is a family_id
    targets <- trna_df[trna_df$family_id == target_ids, ]
    if (nrow(targets) == 0) {
      stop("No tRNAs found for family: ", target_ids)
    }
    # Use first sequence as representative (family members often identical)
    target_seq <- targets$sequence[1]
    target_id <- targets$id[1]
  } else {
    # Single target mode
    target_id <- target_ids[1]
    target_row <- trna_df[trna_df$id == target_id, ]
    if (nrow(target_row) == 0) {
      stop("Target ID not found: ", target_id)
    }
    target_seq <- target_row$sequence
  }

  # If similarity data provided and no region specified, get recommendations
  recommended_regions <- NULL
  if (!is.null(similarity_data) && is.null(region)) {
    recommended_regions <- recommend_probe_regions(
      target_id, trna_df, similarity_data,
      min_length = min_length, max_length = max_length
    )
  }

  # Generate candidate probes
  probes <- generate_candidate_probes(
    target_seq,
    min_length = min_length,
    max_length = max_length,
    region = region
  )

  if (nrow(probes) == 0) {
    warning("No valid probes generated for target: ", target_id)
    return(probes)
  }

  # Add target info
  probes$target_id <- target_id
  if (target_mode == "family") {
    probes$family_id <- target_ids
  }

  # If we have recommended regions, add divergence scores
  if (!is.null(recommended_regions)) {
    # Match probes to recommended regions by position
    probes$divergence_score <- sapply(1:nrow(probes), function(i) {
      # Find matching region in recommendations
      match_idx <- which(
        recommended_regions$start == probes$start[i] &
          recommended_regions$end == probes$end[i]
      )
      if (length(match_idx) > 0) {
        recommended_regions$score[match_idx[1]]
      } else {
        0
      }
    })
  } else {
    probes$divergence_score <- NA_real_
  }

  # Score probes
  probes$score <- sapply(1:nrow(probes), function(i) {
    score_probe(
      probes[i, ],
      specificity_score = if (!is.na(probes$divergence_score[i])) probes$divergence_score[i] else NULL
    )
  })

  # Sort by score and return top N
  probes <- probes[order(-probes$score), ]
  probes$rank <- 1:nrow(probes)

  if (nrow(probes) > top_n) {
    probes <- probes[1:top_n, ]
  }

  probes
}

#' Design probes using hierarchical target selection
#'
#' Uses a target selection (desired/avoid groups) to find optimal probe regions
#' that hit all desired targets while avoiding off-targets.
#'
#' @param selection Target selection from create_target_selection
#' @param trna_df Full data frame from load_organism_data
#' @param min_length Minimum probe length (default 20)
#' @param max_length Maximum probe length (default 25)
#' @param min_conservation Minimum conservation in desired group (default 90%)
#' @param min_divergence Minimum divergence from avoid group (default 30%)
#' @param top_n Return top N probes (default 20)
#' @return List with: probes (data frame), selection_analysis, region_analysis
#' @export
design_probes_selective <- function(selection,
                                    trna_df,
                                    min_length = 20,
                                    max_length = 25,
                                    min_conservation = 90,
                                    min_divergence = 30,
                                    top_n = 20,
                                    region_pref = "any",
                                    avoid_anticodon = TRUE) {

  # Analyze conservation within desired group
  conservation <- analyze_group_conservation(selection$desired)

  # Get reference sequence length for region filtering
  ref_seq <- selection$desired$sequence[1]
  ref_len <- nchar(ref_seq)
  midpoint <- ref_len / 2

  # Get anticodon position from reference (if available)
  ref_row <- selection$desired[1, ]
  anticodon_start <- ref_row$anticodon_start
  anticodon_end <- ref_row$anticodon_end
  # Expand anticodon zone to include 1 nt on each side (heavily modified region)
  if (!is.na(anticodon_start) && !is.na(anticodon_end)) {
    anticodon_zone_start <- max(1, anticodon_start - 1)
    anticodon_zone_end <- min(ref_len, anticodon_end + 1)
  } else {
    # Default anticodon position if not available (~34-36 in standard numbering)
    anticodon_zone_start <- 33
    anticodon_zone_end <- 37
  }

  # Find selective regions
  regions <- find_selective_regions(
    selection,
    min_length = min_length,
    max_length = max_length,
    min_conservation = min_conservation,
    min_divergence = min_divergence
  )

  if (nrow(regions) == 0) {
    warning("No regions found meeting conservation/divergence criteria. ",
            "Try relaxing min_conservation or min_divergence.")
    return(list(
      probes = data.frame(),
      conservation_analysis = conservation,
      region_analysis = regions
    ))
  }

  # Use first desired target as reference
  ref_id <- selection$desired$id[1]
  ref_seq <- selection$desired$sequence[1]

  # Filter regions by preference
  if (region_pref == "5prime") {
    regions <- regions[regions$end <= midpoint + 5, ]  # Allow slight overlap
  } else if (region_pref == "3prime") {
    regions <- regions[regions$start >= midpoint - 5, ]
  }

  if (nrow(regions) == 0) {
    warning("No regions found in preferred region. Try 'any' region preference.")
    return(list(
      probes = data.frame(),
      conservation_analysis = conservation,
      region_analysis = data.frame()
    ))
  }

  # Filter out anticodon-overlapping regions FIRST if avoidance is enabled
  # This ensures we prioritize non-anticodon regions
  if (avoid_anticodon) {
    # Mark which regions overlap anticodon
    regions$overlaps_ac <- (regions$start <= anticodon_zone_end & regions$end >= anticodon_zone_start)

    # Try to get non-overlapping regions first
    non_ac_regions <- regions[!regions$overlaps_ac, ]

    if (nrow(non_ac_regions) > 0) {
      # Use non-anticodon regions preferentially
      best_regions <- head(non_ac_regions[non_ac_regions$quality %in% c("EXCELLENT", "GOOD", "FAIR"), ], 10)
      if (nrow(best_regions) < 5) {
        # Add more non-anticodon regions even if lower quality
        best_regions <- head(non_ac_regions, 10)
      }
      # Also include some anticodon regions for comparison (but fewer)
      ac_regions <- head(regions[regions$overlaps_ac & regions$quality %in% c("EXCELLENT", "GOOD"), ], 5)
      best_regions <- rbind(best_regions, ac_regions)
    } else {
      # No non-anticodon regions available - use all with warning
      warning("All viable probe regions overlap the anticodon. Consider this when interpreting results.")
      best_regions <- head(regions[regions$quality %in% c("EXCELLENT", "GOOD", "FAIR"), ], 15)
    }
  } else {
    # No anticodon avoidance - just take top regions
    best_regions <- head(regions[regions$quality %in% c("EXCELLENT", "GOOD", "FAIR"), ], 15)
  }

  if (nrow(best_regions) == 0) {
    best_regions <- head(regions, 10)  # Fall back to top 10 regardless of quality
  }

  all_probes <- list()

  for (i in 1:nrow(best_regions)) {
    reg <- best_regions[i, ]
    target_region <- substr(ref_seq, reg$start, reg$end)
    probe_seq <- reverse_complement(target_region)

    # Determine which region of tRNA this probe targets
    probe_midpoint <- (reg$start + reg$end) / 2
    if (probe_midpoint <= midpoint) {
      trna_region <- "5' half"
    } else {
      trna_region <- "3' half"
    }

    # Check anticodon overlap
    overlaps_anticodon <- (reg$start <= anticodon_zone_end && reg$end >= anticodon_zone_start)

    all_probes[[i]] <- data.frame(
      start = reg$start,
      end = reg$end,
      length = reg$length,
      target_region = target_region,
      probe_sequence = probe_seq,
      gc_content = calculate_gc_content(probe_seq),
      tm_basic = calculate_tm_basic(probe_seq),
      tm_nn = calculate_tm_nn(probe_seq),
      desired_conservation = reg$mean_conservation,
      avoid_divergence = reg$mean_divergence,
      selectivity_score = reg$selectivity_score,
      quality = reg$quality,
      trna_region = trna_region,
      overlaps_anticodon = overlaps_anticodon,
      stringsAsFactors = FALSE
    )
  }

  probes <- bind_rows(all_probes)

  # Calculate modification penalty for each probe
  probes$modification_penalty <- sapply(1:nrow(probes), function(i) {
    calculate_modification_penalty(probes$start[i], probes$end[i], ref_len)
  })

  # Score probes with modification penalty
  probes$score <- sapply(1:nrow(probes), function(i) {
    base_score <- score_probe(probes[i, ])
    # Bonus for selectivity
    selectivity_bonus <- probes$selectivity_score[i] / 2

    # Modification penalty - probes overlapping modified regions score lower
    # This replaces the simple anticodon penalty with position-aware scoring
    mod_penalty <- probes$modification_penalty[i]

    base_score + selectivity_bonus - mod_penalty
  })

  # Sort and rank
  probes <- probes[order(-probes$score), ]
  probes$rank <- 1:nrow(probes)

  if (nrow(probes) > top_n) {
    probes <- probes[1:top_n, ]
  }

  # Add reference info
  probes$reference_id <- ref_id

  list(
    probes = probes,
    conservation_analysis = conservation,
    region_analysis = regions,
    selection = selection
  )
}

# =============================================================================
# Coverage Analysis Functions
# =============================================================================

#' Calculate which targets a probe covers
#'
#' Checks binding of a probe sequence against multiple target sequences.
#' A probe "covers" a target if the number of mismatches is <= max_mismatches.
#'
#' @param probe_sequence The probe sequence (5'->3', what you order)
#' @param probe_start Start position of probe binding site
#' @param probe_end End position of probe binding site
#' @param targets_df Data frame with id and sequence columns
#' @param max_mismatches Maximum mismatches to consider "covered" (default 3)
#' @return Data frame with target_id, binds (logical), mismatches (count)
#' @export
calculate_probe_coverage <- function(probe_sequence, probe_start, probe_end,
                                     targets_df, max_mismatches = 3) {
  # Reverse the probe to get 3'->5' for alignment
  probe_reversed <- paste(rev(strsplit(probe_sequence, "")[[1]]), collapse = "")
  probe_chars <- strsplit(probe_reversed, "")[[1]]
  probe_len <- length(probe_chars)

  results <- lapply(seq_len(nrow(targets_df)), function(i) {
    target_id <- targets_df$id[i]
    target_seq <- targets_df$sequence[i]
    target_len <- nchar(target_seq)

    # Extract binding region from target
    # Handle case where probe extends beyond target
    actual_start <- max(1, probe_start)
    actual_end <- min(target_len, probe_end)

    if (actual_end < actual_start || (actual_end - actual_start + 1) < probe_len * 0.5) {
      # Target too short or binding region invalid
      return(data.frame(
        target_id = target_id,
        binds = FALSE,
        mismatches = NA_integer_,
        stringsAsFactors = FALSE
      ))
    }

    binding_region <- substr(target_seq, actual_start, actual_end)
    region_chars <- strsplit(binding_region, "")[[1]]

    # Count mismatches (Watson-Crick complement check)
    mismatches <- 0
    for (j in seq_along(region_chars)) {
      if (j > probe_len) break

      target_base <- toupper(region_chars[j])
      probe_base <- toupper(probe_chars[j])

      # Check Watson-Crick pairing (probe binds to complement)
      is_match <- (target_base == "A" && probe_base == "T") ||
                  (target_base == "T" && probe_base == "A") ||
                  (target_base == "U" && probe_base == "A") ||
                  (target_base == "G" && probe_base == "C") ||
                  (target_base == "C" && probe_base == "G")

      if (!is_match) mismatches <- mismatches + 1
    }

    # Account for length differences as mismatches
    len_diff <- abs(length(region_chars) - probe_len)
    mismatches <- mismatches + len_diff

    data.frame(
      target_id = target_id,
      binds = mismatches <= max_mismatches,
      mismatches = mismatches,
      stringsAsFactors = FALSE
    )
  })

  bind_rows(results)
}

#' Build coverage matrix for all probes against all targets
#'
#' Creates a matrix showing which probes cover which targets, plus summary stats.
#'
#' @param probes Data frame of probes (from design_probes_selective)
#' @param targets_df Data frame of target tRNAs
#' @param max_mismatches Maximum mismatches for coverage (default 3)
#' @return List with:
#'   - matrix: logical matrix [probe_rank, target_id]
#'   - probes: updated probes df with coverage columns
#'   - targets_covered_by: list mapping target_id -> probe_ranks that cover it
#' @export
build_coverage_matrix <- function(probes, targets_df, max_mismatches = 3) {
  n_probes <- nrow(probes)
  n_targets <- nrow(targets_df)
  target_ids <- targets_df$id

  # Initialize coverage matrix
  cov_matrix <- matrix(FALSE, nrow = n_probes, ncol = n_targets)
  rownames(cov_matrix) <- probes$rank

colnames(cov_matrix) <- target_ids

  # Also track mismatches for each probe-target pair
  mismatch_matrix <- matrix(NA_integer_, nrow = n_probes, ncol = n_targets)
  rownames(mismatch_matrix) <- probes$rank
  colnames(mismatch_matrix) <- target_ids

  # Calculate coverage for each probe
  for (i in seq_len(n_probes)) {
    coverage <- calculate_probe_coverage(
      probe_sequence = probes$probe_sequence[i],
      probe_start = probes$start[i],
      probe_end = probes$end[i],
      targets_df = targets_df,
      max_mismatches = max_mismatches
    )

    for (j in seq_len(nrow(coverage))) {
      target_idx <- which(target_ids == coverage$target_id[j])
      if (length(target_idx) > 0) {
        cov_matrix[i, target_idx] <- coverage$binds[j]
        mismatch_matrix[i, target_idx] <- coverage$mismatches[j]
      }
    }
  }

  # Add coverage columns to probes dataframe
  probes$n_targets_hit <- rowSums(cov_matrix)
  probes$coverage_pct <- round(probes$n_targets_hit / n_targets * 100, 1)
  probes$targets_covered <- sapply(seq_len(n_probes), function(i) {
    covered_ids <- target_ids[cov_matrix[i, ]]
    paste(covered_ids, collapse = ",")
  })

  # Calculate "new coverage" - targets not covered by higher-ranked probes
  probes$new_coverage <- 0L
  covered_so_far <- logical(n_targets)
  for (i in seq_len(n_probes)) {
    new_hits <- cov_matrix[i, ] & !covered_so_far
    probes$new_coverage[i] <- sum(new_hits)
    covered_so_far <- covered_so_far | cov_matrix[i, ]
  }

  # Build reverse mapping: which probes cover each target
  targets_covered_by <- lapply(target_ids, function(tid) {
    idx <- which(target_ids == tid)
    probe_ranks <- which(cov_matrix[, idx])
    as.integer(probe_ranks)
  })
  names(targets_covered_by) <- target_ids

  list(
    matrix = cov_matrix,
    mismatch_matrix = mismatch_matrix,
    probes = probes,
    targets_covered_by = targets_covered_by,
    n_targets = n_targets
  )
}

#' Estimate minimum probes needed for various coverage levels
#'
#' Uses greedy set cover algorithm to find probe combinations that
#' achieve different coverage levels.
#'
#' @param coverage_matrix Coverage matrix from build_coverage_matrix()
#' @param target_ids Vector of target IDs
#' @return Data frame with n_probes, coverage_count, coverage_pct, probe_ranks
#' @export
estimate_probes_needed <- function(coverage_matrix, target_ids) {
  n_targets <- length(target_ids)
  n_probes <- nrow(coverage_matrix)

  if (n_probes == 0 || n_targets == 0) {
    return(data.frame(
      n_probes = integer(),
      coverage_count = integer(),
      coverage_pct = numeric(),
      probe_ranks = character(),
      stringsAsFactors = FALSE
    ))
  }

  # Greedy set cover
  covered <- logical(n_targets)
  selected_probes <- integer()
  results <- list()

  remaining_probes <- seq_len(n_probes)

  while (sum(covered) < n_targets && length(remaining_probes) > 0) {
    # Find probe that covers most uncovered targets
    best_probe <- NULL
    best_new_coverage <- 0

    for (p in remaining_probes) {
      new_coverage <- sum(coverage_matrix[p, ] & !covered)
      if (new_coverage > best_new_coverage) {
        best_new_coverage <- new_coverage
        best_probe <- p
      }
    }

    if (is.null(best_probe) || best_new_coverage == 0) break

    # Add this probe
    selected_probes <- c(selected_probes, best_probe)
    covered <- covered | coverage_matrix[best_probe, ]
    remaining_probes <- setdiff(remaining_probes, best_probe)

    # Record this step
    results[[length(results) + 1]] <- data.frame(
      n_probes = length(selected_probes),
      coverage_count = sum(covered),
      coverage_pct = round(sum(covered) / n_targets * 100, 1),
      probe_ranks = paste(selected_probes, collapse = ","),
      stringsAsFactors = FALSE
    )
  }

  if (length(results) == 0) {
    return(data.frame(
      n_probes = 0L,
      coverage_count = 0L,
      coverage_pct = 0,
      probe_ranks = "",
      stringsAsFactors = FALSE
    ))
  }

  bind_rows(results)
}

#' Get cumulative coverage for a set of selected probes
#'
#' @param selected_ranks Vector of probe ranks that user has selected
#' @param coverage_matrix Coverage matrix from build_coverage_matrix()
#' @param target_ids Vector of target IDs
#' @return List with coverage stats and lists of covered/uncovered targets
#' @export
get_cumulative_coverage <- function(selected_ranks, coverage_matrix, target_ids) {
  n_targets <- length(target_ids)

  if (length(selected_ranks) == 0) {
    return(list(
      n_selected = 0,
      n_covered = 0,
      n_uncovered = n_targets,
      coverage_pct = 0,
      covered_ids = character(),
      uncovered_ids = target_ids
    ))
  }

  # Combine coverage from all selected probes
  combined_coverage <- logical(n_targets)
  for (rank in selected_ranks) {
    if (rank <= nrow(coverage_matrix)) {
      combined_coverage <- combined_coverage | coverage_matrix[rank, ]
    }
  }

  covered_ids <- target_ids[combined_coverage]
  uncovered_ids <- target_ids[!combined_coverage]

  list(
    n_selected = length(selected_ranks),
    n_covered = length(covered_ids),
    n_uncovered = length(uncovered_ids),
    coverage_pct = round(length(covered_ids) / n_targets * 100, 1),
    covered_ids = covered_ids,
    uncovered_ids = uncovered_ids
  )
}

#' Re-rank probes for optimal coverage using greedy set cover
#'
#' Takes probes ranked by individual quality and re-ranks them so that
#' each successive probe maximizes NEW target coverage. This produces
#' a complementary probe set rather than variations of the same best probe.
#'
#' @param probes Data frame of probes with original quality-based ranking
#' @param coverage_matrix Coverage matrix from build_coverage_matrix()
#' @param target_ids Vector of target IDs
#' @return Data frame of probes with new coverage-optimized ranking
#' @export
rerank_probes_for_coverage <- function(probes, coverage_matrix, target_ids) {
  n_probes <- nrow(probes)
  n_targets <- length(target_ids)

  if (n_probes == 0 || n_targets == 0) {
    return(probes)
  }

  # Track which targets are covered and which probes are assigned
  covered <- logical(n_targets)
  remaining_probes <- seq_len(n_probes)
  new_order <- integer()

  # Greedy set cover: pick probe covering most uncovered targets
  # Use original quality score as tiebreaker
  while (length(remaining_probes) > 0) {
    best_probe <- NULL
    best_new_coverage <- -1
    best_quality <- -Inf

    for (probe_idx in remaining_probes) {
      # Count new targets this probe would cover
      probe_covers <- coverage_matrix[probe_idx, ]
      new_coverage <- sum(probe_covers & !covered)

      # Get quality score for tiebreaker (higher is better)
      quality <- probes$selectivity_score[probe_idx]

      # Select if covers more new targets, or same but higher quality
      if (new_coverage > best_new_coverage ||
          (new_coverage == best_new_coverage && quality > best_quality)) {
        best_probe <- probe_idx
        best_new_coverage <- new_coverage
        best_quality <- quality
      }
    }

    if (is.null(best_probe)) break

    # Add this probe to the new order
    new_order <- c(new_order, best_probe)
    covered <- covered | coverage_matrix[best_probe, ]
    remaining_probes <- setdiff(remaining_probes, best_probe)

    # If all targets covered, remaining probes get appended by quality
    if (all(covered)) {
      # Sort remaining by original quality score (descending)
      if (length(remaining_probes) > 0) {
        remaining_by_quality <- remaining_probes[order(
          probes$selectivity_score[remaining_probes], decreasing = TRUE
        )]
        new_order <- c(new_order, remaining_by_quality)
      }
      break
    }
  }

  # Reorder probes and assign new ranks
  probes_reordered <- probes[new_order, ]
  probes_reordered$original_rank <- probes_reordered$rank
  probes_reordered$rank <- seq_len(nrow(probes_reordered))

  probes_reordered
}
