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
                                    top_n = 20) {

  # Analyze conservation within desired group
  conservation <- analyze_group_conservation(selection$desired)

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

  # Generate probes from top regions
  best_regions <- head(regions[regions$quality %in% c("EXCELLENT", "GOOD", "FAIR"), ], 10)

  if (nrow(best_regions) == 0) {
    best_regions <- head(regions, 5)  # Fall back to top 5 regardless of quality
  }

  all_probes <- list()

  for (i in 1:nrow(best_regions)) {
    reg <- best_regions[i, ]
    target_region <- substr(ref_seq, reg$start, reg$end)
    probe_seq <- reverse_complement(target_region)

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
      stringsAsFactors = FALSE
    )
  }

  probes <- bind_rows(all_probes)

  # Score probes
  probes$score <- sapply(1:nrow(probes), function(i) {
    base_score <- score_probe(probes[i, ])
    # Bonus for selectivity
    selectivity_bonus <- probes$selectivity_score[i] / 2
    base_score + selectivity_bonus
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
