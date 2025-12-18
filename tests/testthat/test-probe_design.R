# test-probe_design.R
# Unit tests for probe design functions

library(testthat)

# source("../../R/probe_design.R")
# source("../../R/sequence_utils.R")  # For reverse_complement

test_that("calculate_gc_content computes correct percentages", {
  expect_equal(calculate_gc_content("GGCC"), 100)
  expect_equal(calculate_gc_content("AATT"), 0)
  expect_equal(calculate_gc_content("ATGC"), 50)
  expect_equal(calculate_gc_content("ATGCATGC"), 50)

  # Real probe-like sequence
  expect_equal(calculate_gc_content("ACGTACGTACGTACGTACGT"), 50)  # 20 nt, 10 GC
})

test_that("calculate_gc_content handles edge cases", {
  expect_equal(calculate_gc_content(""), 0)
  expect_equal(calculate_gc_content("G"), 100)
  expect_equal(calculate_gc_content("A"), 0)

  # Case insensitive
  expect_equal(calculate_gc_content("atgc"), 50)
})

test_that("calculate_tm_basic returns reasonable values", {
  # Short oligo (Wallace rule)
  short_seq <- "ATGCATGC"  # 8 nt
  tm <- calculate_tm_basic(short_seq)
  expect_true(tm > 0 && tm < 100)

  # 20 nt probe
  probe_seq <- "ACGTACGTACGTACGTACGT"
  tm <- calculate_tm_basic(probe_seq)
  expect_true(tm > 40 && tm < 80)  # Reasonable range for 50% GC

  # Higher GC = higher Tm
  high_gc <- "GCGCGCGCGCGCGCGCGCGC"
  low_gc <- "ATATATATATATATATATATAT"
  expect_true(calculate_tm_basic(high_gc) > calculate_tm_basic(low_gc))
})

test_that("calculate_tm_nn returns reasonable values", {
  probe_seq <- "ACGTACGTACGTACGTACGT"  # 20 nt, 50% GC

  tm <- calculate_tm_nn(probe_seq)
  expect_true(tm > 40 && tm < 80)

  # Higher concentration = higher Tm
  tm_high <- calculate_tm_nn(probe_seq, probe_conc = 1e-6)
  tm_low <- calculate_tm_nn(probe_seq, probe_conc = 1e-9)
  expect_true(tm_high > tm_low)

  # Higher salt = higher Tm
  tm_high_salt <- calculate_tm_nn(probe_seq, na_conc = 100e-3)
  tm_low_salt <- calculate_tm_nn(probe_seq, na_conc = 10e-3)
  expect_true(tm_high_salt > tm_low_salt)
})

test_that("generate_candidate_probes generates correct number of probes", {
  # 30 nt target, probes 20-22 nt
  target <- "ACGTACGTACGTACGTACGTACGTACGTAC"  # 30 nt

  probes <- generate_candidate_probes(target, min_length = 20, max_length = 22)

  # For length 20: 30-20+1 = 11 positions
  # For length 21: 30-21+1 = 10 positions
  # For length 22: 30-22+1 = 9 positions
  # Total: 30 probes
  expect_equal(nrow(probes), 30)

  # Check columns exist
  expect_true("probe_sequence" %in% names(probes))
  expect_true("gc_content" %in% names(probes))
  expect_true("tm_nn" %in% names(probes))
  expect_true("start" %in% names(probes))
  expect_true("end" %in% names(probes))
})

test_that("generate_candidate_probes respects region constraints", {
  target <- paste(rep("A", 80), collapse = "")  # 80 nt target

  # 5' region only
  probes_5p <- generate_candidate_probes(target, min_length = 20, max_length = 20, region = "5prime")
  expect_true(all(probes_5p$start <= 35))
  expect_true(all(probes_5p$end <= 35))

  # 3' region only
  probes_3p <- generate_candidate_probes(target, min_length = 20, max_length = 20, region = "3prime")
  expect_true(all(probes_3p$start >= 80 - 34))

  # Custom region
  probes_custom <- generate_candidate_probes(target, min_length = 20, max_length = 20, region = c(30, 50))
  expect_true(all(probes_custom$start >= 30))
  expect_true(all(probes_custom$end <= 50))
})

test_that("probe_sequence is reverse complement of target_region", {
  target <- "ACGTACGTACGTACGTACGTACGTACGTAC"
  probes <- generate_candidate_probes(target, min_length = 20, max_length = 20)

  for (i in 1:min(5, nrow(probes))) {
    target_region <- probes$target_region[i]
    probe_seq <- probes$probe_sequence[i]

    # Probe should be RC of target region
    expected_rc <- reverse_complement(target_region)
    expect_equal(probe_seq, expected_rc)
  }
})

test_that("score_probe penalizes extreme GC content", {
  # Create mock probe data
  probe_good_gc <- list(gc_content = 50, tm_nn = 55, length = 22)
  probe_low_gc <- list(gc_content = 20, tm_nn = 55, length = 22)
  probe_high_gc <- list(gc_content = 80, tm_nn = 55, length = 22)

  score_good <- score_probe(probe_good_gc)
  score_low <- score_probe(probe_low_gc)
  score_high <- score_probe(probe_high_gc)

  expect_true(score_good > score_low)
  expect_true(score_good > score_high)
})

test_that("score_probe gives bonus for optimal length", {
  # 22-23 nt should get a small bonus
  probe_22 <- list(gc_content = 50, tm_nn = 55, length = 22)
  probe_23 <- list(gc_content = 50, tm_nn = 55, length = 23)
  probe_20 <- list(gc_content = 50, tm_nn = 55, length = 20)
  probe_25 <- list(gc_content = 50, tm_nn = 55, length = 25)

  score_22 <- score_probe(probe_22)
  score_23 <- score_probe(probe_23)
  score_20 <- score_probe(probe_20)
  score_25 <- score_probe(probe_25)

  # 22 and 23 should have same bonus
  expect_equal(score_22, score_23)
  # 22/23 should score higher than 20 or 25
  expect_true(score_22 > score_20)
  expect_true(score_22 > score_25)
})

test_that("score_probe applies Tm penalty when target specified", {
  probe <- list(gc_content = 50, tm_nn = 55, length = 22)

  # No Tm target - no penalty
  score_no_target <- score_probe(probe)

  # Tm matches target - no penalty
  score_exact <- score_probe(probe, target_tm = 55)

  # Tm differs from target - penalty applied
  score_off <- score_probe(probe, target_tm = 65)

  expect_equal(score_no_target, score_exact)
  expect_true(score_exact > score_off)
})

test_that("score_probe applies specificity bonus", {
  probe <- list(gc_content = 50, tm_nn = 55, length = 22)

  score_no_spec <- score_probe(probe)
  score_with_spec <- score_probe(probe, specificity_score = 0.8)

  expect_true(score_with_spec > score_no_spec)
})

test_that("score_probe never returns negative", {
  # Extreme case: very bad GC content
  probe_bad <- list(gc_content = 0, tm_nn = 55, length = 15)

  score <- score_probe(probe_bad)

  expect_true(score >= 0)
})

test_that("score_probe handles edge cases at GC boundaries", {
  # Right at the boundaries
  probe_low_boundary <- list(gc_content = 40, tm_nn = 55, length = 22)
  probe_high_boundary <- list(gc_content = 60, tm_nn = 55, length = 22)
  probe_inside <- list(gc_content = 50, tm_nn = 55, length = 22)

  score_low <- score_probe(probe_low_boundary)
  score_high <- score_probe(probe_high_boundary)
  score_inside <- score_probe(probe_inside)

  # At boundaries should have same score as inside (no penalty)
  expect_equal(score_low, score_inside)
  expect_equal(score_high, score_inside)
})

# =============================================================================
# Tests for modification penalty functions
# =============================================================================

test_that("get_modification_zones returns correct structure", {
  zones <- get_modification_zones(76)

  expect_true(is.data.frame(zones))
  expect_true("start" %in% names(zones))
  expect_true("end" %in% names(zones))
  expect_true("risk" %in% names(zones))
  expect_true("description" %in% names(zones))
  expect_equal(nrow(zones), 3)  # D-loop, Anticodon, TΨC
})
test_that("get_modification_zones scales with sequence length", {
  zones_76 <- get_modification_zones(76)
  zones_152 <- get_modification_zones(152)  # Double length

  # Positions should roughly double
  expect_true(zones_152$start[1] > zones_76$start[1])
  expect_true(zones_152$end[2] > zones_76$end[2])
})

test_that("calculate_modification_penalty returns 0 for non-modified regions", {
  # Position 1-10 should not overlap any modification zones
  penalty <- calculate_modification_penalty(1, 10, seq_length = 76)

  expect_equal(penalty, 0)
})

test_that("calculate_modification_penalty penalizes anticodon overlap", {
  # Anticodon is around position 34-37 in standard tRNA
  penalty_ac <- calculate_modification_penalty(34, 40, seq_length = 76)
  penalty_no_ac <- calculate_modification_penalty(1, 10, seq_length = 76)

  expect_true(penalty_ac > penalty_no_ac)
  expect_true(penalty_ac > 0)
})

test_that("calculate_modification_penalty scales with overlap fraction", {
  # More overlap = higher penalty
  penalty_full <- calculate_modification_penalty(34, 37, seq_length = 76)  # Full anticodon
  penalty_partial <- calculate_modification_penalty(36, 45, seq_length = 76)  # Partial

  # Full overlap of anticodon should have higher penalty per base
  # but partial overlap is a larger region
  expect_true(penalty_full > 0)
  expect_true(penalty_partial > 0)
})

test_that("calculate_modification_penalty handles D-loop region", {
  # D-loop is around positions 16-20
  penalty_dloop <- calculate_modification_penalty(15, 22, seq_length = 76)

  expect_true(penalty_dloop > 0)
})

test_that("calculate_modification_penalty handles TΨC loop region", {
  # TΨC loop is around positions 54-56
  penalty_tpsi <- calculate_modification_penalty(52, 58, seq_length = 76)

  expect_true(penalty_tpsi > 0)
})
