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
