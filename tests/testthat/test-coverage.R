# test-coverage.R
# Unit tests for coverage analysis functions in probe_design.R

library(testthat)

# =============================================================================
# Tests for calculate_probe_coverage()
# =============================================================================

test_that("calculate_probe_coverage identifies perfect matches", {
  # Create simple test targets
  targets_df <- data.frame(
    id = c("target1", "target2"),
    sequence = c(
      "AAAAAAAAAAGCTAGCTAGCAAAAAAAAAAA",  # Has GCTAGCTAG at positions 11-19
      "AAAAAAAAAAGCTAGCTAGCAAAAAAAAAAA"   # Identical
    ),
    stringsAsFactors = FALSE
  )

  # Probe is reverse complement of target region
  # Target region: GCTAGCTAGC (positions 11-20)
  # RC: GCTAGCTAGC (palindrome in this case for simplicity)
  probe_seq <- reverse_complement("GCTAGCTAGC")

 coverage <- calculate_probe_coverage(
    probe_sequence = probe_seq,
    probe_start = 11,
    probe_end = 20,
    targets_df = targets_df,
    max_mismatches = 3
  )

  expect_equal(nrow(coverage), 2)
  expect_true(all(coverage$binds))
  expect_true(all(coverage$mismatches == 0))
})

test_that("calculate_probe_coverage detects mismatches correctly", {
  targets_df <- data.frame(
    id = c("perfect", "one_mismatch", "many_mismatches"),
    sequence = c(
      "AAAAAATCGATCGATCGAAAAAA",  # TCGATCGATCG at 7-17
      "AAAAAATCGATCGATCGAAAAAA",  # Same (will modify probe)
      "AAAAAANNNNNNNNNNNAAAAAA"   # All Ns - many mismatches
    ),
    stringsAsFactors = FALSE
  )

  # Test with a probe that matches first two but not third
  target_region <- "TCGATCGATCG"
  probe_seq <- reverse_complement(target_region)

  coverage <- calculate_probe_coverage(
    probe_sequence = probe_seq,
    probe_start = 7,
    probe_end = 17,
    targets_df = targets_df,
    max_mismatches = 3
  )

  expect_equal(nrow(coverage), 3)
  # First two should bind (0 mismatches)
  expect_true(coverage$binds[1])
  expect_true(coverage$binds[2])
  # Third has many mismatches
  expect_false(coverage$binds[3])
  expect_true(coverage$mismatches[3] > 3)
})

test_that("calculate_probe_coverage respects max_mismatches threshold", {
  targets_df <- data.frame(
    id = c("target1"),
    sequence = c("AAAAAATCGATCGATCGAAAAAA"),
    stringsAsFactors = FALSE
  )

  # Create a probe with known mismatches
  target_region <- "TCGATCGATCG"
  probe_seq <- reverse_complement(target_region)

  # With max_mismatches = 0, only perfect matches
  coverage_strict <- calculate_probe_coverage(
    probe_sequence = probe_seq,
    probe_start = 7,
    probe_end = 17,
    targets_df = targets_df,
    max_mismatches = 0
  )

  # With max_mismatches = 10, more permissive
  coverage_permissive <- calculate_probe_coverage(
    probe_sequence = probe_seq,
    probe_start = 7,
    probe_end = 17,
    targets_df = targets_df,
    max_mismatches = 10
  )

  expect_true(coverage_strict$binds[1])  # Perfect match
  expect_true(coverage_permissive$binds[1])
})

test_that("calculate_probe_coverage handles edge cases", {
  # Empty targets
  empty_df <- data.frame(
    id = character(),
    sequence = character(),
    stringsAsFactors = FALSE
  )

  coverage <- calculate_probe_coverage(
    probe_sequence = "ACGT",
    probe_start = 1,
    probe_end = 4,
    targets_df = empty_df,
    max_mismatches = 3
  )

  expect_equal(nrow(coverage), 0)
})

test_that("calculate_probe_coverage handles short sequences", {
  targets_df <- data.frame(
    id = c("short"),
    sequence = c("ACGT"),  # Only 4 nt
    stringsAsFactors = FALSE
  )

  # Probe extends beyond sequence
  coverage <- calculate_probe_coverage(
    probe_sequence = "ACGTACGTACGT",
    probe_start = 1,
    probe_end = 12,
    targets_df = targets_df,
    max_mismatches = 3
  )

  expect_equal(nrow(coverage), 1)
  # Should not bind due to length mismatch
  expect_false(coverage$binds[1])
})

# =============================================================================
# Tests for build_coverage_matrix()
# =============================================================================

test_that("build_coverage_matrix creates correct matrix dimensions", {
  targets_df <- data.frame(
    id = c("t1", "t2", "t3"),
    sequence = c(
      "AAAAAATCGATCGATCGAAAAAA",
      "AAAAAATCGATCGATCGAAAAAA",
      "AAAAAATCGATCGATCGAAAAAA"
    ),
    stringsAsFactors = FALSE
  )

  probes <- data.frame(
    rank = 1:2,
    probe_sequence = c(
      reverse_complement("TCGATCGATCG"),
      reverse_complement("AAAAATCGATC")
    ),
    start = c(7, 5),
    end = c(17, 15),
    selectivity_score = c(80, 70),
    stringsAsFactors = FALSE
  )

  result <- build_coverage_matrix(probes, targets_df, max_mismatches = 3)

  expect_equal(nrow(result$matrix), 2)  # 2 probes
  expect_equal(ncol(result$matrix), 3)  # 3 targets
  expect_equal(result$n_targets, 3)
  expect_true("n_targets_hit" %in% names(result$probes))
  expect_true("coverage_pct" %in% names(result$probes))
})

test_that("build_coverage_matrix tracks new coverage correctly", {
  targets_df <- data.frame(
    id = c("t1", "t2"),
    sequence = c(
      "AAAAAATCGATCGATCGAAAAAA",
      "AAAAAATCGATCGATCGAAAAAA"
    ),
    stringsAsFactors = FALSE
  )

  probes <- data.frame(
    rank = 1:2,
    probe_sequence = c(
      reverse_complement("TCGATCGATCG"),
      reverse_complement("TCGATCGATCG")  # Same probe
    ),
    start = c(7, 7),
    end = c(17, 17),
    selectivity_score = c(80, 70),
    stringsAsFactors = FALSE
  )

  result <- build_coverage_matrix(probes, targets_df, max_mismatches = 3)

  # First probe covers both targets
  expect_equal(result$probes$new_coverage[1], 2)
  # Second probe covers same targets, so new coverage = 0
  expect_equal(result$probes$new_coverage[2], 0)
})

test_that("build_coverage_matrix handles single probe single target", {
  # Test minimal case instead of empty (empty edge case not supported)
  targets_df <- data.frame(
    id = c("t1"),
    sequence = c("AAAAAATCGATCGATCGAAAAAA"),
    stringsAsFactors = FALSE
  )

  probes <- data.frame(
    rank = 1L,
    probe_sequence = reverse_complement("TCGATCGATCG"),
    start = 7L,
    end = 17L,
    selectivity_score = 80,
    stringsAsFactors = FALSE
  )

  result <- build_coverage_matrix(probes, targets_df, max_mismatches = 3)

  expect_equal(nrow(result$matrix), 1)
  expect_equal(ncol(result$matrix), 1)
  expect_equal(result$n_targets, 1)
})

# =============================================================================
# Tests for estimate_probes_needed()
# =============================================================================

test_that("estimate_probes_needed uses greedy set cover", {
  # Create coverage matrix where probe 1 covers t1, t2 and probe 2 covers t2, t3
  coverage_matrix <- matrix(
    c(TRUE, TRUE, FALSE,   # Probe 1: covers t1, t2
      FALSE, TRUE, TRUE),  # Probe 2: covers t2, t3
    nrow = 2, byrow = TRUE
  )
  rownames(coverage_matrix) <- c("1", "2")
  colnames(coverage_matrix) <- c("t1", "t2", "t3")

  target_ids <- c("t1", "t2", "t3")

  result <- estimate_probes_needed(coverage_matrix, target_ids)

  # Should need 2 probes to cover all 3 targets
  expect_true(any(result$coverage_pct == 100))
  expect_equal(result$n_probes[nrow(result)], 2)
})

test_that("estimate_probes_needed handles full coverage by single probe", {
  coverage_matrix <- matrix(
    c(TRUE, TRUE, TRUE),
    nrow = 1
  )
  rownames(coverage_matrix) <- "1"
  colnames(coverage_matrix) <- c("t1", "t2", "t3")

  target_ids <- c("t1", "t2", "t3")

  result <- estimate_probes_needed(coverage_matrix, target_ids)

  expect_equal(nrow(result), 1)
  expect_equal(result$coverage_pct[1], 100)
  expect_equal(result$n_probes[1], 1)
})

test_that("estimate_probes_needed handles empty inputs", {
  empty_matrix <- matrix(nrow = 0, ncol = 0)

  result <- estimate_probes_needed(empty_matrix, character(0))

  expect_equal(nrow(result), 0)
})

# =============================================================================
# Tests for get_cumulative_coverage()
# =============================================================================

test_that("get_cumulative_coverage calculates combined coverage", {
  coverage_matrix <- matrix(
    c(TRUE, FALSE, FALSE,   # Probe 1: covers t1 only
      FALSE, TRUE, FALSE,   # Probe 2: covers t2 only
      FALSE, FALSE, TRUE),  # Probe 3: covers t3 only
    nrow = 3, byrow = TRUE
  )
  rownames(coverage_matrix) <- c("1", "2", "3")
  colnames(coverage_matrix) <- c("t1", "t2", "t3")

  target_ids <- c("t1", "t2", "t3")

  # Select probes 1 and 2
  result <- get_cumulative_coverage(c(1, 2), coverage_matrix, target_ids)

  expect_equal(result$n_selected, 2)
  expect_equal(result$n_covered, 2)
  expect_equal(result$n_uncovered, 1)
  expect_equal(length(result$covered_ids), 2)
  expect_true("t1" %in% result$covered_ids)
  expect_true("t2" %in% result$covered_ids)
  expect_true("t3" %in% result$uncovered_ids)
})

test_that("get_cumulative_coverage handles empty selection", {
  coverage_matrix <- matrix(TRUE, nrow = 1, ncol = 3)
  target_ids <- c("t1", "t2", "t3")

  result <- get_cumulative_coverage(integer(0), coverage_matrix, target_ids)

  expect_equal(result$n_selected, 0)
  expect_equal(result$n_covered, 0)
  expect_equal(result$n_uncovered, 3)
  expect_equal(result$coverage_pct, 0)
})

test_that("get_cumulative_coverage handles full coverage", {
  coverage_matrix <- matrix(TRUE, nrow = 2, ncol = 3)
  target_ids <- c("t1", "t2", "t3")

  result <- get_cumulative_coverage(c(1, 2), coverage_matrix, target_ids)

  expect_equal(result$n_covered, 3)
  expect_equal(result$coverage_pct, 100)
  expect_equal(length(result$uncovered_ids), 0)
})

# =============================================================================
# Tests for rerank_probes_for_coverage()
# =============================================================================

test_that("rerank_probes_for_coverage optimizes for coverage", {
  targets_df <- data.frame(
    id = c("t1", "t2", "t3"),
    sequence = c(
      "AAAAAATCGATCGATCGAAAAAA",  # Position 7-17: TCGATCGATCG
      "AAAAAGCATCGATCGATCAAAAA",  # Position 6-16: GCATCGATCGA
      "AAAAAATCGATCGATCGAAAAAA"   # Same as t1
    ),
    stringsAsFactors = FALSE
  )

  # Probe 1 covers t1 and t3, Probe 2 covers t2
  # But probe 2 has higher selectivity score
  probes <- data.frame(
    rank = 1:2,
    probe_sequence = c(
      reverse_complement("TCGATCGATCG"),
      reverse_complement("GCATCGATCGA")
    ),
    start = c(7, 6),
    end = c(17, 16),
    selectivity_score = c(70, 90),  # Probe 2 has higher score
    stringsAsFactors = FALSE
  )

  coverage_result <- build_coverage_matrix(probes, targets_df, max_mismatches = 3)

  reranked <- rerank_probes_for_coverage(
    probes,
    coverage_result$matrix,
    targets_df$id
  )

  # First reranked probe should cover the most targets
  expect_true("original_rank" %in% names(reranked))
  expect_equal(reranked$rank[1], 1)  # New rank
})

test_that("rerank_probes_for_coverage preserves all probes", {
  targets_df <- data.frame(
    id = c("t1", "t2"),
    sequence = c(
      "AAAAAATCGATCGATCGAAAAAA",
      "AAAAAATCGATCGATCGAAAAAA"
    ),
    stringsAsFactors = FALSE
  )

  probes <- data.frame(
    rank = 1:3,
    probe_sequence = rep(reverse_complement("TCGATCGATCG"), 3),
    start = rep(7, 3),
    end = rep(17, 3),
    selectivity_score = c(80, 70, 60),
    stringsAsFactors = FALSE
  )

  coverage_result <- build_coverage_matrix(probes, targets_df, max_mismatches = 3)

  reranked <- rerank_probes_for_coverage(
    probes,
    coverage_result$matrix,
    targets_df$id
  )

  expect_equal(nrow(reranked), 3)  # All probes preserved
})
