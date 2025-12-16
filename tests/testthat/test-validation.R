# test-validation.R
# Unit tests for validation functions

library(testthat)

# source("../../R/validation.R")
# source("../../R/sequence_utils.R")

test_that("count_mismatches_at_position counts correctly", {
  # Perfect match
  probe <- "ACGTACGT"  # Probe binds as RC
  target <- reverse_complement("ACGTACGT")  # Target has exact complement
  expect_equal(count_mismatches_at_position(probe, target, 1), 0)

  # One mismatch
  target_1mm <- paste0(substr(reverse_complement("ACGTACGT"), 1, 3), "A",
                       substr(reverse_complement("ACGTACGT"), 5, 8))
  # Note: this test setup is tricky because of RC
})

test_that("find_best_binding finds optimal position", {
  # Create a target where probe matches perfectly at position 10
  probe <- "GGGGGGGGGGGGGGGGGGGG"  # 20 G's, binds as 20 C's
  target <- paste0(
    paste(rep("A", 9), collapse = ""),   # pos 1-9: A's
    paste(rep("C", 20), collapse = ""),  # pos 10-29: C's (perfect match)
    paste(rep("A", 10), collapse = "")   # pos 30-39: A's
  )

  binding <- find_best_binding(probe, target)

  expect_equal(binding$position, 10)
  expect_equal(binding$mismatches, 0)
  expect_true(binding$binds)
})

test_that("find_best_binding handles no match case", {
  probe <- "GGGGGGGGGGGGGGGGGGGG"  # All G's
  target <- "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"  # All A's (40 nt)

  binding <- find_best_binding(probe, target, max_mismatches = 5)

  expect_equal(binding$mismatches, 20)  # All positions mismatch
  expect_false(binding$binds)  # Exceeds threshold
})

test_that("categorize_off_targets produces correct summary", {
  hits_df <- data.frame(
    id = c("target", "fam1", "fam2", "iso1", "other1"),
    category = c("target", "same_family", "same_family", "same_isoacceptor", "other"),
    mismatches = c(0, 0, 1, 2, 4),
    stringsAsFactors = FALSE
  )

  summary <- categorize_off_targets(hits_df)

  expect_equal(nrow(summary), 4)  # 4 categories
  expect_true("same_family" %in% summary$category)

  family_row <- summary[summary$category == "same_family", ]
  expect_equal(family_row$count, 2)
  expect_equal(family_row$min_mismatches, 0)
  expect_equal(family_row$max_mismatches, 1)
})

test_that("calculate_specificity_score penalizes off-targets", {
  # Perfect specificity - only target hit
  hits_perfect <- data.frame(
    id = c("target"),
    category = c("target"),
    mismatches = c(0),
    stringsAsFactors = FALSE
  )
  expect_equal(calculate_specificity_score(hits_perfect, "target"), 1.0)

  # Some off-targets
  hits_some <- data.frame(
    id = c("target", "off1", "off2"),
    category = c("target", "other", "other"),
    mismatches = c(0, 3, 4),
    stringsAsFactors = FALSE
  )
  score_some <- calculate_specificity_score(hits_some, "target")
  expect_true(score_some > 0 && score_some < 1)

  # Many exact match off-targets = worse score
  hits_bad <- data.frame(
    id = c("target", "off1", "off2", "off3"),
    category = c("target", "other", "other", "other"),
    mismatches = c(0, 0, 0, 1),
    stringsAsFactors = FALSE
  )
  score_bad <- calculate_specificity_score(hits_bad, "target")
  expect_true(score_bad < score_some)
})

test_that("generate_alignment_report produces output", {
  hits_df <- data.frame(
    id = c("target", "off1"),
    category = c("target", "other"),
    mismatches = c(0, 2),
    best_position = c(10, 15),
    amino_acid = c("Ala", "Leu"),
    anticodon = c("AGC", "CAA"),
    stringsAsFactors = FALSE
  )

  trna_df <- data.frame(
    id = c("target", "off1"),
    sequence = c(
      "AAAAAAAAAAACGTACGTACGTACGTACGTAAAAAAAAAA",
      "AAAAAAAAAAAAAAACGTACGTACGTACGTACGTAAAAAA"
    ),
    stringsAsFactors = FALSE
  )

  probe <- "ACGTACGTACGTACGTACGT"

  report <- generate_alignment_report(probe, hits_df, trna_df)

  expect_true(is.character(report))
  expect_true(nchar(report) > 0)
  expect_true(grepl("Probe:", report))
  expect_true(grepl("target", report))
})
