# test-sequence_utils.R
# Unit tests for sequence utilities

library(testthat)

# Source the module (adjust path as needed when running)
# source("../../R/sequence_utils.R")

test_that("parse_trna_header parses E. coli headers correctly", {
  result <- parse_trna_header("tRNA-Ala-GGC-1-2", "ecoli")

  expect_true(is.na(result$compartment))
  expect_equal(result$amino_acid, "Ala")
  expect_equal(result$anticodon, "GGC")
  expect_equal(result$gene_family, 1L)
  expect_equal(result$copy_number, 2L)
})

test_that("parse_trna_header parses human nuclear headers correctly", {
  result <- parse_trna_header("nuc-tRNA-Leu-CAA-12-3", "human")

  expect_equal(result$compartment, "nuclear")
  expect_equal(result$amino_acid, "Leu")
  expect_equal(result$anticodon, "CAA")
  expect_equal(result$gene_family, 12L)
  expect_equal(result$copy_number, 3L)
})

test_that("parse_trna_header parses human mito headers correctly", {
  result <- parse_trna_header("mito-tRNA-Phe-GAA-1-1", "human")

  expect_equal(result$compartment, "mitochondrial")
  expect_equal(result$amino_acid, "Phe")
  expect_equal(result$anticodon, "GAA")
  expect_equal(result$gene_family, 1L)
  expect_equal(result$copy_number, 1L)
})

test_that("parse_trna_header handles simple yeast mito format", {
  result <- parse_trna_header("mito-tRNA-Pro-UGG", "yeast")

  expect_equal(result$compartment, "mitochondrial")
  expect_equal(result$amino_acid, "Pro")
  expect_equal(result$anticodon, "UGG")
  expect_equal(result$gene_family, 1L)  # Default
  expect_equal(result$copy_number, 1L)  # Default
})

test_that("reverse_complement works correctly", {
  # Simple case
  expect_equal(reverse_complement("ATCG"), "CGAT")

  # All bases
  expect_equal(reverse_complement("AATTCCGG"), "CCGGAATT")

  # Longer sequence
  seq <- "GGGGGTATAGCTCAGT"
  rc <- reverse_complement(seq)
  # Verify double RC returns original
  expect_equal(reverse_complement(rc), seq)
})

test_that("filter_trnas filters by compartment", {
  # Create mock data
  mock_df <- data.frame(
    id = c("nuc-1", "nuc-2", "mito-1"),
    compartment = c("nuclear", "nuclear", "mitochondrial"),
    amino_acid = c("Ala", "Leu", "Phe"),
    anticodon = c("AGC", "CAA", "GAA"),
    stringsAsFactors = FALSE
  )

  result <- filter_trnas(mock_df, compartment = "nuclear")
  expect_equal(nrow(result), 2)
  expect_true(all(result$compartment == "nuclear"))

  result <- filter_trnas(mock_df, compartment = "mitochondrial")
  expect_equal(nrow(result), 1)
})

test_that("filter_trnas filters by amino acid", {
  mock_df <- data.frame(
    id = c("1", "2", "3"),
    compartment = c("nuclear", "nuclear", "nuclear"),
    amino_acid = c("Ala", "Ala", "Leu"),
    anticodon = c("AGC", "TGC", "CAA"),
    stringsAsFactors = FALSE
  )

  result <- filter_trnas(mock_df, amino_acid = "Ala")
  expect_equal(nrow(result), 2)
  expect_true(all(result$amino_acid == "Ala"))
})

test_that("extract_region extracts correct subsequence", {
  seq <- "ABCDEFGHIJ"

  expect_equal(extract_region(seq, 1, 3), "ABC")
  expect_equal(extract_region(seq, 5, 7), "EFG")
  expect_equal(extract_region(seq, 1, 10), seq)

  # Error cases
expect_error(extract_region(seq, 0, 5))
  expect_error(extract_region(seq, 1, 15))
  expect_error(extract_region(seq, 5, 3))
})
