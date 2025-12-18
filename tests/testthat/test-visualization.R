# test-visualization.R
# Unit tests for visualization functions

library(testthat)
library(htmltools)

# =============================================================================
# Tests for find_anticodon_position()
# =============================================================================

test_that("find_anticodon_position finds anticodon in expected range", {
  # Standard tRNA-like sequence with AGC at position 35-37
  # Anticodon typically at position ~34-36
  test_seq <- paste0(
    "GGGGGTATAGCTCAGTGGTAGAGCGCGTGCTT",  # 32 nt
    "AGC",                                  # anticodon at 33-35
    "ATGCACGAGGTCCTGGGTTCGATCCCCAGTAC"    # remaining
  )

  result <- find_anticodon_position(test_seq, "AGC", min_pos = 30, max_pos = 45)

  expect_false(is.null(result))
  expect_equal(result$start, 33)
  expect_equal(result$end, 35)
  expect_true(result$validated)
})

test_that("find_anticodon_position converts U to T", {
  test_seq <- "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGCAAAAAAAAAA"  # AGC at position 33-35

  # Should find AGC when searching for UGC
  result_u <- find_anticodon_position(test_seq, "UGC", min_pos = 30, max_pos = 45)

  # UGC converts to TGC, not AGC, so should not find
  expect_null(result_u)

  # But AGC should find it
  result_a <- find_anticodon_position(test_seq, "AGC", min_pos = 30, max_pos = 45)
  expect_false(is.null(result_a))
})

test_that("find_anticodon_position returns NULL for NA inputs", {
  expect_null(find_anticodon_position(NA, "AGC"))
  expect_null(find_anticodon_position("ACGT", NA))
  expect_null(find_anticodon_position(NA, NA))
})

test_that("find_anticodon_position handles anticodon not found", {
  test_seq <- "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"

  result <- find_anticodon_position(test_seq, "GGG", min_pos = 30, max_pos = 45)

  expect_null(result)
})

test_that("find_anticodon_position finds anticodon outside range with warning", {
  # Put anticodon at position 10 (outside expected 30-45)
  test_seq <- "AAAAAAAAAGCGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"

  result <- find_anticodon_position(test_seq, "GCG", min_pos = 30, max_pos = 45)

  # Should find it but mark as not validated
  expect_false(is.null(result))
  expect_false(result$validated)
  expect_false(is.null(result$note))
})

test_that("find_anticodon_position selects closest to position 35", {
  # Put multiple GCG occurrences
  test_seq <- paste0(
    "GCGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",  # GCG at 1 (outside range)
    "GCG",                                  # GCG at 34-36 (closest to 35)
    "AAAAGCGAAAAAAAAA"                      # GCG at 39 (also in range)
  )

  result <- find_anticodon_position(test_seq, "GCG", min_pos = 30, max_pos = 45)

  # Should pick the one closest to position 35
  expect_equal(result$start, 34)
})

# =============================================================================
# Tests for add_anticodon_positions()
# =============================================================================

test_that("add_anticodon_positions adds correct columns", {
  trna_df <- data.frame(
    id = c("test1", "test2"),
    sequence = c(
      "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGCAAAAAAAAAA",
      "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGGGAAAAAAAAAA"
    ),
    anticodon = c("AGC", "GGG"),
    stringsAsFactors = FALSE
  )

  result <- add_anticodon_positions(trna_df)

  expect_true("anticodon_start" %in% names(result))
  expect_true("anticodon_end" %in% names(result))
  expect_true("anticodon_validated" %in% names(result))
})

test_that("add_anticodon_positions handles single row", {
  # Test minimal case (empty not supported by mapply)
  single_df <- data.frame(
    id = "test1",
    sequence = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGCAAAAAAAAAA",
    anticodon = "AGC",
    stringsAsFactors = FALSE
  )

  result <- add_anticodon_positions(single_df)

  expect_equal(nrow(result), 1)
  expect_true("anticodon_start" %in% names(result))
  expect_true("anticodon_end" %in% names(result))
  expect_true("anticodon_validated" %in% names(result))
})

# =============================================================================
# Tests for group_by_sequence()
# =============================================================================

test_that("group_by_sequence groups identical sequences", {
  trna_df <- data.frame(
    id = c("a", "b", "c"),
    sequence = c("ACGT", "ACGT", "GGGG"),
    stringsAsFactors = FALSE
  )

  result <- group_by_sequence(trna_df)

  expect_true("sequence_group" %in% names(result))
  expect_true("group_size" %in% names(result))

  # a and b should be in same group
  expect_equal(result$sequence_group[1], result$sequence_group[2])
  # c should be in different group
  expect_false(result$sequence_group[1] == result$sequence_group[3])

  # Group sizes
  expect_equal(result$group_size[1], 2)
  expect_equal(result$group_size[3], 1)
})

test_that("group_by_sequence handles single sequence", {
  trna_df <- data.frame(
    id = "single",
    sequence = "ACGT",
    stringsAsFactors = FALSE
  )

  result <- group_by_sequence(trna_df)

  expect_equal(nrow(result), 1)
  expect_equal(result$group_size[1], 1)
})

test_that("group_by_sequence handles empty data frame", {
  empty_df <- data.frame(
    id = character(),
    sequence = character(),
    stringsAsFactors = FALSE
  )

  result <- group_by_sequence(empty_df)

  expect_equal(nrow(result), 0)
})

test_that("group_by_sequence marks representatives correctly", {
  trna_df <- data.frame(
    id = c("a", "b", "c"),
    sequence = c("ACGT", "ACGT", "ACGT"),
    stringsAsFactors = FALSE
  )

  result <- group_by_sequence(trna_df)

  expect_true("is_representative" %in% names(result))
  expect_equal(sum(result$is_representative), 1)
  expect_true(result$is_representative[1])
})

# =============================================================================
# Tests for count_differences() and find_difference_positions()
# =============================================================================

test_that("count_differences counts correctly", {
  expect_equal(count_differences("ACGT", "ACGT"), 0)
  expect_equal(count_differences("ACGT", "ACGG"), 1)
  expect_equal(count_differences("ACGT", "GGGG"), 3)
  expect_equal(count_differences("AAAA", "TTTT"), 4)
})

test_that("count_differences handles length differences", {
  # Different lengths should count extra as differences
  result <- count_differences("ACGT", "ACGTAA")

  expect_true(result >= 2)  # At least the extra bases
})

test_that("find_difference_positions returns correct positions", {
  positions <- find_difference_positions("ACGT", "ACGG")

  expect_true(length(positions) > 0)
  expect_true(any(grepl("4:", positions)))  # Difference at position 4
})

test_that("find_difference_positions handles identical sequences", {
  positions <- find_difference_positions("ACGT", "ACGT")

  expect_equal(length(positions), 0)
})

# =============================================================================
# Tests for collapse_identical_sequences()
# =============================================================================

test_that("collapse_identical_sequences collapses correctly", {
  group_df <- data.frame(
    id = c("a", "b", "c"),
    sequence = c("ACGT", "ACGT", "GGGG"),
    stringsAsFactors = FALSE
  )

  result <- collapse_identical_sequences(group_df)

  expect_equal(nrow(result$unique_df), 2)  # 2 unique sequences
  expect_equal(length(result$group_info), 2)
})

test_that("collapse_identical_sequences preserves member IDs", {
  group_df <- data.frame(
    id = c("a", "b", "c"),
    sequence = c("ACGT", "ACGT", "GGGG"),
    stringsAsFactors = FALSE
  )

  result <- collapse_identical_sequences(group_df)

  # Find group for ACGT
  acgt_members <- result$group_info[["ACGT"]]
  expect_true("a" %in% acgt_members)
  expect_true("b" %in% acgt_members)
})

# =============================================================================
# Tests for compute_consensus()
# =============================================================================

test_that("compute_consensus returns most common base", {
  sequences <- c("ACGT", "ACGT", "AGGT")

  consensus <- compute_consensus(sequences)

  # Position 2: C appears twice, G once -> should be C
  expect_equal(nchar(consensus), 4)
  expect_equal(substr(consensus, 2, 2), "C")
})

test_that("compute_consensus handles single sequence", {
  consensus <- compute_consensus("ACGT")

  expect_equal(consensus, "ACGT")
})

test_that("compute_consensus pads shorter sequences", {
  sequences <- c("ACGT", "AC")

  consensus <- compute_consensus(sequences)

  expect_equal(nchar(consensus), 4)
})

# =============================================================================
# Tests for format_differences()
# =============================================================================

test_that("format_differences returns empty for no differences", {
  result <- format_differences(0)

  expect_equal(result, "")
})

test_that("format_differences shows positions for few differences", {
  result <- format_differences(2, c("10:A>G", "20:C>T"), seq_length = 100)

  expect_true(grepl("10:A>G", result))
  expect_true(grepl("identity", result))
})

test_that("format_differences shows percentage for many differences", {
  result <- format_differences(10, NULL, seq_length = 100)

  expect_true(grepl("90%", result))  # 90% identity
  expect_true(grepl("10 nt differ", result))
})

# =============================================================================
# Tests for get_base_color()
# =============================================================================

test_that("get_base_color returns valid colors", {
  expect_true(nchar(get_base_color("A")) > 0)
  expect_true(nchar(get_base_color("T")) > 0)
  expect_true(nchar(get_base_color("G")) > 0)
  expect_true(nchar(get_base_color("C")) > 0)
})

test_that("get_base_color handles lowercase", {
  expect_equal(get_base_color("a"), get_base_color("A"))
  expect_equal(get_base_color("t"), get_base_color("T"))
})

test_that("get_base_color treats U same as T", {
  expect_equal(get_base_color("U"), get_base_color("T"))
})

test_that("get_base_color returns gray for unknown", {
  result <- get_base_color("N")

  expect_true(grepl("999", result))  # Gray color
})

# =============================================================================
# Tests for format_trna_id()
# =============================================================================

test_that("format_trna_id removes nuc-tRNA prefix", {
  expect_equal(format_trna_id("nuc-tRNA-Ala-AGC-1-1"), "Ala-AGC-1-1")
})

test_that("format_trna_id keeps mito prefix simplified", {
  expect_equal(format_trna_id("mito-tRNA-Ala-AGC"), "mito-Ala-AGC")
})

test_that("format_trna_id leaves other IDs unchanged", {
  expect_equal(format_trna_id("Ala-AGC-1-1"), "Ala-AGC-1-1")
})

# =============================================================================
# Tests for find_variable_positions()
# =============================================================================

test_that("find_variable_positions finds differences", {
  sequences <- c("ACGT", "AGGT")

  result <- find_variable_positions(sequences)

  expect_true(2 %in% result)  # Position 2 differs
})

test_that("find_variable_positions returns empty for identical", {
  sequences <- c("ACGT", "ACGT", "ACGT")

  result <- find_variable_positions(sequences)

  expect_equal(length(result), 0)
})

test_that("find_variable_positions handles single sequence", {
  result <- find_variable_positions("ACGT")

  expect_equal(length(result), 0)
})

# =============================================================================
# Tests for HTML rendering functions
# =============================================================================

test_that("format_sequence_html returns valid HTML", {
  html <- format_sequence_html("ACGTACGT", 3, 5)

  expect_s3_class(html, "shiny.tag")
})

test_that("format_sequence_html handles missing anticodon positions", {
  html <- format_sequence_html("ACGTACGT", NA, NA)

  expect_s3_class(html, "shiny.tag")
})

test_that("get_visualization_css returns style tag", {
  css <- get_visualization_css()

  expect_s3_class(css, "shiny.tag")
  expect_equal(css$name, "style")
})

test_that("create_terminology_html returns valid structure", {
  html <- create_terminology_html()

  expect_s3_class(html, "shiny.tag")
})
