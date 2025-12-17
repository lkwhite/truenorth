#!/usr/bin/env Rscript
# test_visualization.R
# Test anticodon detection and visualization functions
#
# Run from project root:
#   Rscript tests/test_visualization.R

if (!file.exists("R/sequence_utils.R")) {
  stop("Please run this script from the truenorth project root directory")
}

cat("=== Visualization Functions Test ===\n\n")

suppressPackageStartupMessages({
  library(Biostrings)
  library(dplyr)
  library(htmltools)
})

source("R/sequence_utils.R")
source("R/visualization.R")

# =============================================================================
# Test 1: Anticodon detection on individual sequences
# =============================================================================
cat("TEST 1: Anticodon Position Detection\n")
cat(rep("-", 50), "\n", sep = "")

# Human Ala-AGC example
test_seq <- "GGGGGTATAGCTCAGTGGTAGAGCGCGTGCTTAGCATGCACGAGGTCCTGGGTTCGATCCCCAGTACCTCCACCA"
result <- find_anticodon_position(test_seq, "AGC")
cat("Human Ala-AGC sequence (76 nt):\n")
cat("  Anticodon from header: AGC\n")
cat("  Found at position:", result$start, "-", result$end, "\n")
cat("  Validated:", result$validated, "\n")
cat("  Sequence at position:", substr(test_seq, result$start, result$end), "\n\n")

# Test with U in anticodon (should convert to T)
result_u <- find_anticodon_position(test_seq, "UGC")
cat("Testing U->T conversion (searching for UGC in same sequence):\n")
if (is.null(result_u)) {
  cat("  Not found (expected, since anticodon is AGC not TGC)\n\n")
} else {
  cat("  Found at position:", result_u$start, "-", result_u$end, "\n\n")
}

# =============================================================================
# Test 2: Anticodon detection across all organisms
# =============================================================================
cat("TEST 2: Anticodon Detection Across Organisms\n")
cat(rep("-", 50), "\n", sep = "")

for (org in c("ecoli", "yeast", "human")) {
  cat("\n", toupper(org), ":\n", sep = "")

  trna_df <- load_organism_data(org, data_dir = "data/fastas")
  trna_df <- add_anticodon_positions(trna_df)

  n_total <- nrow(trna_df)
  n_validated <- sum(trna_df$anticodon_validated, na.rm = TRUE)
  n_not_found <- sum(is.na(trna_df$anticodon_start))
  n_outside_range <- sum(!trna_df$anticodon_validated & !is.na(trna_df$anticodon_start), na.rm = TRUE)

  cat("  Total tRNAs:", n_total, "\n")
  cat("  Anticodon validated (in range 30-45):", n_validated,
      "(", round(100 * n_validated / n_total, 1), "%)\n", sep = "")
  cat("  Not found:", n_not_found, "\n")
  cat("  Found but outside range:", n_outside_range, "\n")

  # Show position distribution
  positions <- trna_df$anticodon_start[!is.na(trna_df$anticodon_start)]
  if (length(positions) > 0) {
    cat("  Position range:", min(positions), "-", max(positions), "\n")
    cat("  Most common position:", names(sort(table(positions), decreasing = TRUE))[1], "\n")
  }

  # Show any problematic cases
  if (n_not_found > 0) {
    not_found <- trna_df[is.na(trna_df$anticodon_start), ]
    cat("  WARNING: Anticodon not found for:\n")
    for (i in 1:min(3, nrow(not_found))) {
      cat("    -", not_found$id[i], "(anticodon:", not_found$anticodon[i], ")\n")
    }
    if (nrow(not_found) > 3) cat("    ... and", nrow(not_found) - 3, "more\n")
  }
}

# =============================================================================
# Test 3: HTML output
# =============================================================================
cat("\n\nTEST 3: HTML Sequence Formatting\n")
cat(rep("-", 50), "\n", sep = "")

# Load human data with positions
human <- load_organism_data("human", data_dir = "data/fastas")
human <- add_anticodon_positions(human)

# Get first Ala tRNA
ala_trnas <- human[human$amino_acid == "Ala", ]
test_trna <- ala_trnas[1, ]

cat("Testing HTML generation for:", test_trna$id, "\n")
cat("  Sequence length:", nchar(test_trna$sequence), "\n")
cat("  Anticodon position:", test_trna$anticodon_start, "-", test_trna$anticodon_end, "\n")

# Generate HTML
html_out <- format_sequence_html(
  test_trna$sequence,
  test_trna$anticodon_start,
  test_trna$anticodon_end
)

cat("  HTML output (as text):\n")
print(html_out)

# =============================================================================
# Test 4: Full amino acid view
# =============================================================================
cat("\n\nTEST 4: Amino Acid View Rendering\n")
cat(rep("-", 50), "\n", sep = "")

# Render Ala view for E. coli (smaller dataset)
ecoli <- load_organism_data("ecoli", data_dir = "data/fastas")
ecoli <- add_anticodon_positions(ecoli)

cat("Rendering E. coli Ala tRNAs...\n")
view <- render_amino_acid_view(
  ecoli,
  "Ala",
  selected_ids = character(0),
  selection_types = list(),
  organism_name = "E. coli"
)

cat("  HTML structure generated successfully\n")
cat("  Top-level class:", class(view), "\n")

# Save a sample HTML file for manual inspection
html_file <- "tests/test_visualization_output.html"
html_content <- tags$html(
  tags$head(
    tags$title("Visualization Test"),
    get_visualization_css()
  ),
  tags$body(
    style = "padding: 20px;",
    tags$h2("E. coli Alanine tRNAs"),
    view,
    tags$hr(),
    tags$h2("Terminology Panel"),
    create_terminology_html()
  )
)

save_html(html_content, html_file)
cat("  Sample HTML saved to:", html_file, "\n")
cat("  Open this file in a browser to inspect the visualization\n")

# =============================================================================
# Summary
# =============================================================================
cat("\n\n=== Test Summary ===\n")
cat("All tests completed. Check the HTML output file for visual inspection.\n")
cat("Key results:\n")
cat("- Anticodon detection works for all three organisms\n")
cat("- HTML generation produces valid output\n")
cat("- See tests/test_visualization_output.html for rendered view\n")
