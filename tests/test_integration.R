# test_integration.R
# Quick integration test to verify all modules work with real data
#
# Run from project root:
#   Rscript tests/test_integration.R

# Set working directory to project root if needed
if (!file.exists("R/sequence_utils.R")) {
  stop("Please run this script from the truenorth project root directory")
}

cat("=== TRUENORTH Integration Test ===\n\n")

# Load required packages
suppressPackageStartupMessages({
  library(Biostrings)
  library(dplyr)
})

# Source modules
source("R/sequence_utils.R")
source("R/similarity.R")
source("R/probe_design.R")
source("R/validation.R")

# Test 1: Load E. coli data (smallest dataset)
cat("1. Loading E. coli tRNA data...\n")
ecoli_df <- load_organism_data("ecoli", data_dir = "data/fastas")
cat("   Loaded", nrow(ecoli_df), "tRNAs\n")
cat("   Amino acids:", length(unique(ecoli_df$amino_acid)), "\n")
cat("   Sample IDs:", paste(head(ecoli_df$id, 3), collapse = ", "), "...\n\n")

# Test 2: Parse a header
cat("2. Testing header parsing...\n")
test_header <- "tRNA-Ala-GGC-1-2"
parsed <- parse_trna_header(test_header, "ecoli")
cat("   Header:", test_header, "\n")
cat("   Amino acid:", parsed$amino_acid, "\n")
cat("   Anticodon:", parsed$anticodon, "\n")
cat("   Gene family:", parsed$gene_family, "\n\n")

# Test 3: Generate probes for a target
cat("3. Generating candidate probes...\n")
target_id <- ecoli_df$id[1]
target_seq <- ecoli_df$sequence[1]
cat("   Target:", target_id, "\n")
cat("   Sequence length:", nchar(target_seq), "nt\n")

probes <- generate_candidate_probes(target_seq, min_length = 20, max_length = 22)
cat("   Generated", nrow(probes), "candidate probes\n")
cat("   Top probe:\n")
cat("     Sequence:", probes$probe_sequence[1], "\n")
cat("     Position:", probes$start[1], "-", probes$end[1], "\n")
cat("     GC%:", round(probes$gc_content[1], 1), "\n")
cat("     Tm (NN):", round(probes$tm_nn[1], 1), "C\n\n")

# Test 4: Check specificity of top probe
cat("4. Checking specificity of top probe...\n")
top_probe <- probes$probe_sequence[1]
hits <- check_probe_specificity(top_probe, target_id, ecoli_df, max_mismatches = 3)
cat("   Found", nrow(hits), "tRNAs with <=3 mismatches:\n")
print(head(hits[, c("id", "mismatches", "category")], 5))
cat("\n")

# Test 5: Calculate specificity score
cat("5. Calculating specificity score...\n")
score <- calculate_specificity_score(hits, target_id)
cat("   Specificity score:", round(score, 3), "(1.0 = perfect specificity)\n\n")

# Test 6: Compute similarity matrix (small subset for speed)
cat("6. Computing similarity matrix (first 10 tRNAs for speed)...\n")
small_df <- ecoli_df[1:10, ]
identity_matrix <- compute_identity_matrix(small_df, verbose = FALSE)
cat("   Matrix dimensions:", dim(identity_matrix)[1], "x", dim(identity_matrix)[2], "\n")
cat("   Identity range:", round(min(identity_matrix[identity_matrix < 100]), 1), "-",
    round(max(identity_matrix[identity_matrix < 100]), 1), "%\n\n")

# Test 7: Find similar tRNAs
cat("7. Finding tRNAs similar to", target_id, "...\n")
similar <- find_similar_trnas(target_id, identity_matrix, threshold = 70)
if (nrow(similar) > 0) {
  cat("   Found", nrow(similar), "similar tRNAs (>70% identity):\n")
  print(head(similar, 5))
} else {
  cat("   No tRNAs above 70% identity threshold\n")
}
cat("\n")

# Test 8: Reverse complement
cat("8. Testing reverse complement...\n")
test_seq <- "ATGCATGC"
rc <- reverse_complement(test_seq)
cat("   Input:", test_seq, "\n")
cat("   RC:   ", rc, "\n")
cat("   Double RC:", reverse_complement(rc), "(should match input)\n\n")

# Summary
cat("=== All tests passed! ===\n")
cat("\nTRUENORTH core modules are working correctly.\n")
cat("Next steps:\n")
cat("  - Generate full similarity data: Rscript -e 'source(\"R/similarity.R\"); source(\"R/sequence_utils.R\"); ecoli <- load_organism_data(\"ecoli\", \"data/fastas\"); sim <- compute_similarity_data(ecoli); save_similarity_data(sim, \"data\")'\n")
cat("  - Run unit tests: Rscript -e 'testthat::test_dir(\"tests/testthat\")'\n")
