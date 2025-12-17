#!/usr/bin/env Rscript
# test_hierarchical_selection.R
# Demonstrate hierarchical target selection workflow
#
# Run from project root:
#   Rscript tests/test_hierarchical_selection.R

if (!file.exists("R/sequence_utils.R")) {
  stop("Please run this script from the truenorth project root directory")
}

cat("=== Hierarchical Target Selection Demo ===\n\n")

suppressPackageStartupMessages({
  library(Biostrings)
  library(dplyr)
})

source("R/sequence_utils.R")
source("R/similarity.R")
source("R/probe_design.R")
source("R/validation.R")
source("R/target_selection.R")

# Load E. coli data
cat("Loading E. coli tRNA data...\n")
ecoli <- load_organism_data("ecoli", data_dir = "data/fastas")
sim <- load_similarity_data("ecoli", data_dir = "data")
cat("Loaded", nrow(ecoli), "tRNAs\n\n")

# ============================================================
# Example 1: View the hierarchy
# ============================================================
cat("=" , rep("=", 50), "\n", sep = "")
cat("EXAMPLE 1: View tRNA Hierarchy\n")
cat(rep("=", 51), "\n\n", sep = "")

hierarchy <- build_trna_hierarchy(ecoli)
print(hierarchy)

# ============================================================
# Example 2: Select all Ala tRNAs, avoid all others
# ============================================================
cat("\n", rep("=", 51), "\n", sep = "")
cat("EXAMPLE 2: Target all Ala tRNAs\n")
cat(rep("=", 51), "\n\n", sep = "")

# Select using amino acid filter
ala_trnas <- select_trnas(ecoli, amino_acids = "Ala")
cat("Found", nrow(ala_trnas), "Ala tRNAs:\n")
print(ala_trnas[, c("id", "anticodon", "gene_family")])

# Create target selection
selection_ala <- create_target_selection(ecoli, ala_trnas$id)
print(selection_ala)

# Analyze conservation within Ala tRNAs
cat("\nConservation within desired group:\n")
conservation <- analyze_group_conservation(selection_ala$desired, sim)
cat("  Sequences:", conservation$n_sequences, "\n")
cat("  Unique:", conservation$n_unique_sequences, "\n")
cat("  Identity range:", conservation$min_identity, "-", conservation$max_identity, "%\n")
cat("  Mean identity:", conservation$mean_identity, "%\n")

# Analyze divergence from non-Ala tRNAs
cat("\nDivergence from avoid group:\n")
divergence <- analyze_group_divergence(selection_ala, sim)
cat("  Avoid targets:", divergence$n_avoid, "\n")
cat("  Cross-identity range:", divergence$min_cross_identity, "-", divergence$max_cross_identity, "%\n")
cat("  Specificity challenge:", divergence$specificity_challenge, "\n")

# ============================================================
# Example 3: Target Ala-GGC but NOT Ala-TGC
# ============================================================
cat("\n", rep("=", 51), "\n", sep = "")
cat("EXAMPLE 3: Target Ala-GGC, avoid Ala-TGC\n")
cat(rep("=", 51), "\n\n", sep = "")

# Select Ala-GGC
ala_ggc <- select_trnas(ecoli, amino_acids = "Ala", anticodons = "GGC")
# Select Ala-TGC as avoid
ala_tgc <- select_trnas(ecoli, amino_acids = "Ala", anticodons = "TGC")

cat("Desired (Ala-GGC):", nrow(ala_ggc), "tRNAs\n")
cat("Avoid (Ala-TGC):", nrow(ala_tgc), "tRNAs\n\n")

selection_ggc <- create_target_selection(ecoli, ala_ggc$id, ala_tgc$id)

# Analyze position-by-position divergence
cat("Position-by-position analysis:\n")
pos_analysis <- analyze_region_divergence(selection_ggc)
cat("  Reference:", attr(pos_analysis, "reference_id"), "\n")
cat("  Sequence length:", nrow(pos_analysis), "nt\n\n")

# Show positions with best selectivity
cat("Top 10 most selective positions:\n")
top_pos <- head(pos_analysis[order(-pos_analysis$selectivity_score), ], 10)
print(top_pos)

# Find selective regions
cat("\nFinding selective probe regions (20-25 nt)...\n")
regions <- find_selective_regions(selection_ggc, min_length = 20, max_length = 25)
cat("Found", nrow(regions), "possible regions\n\n")

cat("Top 5 regions by selectivity:\n")
print(head(regions[, c("start", "end", "length", "mean_conservation",
                       "mean_divergence", "selectivity_score", "quality")], 5))

# ============================================================
# Example 4: Design selective probes
# ============================================================
cat("\n", rep("=", 51), "\n", sep = "")
cat("EXAMPLE 4: Design Selective Probes\n")
cat(rep("=", 51), "\n\n", sep = "")

result <- design_probes_selective(selection_ggc, ecoli)

cat("Conservation analysis:\n")
cat("  All", result$conservation_analysis$n_sequences, "Ala-GGC sequences are",
    ifelse(result$conservation_analysis$all_identical, "IDENTICAL", "similar"), "\n")
cat("  Mean identity:", result$conservation_analysis$mean_identity, "%\n\n")

cat("Top 5 probes:\n")
print(result$probes[1:5, c("rank", "probe_sequence", "start", "end",
                           "gc_content", "tm_nn", "desired_conservation",
                           "avoid_divergence", "quality")])

# Validate best probe
cat("\nValidating best probe against ALL E. coli tRNAs:\n")
best_probe <- result$probes$probe_sequence[1]
validation <- validate_probe(best_probe, ala_ggc$id[1], ecoli, sim)

cat("Specificity score:", round(validation$specificity_score, 3), "\n\n")
cat("Off-target summary:\n")
print(validation$summary)

cat("\n=== Demo Complete ===\n")
