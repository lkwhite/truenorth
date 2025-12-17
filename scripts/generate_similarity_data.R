#!/usr/bin/env Rscript
# generate_similarity_data.R
# Pre-compute similarity matrices for all organisms
#
# Run from project root:
#   Rscript scripts/generate_similarity_data.R
#
# This may take a few minutes for human (454 tRNAs = ~103k pairs)

# Set working directory to project root if needed
if (!file.exists("R/sequence_utils.R")) {
  stop("Please run this script from the truenorth project root directory")
}

cat("=== Generating Similarity Data for TRUENORTH ===\n\n")

# Load required packages
suppressPackageStartupMessages({
  library(Biostrings)
  library(dplyr)
})

# Source modules
source("R/sequence_utils.R")
source("R/similarity.R")

# Process each organism
organisms <- c("ecoli", "yeast", "human")

for (org in organisms) {
  cat("\n", paste(rep("=", 50), collapse = ""), "\n")
  cat("Processing:", org, "\n")
  cat(paste(rep("=", 50), collapse = ""), "\n\n")

  # Load data
  trna_df <- load_organism_data(org, data_dir = "data/fastas")
  cat("Loaded", nrow(trna_df), "tRNAs\n\n")

  # Compute similarity
  start_time <- Sys.time()
  similarity_data <- compute_similarity_data(trna_df, verbose = TRUE)
  elapsed <- difftime(Sys.time(), start_time, units = "secs")

  cat("\nComputation time:", round(elapsed, 1), "seconds\n")

  # Save
  save_similarity_data(similarity_data, data_dir = "data")

  cat("\n")
}

cat("\n=== Done! ===\n")
cat("Similarity data saved to data/similarity/\n")
cat("Files created:\n")
list.files("data/similarity", full.names = TRUE) |> cat(sep = "\n")
