# global.R
# Load packages and source R modules for TRUENORTH Shiny app

# =============================================================================
# Packages
# =============================================================================

library(shiny)
library(bslib)
library(DT)
library(htmltools)
library(dplyr)
library(Biostrings)

# =============================================================================
# Source R modules (relative to app directory)
# =============================================================================

source("../R/sequence_utils.R")
source("../R/similarity.R")
source("../R/probe_design.R")
source("../R/validation.R")
source("../R/target_selection.R")
source("../R/visualization.R")

# =============================================================================
# Source Shiny modules - Wizard-based UI
# =============================================================================

source("modules/help_content.R")  # Help content and tooltips
source("modules/wizard_step1.R")
source("modules/wizard_step2.R")
source("modules/wizard_step3.R")
source("modules/wizard_step4.R")

# Old modules (archived)
# source("modules/browse_module.R")
# source("modules/design_module.R")
# source("modules/results_module.R")

# =============================================================================
# Pre-load data for faster startup
# =============================================================================

# Available organisms
ORGANISMS <- c(
  "Human" = "human",
  "Yeast" = "yeast",
  "E. coli" = "ecoli"
)

# Amino acid display names
AMINO_ACID_NAMES <- c(
  Ala = "Alanine", Arg = "Arginine", Asn = "Asparagine", Asp = "Aspartic acid",
  Cys = "Cysteine", Gln = "Glutamine", Glu = "Glutamic acid", Gly = "Glycine",
  His = "Histidine", Ile = "Isoleucine", Leu = "Leucine", Lys = "Lysine",
  Met = "Methionine", Phe = "Phenylalanine", Pro = "Proline", Ser = "Serine",
  Thr = "Threonine", Trp = "Tryptophan", Tyr = "Tyrosine", Val = "Valine",
  iMet = "Initiator Met", fMet = "Formyl-Met", Ile2 = "Isoleucine-2",
  SeC = "Selenocysteine", Sup = "Suppressor", Undet = "Undetermined"
)

# Data directory (relative to app)
DATA_DIR <- "../data"

# =============================================================================
# Helper functions
# =============================================================================

#' Load organism data with anticodon positions
#'
#' @param organism One of "human", "yeast", "ecoli"
#' @return Data frame with tRNA data and anticodon positions
load_trna_data <- function(organism) {
  trna_df <- load_organism_data(organism, data_dir = file.path(DATA_DIR, "fastas"))
  trna_df <- add_anticodon_positions(trna_df)
  trna_df
}

#' Load similarity data for an organism
#'
#' @param organism One of "human", "yeast", "ecoli"
#' @return Similarity data list
load_sim_data <- function(organism) {
  load_similarity_data(organism, data_dir = DATA_DIR)
}

#' Format amino acid for display
#'
#' @param aa Three-letter amino acid code
#' @return Display string like "Alanine (Ala)"
format_amino_acid <- function(aa) {
  full_name <- AMINO_ACID_NAMES[aa]
  if (is.na(full_name)) full_name <- aa

  paste0(full_name, " (", aa, ")")
}
