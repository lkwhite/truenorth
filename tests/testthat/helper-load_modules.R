# helper-load_modules.R
# This file is automatically sourced by testthat before running tests
# It loads all R modules needed for testing

library(Biostrings)
library(dplyr)
library(htmltools)

# Source all R modules from project root
project_root <- normalizePath(file.path(dirname(getwd()), ".."))

source(file.path(project_root, "R/sequence_utils.R"))
source(file.path(project_root, "R/similarity.R"))
source(file.path(project_root, "R/probe_design.R"))
source(file.path(project_root, "R/validation.R"))
source(file.path(project_root, "R/target_selection.R"))
source(file.path(project_root, "R/visualization.R"))
