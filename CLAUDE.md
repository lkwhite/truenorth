# Compass - Northern Probe Design Tool

## Project Overview

Compass is a standalone GUI application for molecular biologists to design northern blot probes, starting with the use case of tRNA probing.

**Target Users:** Molecular biologists designing northern probes
**Initial Scope:** tRNA probe design for human, yeast (*S. cerevisiae*), and *E. coli*
**Tech Stack:** R Shiny application with web hosting capability

## Architecture Decisions

- **GUI Framework:** R Shiny (enables rapid prototyping and web deployment)
- **Backend:** R for sequence analysis and probe design algorithms
- **Deployment:** Shinyapps.io or self-hosted Shiny Server

## Project Structure

```
compass/
├── CLAUDE.md                    # This file - AI coding instructions
├── PROJECT_TRACKER.md           # Current tasks and progress
├── README.md                    # User-facing documentation
├── app/                         # Shiny application
│   ├── ui.R                     # User interface
│   ├── server.R                 # Server logic
│   ├── global.R                 # Global variables and setup
│   └── modules/                 # Shiny modules
├── R/                           # Core R functions (non-Shiny)
│   ├── probe_design.R           # Probe design algorithms
│   ├── sequence_utils.R         # Sequence manipulation
│   ├── thermodynamics.R         # Tm calculations, etc.
│   └── validation.R             # Input validation
├── data/                        # Reference data
│   └── fastas/                  # Reference tRNA sequences
│       ├── human/
│       ├── yeast/
│       └── ecoli/
├── tests/                       # Unit tests
│   └── testthat/
├── www/                         # Static web assets
│   ├── css/
│   └── js/
└── renv/                        # R environment management
```

## File Placement Rules

| File Type | Location | Notes |
|-----------|----------|-------|
| Shiny UI components | `app/ui.R` or `app/modules/` | Use modules for reusable components |
| Shiny server logic | `app/server.R` or `app/modules/` | Keep logic modular |
| Core algorithms | `R/` | Non-Shiny-dependent code |
| Reference sequences | `data/fastas/` | Organized by organism |
| Unit tests | `tests/testthat/` | Match source file names |
| CSS/JavaScript | `www/` | Static assets |

## Coding Standards

### R Code
- Use tidyverse style guide
- Document functions with roxygen2 comments
- Prefer `snake_case` for function and variable names
- Use Shiny modules for complex UI components
- Keep reactive expressions minimal and focused

### Testing
- Write unit tests for all probe design algorithms
- Test edge cases (short sequences, ambiguous bases, etc.)
- Use testthat framework

### Git Commits
- Use conventional commits: `feat:`, `fix:`, `docs:`, `test:`, `refactor:`
- Reference GitHub issues when applicable

## Development Workflow

1. **Before starting work:** Check PROJECT_TRACKER.md for current status
2. **Feature branches:** Use `feature/description` naming
3. **Testing:** Run `testthat::test_dir("tests/testthat")` before commits
4. **Documentation:** Update README.md for user-facing changes

## Key Domain Concepts

### Northern Blot Probes
- Short oligonucleotides (typically 20-25 nt for tRNAs) complementary to target RNA
- Must have appropriate melting temperature (Tm) for hybridization
- Should be specific to target (avoid cross-hybridization)

### tRNA-Specific Considerations
- tRNAs are ~70-90 nucleotides long
- High sequence similarity between isoacceptors
- Modified bases may affect hybridization
- Consider secondary structure when designing probes

### Probe Design Parameters
- **Length:** 20-25 nucleotides (for tRNA probes)
- **Tm:** Target hybridization temperature ± buffer
- **GC content:** Typically 40-60%
- **Specificity:** Check for off-target matches

## Reference Data

Source: `/Users/laurawhite/Projects/tRNAs-in-space/fastas/`

| Organism | File | Notes |
|----------|------|-------|
| Human | `hg38-mito-and-nuclear-tRNAs.fa` | Nuclear and mitochondrial |
| Yeast | `sacCer-mito-and-nuclear-tRNAs.fa` | *S. cerevisiae* |
| E. coli | `ecoliK12MG1655-tRNAs.fa` | K-12 MG1655 strain |

## Privacy & Deployment

- Local development: All processing client-side or on localhost
- Web deployment: Consider data privacy for uploaded sequences
- No external API dependencies for core functionality

## Session Workflow

### Starting a Session
1. Read PROJECT_TRACKER.md for current status
2. Check for any blockers or pending decisions
3. Review recent commits if needed

### Ending a Session
1. Update PROJECT_TRACKER.md with progress
2. Commit all changes with descriptive messages
3. Note any blockers or next steps
