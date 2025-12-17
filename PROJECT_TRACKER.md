# Compass Project Tracker

**Last Updated:** 2024-12-16

## Project Status: Core Algorithms Complete

## Current Phase: Phase 1 MVP (Algorithms Done, UI Pending)

Core probe design algorithms implemented. Pre-computed similarity approach for specificity-aware design.

---

## Active Tasks

| Task | Status | Notes |
|------|--------|-------|
| Shiny app skeleton | Not Started | Basic UI layout |
| Visualization layer | Not Started | Heatmaps, alignment views |

---

## Completed Tasks

| Task | Date | Notes |
|------|------|-------|
| Create CLAUDE.md | 2024-12-16 | Project instructions |
| Create PROJECT_TRACKER.md | 2024-12-16 | This file |
| R/sequence_utils.R | 2024-12-16 | FASTA parsing, header parsing, filtering |
| R/similarity.R | 2024-12-16 | Pre-computed pairwise alignment, divergent region detection |
| R/probe_design.R | 2024-12-16 | Probe generation, Tm (NN method), GC calculation |
| R/validation.R | 2024-12-16 | Specificity checking, off-target analysis |
| Unit tests | 2024-12-16 | tests/testthat/ for all modules |
| Similarity data cache | 2024-12-16 | Pre-computed for human, yeast, ecoli |
| R/target_selection.R | 2024-12-16 | Hierarchical target selection, divergence analysis |
| Integration tests | 2024-12-16 | Full workflow verified |

---

## Backlog

### Phase 1: MVP
- [x] Load and parse FASTA sequences
- [x] Basic probe design (user-specified region)
- [x] Tm calculation (nearest-neighbor method)
- [x] Specificity check with pre-computed similarity
- [ ] Display probe sequences with properties (Shiny UI)
- [ ] Export probe list (CSV)

### Phase 2: Enhanced Specificity
- [x] Isoacceptor discrimination analysis
- [x] Visual alignment of probes to targets (text-based)
- [x] Highlight unique/divergent regions
- [x] Hierarchical target selection (desired/avoid groups)
- [x] Position-by-position divergence analysis
- [ ] Cross-organism specificity checking

### Phase 3: Advanced Features
- [ ] Secondary structure consideration
- [ ] Batch probe design for multiple targets
- [ ] Probe optimization suggestions
- [ ] Tm normalization across probe set

### Phase 4: Polish & Deploy
- [ ] User documentation
- [ ] Example workflows
- [ ] Deploy to shinyapps.io
- [ ] Performance optimization

---

## Blockers

*None currently*

---

## Decisions Log

| Date | Decision | Rationale |
|------|----------|-----------|
| 2024-12-16 | Use R Shiny for GUI | Rapid prototyping, easy web deployment, good for bioinformatics |
| 2024-12-16 | Start with tRNA probes | Specific use case with known requirements |
| 2024-12-16 | Probe length 20-25 nt | Appropriate for tRNA northern blots |

---

## Key Files

| File | Purpose |
|------|---------|
| `CLAUDE.md` | AI coding instructions |
| `PROJECT_TRACKER.md` | This file - task tracking |
| `R/sequence_utils.R` | FASTA parsing, header parsing, filtering |
| `R/similarity.R` | Pre-computed similarity, divergent region detection |
| `R/probe_design.R` | Probe generation, Tm/GC calculations, selective design |
| `R/validation.R` | Specificity checking, off-target analysis |
| `R/target_selection.R` | Hierarchical selection, divergence analysis |
| `scripts/generate_similarity_data.R` | Generate cached similarity matrices |
| `tests/test_integration.R` | Quick integration test |
| `tests/test_hierarchical_selection.R` | Hierarchical selection demo |
| `data/similarity/*.rds` | Pre-computed similarity matrices |
| `app/ui.R` | Shiny user interface (to create) |
| `app/server.R` | Shiny server logic (to create) |

---

## Reference Links

- [Shiny Documentation](https://shiny.posit.co/)
- [Biostrings Package](https://bioconductor.org/packages/Biostrings/)
- Nearest-neighbor Tm calculation methods

---

## Session Notes

### 2024-12-16 - Hierarchical Target Selection
- Added `R/target_selection.R` with:
  - `build_trna_hierarchy()`: View AA → anticodon → family → copy tree
  - `select_trnas()`: Filter by include/exclude at any level
  - `create_target_selection()`: Define desired and avoid groups
  - `analyze_group_conservation()`: How similar are desired targets?
  - `analyze_group_divergence()`: How different are desired vs avoid?
  - `analyze_region_divergence()`: Position-by-position analysis
  - `find_selective_regions()`: Find optimal probe regions
- Added `design_probes_selective()` to probe_design.R
- Generated and cached similarity data for all 3 organisms
- All integration tests passing

### 2024-12-16 - Core Algorithms Implementation
- Implemented all 4 R modules:
  - `sequence_utils.R`: FASTA parsing, header parsing for all 3 organisms
  - `similarity.R`: Pre-computed pairwise alignment for specificity-aware design
  - `probe_design.R`: Probe generation with nearest-neighbor Tm calculation
  - `validation.R`: Specificity checking, off-target categorization, alignment reports
- Created unit tests for all modules
- Created integration test script
- Key design: Pre-compute similarity matrices to instantly find divergent regions

### 2024-12-16 - Initial Setup
- Created project structure and documentation
- Set up GitHub repository
- Copied reference FASTA files from tRNAs-in-space project
