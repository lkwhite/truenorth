# Compass Project Tracker

**Last Updated:** 2024-12-16

## Project Status: Initial Setup

## Current Phase: Foundation

Building the core infrastructure and basic probe design functionality.

---

## Active Tasks

| Task | Status | Notes |
|------|--------|-------|
| Project scaffolding | In Progress | Setting up directories, configs |
| Core probe design algorithm | Not Started | Basic Tm calculation, specificity |
| Shiny app skeleton | Not Started | Basic UI layout |

---

## Completed Tasks

| Task | Date | Notes |
|------|------|-------|
| Create CLAUDE.md | 2024-12-16 | Project instructions |
| Create PROJECT_TRACKER.md | 2024-12-16 | This file |

---

## Backlog

### Phase 1: MVP
- [ ] Load and parse FASTA sequences
- [ ] Basic probe design (user-specified region)
- [ ] Tm calculation (nearest-neighbor method)
- [ ] Simple specificity check within organism
- [ ] Display probe sequences with properties
- [ ] Export probe list (CSV)

### Phase 2: Enhanced Specificity
- [ ] Cross-organism specificity checking
- [ ] Isoacceptor discrimination analysis
- [ ] Visual alignment of probes to targets
- [ ] Highlight unique regions

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
| `app/ui.R` | Shiny user interface (to create) |
| `app/server.R` | Shiny server logic (to create) |
| `R/probe_design.R` | Core algorithms (to create) |

---

## Reference Links

- [Shiny Documentation](https://shiny.posit.co/)
- [Biostrings Package](https://bioconductor.org/packages/Biostrings/)
- Nearest-neighbor Tm calculation methods

---

## Session Notes

### 2024-12-16 - Initial Setup
- Created project structure and documentation
- Set up GitHub repository
- Copied reference FASTA files from tRNAs-in-space project
