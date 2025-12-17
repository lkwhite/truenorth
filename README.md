# TRUENORTH

**tRNA-Resolved Utility for Evaluating Northern Oligo Reference Targets via Hybridization**

A Shiny web application for designing northern blot probes targeting tRNAs with high specificity.

![TRUENORTH Screenshot](www/screenshot.png)

## Features

- **Intelligent Probe Design**: Automatically identifies optimal probe sequences that distinguish your target tRNAs from related sequences
- **Coverage Optimization**: For multi-target experiments, ranks probes to maximize coverage with minimal probe sets
- **Three Detection Goals**:
  - **Specific Isodecoders**: Target individual tRNA genes
  - **Isoacceptor Families**: Detect all tRNAs with the same anticodon
  - **Amino Acid Groups**: Detect all tRNAs charging the same amino acid
- **Built-in Specificity Analysis**: Visualizes probe binding across all targets and potential off-targets
- **Modification Awareness**: Warns about probe regions overlapping heavily modified tRNA positions

## Supported Organisms

- Human (*Homo sapiens*) - nuclear and mitochondrial tRNAs
- Yeast (*Saccharomyces cerevisiae*) - nuclear and mitochondrial tRNAs
- *E. coli* K-12 MG1655

## Quick Start

### Run Locally

```r
# Install dependencies
install.packages(c("shiny", "bslib", "DT", "htmltools", "dplyr"))

# Install Bioconductor package
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("Biostrings")

# Clone and run
git clone https://github.com/lkwhite/truenorth.git
cd truenorth
Rscript -e "shiny::runApp('app')"
```

Then open http://127.0.0.1:7482 in your browser.

### Deploy to shinyapps.io

```r
# Install rsconnect
install.packages("rsconnect")

# Configure your shinyapps.io account
rsconnect::setAccountInfo(name='YOUR_ACCOUNT',
                          token='YOUR_TOKEN',
                          secret='YOUR_SECRET')

# Deploy
rsconnect::deployApp('app')
```

## Usage Guide

### 1. Select Your Goal

Choose what you want to detect:
- **Specific tRNAs**: Pick individual isodecoders from a filterable list
- **Isoacceptor**: Select an anticodon family (e.g., all tRNA-Gly-GCC genes)
- **Amino Acid**: Select all tRNAs for an amino acid (e.g., all Glycine tRNAs)

### 2. Review Feasibility

TRUENORTH analyzes your selection and reports:
- Target conservation (how similar your targets are to each other)
- Specificity potential (how different from non-targets)
- Estimated probes needed for full coverage

### 3. Design & Select Probes

- Probes are ranked by coverage optimization (each successive probe covers new targets)
- Click rows to select probes for your set
- Watch the coverage tracker update in real-time
- Export your selected probes with coverage summary

## Project Structure

```
truenorth/
├── app/                    # Shiny application
│   ├── app.R              # Entry point
│   ├── global.R           # Package loading & setup
│   ├── ui.R               # User interface
│   ├── server.R           # Server logic
│   └── modules/           # Shiny modules for wizard steps
├── R/                      # Core algorithms
│   ├── sequence_utils.R   # FASTA parsing
│   ├── similarity.R       # Sequence alignment
│   ├── probe_design.R     # Probe generation & ranking
│   ├── validation.R       # Specificity checking
│   ├── target_selection.R # Target group handling
│   └── visualization.R    # Alignment diagrams
├── data/                   # Reference data
│   ├── fastas/            # tRNA sequences by organism
│   └── similarity/        # Pre-computed similarity matrices
├── tests/                  # Unit tests
└── www/                    # Static assets (CSS)
```

## Probe Design Algorithm

TRUENORTH uses a multi-step algorithm:

1. **Candidate Generation**: Slides a window across aligned target sequences to find conserved regions
2. **Thermodynamic Filtering**: Calculates nearest-neighbor Tm and filters by GC content
3. **Specificity Scoring**: Penalizes probes that bind well to non-target tRNAs
4. **Modification Awareness**: Penalizes regions overlapping anticodon loop and other heavily modified positions
5. **Coverage Ranking**: Re-ranks probes using greedy set cover to maximize target coverage

## Requirements

- R >= 4.0
- Shiny >= 1.7
- Biostrings (Bioconductor)
- bslib, DT, htmltools, dplyr

## Citation

If you use TRUENORTH in your research, please cite:

> TRUENORTH: A web tool for designing specific northern blot probes for tRNA detection. [URL]

## License

MIT License - see [LICENSE](LICENSE) file.

## Contributing

Contributions welcome! Please open an issue or pull request.

## Acknowledgments

tRNA sequence data from:
- GtRNAdb (human, yeast)
- EcoCyc (E. coli)
