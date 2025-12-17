# TRUENORTH

**tRNA-Resolved Utility for Evaluating Northern Oligo Reference Targets via Hybridization**

> ⚠️ **In Development** - This tool is under active development. Features may change, and results should be validated experimentally before use. Use at your own risk.

🔗 **Try it: https://lkwhite.shinyapps.io/truenorth/**

## What it does

TRUENORTH designs northern blot probes for tRNA detection with high specificity. It helps you find probe sequences that will detect your target tRNAs while avoiding cross-hybridization with related sequences.

## Features

- **Three detection goals**: Target specific isodecoders, entire anticodon families, or all tRNAs for an amino acid
- **Coverage optimization**: Ranks probes to maximize target coverage with minimal probe sets
- **Specificity analysis**: Shows probe binding across targets and potential off-targets
- **Modification awareness**: Warns about regions overlapping heavily modified tRNA positions

## Supported organisms

- Human (nuclear and mitochondrial)
- Yeast (*S. cerevisiae*)
- *E. coli* K-12

## How to use

1. **Select your goal** - What do you want to detect?
2. **Pick your targets** - Select specific tRNAs, an anticodon, or amino acid
3. **Review feasibility** - See conservation and specificity analysis
4. **Select probes** - Click to select probes; watch coverage update in real-time
5. **Export** - Download your selected probes with coverage summary

## License

MIT
