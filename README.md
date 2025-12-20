# TRUENORTH

**tRNA-Resolved Utility for Evaluating Northern Oligo Reference Targets via Hybridization**

> **In Development** - Features may change. Validate results experimentally before use.

**Use it here:** https://lkwhite.shinyapps.io/truenorth/

---

## What It Does

TRUENORTH designs northern blot probes for tRNA detection with high specificity. It helps you find probe sequences that detect your target tRNAs while avoiding cross-hybridization.

### Features

- **Goal-oriented workflow** - Target specific isodecoders, isoacceptor families, or all tRNAs for an amino acid
- **Distinguish mode** - Design probes that differentiate between similar isodecoders with cross-hybridization analysis
- **Coverage optimization** - Probes ranked to maximize coverage with minimal probe sets
- **Specificity analysis** - Binding diagrams show hybridization across targets and off-targets
- **Modification warnings** - Flags probes overlapping heavily modified tRNA regions

### Supported Organisms

- Human (nuclear + mitochondrial)
- Yeast (*S. cerevisiae*)
- *E. coli* K-12

---

## How to Use

1. **Select your goal** - What do you want to detect?
2. **Pick targets** - Browse tRNAs by amino acid, anticodon, or gene family
3. **Review feasibility** - See conservation and specificity analysis
4. **Select probes** - Build your probe set with real-time coverage tracking
5. **Export** - Download selected probes as CSV

---

## Terminology

| Term | Definition |
|------|------------|
| **Isoacceptor** | tRNAs sharing the same anticodon |
| **Isodecoder** | Same anticodon, different body sequence |
| **Gene family** | Copies from same ancestral gene |

---

## License

MIT
