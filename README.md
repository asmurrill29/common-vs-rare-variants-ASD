# Genetic Architecture and Pathway Convergence Across the Allele Frequency Spectrum in ASD
---

## Overview

> Brief Description

## Project Structure

> Use draw folder structure from VSCode

---

## Data Sources

| Dataset | Description | Access |
|--------|-------------|--------|
| Grove et al. 2019 | Common variant GWAS summary statistics | [PGC Download](https://pgc.unc.edu/for-researchers/download-results/) |
| Satterstrom et al. 2020 | Rare-variant WES-derived gene burden statistics | [ASC Download](https://asc.broadinstitute.org/downloads) |
| gnomAD v2.1.1 | LOEUF constraint scores | [gnomAD Browser](https://gnomad.broadinstitute.org/downloads)|
| MSigDB | Curated gene sets for pathway enrichment (GO Biological Process, KEGG) | See dependencies for download options. |

> ⚠️ Raw data files will not be committed to this repository. 

---

## Dependencies

### R
```r
 install.packages("tidyverse")
 install.packages("ggplot2")
 BiocManager::install("fgsea")
```
**MSigDB Download**

Option 1: within R (recommended) 

```r
# install and load msigdbr
install.packages("msigdbr")
library(msigdbr)
```

Option 2: Direct download

Download gene sets manually from [MSigDB](https://www.gsea-msigdb.org/gsea/msigdb/)

> All R session dependencies summarized in "session_info.txt".
### Python
- `subprocess` (standard library, no installation required) 


### Command Line Tools
- [MAGMA v1.10](https://cncr.nl/research/magma/)

For those working on the Boston University SCC:
```bash
module load magma-bio/1.10
```

### Environment
- VSCode with Remote SSH on BU SCC
- R version: 4.5.2
- Python version: 3.6.8

---

## Pipeline

### Step 1 — Data Preparation
```bash
Rscript scripts/magma_prep.R
```

### Step 2 — MAGMA Gene-Level Analysis
```bash
python scripts/run_magma.py
```

### Step 3 — Pathway Enrichment
```bash
Rscript scripts/fgsea.R
```

### Step 4 — LOEUF Annotation
```bash
Rscript scripts/loeuf_annotation.R
```

---

## Aims

1. Perform gene-level association analysis on common-variant GWAS summary statistics using MAGMA.
2. Perform gene-level association analysis on rare-variant WES gene burden scores using MAGMA.
3. Implement pathway enrichment analysis on both gene-level outputs using fgsea, restricting to canonical ASD-relevant gene sets. Cross-reference enriched genes against TADA+ genes and SFARI-curated risk genes.
4. Compare results to identify convergence or divergence across the allele frequency spectrum.
5. Annotate genes contributing to enriched pathways with LOEUF scores from gnomAD to determine whether convergent pathways are disproportionately driven by constrained genes.

---
## Limitations
> Fill in. 
---


## References

1. Grove, J. et al. (2019). Identification of common genetic risk variants for autism spectrum disorder. *Nat Genet* 51, 431–444.
2. Nóbrega, I. d. S. et al. (2024). The Importance of Large-Scale Genomic Studies to Unravel Genetic Risk Factors for Autism. *IJMS* 25(11), 5816.
3. Satterstrom, F. K. et al. (2020). Large-Scale Exome Sequencing Study Implicates Both Developmental and Functional Changes in the Neurobiology of Autism. *Cell* 180(3), 568–584.
4. Trost, B. et al. (2022). Genomic architecture of autism from comprehensive whole-genome sequence annotation. *Cell* 185(23), 4409–4427.

---

## Author

**Amalya S. Murrill**

BS859 Final Project — Spring 2026

Boston University