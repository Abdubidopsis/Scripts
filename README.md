# 🧬 Arabis Genomic and Transcriptomic Scripts

This repository contains a set of scripts for processing genomic and transcriptomic data, specifically tailored for Arabis species research involving QTL mapping, Stacks data processing, and differential expression analysis.

---

## 📄 Scripts Included

| File Name                    | Description |
|-----------------------------|-------------|
| `QTLmapping.R`              | R script for performing QTL analysis and visualizing genotype–phenotype associations. |
| `gstacks_distribs.R`        | R script to summarize and visualize read depth and SNP distributions from Stacks' `gstacks` output. |
| `stacks_raw_process.py`     | Python script for initial raw data processing in Stacks pipeline (e.g., cleaning and demultiplexing). |
| `deseq_drought_transcriptome.py` | Python script for differential expression analysis using DESeq2 via rpy2 or prepared count data. Focuses on drought response in transcriptome datasets. |

---

## 🔧 Requirements

### R
- `tidyverse`
- `qtl`
- `ggplot2`

### Python 3
- `pandas`
- `argparse`
- `os`
- *(optional)* `rpy2` if interfacing with R from Python

---

## ▶️ Usage

Clone the repository:

```bash
git clone https://github.com/Abdubidopsis/Scripts.git
cd Scripts
