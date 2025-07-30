# 🧬 Arabis Genomic and Transcriptomic Scripts

This repository contains a set of scripts and pipelines for genomic and transcriptomic analysis, focused on Arabis species. It includes tools for QTL mapping, RNA-seq analysis, Stacks processing, and genome assembly quality checks.

---

## 📄 Files Included

| File Name                               | Description |
|----------------------------------------|-------------|
| `QTLmapping.R`                          | R script for performing QTL analysis and visualizing genotype–phenotype associations. |
| `gstacks_distribs.R`                   | R script to summarize and visualize read depth and SNP distributions from Stacks `gstacks` output. |
| `stacks_raw_process.py`                | Python script for initial processing of raw reads using the Stacks pipeline. |
| `deseq_drought_transcriptome.py`       | Python script for differential expression analysis using DESeq2 (via count data or R integration). |
| `Genome_assmebly_Quality_Curation_pipeline.txt` | A step-by-step protocol to assess and curate genome assembly quality. Developed using non-model plant genomes. |

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
- *(optional)* `rpy2` if integrating with R

---

## ▶️ Usage

Clone the repository:

```bash
git clone https://github.com/Abdubidopsis/Scripts.git
cd Scripts
