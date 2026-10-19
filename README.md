# Paired RNA-seq Differential Expression Analysis

## Overview figure

![RNA-seq differential-expression results overview](figures/03_overview.svg)

Differential-expression analysis of paired human CD4 T-cell bulk RNA-seq samples from oral mucosa and blood.

## What this project does

Use a public GEO count matrix to identify genes that differ between **oral mucosa and blood** while accounting for the matched subject design.

The project demonstrates:

- public GEO count-matrix processing
- paired sample design
- low-count filtering
- DESeq2 normalization and Wald testing
- multiple-testing correction
- PCA-based sample structure assessment
- differential-expression reporting
- reproducible analysis with R and Python metadata handling

## Dataset

**NCBI GEO:** GSE116139 — *A human TH17 population with a tissue-resident signature in healthy and inflamed oral mucosal tissues [Bulk RNA-seq]*

- Organism: *Homo sapiens*
- Platform: Illumina HiSeq 2500
- Input: processed bulk RNA-seq count matrix
- Analysis subset: matched CD69− samples
- Samples analyzed: 18
- Subjects: 9
- Comparison: oral mucosa vs blood

The full GEO series contains additional sample types and experimental components. This project deliberately uses the matched CD69− subset so that the comparison remains focused and paired.

## Analysis design

The analysis uses a subject-blocked design: `~ subject + tissue`. This accounts for subject-to-subject variation while testing the tissue effect.

Significance is defined as adjusted p-value < 0.05 and |log2 fold change| ≥ 1. Benjamini–Hochberg correction is used for multiple testing.

## Workflow

![Project 06 RNA-seq workflow](figures/01_workflow.svg)

GEO count matrix → matched CD69− sample selection → subject/tissue metadata → low-count filtering → DESeq2 normalization and Wald test → PCA and differential-expression results.

## Results

The committed summary reports 18 matched samples from 9 subjects, 11,696 genes after filtering, and 519 significant genes (294 upregulated and 225 downregulated).

| Result | Value |
|---|---:|
| Samples | 18 |
| Subjects | 9 |
| Genes after filtering | 11,696 |
| Significant genes | 519 |
| Upregulated | 294 |
| Downregulated | 225 |

[Verified results report](results/REPORT.md) · [Metrics table](results/metrics.tsv)

![Differential-expression evidence](results/figures/03_differential_expression_evidence.svg)

## Output map

![Project 06 output map](figures/02_output_file_map.svg)

The workflow also generates the full differential-expression table, normalized counts, PCA, volcano plot and sample/library mapping. Larger generated outputs are uploaded as a workflow artifact rather than committed to the repository.

## Repository structure

```text
.
├── README.md
├── LICENSE
├── metadata/DATASET.md
├── scripts/run_deseq2.R
├── scripts/resolve_library_map.py
├── figures/
└── results/
```

## Limitations

This project analyzes only the selected CD69− paired subset, does not reproduce every comparison in the original study, and does not perform pathway enrichment or independent biological validation. It uses a processed GEO count matrix rather than starting from raw FASTQ files. Results describe this dataset and are not presented as clinical biomarkers or proof of tissue-specific causality.

## Reproducibility

The analysis records the GEO accession, sample subset, paired statistical design, filtering threshold and outputs. Large raw/generated files are kept outside the repository so the project remains lightweight and easy to inspect.
