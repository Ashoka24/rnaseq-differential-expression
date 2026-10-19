# Dataset — GSE116139

## Source

NCBI Gene Expression Omnibus (GEO): **GSE116139**

The study is a human bulk RNA-seq dataset investigating CD4 T-cell transcriptional profiles in oral mucosa and blood. GEO reports 28 samples in the series and describes paired oral mucosa and blood CD4 T-cell subsets from 10 subjects.

## Selected analysis subset

This project uses the **CD69−** samples that have a matched mucosa and blood sample from the same subject.

Subjects: 19, 34, 40, 42, 43, 47, 49, 53, 56

Total: **18 samples / 9 subjects**

## Sample design

| Subject | Mucosa CD69− | Blood CD69− |
|---|---|---|
| 19 | GSM3211196 | GSM3211194 |
| 34 | GSM3211169 | GSM3211171 |
| 40 | GSM3211174 | GSM3211175 |
| 42 | GSM3211177 | GSM3211178 |
| 43 | GSM3211186 | GSM3211185 |
| 47 | GSM3211180 | GSM3211179 |
| 49 | GSM3211183 | GSM3211182 |
| 53 | GSM3211189 | GSM3211190 |
| 56 | GSM3211193 | GSM3211191 |

## Input

GEO supplementary processed count matrix: `GSE116139_BulkRNAseqCounts.txt.gz`

Raw sequencing reads remain in SRA and are not stored in this repository.

## Statistical model

The paired structure is preserved using `~ subject + tissue`, where subject accounts for between-person variation and tissue is the biological contrast of interest.

## Scope

This repository does not attempt to reproduce the entire original study. It is a focused portfolio analysis built around one well-defined paired comparison. It does not reproduce the full study, single-cell component, or every comparison described by the authors.
