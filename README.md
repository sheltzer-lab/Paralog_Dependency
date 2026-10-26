# Paralog dependency code

Code accompanying *Systematic analysis of genetic vulnerabilities caused by chromosome loss events in cancer*.

Analysis of aneuploidy-associated genetic vulnerabilities in cancer cell lines, covering paralog insufficiency, CYCLOPS genes, and interaction-network dependencies.

## Repository contents

```
Codes/                    analysis scripts, one per figure
R_package_versions.csv    package versions from the analysis environment
```

## Scripts by figure

| Figure | Script | Analysis |
|---|---|---|
| 1, S1 | [Figure 1 and Figure S1.r](Codes/Fig1%20%26%20Fig%20S1/Figure%201%20and%20Figure%20S1.r) | Chromosome-arm loss vs neutral dependency screen (DepMap CRISPR); aneuploidy-score filtering |
| 2 | [Figure 2.r](Codes/Fig2/Figure%202.r) | Paralog pair identification and characterisation: identity, family size, interactors, co-localisation, conservation |
| 3 | [Figure 3.r](Codes/Fig3/Figure%203.r) | Dosage compensation (panels B-H): TCGA ATAC-seq and H3K27ac HiChIP TSS profiles, stabilised log2FC and permutation tests |
| 4 | [Figure 4.r](Codes/Fig4/Figure%204.r) | CYCLOPS genes and interaction-network dependencies; GO enrichment |
| 6 | [Hotspot mutations Fixed.qmd](Codes/Fig6/Hotspot%20mutations%20Fixed.qmd) | Benchmarking against published biomarker-associated vulnerabilities |
| S2 | [Figure S2.r](Codes/FigS2/Figure%20S2.r) | CRISPR off-target assessment: sgRNA identity to partner paralog, VAF from targeted sequencing, TMB |
| S3 | [Figure S3.r](Codes/FigS3/Figure%20S3.r) | RNAi replication of top paralog dependencies (DEMETER2) |


## Data

All primary datasets are public and can be downloaded separately. File names in the scripts match the original source file names; dataset versions are given in the Methods.

| Source | Data |
|---|---|
| DepMap 24Q2 | CRISPR gene effect, expression (TPM, log2+1), copy number, model annotations |
| DepMap | DEMETER2 combined RNAi gene dependency scores |
| CCLE | Cell line annotations (`Cell_lines_annotations_20181226.txt`) |
| Published arm-level calls | Aneuploidy scores and chromosome-arm call scores |
| Ensembl (biomaRt) | Gene annotation and human paralog pairs |
| TCGA | ATAC-seq and H3K27ac HiChIP |
| Human Protein Atlas | Subcellular localisation |
| STRING (via `rbioapi`) | Protein-protein interactions |
| CORUM | Protein complex membership |
| UCSC | phyloP conservation (hg38) |

Sequencing data generated in this study are deposited at [ACCESSION]. Significant dependency results are provided as Table S1.

Genome build: hg38.

## Environment

R 4.4.2. Versions of all packages present in the analysis environment are listed in [R_package_versions.csv](R_package_versions.csv).

## Contact

[Yi Bei, yibei@stanford.edu]