# Paralog dependency code

Code accompanying Systematic analysis of genetic vulnerabilities caused by chromosome loss events in cancer.

Analysis of aneuploidy-associated genetic vulnerabilities in cancer cell lines, including paralog insufficiency, CYCLOPS genes, and interaction-network dependencies.

## Repository structure

```
Codes/      analysis scripts, one or more per figure
Data/       input files (see Data section)
Outputs/    generated tables and figure panels
```

## Scripts by figure

| Figure | Script |
|---|---|
| 1, S1 | `chr8p_dep.R`, `chrY_loss_dep.R` |
| 2 | `significant_37_paralog_pairs_characteristic_analysis.r`, `paralog_pair_protein_colocalisation.r` |
| S2 | `sgRNA_identity_to_partner_paralog.R`, `VAF_summary.R`, `VAF_summay_ontarget.R` |
| 3 | `ATAC_seq_TSS_new_combo_plot.r`, `ATACseq_TSS_stabilized_log2FC_analysis.R`, `HiChIP_H3K27ac_peakcalling_allsample.r` |
| S3 | `RNAi_dependency.R` |
| 4 | `CYCLOPS_singleton.r`, `sig_72_interactor_tendency.r`, `GSEA_analysis_CYCLOPS.R`, `GSEA_analysis_interaction.R` |
| 5, 6, S4-S6 | [TO ADD] |

## Data

Genome build: hg38.

## Environment

Package versions listed in R_package_versions.csv

## Running

Scripts are standalone and can be run in any order. Each reads from `Data/` and writes to `Outputs/`.

Random seeds are set where applicable. Permutation and resampling results are reproducible given the same seed.

## Contact

[NAME, EMAIL]

## Citation

[CITATION]# Paralog_Dependency
