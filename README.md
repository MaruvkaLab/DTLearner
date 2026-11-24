# Predicting tumour doubling time from gene expression

<!-- TODO: replace with the paper title, authors and a link once published -->

Analysis code for **[A machine learning model for predicting tumor doubling time from RNA expression]** ([JOURNAL / preprint], [YEAR]).

This repository contains the R scripts that reproduce every main and supplementary figure and supplementary table of the paper. A 100-gene expression signature associated with doubling time (DT) is identified in cancer cell lines (CCLE), patient-derived xenografts (PDX) and hepatocellular carcinoma (HCC). An ensemble of three models (Ridge regression, Random Forest and radial SVR) trained on PDX models with measured DT is then used to predict doubling time (pDT) in independent patient cohorts (HCC, lower-grade glioma, TCGA and a drug-response cohort) and in CCLE cell lines.

## Repository contents

| Script | What it does | Output folder | Paper items |
|---|---|---|---|
| `plot_figure_1.R` | Expression vs doubling time in CCLE, PDX and HCC; GSEA (GO BP) | `figure_1/` | Fig. 1, Supp. Fig. 1, Supp. Table 2 |
| `plot_figure_2.R` | Cross-validated DT prediction models on PDX and CCLE, CCLE ↔ PDX domain swap | `figure_2/` | Fig. 2, Supp. Fig. 2, Supp. Table 3 |
| `predict_hcc.R` | Train on PDX, predict pDT for HCC (GSE54236) | `hcc_pred/` | input to Fig. 3 |
| `predict_lgg.R` | Train on PDX, predict pDT for LGG | `lgg_pred/` | input to Fig. 3 |
| `predict_hcc_LOOCV.R` | Leave-one-out CV within the HCC cohort only | `loocv_hcc/` | input to Supp. Fig. 3 |
| `predict_tcga.R` | Train on PDX, predict pDT for TCGA tumours | `tcga_pred/` | input to Fig. 4 |
| `predict_drug_cohort.R` | Train on PDX, predict pDT for the patient drug-response cohort | `drug_cohort_pred/` | input to Fig. 5 |
| `predict_ccle.R` | Train on CCLE lines with measured DT, predict the rest | `ccle_pred/` | input to Fig. 5 |
| `plot_figure_3.R` | pDT vs measured DT, MKI67, survival (HCC); subtype and TVDT (LGG) | `figure_3/` | Fig. 3, Supp. Fig. 3, Supp. Table 4 |
| `plot_figure_4.R` | pDT across TCGA: mutation count, stage, survival | `figure_4/` | Fig. 4, Supp. Figs. 4–5, Supp. Table 5 |
| `plot_figure_5.R` | pDT and drug response (patients and DepMap PRISM) | `figure_5/` | Fig. 5, Supp. Fig. 6, Supp. Tables 6–7 |
| `run_all.sh` | Runs all of the above in the correct order | `logs/` | — |

## Requirements

- **R** (version 4.4.0 (2024-04-24) -- "Puppy Cup" was used; see *Session info* below) and `Rscript` on the PATH
- **bash** (Linux, macOS, or Git Bash / WSL on Windows) to use `run_all.sh`

Each script installs any missing packages automatically the first time it runs. The packages used are:

- **CRAN:** dplyr, e1071, future, future.apply, ggplot2, ggpubr, ggraph, ggrepel, glmnet, HGNChelper, igraph, latex2exp, msigdbr, openxlsx, parallelly, patchwork, pls, png, ranger, readxl, scales, stringr, survival, survminer, tibble, tidyr, xgboost, data.table
- **Bioconductor:** clusterProfiler, GOSemSim, GO.db, org.Hs.eg.db

## Input data

The input data are **not included** in this repository. Place them in an `input_files/` folder next to the scripts, with this layout:

```
input_files/
├── hgnc_complete_set.txt
├── ccle/
│   ├── DT_Expression_CCLE_Data.csv
│   ├── Models_DT.csv
│   ├── Repurposing_Public_24Q2_Extended_Primary_Data_Matrix.csv
│   ├── Repurposing_Public_24Q2_Extended_Primary_Compound_List.csv
│   └── PortalCompounds.csv
├── pdx/
│   ├── logTPM_matrix_aver.csv
│   ├── PDM_Samples_com_gc.csv
│   ├── gene_symbol_map_pdx.csv
│   └── growth_curve_points/   (all_markers.csv, gc.png, points.csv)
├── hcc/
│   ├── GSE54236_qnorm.tsv
│   ├── GSE54236_clinical_metadata_compact.tsv
│   └── GSE54236_clinical_meta_compact_EV.xlsx
├── lgg/
│   ├── rna_input_lgg.csv
│   └── TableS1_GRvolume.xlsx
├── tcga/
│   ├── TcgaTargetGtex_rsem_gene_tpm.gz
│   └── TCGA_clinical_data.tsv
└── patient_drug_cohort/
    ├── rna_tab_pat_drag_101_genes.csv
    ├── clinical.tsv
    ├── regimen_classification.csv
    ├── regimens_deduplicated.csv
    ├── regimen_short_names.csv
    └── harmonized/
```

<!-- TODO: fill in where each dataset can be obtained (and any processing applied) -->

| Dataset | Source |
|---|---|
| CCLE expression and doubling times | DepMap portal, [release / link] |
| PRISM drug response | DepMap PRISM Repurposing Public 24Q2 |
| PDX expression and growth curves | [source / link] |
| HCC | GEO accession GSE54236 |
| LGG expression and tumour-volume doubling time | [source / link] |
| TCGA expression | UCSC Xena, TCGA/TARGET/GTEx `TcgaTargetGtex_rsem_gene_tpm` |
| TCGA clinical data | [source / link] |
| Patient drug-response cohort | [source / link] |
| HGNC gene symbols | HGNC complete set, [download date] |

## How to run

Clone or download the repository, add the `input_files/` folder, then from the repository folder run:

```bash
bash run_all.sh
```

Each script runs in its own R session and writes its log to `logs/<script>.log`; a summary is written to `logs/run_all_summary.csv`. Useful options:

```bash
bash run_all.sh --keep-going              # don't stop at the first failure
bash run_all.sh --from plot_figure_3.R    # resume from a given script
RSCRIPT=/path/to/Rscript bash run_all.sh  # if Rscript is not on the PATH
```

The scripts can also be run one at a time (e.g. with **Source** in RStudio), as long as the order in the table above is kept: `plot_figure_2.R` needs `plot_figure_1.R`, and each figure 3–5 script needs its `predict_*.R` scripts. Every script uses the folder it is saved in as its working directory; to use another location, set `WORK_DIR` at the top of the script.

Random seeds are fixed in the scripts, so results should be reproducible on the same package versions.

## Session info

<!-- TODO: paste the output of sessionInfo() from the R session used for the paper -->

## Citation

If you use this code, please cite:

> [AUTHORS]. [TITLE]. [JOURNAL] [YEAR]. doi: [DOI]

## License

[LICENSE — e.g. MIT, see `LICENSE`]

## Contact

[NAME], [LAB], [INSTITUTION] — [EMAIL]
