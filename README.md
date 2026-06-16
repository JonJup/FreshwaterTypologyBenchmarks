**Grading on a Curve: Generating context-specific Benchmarks for the Biological Evaluation of Freshwater Typology Systems**

<!-- Optional badges -->
<!-- [![DOI](https://zenodo.org/badge/DOI/xx.xxxx/zenodo.xxxxxxx.svg)](https://doi.org/xx.xxxx/zenodo.xxxxxxx) -->
<!-- [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT) -->
<!-- [![Shiny app](https://img.shields.io/badge/Shiny-app-blue)](https://<your-shinyapps-url>) -->

This repository accompanies the manuscript:    

Grading on a Curve: Simulation-Based Benchmarks for the Biological Evaluation of Freshwater Typology Performance
<!-- > **[Full title]** (202x). [Authors]. *[Journal]*. DOI: [xx.xxxx/xxxxxx]. Preprint: [link]. --> 

It is an annonymized version of the repository which will be linked and made available upon acceptance. 
Annonymiztaion was achieved through https://anonymous.4open.science/dashboard

It contains the code to reproduce the results reported in the paper.
The data can be found [here](https://zenodo.org/records/20701841?preview=1&token=eyJhbGciOiJIUzUxMiJ9.eyJpZCI6ImEyZmUzNDQxLWQ5OTMtNDE4Zi05ZTNiLTEwNWM2YThlYmI2ZiIsImRhdGEiOnt9LCJyYW5kb20iOiJhNzhiNjRjMGI1NzM3YmFmZWM2ZDBjYTcxZWZlMTE3NiJ9.FdzRNr2_NgbB6VD6KOOAgEaH1Mp1WeDVWO_6lzwYJtlIiO-CGFY_3d8xvZT0OTiMdgaGdiePCZ4r75XMIegH8A).
The linked Zenodo contains two kinds of data. 
- Raw data which will be made available with extensive meta data upon acceptance of the manuscript
- Intermediate Products from the analysis pipleine provided to ease the review process. The data will not be made broadly available upon acceptance. 

The raw environmental data are openly available and not provided again here, as they contain some large files. All links are provided in the respective R scirpts (see code/prepare_environment).

We provide the following intermediate products:
- biological data (i.e., the raw data) with environmental variables from respective catchments (biota/*taxon*_w_environment.rds).
- The EU HydroDEM catchments with an assigned Environmental Zone (eu_hydro_dem_w_enz.rds).
- 

Please note that after the creation of schemes in [04_define_schemes.R](code/04_define_schemes.R) all steps were run seperately for each taxonomic group. Many scripts thus assume a prallel folder structure where data/ lives inside a folder called e.g., fish/. In the code these are marked with placeholders: e.g., fish_folder.

---

## Overview

We have developed **pan-European ecological benchmarks** for freshwater typology systems using a large database of diatoms, fish, invertebrates, and macrophytes, [joint species distribution models](https://github.com/hmsc-r/hmsc-hpc/tree/main) fitted independently to subsets of this database, a novel [modification algorithm](code/09_simulate_data.R) for environments  

The pipeline is deliberately modular so that each component — data preparation, HMSC fitting, QRF benchmarking, and evaluation — can be re-run independently against updated inputs.

## Scope of the paper

- **Taxonomic groups.** Diatoms, fish, macroinvertebrates, macrophytes.
- **Geographic scope.** Continental Europe.
- **Dataset.** ~400,000 biological samples harmonized from ~90 national and regional datasets, paired with catchment-derived environmental descriptors.

## Repository structure

```
.
├── R/                    # Reusable functions sourced by the pipeline
├── code/                 # Top-level, numbered scripts (one per pipeline stage)
│   ├── 01_add_eu_hydro_to_biota.R
│   ├── 02_combine_env.R
│   ├── 03_add_env_biota.R
│   ├── 04_define_schemes.R
│   ├── ...
│   ├── results/
├── docs/
│   ├── 01_preparation.md
│   ├── 02_hmsc.md
│   ├── 03_simulation_evaluation_qrf.md
│   ├── 04_results_to_paper.md
├── hpc/                  # SLURM / Singularity definition files
└── README.md
```

## Data

Raw biological and environmental data are **not redistributed** in this repository 
but are available on [Zenodo](https://zenodo.org/records/20701841?preview=1&token=eyJhbGciOiJIUzUxMiJ9.eyJpZCI6ImVmMTZmMjg4LThlNzEtNDlmZS1hNmUwLTBjN2RkY2NjNTU2MiIsImRhdGEiOnt9LCJyYW5kb20iOiJlMTU0NTI0MDM1Y2MzNTQzMTFhMTNiZTAzYjgyMmIxNSJ9.ESQO7owKnBSoYGuFdBs73TONFFTaGuvP5BhJ--YTgzEssIlaYLgx25gINDf__tttRC9jSPkUTXj7Q06CAaEPNQ)
Data from different processing steps will be shared upon request.

## Reproducing the analyses

The full pipeline can be reproduced by running the numbered scripts in order.
Dependencies between scripts and files is shown in the 'docs/' folder.  

Parts of the analysis were run on HPC at ANNONYMIZED FOR PEER REVIEW. 
These used singularity containers. The code to compile the container is 
available under 'hpc/'

### HPC workflow

HPC jobs were managed via SLURM and the shell scripts used to submit jobs are 
available under 'hpc/'.

## Shiny application

The Shiny app provides an interactive view of the benchmarks:

- QRF prediction and density visualisation for a user-supplied site or catchment
- Batch prediction from uploaded CSV/Excel tables
- KNN-based imputation of missing environmental predictors

A hosted version is available **[HERE](https://pulse-shiny2-typecoherenceprediction.2.rahtiapp.fi/)**.

## Funding and acknowledgments

This work was supported by ANNONYMIZED FOR PEER REVIEW. 

<!-- the **DFG Walter Benjamin Fellowship** (grant no. 557888845). We thank the data providers listed in the manuscript supplement for sharing harmonised biological and environmental records. --> 
<!-- We also wish to thank the Finnish Computing Competence Infrastructure (FCCI) for supporting this project with computational and data storage resources. --> 


## Contact

Questions, bug reports, and collaboration enquiries are welcome via GitHub issues, or by email to ANNONYMIZED FOR PEER REVIEW.


- **Code:** GPL-3.0
- **Derived data and figures:** CC BY 4.0, unless otherwise stated in the relevant Zenodo record.
- Source datasets retain their original licences; contact the listed providers for redistribution terms.
