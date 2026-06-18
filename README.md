**Grading on a curve: generating context-specific benchmarks for the biological evaluation of freshwater typology systems**

<!-- Optional badges -->
<!-- [![DOI](https://zenodo.org/badge/DOI/xx.xxxx/zenodo.xxxxxxx.svg)](https://doi.org/xx.xxxx/zenodo.xxxxxxx) -->
<!-- [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT) -->
<!-- [![Shiny app](https://img.shields.io/badge/Shiny-app-blue)](https://<your-shinyapps-url>) -->

This repository accompanies the manuscript:    

Grading on a curve: generating context-specific benchmarks for the biological evaluation of freshwater typology systems
<!-- > **[Full title]** (202x). [Authors]. *[Journal]*. DOI: [xx.xxxx/xxxxxx]. Preprint: [link]. --> 

It is an anonymized version of the repository which will be linked and made available upon acceptance. 
Anonymiztaion was achieved through https://anonymous.4open.science/dashboard

It contains the code to reproduce the results reported in the paper.
The data can be found [here](https://zenodo.org/records/20701841?preview=1&token=eyJhbGciOiJIUzUxMiJ9.eyJpZCI6ImEyZmUzNDQxLWQ5OTMtNDE4Zi05ZTNiLTEwNWM2YThlYmI2ZiIsImRhdGEiOnt9LCJyYW5kb20iOiJhNzhiNjRjMGI1NzM3YmFmZWM2ZDBjYTcxZWZlMTE3NiJ9.FdzRNr2_NgbB6VD6KOOAgEaH1Mp1WeDVWO_6lzwYJtlIiO-CGFY_3d8xvZT0OTiMdgaGdiePCZ4r75XMIegH8A).
The linked Zenodo contains two kinds of data. 
- Raw data which will be made available with extensive meta data upon acceptance of the manuscript
- Intermediate Products from the analysis pipleine provided to ease the review process. The data will not be made broadly available upon acceptance. 

The raw environmental data are openly available and not provided again here, as they contain some large files. All links are provided in the respective R scirpts (see code/prepare_environment).

We provide the following intermediate products:
- biological data (i.e., the raw data) with environmental variables from respective catchments (biota/*taxon*_w_environment.rds).
- The EU HydroDEM catchments with an assigned Environmental Zone (eu_hydro_dem_w_enz.rds).
- For each of the taxon-group specific steps, we provide some example files for the invertebrates. These examples enable reviews to run the code and confirm that it is working. Uploading all data would exceed Zenodo storage limits. 

---

## Overview

We have developed **pan-European ecological benchmarks** for freshwater typology systems using a large database of diatoms, fish, invertebrates, and macrophytes, [joint species distribution models](https://github.com/hmsc-r/hmsc-hpc/tree/main) fitted independently to subsets of this database, a novel [modification algorithm](code/09_simulate_data.R) for environments  

The pipeline is deliberately modular so that each component — data preparation, HMSC fitting, QRF benchmarking, and evaluation — can be re-run independently against updated inputs.
Please note that after the creation of schemes in [04_define_schemes.R](parent/code/04_define_schemes.R) all steps were run seperately for each taxonomic group. Many scripts thus assume a prallel folder structure where data/ lives inside a folder called e.g., fish/. In the code these are marked with placeholders: e.g., fish_folder.

## Scope of the paper

- **Taxonomic groups.** Diatoms, fish, macroinvertebrates, macrophytes.
- **Geographic scope.** Continental Europe.
- **Dataset.** ~400,000 biological samples harmonized from ~90 national and regional datasets, paired with catchment-derived environmental descriptors.

## Repository structure

```
.
├── R/                    # Reusable functions sourced by the pipeline
├── parent/code/                 # Top-level, numbered scripts (one per pipeline stage)
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
├── shell/                  # SLURM / Singularity definition files
└── README.md
```

## Data

Raw biological and environmental data are **not redistributed** in this repository 
but biological data are available on [Zenodo](https://zenodo.org/records/20701841?preview=1&token=eyJhbGciOiJIUzUxMiJ9.eyJpZCI6ImVmMTZmMjg4LThlNzEtNDlmZS1hNmUwLTBjN2RkY2NjNTU2MiIsImRhdGEiOnt9LCJyYW5kb20iOiJlMTU0NTI0MDM1Y2MzNTQzMTFhMTNiZTAzYjgyMmIxNSJ9.ESQO7owKnBSoYGuFdBs73TONFFTaGuvP5BhJ--YTgzEssIlaYLgx25gINDf__tttRC9jSPkUTXj7Q06CAaEPNQ)
Some data from intermediate processing are made available alongside them.

## Reproducing the analyses

The full pipeline can be reproduced by running the numbered scripts in order.
Dependencies between scripts and files is shown in the ['docs/'](docs/) folder.  

Parts of the analysis were run on HPC at ANONYMIZED FOR PEER REVIEW. 

You can start by running the scripts numbered [01](parent/code/01_add_eu_hydro_to_biota), [03](parent/code/03_add_env_biota.R), and [04](parent/code/04_define_schemes.R) ([02](parent/code/02_combine_env.R) requires an external download; see below). All required data is provided in the complimentary [Zenodo repository](https://zenodo.org/records/20701841?preview=1&token=eyJhbGciOiJIUzUxMiJ9.eyJpZCI6ImVmMTZmMjg4LThlNzEtNDlmZS1hNmUwLTBjN2RkY2NjNTU2MiIsImRhdGEiOnt9LCJyYW5kb20iOiJlMTU0NTI0MDM1Y2MzNTQzMTFhMTNiZTAzYjgyMmIxNSJ9.ESQO7owKnBSoYGuFdBs73TONFFTaGuvP5BhJ--YTgzEssIlaYLgx25gINDf__tttRC9jSPkUTXj7Q06CAaEPNQ). Script two requires the download of openly available data which we do not reprovide here. The products of script [02_combine_env.R](parent/code/02_combine_env.R) are provided in the Zenodo repositorty. 

The scripts [05_build_hmsc_hpc_models.R](parent/code/05_build_hmsc_hpc_models.R) to [11_fit_qrf.R](parent/code/11_fit_qrf.R] were run on servers.
Each is associated with a shell script. These can be found in [`shell`](shell/)
We fit the models running the shell script [02_hmsc_array.sh](shell/02_hmsc_array.sh). This script needs to be called for each taxonomic separately and from inside the respective folder (i.e., diatom_folder/ for diatoms). 
On servers we ran two singulariy containers. One for the HMSC-HPC model which is running Tensor Flow in Python, and one for R scripts. Both container (.sif) files are in the Zenodo repository and definitions (.def) are included in this github respositroy ([here](fit_hmsc.def) and [here](r_v1-4.def)).
The final scirpts ([`parent/code/results/`](parent/code/results/)), which compile the results and create figures where run locally. 


## Shiny application

The Shiny app provides an interactive view of the benchmarks:

- QRF prediction and density visualisation for a user-supplied site or catchment
- Batch prediction from uploaded CSV/Excel tables
- KNN-based imputation of missing environmental predictors

A hosted version is available **[HERE](https://pulse-shiny2-typecoherenceprediction.2.rahtiapp.fi/)**.

## Funding and acknowledgments

This work was supported by ANONYMIZED FOR PEER REVIEW. 

<!-- the **DFG Walter Benjamin Fellowship** (grant no. 557888845). We thank the data providers listed in the manuscript supplement for sharing harmonised biological and environmental records. --> 
<!-- We also wish to thank the Finnish Computing Competence Infrastructure (FCCI) for supporting this project with computational and data storage resources. --> 


## Contact

Questions, bug reports, and collaboration enquiries are welcome via GitHub issues, or by email to ANONYMIZED FOR PEER REVIEW.


- **Code:** GPL-3.0
- **Derived data and figures:** CC BY 4.0, unless otherwise stated in the relevant Zenodo record.
- Source datasets retain their original licences; contact the listed providers for redistribution terms.
