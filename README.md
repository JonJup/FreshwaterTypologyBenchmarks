**Pan-European freshwater ecotypologies and ecological benchmarks from data-driven typologies and joint species distribution models**

<!-- Optional badges -->
<!-- [![DOI](https://zenodo.org/badge/DOI/xx.xxxx/zenodo.xxxxxxx.svg)](https://doi.org/xx.xxxx/zenodo.xxxxxxx) -->
<!-- [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT) -->
<!-- [![Shiny app](https://img.shields.io/badge/Shiny-app-blue)](https://<your-shinyapps-url>) -->

This repository accompanies the manuscript:    

Jupke et al (*in preparation*) Grading on a Curve: Simulation-Based Benchmarks for the Biological Evaluation of Freshwater Typology Performance


<!-- > **[Full title]** (202x). [Authors]. *[Journal]*. DOI: [xx.xxxx/xxxxxx]. Preprint: [link]. --> 

It contains the code to reproduce the results reported in the paper.

---

## Overview

We have developed **pan-European ecological benchmarks** for freshwater typology systems using a large [database](MIDFIRELINE) of diatoms, fish, invertebrates, and macrophytes, [joint species distribution models](https://github.com/hmsc-r/hmsc-hpc/tree/main) fitted independently to subsets of this database, a novel [modification algorithm](code/09_simulate_data.R) for environments  

The pipeline is deliberately modular so that each component — data preparation, HMSC fitting, QRF benchmarking, and evaluation — can be re-run independently against updated inputs.

## Scope of the paper

- **Taxonomic groups.** Diatoms, fish, macroinvertebrates, macrophytes.
- **Geographic scope.** Continental Europe.
- **Dataset.** ~400,000 biological samples harmonised from ~90 national and regional datasets, paired with catchment-derived environmental descriptors.

## Repository structure

```
.
├── R/                    # Reusable functions sourced by the pipeline
├── code/                # Top-level, numbered scripts (one per pipeline stage)
│   ├── 01_prepare_data.R
│   ├── 02_build_typology.R
│   ├── 03_fit_hmsc.R
│   ├── 04_simulate_counterfactuals.R
│   ├── 05_fit_qrf_benchmarks.R
│   └── 06_evaluate.R
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
but are available on [Zenodo](ADD MIDFIRE).

Data from different processing steps will be shared upon request (jonjup@protonmail.com).

The Shiny app loads fitted models remotely from the [Zenodo archive](https://zenodo.org/records/19554225) by default, so end users do not need to refit the models to explore the benchmarks.

## Reproducing the analyses

The full pipeline can be reproduced by running the numbered scripts in order.
Dependencies between scripts and files is shown in the 'docs/' folder.  

Parts of the analysis were run on HPC at the University of Helsinki. 
These used singularity containers. The code to compile the container is 
available under 'hpc/'

### HPC workflow

HPC jobs were managed via SLURM and the shell scripts used to submit jobs are 
available under 'hpc/'.

## Shiny application

The PULSE Shiny app provides an interactive view of the benchmarks:

- QRF prediction and density visualisation for a user-supplied site or catchment
- Batch prediction from uploaded CSV/Excel tables
- KNN-based imputation of missing environmental predictors

A hosted version is available at: **[URL placeholder]**.
The GitHub repository is [here](https://github.com/JonJup/TypologyBenchmarkPrediction).

## Funding and acknowledgments

This work was supported by the **DFG Walter Benjamin Fellowship** (grant no. 557888845). We thank the data providers listed in the manuscript supplement for sharing harmonised biological and environmental records.
We also wish to thank the Finnish Computing Competence Infrastructure (FCCI) for supporting this project with computational and data storage resources.


## Contact

Questions, bug reports, and collaboration enquiries are welcome via GitHub issues, or by email to jonjup[at]protonmail.com
=======
- **Code:** [GPL-3.0 ]
- **Derived data and figures:** CC BY 4.0, unless otherwise stated in the relevant Zenodo record.
- Source datasets retain their original licences; contact the listed providers for redistribution terms.

## Contact

Questions, bug reports, and collaboration enquiries are welcome via GitHub issues, or by email to **[corresponding author]**.
>>>>>>> b95c47e54e6f9134c6d1ad6cc3dab63772d84931
