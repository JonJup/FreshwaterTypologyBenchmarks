**Pan-European freshwater ecotypologies and ecological benchmarks from data-driven typologies and joint species distribution models**

<!-- Optional badges -->
<!-- [![DOI](https://zenodo.org/badge/DOI/xx.xxxx/zenodo.xxxxxxx.svg)](https://doi.org/xx.xxxx/zenodo.xxxxxxx) -->
<!-- [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT) -->
<!-- [![Shiny app](https://img.shields.io/badge/Shiny-app-blue)](https://<your-shinyapps-url>) -->

This repository accompanies the manuscript:    

Jupke et al (*in preparation*) Grading on a Curve: Simulation-Based Benchmarks for the Biological Evaluation of Freshwater Typology Performance


<!-- > **[Full title]** (202x). [Authors]. *[Journal]*. DOI: [xx.xxxx/xxxxxx]. Preprint: [link]. --> 

It contains the code to reproduce the results reported in the paper. As well as the document itself. 

---

## Overview

We have developed **pan-European ecological benchmarks** for freshwater typology systems using a large [database](MIDFIRELINE) of diatoms, fish, invertebrates, and macrophytes, [joint species distribution models](https://github.com/hmsc-r/hmsc-hpc/tree/main) fitted independently to subsets of this database, a novel [modification algorithm](code/02_models_simulation_evaluation/09_simulate_data.R) for environments  



The pipeline is deliberately modular so that each component — typology construction, HMSC fitting, QRF benchmarking, and evaluation — can be re-run independently against updated inputs.

## Scope of the paper

- **Taxonomic groups.** Diatoms, fish, macroinvertebrates, macrophytes.
- **Geographic scope.** Continental Europe.
- **Dataset.** ~400,000 biological samples harmonised from ~90 national and regional datasets, paired with catchment-derived environmental descriptors.
- **Key outputs.** Ecotypology assignments, HMSC posterior summaries, QRF benchmark surfaces, and a set of reproducible figures and tables.

## Repository structure

```
.
├── R/                    # Reusable functions sourced by the pipeline
├── scripts/              # Top-level, numbered scripts (one per pipeline stage)
│   ├── 01_prepare_data.R
│   ├── 02_build_typology.R
│   ├── 03_fit_hmsc.R
│   ├── 04_simulate_counterfactuals.R
│   ├── 05_fit_qrf_benchmarks.R
│   └── 06_evaluate.R
├── docs/
│   ├── 01_preparation.md
│   ├── Hmsc.md
│   ├── simulation_and_qrf.md
│   ├── results_to_paper.md
├── hpc/                  # SLURM / Singularity templates for HMSC fitting
├── data/
│   ├── raw/              # (Gitignored) source data — see Data section
│   └── processed/        # Harmonised inputs produced by 01_prepare_data.R
├── outputs/              # Fitted models, predictions, figures, tables
└── README.md
```

## Data

Raw biological and environmental data are **not redistributed** in this repository because several source datasets are governed by institutional or national licences. Each script in `scripts/` documents the expected input files and their provenance; the harmonisation code in `01_prepare_data.R` shows exactly how the processed inputs are derived.

Where licensing allows, harmonised inputs and fitted model objects are archived on Zenodo:

- **Processed inputs:** [DOI placeholder]
- **Fitted HMSC models and QRF benchmarks:** [DOI placeholder]

The Shiny app loads fitted models remotely from the Zenodo archive by default, so end users do not need to refit the models to explore the benchmarks.

## Installation

The analyses were run with **R ≥ 4.3**. Package versions are pinned with `renv`.

```r
# Clone the repository, then from the project root:
install.packages("renv")
renv::restore()
```

System-level dependencies:

- GDAL, GEOS, PROJ (for the spatial pipeline)
- A C++17 toolchain (for `Hmsc` and `randomForest`)
- Optional: Singularity/Apptainer for reproducing the HPC runs

## Reproducing the analyses

The full pipeline can be reproduced by running the numbered scripts in order. Intermediate outputs are cached, so individual stages can be re-run without repeating the whole workflow.
Parts of the analysis were run on HPC at the University of Helsinki. These used singularity containers. The code to compile the container is available unter 'hpc/'



Repeat with the configuration file for each taxonomic group.

### HPC workflow

HMSC fitting for all the datasets is computationally intensive. Templates in `hpc/` are provided for SLURM clusters running Singularity. A typical invocation:

```bash
sbatch hpc/fit_hmsc.slurm
```

The template exposes knobs for the number of chains, thinning, transient period, and memory per task, matching the settings reported in the manuscript.

## Shiny application

The PULSE Shiny app provides an interactive view of the benchmarks:

- QRF prediction and density visualisation for a user-supplied site or catchment
- Batch prediction from uploaded CSV/Excel tables
- KNN-based imputation of missing environmental predictors
- Diagnostic overlays against the reference distribution of each ecotype

A hosted version is available at: **[URL placeholder]**.

## Funding and acknowledgments

This work was supported by the **DFG Walter Benjamin Fellowship** (grant no. 557888845). We thank the data providers listed in the manuscript supplement for sharing harmonised biological and environmental records.
We also wish to thank the Finnish Computing Competence Infrastructure (FCCI) for supporting this project with computational and data storage resources.

## License

<<<<<<< HEAD
- Source datasets retain their original licences; see [add MIDFIRE link]

## Contact

Questions, bug reports, and collaboration enquiries are welcome via GitHub issues, or by email to jonjup[at]protonmail.com
=======
- **Code:** [MIT / GPL-3.0 / Apache-2.0 — choose one] — see `LICENSE`.
- **Derived data and figures:** CC BY 4.0, unless otherwise stated in the relevant Zenodo record.
- Source datasets retain their original licences; contact the listed providers for redistribution terms.

## Contact

Questions, bug reports, and collaboration enquiries are welcome via GitHub issues, or by email to **[corresponding author]**.
>>>>>>> b95c47e54e6f9134c6d1ad6cc3dab63772d84931
