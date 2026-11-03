#!/bin/bash
#SBATCH --job-name=qrf_dia
#SBATCH -C ukko,c
#SBATCH --time=12:00:00
#SBATCH --mem=256G
#SBATCH --cpus-per-task=16
#SBATCH --array=1
#SBATCH --output=logs/07_qrf/qrf_%a_%A.out
#SBATCH --error=logs/07_qrf/qrf_%a_%A.err


# Define metrics array
METRICS=(
"Classification Strength"
	"ANOSIM R"
	"AucZeta"
	"PERMANOVA R2"
	"PERMANOVA Fuzzy R2"
	"Fuzzy Mantel"  
)

# Get metric for this array task
METRIC=${METRICS[$SLURM_ARRAY_TASK_ID-1]}

echo "metric is $METRIC"

# Run R script in container
singularity exec ../pulse/container/r_v1-4.sif Rscript ../pulse/code/11_fit_qrf.R \
    --metric "$METRIC" \
    --taxon "diatoms" \
    --cv-intervals \
    --output "data/008_qrf/"
