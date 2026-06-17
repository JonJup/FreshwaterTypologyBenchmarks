#!/bin/bash
#SBATCH --job-name=qrf
#SBATCH --partition=long
#SBATCH --time=15:00:00
#SBATCH --mem=256G
#SBATCH --cpus-per-task=16
#SBATCH --array=1-6
#SBATCH --output=logs/qrf_%a_%A.out
#SBATCH --error=logs/qrf_%a_%A.err


# Define metrics array
METRICS=(
	"classification strength"
	"ANOSIM R mean"
	"AucZeta mean"
	"PERMANOVA R2"
	"PERMANOVA Fuzzy R2"
	"fuzzy_mantel"  
)

# Get metric for this array task
METRIC=${METRICS[$SLURM_ARRAY_TASK_ID-1]}

echo "metric is $METRIC"

# Run R script in container
singularity exec ../parent/shell/r_v1-4.sif Rscript ../parent/code/11_fit_qrf.R \
    --metric "$METRIC" \
    --output "data/008_qrf/"
