#!/bin/bash
#SBATCH --job-name=define_model
#SBATCH --output=logs/def_%A_%a.log
#SBATCH --error=logs/def_%A_%a.err
#SBATCH --partition=short
#SBATCH --time=0:40:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --array=1-515


# ============================================
# Model Definition Job Array
# ============================================

# 1. Define variables
TAXON="diatoms"  

echo "=========================================="
echo "Job ID:       $SLURM_JOB_ID"
echo "Array ID:     $SLURM_ARRAY_TASK_ID"
echo "Taxon:        $TAXON"
echo "Date:         $(date)"
echo ""

# 2. Paths
SCRIPT_PATH="../parent/code/05_build_hmsc_hpc_models.R"
CONTAINER_PATH="../parent/shell/r_v1-4.sif"

# 3. Execution
echo "Starting R analysis for row $SLURM_ARRAY_TASK_ID..."
START_TIME=$(date +%s)

# Pass arguments positionally: 1st=Taxon, 2nd=RowIndex
singularity exec "$CONTAINER_PATH" Rscript "$SCRIPT_PATH" "$TAXON" "$SLURM_ARRAY_TASK_ID"

EXIT_CODE=$?
END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))

if [ $EXIT_CODE -eq 0 ]; then
    echo "SUCCESS: Row $SLURM_ARRAY_TASK_ID completed"
else
    echo "ERROR: Row $SLURM_ARRAY_TASK_ID failed (Exit Code: $EXIT_CODE)"
fi

echo "Duration: ${ELAPSED} seconds"
echo "=========================================="

exit $EXIT_CODE