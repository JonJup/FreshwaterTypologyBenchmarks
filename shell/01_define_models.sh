#!/bin/bash
#SBATCH --job-name=DIA_DEF
#SBATCH --output=logs/01_define/def_%A_%a.log
#SBATCH --error=logs/01_define/def_%A_%a.err
#SBATCH -C ukko,c
#SBATCH --time=0:10:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
##SBATCH --array=1-515
#SBATCH --array=158,159,160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,186,187,188,189,190,191

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
SCRIPT_PATH="../pulse/code/05_build_hmsc_hpc_models.R"
CONTAINER_PATH="../pulse/container/r_v1-4.sif"

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
