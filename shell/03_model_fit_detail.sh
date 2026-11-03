#!/bin/bash
#SBATCH --job-name=evd_mod_d
#SBATCH --output=logs/03_evl_dtl/run2/ev_mod_d_%A_%a.log
#SBATCH --error=logs/03_evl_dtl/run2/ev_mod_d_%A_%a.err
#SBATCH -C ukko,c
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=6
#SBATCH --mem=64G


# ============================================
# Typology Evaluation Job Array
# ============================================

echo "=========================================="
echo "Array Job"
echo "=========================================="
echo "Job ID:       $SLURM_JOB_ID"
echo "Array ID:     $SLURM_ARRAY_TASK_ID"
echo "Hostname:     $(hostname)"
echo "Date:         $(date)"
echo "CPUs:         $SLURM_CPUS_PER_TASK"
echo ""

# Force CPU-only execution
export CUDA_VISIBLE_DEVICES=""

# Optimize threading (scale to CPUs allocated)
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK

# Get the model ID
MODEL_ID=$SLURM_ARRAY_TASK_ID

echo "Array index: $SLURM_ARRAY_TASK_ID"
echo "Model ID:    $MODEL_ID"
echo ""

# Construct file paths using the actual model ID
MODEL_ID_PADDED=$(printf "%04d" $MODEL_ID)
INPUT_FILE="data/001_unfitted_hmsc_models/diatoms_${MODEL_ID_PADDED}.rds"
OUTPUT_FILE="data/004_model_fit_detail/diatoms_${MODEL_ID_PADDED}.rds"

# Check if input file exists
if [ ! -f "$INPUT_FILE" ]; then
    echo "ERROR: Input file not found: $INPUT_FILE"
    exit 1
fi

# Check if output already exists (skip if already done)
if [ -f "$OUTPUT_FILE" ]; then
    echo "SKIPPING: Output file already exists: $OUTPUT_FILE"
    exit 0
fi

# FIXED: Changed ITER_ID to MODEL_ID
echo "Processing iteration: $MODEL_ID"
echo "Input:  $INPUT_FILE"
echo "Output: $OUTPUT_FILE"
echo ""

# Run analysis
echo "Starting R analysis..."
START_TIME=$(date +%s)

# RECOMMENDATION: Use absolute paths if possible, or ensure you submit from the correct folder
# Run R script with singularity
singularity exec ../pulse/container/r_v1-4.sif Rscript ../pulse/code/07_evaluate_model_fit_detail.R \
    --iter_id $SLURM_ARRAY_TASK_ID \
    --input "$INPUT_FILE" \
    --output "$OUTPUT_FILE"

EXIT_CODE=$?
END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))

if [ $EXIT_CODE -eq 0 ]; then
    echo ""
    echo "SUCCESS: Analysis completed"
else
    echo ""
    echo "ERROR: Analysis failed with exit code $EXIT_CODE"
fi

# Use awk for calculation (bc is not always installed on all minimal containers/nodes)
DURATION_MINS=$(awk "BEGIN {print $ELAPSED/60}")
echo "Duration: ${ELAPSED} seconds (${DURATION_MINS} minutes)"
echo "=========================================="

exit $EXIT_CODE
