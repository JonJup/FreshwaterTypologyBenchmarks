#!/bin/bash
#SBATCH --job-name=varPartDIA
#SBATCH --output=logs/vp_%A_%a.log
#SBATCH --error=logs/vp_%A_%a.err
#SBATCH --partition=short
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=12
#SBATCH --mem=256G
#SBATCH --array=1-503  


# ============================================
# Variation Partitioning Job Array
# ============================================

# --- Configuration ---
# Taxon group used to build input/output file names.
# Change this to switch e.g. to "macrophytes", "fish", etc.
TAXON="diatoms"

echo "=========================================="
echo "Array Job"
echo "=========================================="
echo "Job ID:     $SLURM_JOB_ID"
echo "Array ID:   $SLURM_ARRAY_TASK_ID"
echo "Hostname:   $(hostname)"
echo "Date:       $(date)"
echo "CPUs:       $SLURM_CPUS_PER_TASK"
echo "Taxon:        $TAXON"
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

# Construct file paths using the 0-padded model ID
MODEL_ID_PADDED=$(printf "%04d" $MODEL_ID)
INPUT_FILE="data/001_unfitted_hmsc_models/${TAXON}_${MODEL_ID_PADDED}.rds"
OUTPUT_FILE="data/005_variation_partitioning/${TAXON}_${MODEL_ID_PADDED}.rds"

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

echo "Processing iteration: $MODEL_ID"
echo "Input:  $INPUT_FILE"
echo "Output: $OUTPUT_FILE"
echo ""

# Run analysis
echo "Starting R analysis..."
START_TIME=$(date +%s)

# Run R script with singularity, passing the array task ID
singularity exec ../parent/shell/r_v1-4.sif Rscript ../parent/code/08_variation_partitioning.R \
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

echo "Duration: ${ELAPSED} seconds ($(echo "scale=1; $ELAPSED/60" | bc) minutes)"
echo "=========================================="

exit $EXIT_CODE
