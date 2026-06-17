#!/bin/bash
#SBATCH --job-name=EvTyD_TK4
#SBATCH --output=logs/06_EvTy/et_%A_%a.log
#SBATCH --error=logs/06_EvTy/et_%A_%a.err
#SBATCH --partition=long
#SBATCH --time=72:00:00
#SBATCH --cpus-per-task=12
#SBATCH --mem=256G
#SBATCH --array=1-515



# ============================================
# Typology Evaluation Job Array
# ============================================

# --- Configuration ---
# Taxon group used to build input/output file names.
# Change this to switch e.g. to "macrophytes", "fish", etc.
TAXON="diatoms"

echo "=========================================="
echo "Typology Evaluation Array Job"
echo "=========================================="
echo "Job ID:     $SLURM_JOB_ID"
echo "Array ID:   $SLURM_ARRAY_TASK_ID"
echo "Hostname:   $(hostname)"
echo "Date:       $(date)"
echo "CPUs:       $SLURM_CPUS_PER_TASK"
echo "Taxon:      $TAXON"
echo ""

# Force CPU-only execution
export CUDA_VISIBLE_DEVICES=""

# Optimize threading (scale to CPUs allocated)
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK


# Get the model ID from the list (line number = array task ID)
MODEL_ID=$SLURM_ARRAY_TASK_ID

echo "Array index: $SLURM_ARRAY_TASK_ID"
echo "Model ID:    $MODEL_ID"
echo ""

# Construct file paths using the model ID
MODEL_ID_PADDED=$(printf "%04d" $MODEL_ID)
INPUT_FILE="data/006_simulated_data/${TAXON}_${MODEL_ID_PADDED}.rds"
OUTPUT_FILE="data/007_evaluations/${TAXON}_${MODEL_ID_PADDED}.rds"

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

echo "Processing iteration: $ITER_ID"
echo "Input:  $INPUT_FILE"
echo "Output: $OUTPUT_FILE"
echo ""

# Run analysis
echo "Starting R analysis..."
START_TIME=$(date +%s)

# Run R script with singularity, passing the array task ID
singularity exec ../parent/shell/r_v1-4.sif Rscript ../parent/code/10_evaluate_typologies.R \
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
