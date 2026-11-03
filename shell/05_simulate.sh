#!/bin/bash
#SBATCH --job-name=D_SIM
#SBATCH --output=logs/05_sim/sim_%A_%a.log
#SBATCH --error=logs/05_sim/sim_%A_%a.err
#SBATCH -C ukko,c
#SBATCH --time=3:30:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
##SBATCH --array=178,179,180,181,182,183,184,185,187,188,189,190,191,192,193,194,195,196,197,198,199,200,201,202,203,204,205,206,207,208,209,210,211,212,213
##SBATCH --array=157,158,159,160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,186,187,188,189,190,191
##SBATCH --array=1

# ============================================
# Simulation Job Array
# One iteration per job
# ============================================

echo "=========================================="
echo "Array Job"
echo "=========================================="
echo "Job ID:     $SLURM_JOB_ID"
echo "Array ID:   $SLURM_ARRAY_TASK_ID"
echo "Hostname:   $(hostname)"
echo "Date:       $(date)"
echo "CPUs:       $SLURM_CPUS_PER_TASK"
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
INPUT_FILE="data/001_unfitted_hmsc_models/diatoms_${MODEL_ID_PADDED}.rds"
OUTPUT_FILE="data/006_simulated_data/diatoms_${MODEL_ID_PADDED}.rds"

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
singularity exec ../pulse/container/r_v1-4.sif Rscript ../pulse/code/09_simulate_data_v6.R \
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
