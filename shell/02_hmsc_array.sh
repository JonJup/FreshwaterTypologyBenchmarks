#!/bin/bash
#SBATCH --job-name=hmsc_*taxon*
#SBATCH --output=logs/hmsc_%A_%a.log
#SBATCH --error=logs/hmsc_%A_%a.err
#SBATCH --partition=short
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=24            
#SBATCH --mem=128G                   
#SBATCH --array=1-515

# ============================================
# HMSC Job Array - Two chains per job (Parallel)
# ============================================

# --- Configuration ---
# Taxon group used to build input/output file names.
# Change this to switch e.g. to "macrophytes", "fish", etc.
TAXON="diatoms"

echo "=========================================="
echo "HMSC-HPC Array Job"
echo "=========================================="
echo "Job ID:     $SLURM_JOB_ID"
echo "Array ID:   $SLURM_ARRAY_TASK_ID"
echo "Hostname:   $(hostname)"
echo "CPUs Total: $SLURM_CPUS_PER_TASK"
echo "Taxon:        $TAXON"
echo ""

# Force CPU-only execution
export CUDA_VISIBLE_DEVICES=""

# --- CRITICAL CHANGE FOR PARALLEL CHAINS ---
# You are running 2 chains at once on 24 cores.
# Assign 12 threads per chain so they don't overlap.
THREADS_PER_CHAIN=12

export TF_CPP_MIN_LOG_LEVEL=2
export TF_NUM_INTRAOP_THREADS=$THREADS_PER_CHAIN
export TF_NUM_INTEROP_THREADS=2
export OMP_NUM_THREADS=$THREADS_PER_CHAIN
export MKL_NUM_THREADS=$THREADS_PER_CHAIN
export OPENBLAS_NUM_THREADS=$THREADS_PER_CHAIN
export NUMEXPR_NUM_THREADS=$THREADS_PER_CHAIN

MODEL_ID=$(printf "%04d" $SLURM_ARRAY_TASK_ID)

echo "Array index: $SLURM_ARRAY_TASK_ID"
echo "Model ID:    $MODEL_ID"
echo ""

# Construct file paths
INPUT_FILE="data/002_initialized_hmsc_models/${TAXON}_${MODEL_ID}.rds"

# Define separate output files for each chain
OUTPUT_FILE_C0="data/003_fitted_hmsc_models/${TAXON}_${MODEL_ID}_c0.rds"
OUTPUT_FILE_C1="data/003_fitted_hmsc_models/${TAXON}_${MODEL_ID}_c1.rds"

# Check if input file exists
if [ ! -f "$INPUT_FILE" ]; then
    echo "ERROR: Input file not found: $INPUT_FILE"
    exit 1
fi

# Check if BOTH outputs already exist
if [ -f "$OUTPUT_FILE_C0" ] && [ -f "$OUTPUT_FILE_C1" ]; then
    echo "SKIPPING: Both chain files already exist."
    exit 0
fi

echo "Processing model: $MODEL_ID"
echo "Input:  $INPUT_FILE"
echo "Output Chain 0: $OUTPUT_FILE_C0"
echo "Output Chain 1: $OUTPUT_FILE_C1"
echo ""

# Run analysis
echo "Starting HMSC analysis (2 Parallel Chains)..."
START_TIME=$(date +%s)

# --- RUN CHAIN 0 IN BACKGROUND ---
# Note the '&' at the end
echo "Launching Chain 0..."
singularity exec ../parent/shell/hmsc.sif python3 -m hmsc.run_gibbs_sampler \
    --input "$INPUT_FILE" \
    --output "$OUTPUT_FILE_C0" \
    --samples 2150 \
    --transient 10000 \
    --thin 150 \
    --fp 32 \
	--verbose 1000 \
    --chains 0 &

PID_C0=$!

# --- RUN CHAIN 1 IN BACKGROUND ---
# Note the '&' at the end
echo "Launching Chain 1..."
singularity exec ../parent/shell/hmsc.sif python3 -m hmsc.run_gibbs_sampler \
    --input "$INPUT_FILE" \
    --output "$OUTPUT_FILE_C1" \
    --samples 2150	 \
    --transient 10000 \
    --thin 150 \
    --fp 32 \
	--verbose 1000 \
    --chains 1 &

PID_C1=$!

# --- WAIT FOR COMPLETION ---
echo "Both chains running. Waiting for completion..."
wait $PID_C0
wait $PID_C1

END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))

echo ""
echo "Duration: ${ELAPSED} seconds ($(echo "scale=1; $ELAPSED/60" | bc) minutes)"
echo "=========================================="
