#!/bin/bash

#SBATCH --job-name=orbweaver_gpu
#SBATCH --output=orbweaver_gpu_%j.out
#SBATCH --error=orbweaver_gpu_%j.err
#SBATCH --partition=day-long-gpu
#SBATCH --gres=gpu:1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1   # Typically 1 task, then use cpus-per-task for multi-threading
#SBATCH --cpus-per-task=4     # Number of CPU cores for Orbweaver's --threads
#SBATCH --mem=32G             # Memory request (e.g., 32GB)
#SBATCH --time=1-00:00:00     # Max walltime D-HH:MM:SS (1 day)

# --- Configuration ---
# Path to the Orbweaver root directory
# This assumes the script is located in a 'scripts' subdirectory of the project root.
SCRIPT_DIR_HPC="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
ROOT_DIR_HPC="$(dirname "$SCRIPT_DIR_HPC")"

# Orbweaver binary (adjust if your structure is different)
ORBWEAVER_BIN="$ROOT_DIR_HPC/target/release/orbweaver"

# Input and Output (IMPORTANT: Adjust these paths for your specific run)
INPUT_FASTA="$ROOT_DIR_HPC/test_data/prepared_benchmark/real_genome_1MB.fasta" # Example: Use a real path to your large input FASTA
OUTPUT_DIR_HPC="$ROOT_DIR_HPC/hpc_gpu_results_$(date +%Y%m%d_%H%M%S)" # Unique output directory
OUTPUT_JSON="$OUTPUT_DIR_HPC/grammar.json"
OUTPUT_TEXT="$OUTPUT_DIR_HPC/grammar.txt"
# Add other output files if needed (e.g., --output-gfa, --visualize, --export-blocks)

# Orbweaver parameters (adjust as needed for your specific dataset and goals)
MIN_RULE_USAGE=5        # Example: Minimum rule usage
MAX_RULE_COUNT=100000   # Example: Max rule count for large sequences
THREADS=${SLURM_CPUS_PER_TASK:-4} # Use Slurm allocated CPUs, default to 4 if not set

# --- Setup ---
echo "Starting Orbweaver GPU job on $(hostname)"
echo "Job ID: $SLURM_JOB_ID"

# Create a job-specific temporary directory
# Using $SLURM_TMPDIR if available, otherwise creating one in /tmp or a scratch space.
# Adjust JOB_TMPDIR path if your HPC has a preferred scratch location (e.g., /scratch/$USER)
if [ -n "$SLURM_TMPDIR" ]; then
    JOB_TMPDIR="$SLURM_TMPDIR/orbweaver_job_$SLURM_JOB_ID"
else
    # Fallback if $SLURM_TMPDIR is not set. You might want to use /scratch or similar.
    JOB_TMPDIR="/tmp/orbweaver_job_$SLURM_JOB_ID"
fi
mkdir -p "$JOB_TMPDIR"
export TMPDIR="$JOB_TMPDIR"

# Function to clean up temporary directory
cleanup_tmpdir() {
    if [ -d "$JOB_TMPDIR" ]; then
        echo "Cleaning up temporary directory: $JOB_TMPDIR"
        rm -rf "$JOB_TMPDIR"
    fi
}

# Trap EXIT signal to ensure cleanup
trap cleanup_tmpdir EXIT

echo "Temporary directory set to: $TMPDIR"
echo "Submission Directory: $SLURM_SUBMIT_DIR"
echo "Running on partition: $SLURM_JOB_PARTITION"
echo "Allocated GRES: $SLURM_GRES"
echo "Allocated CPUs per task: $SLURM_CPUS_PER_TASK"
echo "Requested Memory: 32G (as per SBATCH directive)"

# Load necessary modules (Uncomment and adjust for your HPC environment)
# echo "Loading modules..."
# module load Rust/1.xx # Or your Rust toolchain module
# module load CUDA/xx.x # Or OpenCL drivers if Orbweaver uses OpenCL and they aren't default
# module load ocl-icd/x.x.x # Example for OpenCL ICD loader
# Ensure your environment provides necessary GPU drivers and libraries (e.g., OpenCL ICD)

# Create output directory
mkdir -p "$OUTPUT_DIR_HPC"
echo "Output will be stored in: $OUTPUT_DIR_HPC"

# Check if Orbweaver binary exists
if [ ! -f "$ORBWEAVER_BIN" ]; then
    echo "Error: Orbweaver binary not found at $ORBWEAVER_BIN" >&2
    echo "Please ensure the path is correct or build the project first (e.g., ./scripts/build.sh)" >&2
    exit 1
fi

# Check if input FASTA exists
if [ ! -f "$INPUT_FASTA" ]; then
    echo "Error: Input FASTA file not found at $INPUT_FASTA" >&2
    echo "Please provide a valid input file path in the script." >&2
    exit 1
fi

# --- Orbweaver Command ---
echo "--- Running Orbweaver ---"
echo "Input FASTA: $INPUT_FASTA"
echo "Output JSON: $OUTPUT_JSON"
echo "Threads: $THREADS"
echo "Orbweaver binary: $ORBWEAVER_BIN"

# Construct the Orbweaver command
# GPU usage is enabled by default in Orbweaver if a compatible GPU is found and OpenCL is set up.
# --no-gpu would disable it. The #SBATCH --gres=gpu:1 directive ensures a GPU is allocated.
# Using --streaming is recommended for large files, which is common on HPC.
CMD="$ORBWEAVER_BIN \
    --input-files "$INPUT_FASTA" \
    -j "$OUTPUT_JSON" \
    --output-text "$OUTPUT_TEXT" \
    --stats \
    --streaming \
    --min-rule-usage $MIN_RULE_USAGE \
    --max-rule-count $MAX_RULE_COUNT \
    --threads $THREADS \
    "
    # Optional flags to consider for large datasets / performance:
    # --chunk-size-streaming <SIZE_IN_BYTES>   # e.g., 1048576 for 1MB streaming chunks
    # --adaptive-chunking                       # If your data benefits from it
    # --max-memory-per-chunk-mb <MB>            # If using adaptive chunking
    # --kmer-size <SIZE>                        # If non-default k-mer size needed for analysis steps
    # --reverse-aware <true|false>              # Default is true
    # --skip-ns <true|false>                    # Default is true
    # --output-repeats "$OUTPUT_DIR_HPC/repeats_summary.tsv" # For a repeat summary

echo "Executing command:"
echo "$CMD"
echo "------------------------"

# Execute the command
# Using eval to handle potential complexities in CMD string, ensure variables are safe.
eval $CMD
EXIT_CODE=$?

# --- Completion ---
echo "------------------------"
if [ $EXIT_CODE -eq 0 ]; then
    echo "✅ Orbweaver HPC Job Completed Successfully"
    echo "Output files are in $OUTPUT_DIR_HPC"
else
    echo "❌ Orbweaver HPC Job FAILED (Exit Code: $EXIT_CODE)"
    echo "Check Slurm output files: orbweaver_gpu_$SLURM_JOB_ID.out and orbweaver_gpu_$SLURM_JOB_ID.err"
fi

echo "Job finished at $(date)"

# Explicitly call cleanup just in case trap doesn't fire as expected under all circumstances
# (though it should). This is more of a safeguard.
cleanup_tmpdir

exit $EXIT_CODE 