#!/bin/bash

#SBATCH --job-name=orbweaver_gpu
#SBATCH --output=orbweaver_gpu_%j.out
#SBATCH --error=orbweaver_gpu_%j.err
#SBATCH --partition=day-long-gpu
#SBATCH --gres=gpu:1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1   # Typically 1 task, then use cpus-per-task for multi-threading
#SBATCH --cpus-per-task=4     # Number of CPU cores for Orbweaver's --threads
#SBATCH --mem=100G             # Memory request (e.g., 32GB)
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

# Base directory for all Orbweaver outputs
BASE_OUTPUT_DIR_HPC="$ROOT_DIR_HPC/orbweaver_runs" # Master output directory

# Run Identifiers (REQUIRED for new output structure)
SPECIES_ID="homo_sapiens"      # Example: species identifier
ASSEMBLY_ID="hg38_primary"     # Example: assembly identifier
CUSTOM_RUN_ID_HPC=""           # Optional: custom run ID (e.g., "my_analysis_v1"). If empty, Orbweaver defaults to a timestamp.

# Resumption and Rerun Flags
ENABLE_RESUME_HPC="false"                 # Set to "true" to attempt to resume a run
RESUME_TARGET_RUN_DIR_HPC=""            # Path to the specific run directory to resume (e.g., "$BASE_OUTPUT_DIR_HPC/$SPECIES_ID/$ASSEMBLY_ID/20231027_103000_default")
                                        # Required if ENABLE_RESUME_HPC is "true"
FORCE_RERUN_HPC="false"                 # Set to "true" to force rerun even if a run is marked 'completed'

# OLD OUTPUT_DIR_HPC - run-specific directory will be determined by Orbweaver or resume logic
# OUTPUT_DIR_HPC="$ROOT_DIR_HPC/hpc_gpu_results_$(date +%Y%m%d_%H%M%S)" # Unique output directory

# These paths will be relative to the run-specific directory created by Orbweaver
# and are mostly for reference/checking in this script.
# Orbweaver will manage their creation within its structured output.
OUTPUT_JSON_BASENAME="grammar.json" # Basename, actual path determined by Orbweaver
OUTPUT_TEXT_BASENAME="grammar.txt"  # Basename, actual path determined by Orbweaver
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

# Create output directory (Orbweaver will create the run-specific one, ensure base exists)
mkdir -p "$BASE_OUTPUT_DIR_HPC"
echo "Base output directory set to: $BASE_OUTPUT_DIR_HPC"
# The old echo for OUTPUT_DIR_HPC is no longer directly applicable here,
# as the run-specific dir will be managed by Orbweaver or resume logic.
# echo "Output will be stored in: $OUTPUT_DIR_HPC"

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
echo "Species ID: $SPECIES_ID"
echo "Assembly ID: $ASSEMBLY_ID"
[ -n "$CUSTOM_RUN_ID_HPC" ] && echo "Custom Run ID: $CUSTOM_RUN_ID_HPC"
echo "Output base directory: $BASE_OUTPUT_DIR_HPC"
# The exact output JSON/Text paths depend on the run ID generated by Orbweaver or provided.
# We can show the base names for now.
echo "Output JSON (basename): $OUTPUT_JSON_BASENAME"
echo "Output Text (basename): $OUTPUT_TEXT_BASENAME"
echo "Threads: $THREADS"
echo "Orbweaver binary: $ORBWEAVER_BIN"

# Construct the Orbweaver command
# GPU usage is enabled by default in Orbweaver if a compatible GPU is found and OpenCL is set up.
# --no-gpu would disable it. The #SBATCH --gres=gpu:1 directive ensures a GPU is allocated.
# Using --streaming is recommended for large files, which is common on HPC.
CMD="$ORBWEAVER_BIN \
    --input-files \"$INPUT_FASTA\" \
    --output-dir \"$BASE_OUTPUT_DIR_HPC\" \
    --species-id \"$SPECIES_ID\" \
    --assembly-id \"$ASSEMBLY_ID\" \
    --stats \
    --streaming \
    --min-rule-usage $MIN_RULE_USAGE \
    --max-rule-count $MAX_RULE_COUNT \
    --threads $THREADS \
    "

# Add optional custom run ID
if [ -n "$CUSTOM_RUN_ID_HPC" ]; then
    CMD+="--run-id \"$CUSTOM_RUN_ID_HPC\" "
fi

# Add resumption flags if enabled
if [ "$ENABLE_RESUME_HPC" == "true" ]; then
    CMD+="--resume "
    if [ -n "$RESUME_TARGET_RUN_DIR_HPC" ]; then
        CMD+="--resume-run-dir \"$RESUME_TARGET_RUN_DIR_HPC\" "
    else
        echo "Warning: ENABLE_RESUME_HPC is true, but RESUME_TARGET_RUN_DIR_HPC is not set. Resumption might not work as expected without a target directory." >&2
        # Orbweaver will error if --resume is present without --resume-run-dir, this is a script-level warning.
    fi
fi

# Add force rerun flag if enabled
if [ "$FORCE_RERUN_HPC" == "true" ]; then
    CMD+="--force-rerun "
fi

# Add other output file basenames (Orbweaver will place them in the structured directory)
# The explicit -j and --output-text flags for basenames might be redundant if Orbweaver
# defaults them based on the new output structure, but we can keep them for now for clarity
# or if specific basenames are desired different from defaults.
# If Orbweaver uses standard names like grammar.json, grammar.txt by default in the run dir,
# these explicit settings for basenames might not be strictly needed when --output-dir is used.
# For now, assuming they are still respected for the *basename* within the run-specific dir.
CMD+="-j \"$OUTPUT_JSON_BASENAME\" "
CMD+="--output-text \"$OUTPUT_TEXT_BASENAME\" "

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
# Prepending stdbuf -oL -eL to make stdout and stderr line-buffered for live output
eval stdbuf -oL -eL $CMD
EXIT_CODE=$?

# --- Completion ---
echo "------------------------"
if [ $EXIT_CODE -eq 0 ]; then
    echo "✅ Orbweaver HPC Job Completed Successfully"
    # The actual output directory needs to be determined based on whether resuming or not,
    # and whether a custom run ID was used. For now, just point to the base.
    echo "Check output in a subdirectory of $BASE_OUTPUT_DIR_HPC/$SPECIES_ID/$ASSEMBLY_ID/"
else
    echo "❌ Orbweaver HPC Job FAILED (Exit Code: $EXIT_CODE)"
    echo "Check Slurm output files: orbweaver_gpu_$SLURM_JOB_ID.out and orbweaver_gpu_$SLURM_JOB_ID.err"
fi

echo "Job finished at $(date)"

# Explicitly call cleanup just in case trap doesn't fire as expected under all circumstances
# (though it should). This is more of a safeguard.
cleanup_tmpdir

exit $EXIT_CODE 