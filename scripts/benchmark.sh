#!/bin/bash
# ======================================================
# Orbweaver Performance Benchmark Script
# ======================================================
# Runs a series of tests to benchmark Orbweaver's performance
# and memory usage against various input sizes and configurations.
# Each test is limited to a 5-minute timeout.

set -e # Exit on error
# set -x # Uncomment for debug mode

# --- Configuration ---
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
ROOT_DIR="$(dirname "$SCRIPT_DIR")"
ORBWEAVER_BIN="$ROOT_DIR/target/release/orbweaver"
RESULTS_DIR="$ROOT_DIR/benchmark_results"
SUMMARY_FILE="$RESULTS_DIR/summary.txt"
TIMEOUT_DURATION="300s" # 5 minutes

# --- Source Genome Configuration ---
# Define the path to the original source genome
SOURCE_GENOME_PATH="genomes/GCF_002237135.1_genome/ncbi_dataset/data/GCF_002237135.1/GCF_002237135.1_ASM223713v2_genomic.fna"

# Define the directory to store prepared (cleaned) data
PREPARED_DATA_DIR="$ROOT_DIR/test_data/prepared_benchmark"

# Define the desired sizes for testing (in bases)
SIZE_1MB=$((1 * 1024 * 1024))
SIZE_10MB=$((10 * 1024 * 1024))
SIZE_50MB=$((50 * 1024 * 1024))
# SIZE_100MB=$((100 * 1024 * 1024)) # Uncomment if needed

# Check if the source genome exists
if [ ! -f "$SOURCE_GENOME_PATH" ]; then
    echo "Error: Source genome file not found at $SOURCE_GENOME_PATH" >&2
    echo "Please download or place the genome file correctly." >&2
    exit 1
fi

# Ensure Orbweaver is built in release mode
if [ ! -f "$ORBWEAVER_BIN" ]; then
    echo "Orbweaver release binary not found. Building..." >&2
    (cd "$ROOT_DIR" && ./scripts/build.sh)
    if [ ! -f "$ORBWEAVER_BIN" ]; then
        echo "Error: Orbweaver binary still not found after build attempt." >&2
        exit 1
    fi
fi

# Ensure necessary directories exist
mkdir -p "$PREPARED_DATA_DIR"
mkdir -p "$RESULTS_DIR"

# --- Helper Functions ---

# Function to run a single benchmark test
# Usage: run_benchmark <test_name> <orbweaver_args_string> <input_fasta_path> <profile_this_test (true/false)>
run_benchmark() {
    local test_name=$1
    local orbweaver_args_param=$2 # Renamed to avoid conflict
    local input_fasta=$3
    local profile_this_test=${4:-false} # Default to false if not provided

    local result_file_base="$RESULTS_DIR/${test_name}"
    local result_stdout="${result_file_base}.stdout"
    local result_stderr="${result_file_base}.stderr" # time output goes here
    local output_json="${result_file_base}_grammar.json"
    local output_flamegraph="${result_file_base}_flamegraph.svg"

    # Clean previous output files for this test
    rm -f "$result_stdout" "$result_stderr" "$output_json" "$output_flamegraph"
    # Clean up any previous global profile output
    rm -rf "$ROOT_DIR/profile"


    local current_orbweaver_args="$orbweaver_args_param"
    if [ "$profile_this_test" = true ]; then
        echo "INFO: Profiling enabled for $test_name. Flamegraph will be saved to $output_flamegraph"
        # The --profile flag tells orbweaver to generate the flamegraph.
        # We no longer need separate flamegraph command logic here.
        current_orbweaver_args="$current_orbweaver_args --profile"
    fi

    echo "----------------------------------------------------------------------"
    echo "Running Benchmark: $test_name"
    # Get input file size for display
    input_size=$(stat -c%s "$input_fasta" 2>/dev/null || stat -f%z "$input_fasta")
    input_size_human=$(numfmt --to=iec-i --suffix=B $input_size 2>/dev/null || echo "unknown size")
    echo "Input: $(basename "$input_fasta") ($input_size_human)"
    echo "Arguments: $current_orbweaver_args"
    echo "Timeout: $TIMEOUT_DURATION"
    echo "Stdout Log: $result_stdout"
    echo "Stderr (Time) Log: $result_stderr"
    echo "----------------------------------------------------------------------"

    # Command to execute using /usr/bin/time -v and timeout
    # Correct structure: /usr/bin/time [time_opts] timeout [timeout_opts] command [command_opts]
    # Note: We execute the whole timeout command within a subshell `sh -c '...'` to handle arguments and redirections correctly, especially with profiling.
    # The final orbweaver command with its args needs careful quoting.
    local orbweaver_cmd_part="\"$ORBWEAVER_BIN\" $current_orbweaver_args -i \"$input_fasta\" -j \"$output_json\""
    local timeout_cmd_part="timeout \"$TIMEOUT_DURATION\" $orbweaver_cmd_part"
    
    # Use sh -c to properly handle the command string with spaces and quotes
    local cmd_to_run="/usr/bin/time -v -o \"$result_stderr\" sh -c '${timeout_cmd_part} > \"$result_stdout\"'"

    # Execute and capture the exit code of the `time` command itself.
    # `eval` is needed here to handle the nested quoting and redirections within cmd_to_run correctly.
    eval "$cmd_to_run"
    local exit_code=$?


    # Check exit code
    if [ $exit_code -eq 0 ]; then
        echo "✅ $test_name PASSED"
        # Try to parse memory and time from time's output
        local peak_mem=$(grep "Maximum resident set size" "$result_stderr" | awk '{print $6}')
        local user_time=$(grep "User time" "$result_stderr" | awk '{print $4}')
        local sys_time=$(grep "System time" "$result_stderr" | awk '{print $4}')
        echo "$test_name,PASS,$user_time,$sys_time,$peak_mem,$(basename "$input_fasta")" >> "$SUMMARY_FILE"
    elif [ $exit_code -eq 124 ]; then
        echo "⚠️ $test_name TIMED OUT (after $TIMEOUT_DURATION)"
        echo "$test_name,TIMEOUT,-,-,-,$(basename "$input_fasta")" >> "$SUMMARY_FILE"
    else
        echo "❌ $test_name FAILED (Exit Code: $exit_code)"
        echo "$test_name,FAIL,$exit_code,-,-,$(basename "$input_fasta")" >> "$SUMMARY_FILE"
        # Optionally, display some of the stdout/stderr from the failed command
        echo "--- STDOUT ($result_stdout) ---"
        tail -n 10 "$result_stdout"
        echo "--- STDERR (Time) ($result_stderr) ---"
        tail -n 10 "$result_stderr"
    fi

    if [ "$profile_this_test" = true ]; then
        if [ -f "$ROOT_DIR/profile/flamegraph.svg" ]; then
            mkdir -p "$(dirname "$output_flamegraph")" # Ensure results dir exists
            mv "$ROOT_DIR/profile/flamegraph.svg" "$output_flamegraph"
            echo "Flamegraph moved to $output_flamegraph"
            rm -rf "$ROOT_DIR/profile" # Clean up the profile directory
        else
            echo "Warning: Profiling was enabled for $test_name, but no flamegraph was found at $ROOT_DIR/profile/flamegraph.svg"
        fi
    fi
    echo "" # Newline for readability
}

# --- Prepare Test Data from Source Genome ---

echo "--- Preparing Test Data from Source Genome ($SOURCE_GENOME_PATH) ---"

# Prepare cleaned files of different sizes
PREP_SCRIPT="$SCRIPT_DIR/prepare_benchmark_data.sh"

# Define output file paths
FASTA_REAL_1MB="$PREPARED_DATA_DIR/real_genome_1MB.fasta"
FASTA_REAL_10MB="$PREPARED_DATA_DIR/real_genome_10MB.fasta"
FASTA_REAL_50MB="$PREPARED_DATA_DIR/real_genome_50MB.fasta"
FASTA_REAL_FULL="$PREPARED_DATA_DIR/real_genome_full.fasta"

# Run preparation script for each size
if ! bash "$PREP_SCRIPT" "$SOURCE_GENOME_PATH" "$FASTA_REAL_1MB" $SIZE_1MB; then exit 1; fi
if ! bash "$PREP_SCRIPT" "$SOURCE_GENOME_PATH" "$FASTA_REAL_10MB" $SIZE_10MB; then exit 1; fi
if ! bash "$PREP_SCRIPT" "$SOURCE_GENOME_PATH" "$FASTA_REAL_50MB" $SIZE_50MB; then exit 1; fi
if ! bash "$PREP_SCRIPT" "$SOURCE_GENOME_PATH" "$FASTA_REAL_FULL" 0; then exit 1; fi # 0 means full size

# Optional: Prepare 100MB if needed
# FASTA_REAL_100MB="$PREPARED_DATA_DIR/real_genome_100MB.fasta"
# if ! bash "$PREP_SCRIPT" "$SOURCE_GENOME_PATH" "$FASTA_REAL_100MB" $SIZE_100MB; then exit 1; fi

echo "--- Test Data Preparation Complete ---"
echo ""

# --- Clear previous summary ---
echo "TestName,Status,UserTime(s),SystemTime(s),PeakMemory(KB),Input" > "$SUMMARY_FILE"

# --- Benchmark Tests ---

echo "*** Running Section 1: Standard Mode (Full Load) ***"
echo "### Section 1: Standard Mode (Full Load) ###" >> "$SUMMARY_FILE"
run_benchmark "Std_Real_1MB" "" "$FASTA_REAL_1MB" true
run_benchmark "Std_Real_10MB" "" "$FASTA_REAL_10MB" false
# run_benchmark "Std_Real_50MB" "" "$FASTA_REAL_50MB" false # Likely to timeout
# run_benchmark "Std_Real_Full" "" "$FASTA_REAL_FULL" false # Likely to timeout

echo "*** Running Section 2: Streaming Mode ***"
echo "### Section 2: Streaming Mode ###" >> "$SUMMARY_FILE"
run_benchmark "Stream_Real_10MB" "--streaming" "$FASTA_REAL_10MB" false
run_benchmark "Stream_Real_50MB" "--streaming" "$FASTA_REAL_50MB" false
run_benchmark "Stream_Real_Full" "--streaming" "$FASTA_REAL_FULL" false

echo "*** Running Section 3: Chunked Processing (Orbweaver Internal Parallel) ***"
echo "### Section 3: Chunked Processing (Orbweaver Internal) ###" >> "$SUMMARY_FILE"
run_benchmark "Chunked_Real_10MB" "--chunk-size 1048576 --chunk-overlap 10240" "$FASTA_REAL_10MB" false
run_benchmark "Chunked_Real_50MB" "--chunk-size 5242880 --chunk-overlap 51200" "$FASTA_REAL_50MB" false
run_benchmark "Chunked_Real_Full" "--chunk-size 10485760 --chunk-overlap 102400" "$FASTA_REAL_FULL" false

echo "*** Running Section 4: GPU Acceleration vs CPU ***"
echo "### Section 4: GPU vs CPU ###" >> "$SUMMARY_FILE"
# Check for GPU capability (basic check using clinfo)
if command -v clinfo &> /dev/null && clinfo | grep -q "Device Type.*GPU"; then
    echo "GPU detected via clinfo. Running GPU comparison tests."
    # Standard Mode Comparison (Smaller size)
    run_benchmark "GPU_Std_10MB_Real" "--use-gpu" "$FASTA_REAL_10MB" false
    run_benchmark "CPU_Std_10MB_Real" "--no-gpu" "$FASTA_REAL_10MB" false

    # Streaming Mode Comparison (Larger size)
    run_benchmark "GPU_Stream_50MB_Real" "--use-gpu --streaming" "$FASTA_REAL_50MB" false
    run_benchmark "CPU_Stream_50MB_Real" "--streaming --no-gpu" "$FASTA_REAL_50MB" false

    # Chunked Mode Comparison (Largest size)
    run_benchmark "GPU_Chunked_Full_Real_cs10M" "--use-gpu --chunk-size 10000000 --chunk-overlap 100000" "$FASTA_REAL_FULL" false
    run_benchmark "CPU_Chunked_Full_Real_cs10M" "--chunk-size 10000000 --chunk-overlap 100000 --no-gpu" "$FASTA_REAL_FULL" false
else
    echo "No GPU detected or clinfo not available. Skipping GPU comparison tests."
    echo "(Skipped GPU Tests - No compatible GPU detected or clinfo missing)" >> "$SUMMARY_FILE"
fi

echo "*** Running Section 5: Memory Optimization Features ***"
echo "### Section 5: Memory Optimizations ###" >> "$SUMMARY_FILE"
run_benchmark "Stream_Real_50MB_Evict10k" "--streaming --max-rule-count 10000" "$FASTA_REAL_50MB" false
run_benchmark "Stream_Real_Full_Evict50k" "--streaming --max-rule-count 50000" "$FASTA_REAL_FULL" false
# Compare baseline streaming with eviction on memory and time
run_benchmark "Stream_Real_50MB_MinUse5" "--streaming --min-rule-usage 5" "$FASTA_REAL_50MB" false
run_benchmark "Stream_Real_Full_MinUse10" "--streaming --min-rule-usage 10" "$FASTA_REAL_FULL" false

echo "*** Running Section 6: Adaptive Chunking ***"
echo "### Section 6: Adaptive Chunking ###" >> "$SUMMARY_FILE"
run_benchmark "AdaptChunk_Real_50MB" "--adaptive-chunking" "$FASTA_REAL_50MB" false
run_benchmark "AdaptChunk_Real_Full" "--adaptive-chunking" "$FASTA_REAL_FULL" false
# Test with memory constraint
run_benchmark "AdaptChunk_Real_Full_Mem500M" "--adaptive-chunking --max-memory-per-chunk-mb 500" "$FASTA_REAL_FULL" false


# --- End of Tests ---
echo ""
echo "========================================"
echo "Benchmark script finished."
echo "Summary report written to: $SUMMARY_FILE"
echo "Detailed logs are in: $RESULTS_DIR/"
echo "========================================" 