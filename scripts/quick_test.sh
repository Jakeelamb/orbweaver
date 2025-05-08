#!/bin/bash
# ======================================================
# Orbweaver Quick Test Script
# ======================================================
# Runs Orbweaver directly for faster iteration and debugging,
# defaulting to streaming and GPU mode.

set -e # Exit on error

# --- Configuration ---
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
ROOT_DIR="$(dirname "$SCRIPT_DIR")"
ORBWEAVER_BIN="$ROOT_DIR/target/release/orbweaver"
INPUT_FASTA="$ROOT_DIR/test_data/prepared_benchmark/real_genome_1MB.fasta" # Use 1MB for now
RESULTS_DIR="$ROOT_DIR/quick_test_results"
OUTPUT_JSON="$RESULTS_DIR/quick_test_grammar.json"
OUTPUT_STDOUT="$RESULTS_DIR/quick_test.stdout"
OUTPUT_STDERR="$RESULTS_DIR/quick_test.stderr"

# Ensure Orbweaver is built
if [ ! -f "$ORBWEAVER_BIN" ]; then
    echo "Orbweaver release binary not found. Building..." >&2
    (cd "$ROOT_DIR" && ./scripts/build.sh)
    if [ ! -f "$ORBWEAVER_BIN" ]; then
        echo "Error: Orbweaver binary still not found after build attempt." >&2
        exit 1
    fi
fi

# Ensure input file exists
if [ ! -f "$INPUT_FASTA" ]; then
    echo "Error: Input FASTA file not found: $INPUT_FASTA" >&2
    echo "Please ensure the benchmark data has been prepared (run benchmark.sh once)." >&2
    exit 1
fi

# Clean previous results and create directory
mkdir -p "$RESULTS_DIR"
rm -f "$OUTPUT_JSON" "$OUTPUT_STDOUT" "$OUTPUT_STDERR"

echo "--- Running Quick Test --- "
echo "Input: $INPUT_FASTA"
echo "Output JSON: $OUTPUT_JSON"
echo "Stdout Log: $OUTPUT_STDOUT"
echo "Stderr Log: $OUTPUT_STDERR"
echo "Command: $ORBWEAVER_BIN -i $INPUT_FASTA -j $OUTPUT_JSON --streaming --stats --min-rule-usage 1000"

# Run Orbweaver directly, redirecting output
$ORBWEAVER_BIN -i "$INPUT_FASTA" -j "$OUTPUT_JSON" --streaming --stats --min-rule-usage 1000 --max-rule-count 10000 > "$OUTPUT_STDOUT" 2> "$OUTPUT_STDERR"
EXIT_CODE=$?

echo "------------------------"
if [ $EXIT_CODE -eq 0 ]; then
    echo "✅ Quick Test Completed Successfully"
echo "Check $OUTPUT_STDOUT and $OUTPUT_STDERR for details."
else
    echo "❌ Quick Test FAILED (Exit Code: $EXIT_CODE)"
echo "Check $OUTPUT_STDOUT and $OUTPUT_STDERR for errors."
fi

exit $EXIT_CODE 