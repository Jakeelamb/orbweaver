#!/bin/bash
# ======================================================
# Memory Efficiency Test Script for Orbweaver
# ======================================================
# Tests the memory efficiency improvements by running with
# different optimization levels and measuring memory usage.

set -e  # Exit on error

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
ROOT_DIR="$(dirname "$SCRIPT_DIR")"
BIN="$ROOT_DIR/target/release/orbweaver"

# Input file to use for testing
INPUT_FILE=$1

if [ -z "$INPUT_FILE" ]; then
    echo "Usage: $0 <fasta_file> [output_dir]"
    echo "Example: $0 data/sample.fasta results"
    exit 1
fi

OUTPUT_DIR=${2:-"memory_test_results"}
mkdir -p "$OUTPUT_DIR"

# Check if orbweaver is built in release mode
if [ ! -f "$BIN" ]; then
    echo "Building orbweaver in release mode..."
    (cd "$ROOT_DIR" && cargo build --release)
fi

echo "=== Orbweaver Memory Efficiency Test ==="
echo "Input file: $INPUT_FILE"
echo "Output directory: $OUTPUT_DIR"

# Function to run a test and measure memory usage
run_test() {
    local name=$1
    local args=$2
    local output_file="$OUTPUT_DIR/${name}_results.log"
    
    echo -n "Running test '$name'... "
    
    # Use time command to measure CPU and memory usage
    /usr/bin/time -v $BIN $args > "$output_file" 2>&1
    
    # Extract maximum resident set size from time output
    local mem_usage=$(grep "Maximum resident set size" "$output_file" | awk '{print $6}')
    local runtime=$(grep "User time" "$output_file" | awk '{print $4}')
    
    echo "Done. Memory: $mem_usage KB, Time: $runtime seconds"
    echo "Results saved to $output_file"
    
    # Return memory usage in KB
    echo $mem_usage
}

# Test 1: Legacy mode (load entire sequence at once)
echo "Test 1: Legacy mode (full sequence loading)"
LEGACY_MEM=$(run_test "legacy" "-i $INPUT_FILE -j $OUTPUT_DIR/legacy_grammar.json --output-text $OUTPUT_DIR/legacy_grammar.txt")

# Test 2: Streaming mode (process sequence in chunks)
echo "Test 2: Streaming mode (chunk processing)"
STREAM_MEM=$(run_test "streaming" "-i $INPUT_FILE -j $OUTPUT_DIR/streaming_grammar.json --output-text $OUTPUT_DIR/streaming_grammar.txt --chunk-size 5000000")

# Test 3: Streaming with rule eviction (memory-constrained)
echo "Test 3: Streaming with rule eviction"
EVICT_MEM=$(run_test "eviction" "-i $INPUT_FILE -j $OUTPUT_DIR/eviction_grammar.json --output-text $OUTPUT_DIR/eviction_grammar.txt --chunk-size 5000000 --max-rule-count 500")

# Calculate memory improvement percentages
STREAM_IMPROVE=$(echo "scale=2; 100 * ($LEGACY_MEM - $STREAM_MEM) / $LEGACY_MEM" | bc)
EVICT_IMPROVE=$(echo "scale=2; 100 * ($LEGACY_MEM - $EVICT_MEM) / $LEGACY_MEM" | bc)

# Generate report
REPORT_FILE="$OUTPUT_DIR/memory_efficiency_report.txt"
echo "=== Orbweaver Memory Efficiency Report ===" > "$REPORT_FILE"
echo "Input file: $INPUT_FILE" >> "$REPORT_FILE"
echo "Test date: $(date)" >> "$REPORT_FILE"
echo "" >> "$REPORT_FILE"
echo "Memory Usage:" >> "$REPORT_FILE"
echo "1. Legacy mode:           $LEGACY_MEM KB" >> "$REPORT_FILE"
echo "2. Streaming mode:        $STREAM_MEM KB ($STREAM_IMPROVE% reduction)" >> "$REPORT_FILE"
echo "3. With rule eviction:    $EVICT_MEM KB ($EVICT_IMPROVE% reduction)" >> "$REPORT_FILE"
echo "" >> "$REPORT_FILE"
echo "Recommendations:" >> "$REPORT_FILE"

if [ $(echo "$STREAM_IMPROVE > 30" | bc) -eq 1 ]; then
    echo "- Use streaming mode for significant memory savings" >> "$REPORT_FILE"
fi

if [ $(echo "$EVICT_IMPROVE > 50" | bc) -eq 1 ]; then
    echo "- Enable rule eviction for large files (--max-rule-count 500)" >> "$REPORT_FILE"
fi

echo "Memory efficiency report generated: $REPORT_FILE"
echo "Done!" 