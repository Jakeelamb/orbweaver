#!/bin/bash
# ======================================================
# Orbweaver Benchmark Script
# ======================================================
# This script runs benchmarks and profiling for Orbweaver.
# It tests different thread counts and configurations to find optimal settings.
#
# Usage:
#   ./scripts/benchmark.sh [input_fasta] [options]
#
# Options:
#   --threads=N,M,P  - Test with N, M, and P threads (default: 1,2,4,8)
#   --profile        - Enable profiling (requires compilation with --features=profiling)
#   --build          - Build with profiling before running benchmarks
#   --help           - Show this help message
# ======================================================

set -e  # Exit on error

# --- Configuration ---
INPUT_FASTA=""
THREAD_COUNTS="1,2,4,8"
PROFILE=false
BUILD=false
MIN_RULE_USAGE=2
ITERATIONS=3

# --- Parse arguments ---
for arg in "$@"; do
  case $arg in
    --threads=*)
      THREAD_COUNTS="${arg#*=}"
      ;;
    --profile)
      PROFILE=true
      ;;
    --build)
      BUILD=true
      ;;
    --help)
      echo "Orbweaver Benchmark Script"
      echo "Usage: ./scripts/benchmark.sh [input_fasta] [options]"
      echo ""
      echo "Options:"
      echo "  --threads=N,M,P  - Test with N, M, and P threads (default: 1,2,4,8)"
      echo "  --profile        - Enable profiling (requires compilation with --features=profiling)"
      echo "  --build          - Build with profiling before running benchmarks"
      echo "  --help           - Show this help message"
      echo ""
      echo "Examples:"
      echo "  ./scripts/benchmark.sh test.fasta"
      echo "  ./scripts/benchmark.sh test.fasta --threads=1,4,8,16 --profile"
      exit 0
      ;;
    *)
      if [[ "$arg" == --* ]]; then
        echo "Unknown option: $arg"
        echo "Run ./scripts/benchmark.sh --help for usage information"
        exit 1
      else
        INPUT_FASTA="$arg"
      fi
      ;;
  esac
done

# Check for input file
if [ -z "$INPUT_FASTA" ]; then
  echo "Error: No input FASTA file specified"
  echo "Run ./scripts/benchmark.sh --help for usage information"
  exit 1
fi

if [ ! -f "$INPUT_FASTA" ]; then
  echo "Error: Input file not found: $INPUT_FASTA"
  exit 1
fi

# --- Build with profiling if requested ---
if [ "$BUILD" = true ]; then
  echo "Building Orbweaver with profiling support..."
  if [ "$PROFILE" = true ]; then
    cargo build --release --features profiling
  else
    cargo build --release
  fi
  echo "Build complete."
fi

# --- Locate binary ---
BINARY_PATH="./target/release/orbweaver"
if [ ! -f "$BINARY_PATH" ]; then
    BINARY_PATH="./target/debug/orbweaver"
    if [ ! -f "$BINARY_PATH" ]; then
        echo "Error: Orbweaver binary not found!"
        echo "Please build the project first:"
        echo "  cargo build --release"
        echo "  or use --build option"
        exit 1
    else
        echo "Notice: Using debug build. For accurate benchmarks, use release build."
    fi
fi

# --- Create results directory ---
RESULTS_DIR="benchmark_results"
mkdir -p "$RESULTS_DIR"

# Get file basename for results
BASENAME=$(basename "$INPUT_FASTA" | sed 's/\.[^.]*$//')
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
RESULTS_FILE="$RESULTS_DIR/${BASENAME}_benchmark_${TIMESTAMP}.csv"

# Write CSV header
echo "threads,iteration,runtime_seconds,rules_count,final_sequence_length,compression_ratio,digram_table_time,replacements,inlining_time,max_memory_mb" > "$RESULTS_FILE"

# Print benchmark info
echo "=========================================================="
echo "Running Orbweaver Benchmarks"
echo "Input file: $INPUT_FASTA"
echo "Thread counts: $THREAD_COUNTS"
echo "Iterations per config: $ITERATIONS"
echo "Profiling: $PROFILE"
echo "Results will be written to: $RESULTS_FILE"
echo "=========================================================="

# Convert comma-separated thread counts to array
IFS=',' read -ra THREAD_ARRAY <<< "$THREAD_COUNTS"

# Get initial filesize for calculations
FILESIZE=$(stat -c%s "$INPUT_FASTA" 2>/dev/null || stat -f%z "$INPUT_FASTA")
echo "Input file size: $(numfmt --to=iec-i --suffix=B $FILESIZE 2>/dev/null || echo "$FILESIZE bytes")"

# Run benchmark
for threads in "${THREAD_ARRAY[@]}"; do
  echo ""
  echo "Testing with $threads thread(s)..."
  
  for i in $(seq 1 $ITERATIONS); do
    echo "  Iteration $i/$ITERATIONS..."
    
    # Clear caches for more consistent results (requires sudo)
    # if command -v sudo &> /dev/null; then
    #   sudo sh -c "echo 3 > /proc/sys/vm/drop_caches" 2>/dev/null || true
    # fi
    
    # Base command
    CMD=("$BINARY_PATH" 
         "-i" "$INPUT_FASTA" 
         "-j" "$RESULTS_DIR/${BASENAME}_${threads}threads_${i}.json" 
         "--min-rule-usage" "$MIN_RULE_USAGE" 
         "--threads" "$threads"
         "--use-encoding")
    
    # Add profiling if enabled
    if [ "$PROFILE" = true ]; then
        CMD+=("--profile")
    fi
    
    # Run the command and parse output
    echo "  Running: ${CMD[*]}"
    output=$("${CMD[@]}" 2>&1)
    
    # Extract metrics from the performance summary file
    SUMMARY_FILE="${BASENAME}_perf_summary.txt"
    if [ -f "$SUMMARY_FILE" ]; then
      # Extract runtime
      runtime=$(grep "Total runtime:" "$SUMMARY_FILE" | grep -oE '[0-9]+\.[0-9]+s' | sed 's/s//')
      rules_count=$(grep "Rules generated:" "$SUMMARY_FILE" | grep -oE '[0-9]+')
      final_length=$(grep "Final sequence length:" "$SUMMARY_FILE" | grep -oE '[0-9]+')
      compression=$(grep "Compression ratio:" "$SUMMARY_FILE" | grep -oE '0\.[0-9]+')
      
      # Extract detailed metrics if available (from stdout)
      digram_time=$(echo "$output" | grep "Digram table rebuilding:" | grep -oE '[0-9]+\.[0-9]+s' | head -1 | sed 's/s//')
      replacements=$(echo "$output" | grep "Total replacements:" | grep -oE '[0-9]+' | head -1)
      inlining_time=$(echo "$output" | grep "Rule inlining:" | grep -oE '[0-9]+\.[0-9]+s' | head -1 | sed 's/s//')
      
      # Attempt to get memory usage
      # Note: Since we don't have direct memory measurement, this is approximated
      max_memory_mb=0
      if command -v ps &> /dev/null; then
        # Get the PID before running and monitor memory
        max_memory_mb=$(echo "$output" | grep "Memory usage" | grep -oE 'allocated=[0-9]+MB' | sed 's/allocated=//' | sed 's/MB//' | sort -n | tail -1)
        if [ -z "$max_memory_mb" ]; then
          max_memory_mb="N/A"
        fi
      fi
      
      # Write to CSV
      echo "$threads,$i,$runtime,$rules_count,$final_length,$compression,$digram_time,$replacements,$inlining_time,$max_memory_mb" >> "$RESULTS_FILE"
      
      echo "  Results: time=${runtime}s, rules=${rules_count}, compression=${compression}"
      
      # Move the summary file to avoid conflicts
      mv "$SUMMARY_FILE" "$RESULTS_DIR/${BASENAME}_${threads}threads_${i}_summary.txt"
    else
      echo "  Warning: Could not find performance summary file."
      # Extract what we can from stdout
      runtime=$(echo "$output" | grep "Orbweaver finished in" | grep -oE '[0-9]+\.[0-9]+s' | sed 's/s//')
      echo "$threads,$i,$runtime,N/A,N/A,N/A,N/A,N/A,N/A,N/A" >> "$RESULTS_FILE"
    fi
    
    # Sleep a bit to let system stabilize
    sleep 2
  done
done

echo ""
echo "=========================================================="
echo "Benchmark complete!"
echo "Results saved to: $RESULTS_FILE"
echo ""

# Calculate and show averages
echo "Summary of results (averages):"
echo "Threads,Avg Time(s),Rules,Final Length,Compression Ratio"

for threads in "${THREAD_ARRAY[@]}"; do
  # Calculate averages using awk
  avg=$(awk -F, -v t="$threads" '$1==t {sum+=$3; count++} END {print sum/count}' "$RESULTS_FILE")
  rules=$(awk -F, -v t="$threads" '$1==t {sum+=$4; count++} END {print sum/count}' "$RESULTS_FILE")
  length=$(awk -F, -v t="$threads" '$1==t {sum+=$5; count++} END {print sum/count}' "$RESULTS_FILE")
  comp=$(awk -F, -v t="$threads" '$1==t {sum+=$6; count++} END {print sum/count}' "$RESULTS_FILE")
  
  printf "%s,%0.2f,%0.0f,%0.0f,%0.4f\n" "$threads" "$avg" "$rules" "$length" "$comp"
done

# If profiling was enabled, point to the profiling results
if [ "$PROFILE" = true ]; then
  echo ""
  echo "Profiling results available in the 'profiles' directory."
  if command -v find &> /dev/null; then
    echo "Latest flamegraph: $(find profiles -name "*flamegraph.svg" -type f -printf '%T@ %p\n' 2>/dev/null | sort -n | tail -1 | cut -f2- -d" ")"
  fi
fi

echo "=========================================================="

# Generate simple plot with gnuplot if available
if command -v gnuplot &> /dev/null; then
  echo "Generating performance plot..."
  PLOT_FILE="$RESULTS_DIR/${BASENAME}_benchmark_${TIMESTAMP}.png"
  
  gnuplot << EOF
set terminal pngcairo enhanced size 1000,600
set output "$PLOT_FILE"
set title "Orbweaver Performance by Thread Count"
set xlabel "Threads"
set ylabel "Runtime (seconds)"
set grid
set key top right
set style data histogram
set style histogram cluster gap 1
set style fill solid border rgb "black"
set auto x
plot "$RESULTS_FILE" using 3:xtic(1) title "Runtime" smooth unique
EOF

  echo "Plot saved to: $PLOT_FILE"
fi 