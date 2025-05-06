#!/bin/bash
# ======================================================
# Orbweaver Benchmark Script
# ======================================================
# This script runs benchmarks and profiling for Orbweaver.
# It tests different thread counts, algorithms, and configurations to find optimal settings.
#
# Usage:
#   ./scripts/benchmark.sh [input_fasta] [options]
#
# Options:
#   --threads=N,M,P  - Test with N, M, and P threads (default: 1,2,4,8)
#   --profile        - Enable profiling (requires compilation with --features=profiling)
#   --build          - Build with profiling before running benchmarks
#   --algorithms     - Test different algorithm combinations
#   --memory         - Test memory optimization algorithms
#   --help           - Show this help message
# ======================================================

set -e  # Exit on error

# --- Configuration ---
INPUT_FASTA=""
THREAD_COUNTS="1,2,4,8"
PROFILE=false
BUILD=false
TEST_ALGORITHMS=false
TEST_MEMORY=false
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
    --algorithms)
      TEST_ALGORITHMS=true
      ;;
    --memory)
      TEST_MEMORY=true
      ;;
    --help)
      echo "Orbweaver Benchmark Script"
      echo "Usage: ./scripts/benchmark.sh [input_fasta] [options]"
      echo ""
      echo "Options:"
      echo "  --threads=N,M,P  - Test with N, M, and P threads (default: 1,2,4,8)"
      echo "  --profile        - Enable profiling (requires compilation with --features=profiling)"
      echo "  --build          - Build with profiling before running benchmarks"
      echo "  --algorithms     - Test different algorithm combinations"
      echo "  --memory         - Test memory optimization algorithms"
      echo "  --help           - Show this help message"
      echo ""
      echo "Examples:"
      echo "  ./scripts/benchmark.sh test.fasta"
      echo "  ./scripts/benchmark.sh test.fasta --threads=1,4,8,16 --profile"
      echo "  ./scripts/benchmark.sh test.fasta --algorithms --memory"
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

# Write CSV header with added fields for new metrics
echo "config,threads,iteration,runtime_seconds,rules_count,final_sequence_length,compression_ratio,digram_table_time,replacement_time,replacements,inlining_time,max_memory_mb,algorithm,chunk_size" > "$RESULTS_FILE"

# Get initial filesize for calculations
FILESIZE=$(stat -c%s "$INPUT_FASTA" 2>/dev/null || stat -f%z "$INPUT_FASTA")
echo "Input file size: $(numfmt --to=iec-i --suffix=B $FILESIZE 2>/dev/null || echo "$FILESIZE bytes")"

# Define test configurations
if [ "$TEST_ALGORITHMS" = true ]; then
  # Testing different algorithm combinations
  CONFIGS=(
    "baseline:--min-rule-usage=$MIN_RULE_USAGE"
    "suffix_array:--min-rule-usage=$MIN_RULE_USAGE"
    "parallel_digram:--min-rule-usage=$MIN_RULE_USAGE --chunk-size=100000"
    "optimized:--min-rule-usage=$MIN_RULE_USAGE --chunk-size=100000 --adaptive-chunking"
  )
elif [ "$TEST_MEMORY" = true ]; then
  # Testing memory optimization configurations
  CONFIGS=(
    "no_encoding:--min-rule-usage=$MIN_RULE_USAGE --use-encoding=false"
    "encoding:--min-rule-usage=$MIN_RULE_USAGE --use-encoding=true"
    "rule_eviction:--min-rule-usage=$MIN_RULE_USAGE --use-encoding=true --max-rule-count=5000"
    "streaming:--min-rule-usage=$MIN_RULE_USAGE --use-encoding=true --streaming"
  )
else
  # Default - just test with threads
  CONFIGS=("default:--min-rule-usage=$MIN_RULE_USAGE")
fi

# Print benchmark info
echo "=========================================================="
echo "Running Orbweaver Benchmarks"
echo "Input file: $INPUT_FASTA"
echo "Thread counts: $THREAD_COUNTS"
if [ "$TEST_ALGORITHMS" = true ]; then
  echo "Testing algorithm variants"
elif [ "$TEST_MEMORY" = true ]; then
  echo "Testing memory optimization variants"
fi
echo "Iterations per config: $ITERATIONS"
echo "Profiling: $PROFILE"
echo "Results will be written to: $RESULTS_FILE"
echo "=========================================================="

# Convert comma-separated thread counts to array
IFS=',' read -ra THREAD_ARRAY <<< "$THREAD_COUNTS"

# Run benchmark
for config in "${CONFIGS[@]}"; do
  CONFIG_NAME="${config%%:*}"
  CONFIG_ARGS="${config#*:}"
  
  echo ""
  echo "Testing configuration: $CONFIG_NAME"
  echo "Args: $CONFIG_ARGS"
  
  for threads in "${THREAD_ARRAY[@]}"; do
    echo ""
    echo "  With $threads thread(s)..."
    
    for i in $(seq 1 $ITERATIONS); do
      echo "    Iteration $i/$ITERATIONS..."
      
      # Clear caches for more consistent results (requires sudo)
      # if command -v sudo &> /dev/null; then
      #   sudo sh -c "echo 3 > /proc/sys/vm/drop_caches" 2>/dev/null || true
      # fi
      
      # Build command with configuration
      CMD=("$BINARY_PATH" 
           "-i" "$INPUT_FASTA" 
           "-j" "$RESULTS_DIR/${BASENAME}_${CONFIG_NAME}_${threads}threads_${i}.json" 
           "--threads" "$threads")
      
      # Add configuration arguments
      for arg in $CONFIG_ARGS; do
        CMD+=("$arg")
      done
      
      # Add profiling if enabled
      if [ "$PROFILE" = true ]; then
          CMD+=("--profile")
      fi
      
      # Extract chunk size for reporting
      chunk_size="N/A"
      if [[ "$CONFIG_ARGS" == *"--chunk-size="* ]]; then
        chunk_size=$(echo "$CONFIG_ARGS" | grep -oE -- "--chunk-size=[0-9]+" | cut -d= -f2)
      fi
      
      # Run the command and parse output
      echo "    Running: ${CMD[*]}"
      output=$("${CMD[@]}" 2>&1)
      
      # Extract metrics from the output
      runtime=$(echo "$output" | grep "Grammar construction complete" | grep -oE 'in [0-9]+\.[0-9]+s' | grep -oE '[0-9]+\.[0-9]+' | head -1)
      rules_count=$(echo "$output" | grep "Rules:" | grep -oE '[0-9]+\.$' | grep -oE '[0-9]+' | head -1)
      final_length=$(echo "$output" | grep "Sequence compressed from" | awk '{print $5}' | head -1)
      compression=$(echo "$output" | grep "compression ratio" | grep -oE '[0-9]+\.[0-9]+x' | sed 's/x//' | head -1)
      
      # Extract detailed metrics
      digram_time=$(echo "$output" | grep "Digram table rebuilding:" | grep -oE '[0-9]+\.[0-9]+s' | head -1 | sed 's/s//')
      replacement_time=$(echo "$output" | grep "Replacement time:" | grep -oE '[0-9]+\.[0-9]+s' | head -1 | sed 's/s//')
      replacements=$(echo "$output" | grep "Total replacements:" | grep -oE '[0-9]+' | head -1)
      inlining_time=$(echo "$output" | grep "Rule inlining:" | grep -oE '[0-9]+\.[0-9]+s' | head -1 | sed 's/s//')
      
      # Attempt to get memory usage
      max_memory_mb=0
      if command -v ps &> /dev/null; then
        max_memory_mb=$(echo "$output" | grep "Memory usage" | grep -oE 'allocated=[0-9]+MB' | sed 's/allocated=//' | sed 's/MB//' | sort -n | tail -1)
        if [ -z "$max_memory_mb" ]; then
          max_memory_mb="N/A"
        fi
      fi
      
      # Identify algorithm used (from log messages)
      algorithm="sequential"
      if echo "$output" | grep -q "Using suffix array optimization"; then
        algorithm="suffix_array"
      elif echo "$output" | grep -q "Using parallel hash map for digram table"; then
        algorithm="parallel_hashmap"
      fi
      
      # Write to CSV
      echo "$CONFIG_NAME,$threads,$i,$runtime,$rules_count,$final_length,$compression,$digram_time,$replacement_time,$replacements,$inlining_time,$max_memory_mb,$algorithm,$chunk_size" >> "$RESULTS_FILE"
      
      echo "    Results: time=${runtime}s, rules=${rules_count}, compression=${compression}"
      
      # Sleep a bit to let system stabilize
      sleep 2
    done
  done
done

echo ""
echo "=========================================================="
echo "Benchmark complete!"
echo "Results saved to: $RESULTS_FILE"
echo ""

# Calculate and show averages
echo "Summary of results (averages):"
echo "Config,Threads,Avg Time(s),Rules,Final Length,Compression Ratio,Digram Table Time(s),Replacement Time(s)"

for config in "${CONFIGS[@]}"; do
  CONFIG_NAME="${config%%:*}"
  
  for threads in "${THREAD_ARRAY[@]}"; do
    # Calculate averages using awk
    avg=$(awk -F, -v c="$CONFIG_NAME" -v t="$threads" '$1==c && $2==t {sum+=$4; count++} END {print sum/count}' "$RESULTS_FILE")
    rules=$(awk -F, -v c="$CONFIG_NAME" -v t="$threads" '$1==c && $2==t {sum+=$5; count++} END {print sum/count}' "$RESULTS_FILE")
    length=$(awk -F, -v c="$CONFIG_NAME" -v t="$threads" '$1==c && $2==t {sum+=$6; count++} END {print sum/count}' "$RESULTS_FILE")
    comp=$(awk -F, -v c="$CONFIG_NAME" -v t="$threads" '$1==c && $2==t {sum+=$7; count++} END {print sum/count}' "$RESULTS_FILE")
    digram_time=$(awk -F, -v c="$CONFIG_NAME" -v t="$threads" '$1==c && $2==t {sum+=$8; count++} END {print sum/count}' "$RESULTS_FILE")
    repl_time=$(awk -F, -v c="$CONFIG_NAME" -v t="$threads" '$1==c && $2==t {sum+=$9; count++} END {print sum/count}' "$RESULTS_FILE")
    
    printf "%s,%s,%0.2f,%0.0f,%0.0f,%0.4f,%0.2f,%0.2f\n" "$CONFIG_NAME" "$threads" "$avg" "$rules" "$length" "$comp" "$digram_time" "$repl_time"
  done
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

# Generate plots with gnuplot if available
if command -v gnuplot &> /dev/null; then
  echo "Generating performance plots..."
  
  # Plot 1: Runtime by thread count and config
  PLOT_FILE="$RESULTS_DIR/${BASENAME}_runtime_${TIMESTAMP}.png"
  gnuplot << EOF
set terminal pngcairo enhanced size 1200,800
set output "$PLOT_FILE"
set title "Orbweaver Runtime by Configuration and Thread Count"
set xlabel "Threads"
set ylabel "Runtime (seconds)"
set grid
set key outside right top
set style data linespoints
set pointsize 1.5
set style line 1 lc rgb '#1B9E77' pt 7
set style line 2 lc rgb '#D95F02' pt 9
set style line 3 lc rgb '#7570B3' pt 5
set style line 4 lc rgb '#E7298A' pt 13
plot "$RESULTS_FILE" using 2:4:(strcol(1)) title "Runtime" with linespoints lw 2 pt 7
EOF
  echo "Runtime plot saved to: $PLOT_FILE"
  
  # Plot 2: Digram table time vs replacement time
  if [ "$TEST_ALGORITHMS" = true ]; then
    PLOT_FILE2="$RESULTS_DIR/${BASENAME}_components_${TIMESTAMP}.png"
    gnuplot << EOF
set terminal pngcairo enhanced size 1200,800
set output "$PLOT_FILE2"
set title "Performance Breakdown by Algorithm"
set xlabel "Algorithm"
set ylabel "Time (seconds)"
set grid
set key outside right top
set style data histogram
set style histogram cluster gap 1
set style fill solid 0.7
set boxwidth 0.9
set xtic rotate by -45 scale 0
plot "$RESULTS_FILE" using 8:xtic(1) title "Digram Table Time" lc rgb '#1B9E77', \
     "" using 9 title "Replacement Time" lc rgb '#D95F02'
EOF
    echo "Component time plot saved to: $PLOT_FILE2"
  fi
fi 