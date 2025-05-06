#!/bin/bash
# Test script for memory efficiency improvements

set -e  # Exit on error

mkdir -p memory_test_results

echo "===== Memory Efficiency Test ====="
echo "Testing memory efficiency improvements in Orbweaver"
echo ""

# Generate a test sequence if not exists
if [ ! -f "test.fasta" ]; then
    echo "Generating test sequence..."
    ./scripts/generate_test_sequence.sh 1000000
fi

# Function to measure memory usage
measure_memory() {
    echo "Running with $1..."
    /usr/bin/time -v cargo run --release -- -i test.fasta -j output.json $2 2>&1 | grep "Maximum resident set size" | awk '{print $6}' > "memory_test_results/$1.txt"
    echo "Peak memory: $(cat "memory_test_results/$1.txt") kB"
    echo ""
}

# Baseline (without memory optimizations)
measure_memory "baseline" "--use-encoding false"

# With 2-bit encoding only
measure_memory "with_encoding" "--use-encoding true"

# With streaming input
measure_memory "with_streaming" "--use-encoding true --streaming"

# With rule eviction (limit to 1000 rules)
measure_memory "with_eviction" "--use-encoding true --max-rule-count 1000"

# With all optimizations
measure_memory "all_optimizations" "--use-encoding true --streaming --max-rule-count 1000"

echo "===== Test Complete ====="
echo "Results saved to memory_test_results/" 