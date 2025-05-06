#!/bin/bash
# ======================================================
# Orbweaver Build Script
# ======================================================
# This script builds the Orbweaver binary with optimizations
# for processing large genomes efficiently.
#
# Options:
#   --debug    Build in debug mode (default: release)
#   --help     Show this help message
# ======================================================

set -e  # Exit on error

# Parse command-line arguments
BUILD_TYPE="release"
OPTIMIZED=true

for arg in "$@"; do
  case $arg in
    --debug)
      BUILD_TYPE="debug"
      ;;
    --help)
      echo "Orbweaver Build Script"
      echo ""
      echo "Usage: ./scripts/build.sh [options]"
      echo ""
      echo "Options:"
      echo "  --debug    Build in debug mode (default: release)"
      echo "  --help     Show this help message"
      echo ""
      echo "Note: The binary will be optimized for memory-efficient processing"
      echo "      of large genomes by default using 2-bit encoding."
      exit 0
      ;;
  esac
done

# Display build configuration
echo "Building Orbweaver in ${BUILD_TYPE} mode..."

# Run the build
if [ "$BUILD_TYPE" = "release" ]; then
    # Release build with optimizations
    cargo build --release
    BINARY_PATH="target/release/orbweaver"
else
    # Debug build
    cargo build
    BINARY_PATH="target/debug/orbweaver"
fi

# Check build result
EXIT_CODE=$?
if [ $EXIT_CODE -eq 0 ]; then
    echo "✅ Build successful!"
    echo "Binary located at: ${BINARY_PATH}"
    echo ""
    # Display usage example with encoding
    echo "Example usage with memory-efficient encoding:"
    echo "  ${BINARY_PATH} -i input.fasta -j output.json"
    echo ""
    echo "Or use the run script for predefined workflows:"
    echo "  ./scripts/run.sh memory --input=large_genome.fasta"
    echo ""
    echo "For maximum memory efficiency:"
    echo "  ./scripts/run.sh memory --input=large_genome.fasta"
else
    echo "❌ Build failed with exit code $EXIT_CODE"
fi

exit $EXIT_CODE 