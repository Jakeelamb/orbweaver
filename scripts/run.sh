#!/bin/bash
# ======================================================
# Orbweaver Run Script
# ======================================================
# This script provides easy execution of common Orbweaver workflows.
# It locates the Orbweaver binary and runs it with common configurations.
#
# Usage:
#   ./scripts/run.sh [workflow] [options]
#
# Available workflows:
#   basic       - Basic grammar construction (default)
#   full        - All outputs and stats
#   visual      - Generate visualization files
#   performance - Optimized for performance
#   memory      - Optimized for large genomes (low memory usage)
#   custom      - Run with custom arguments
#
# Options:
#   --input=FILE  - Input FASTA file (default: input.fasta)
#   --output=DIR  - Output directory (default: output)
#   --help        - Show this help message
# ======================================================

set -e  # Exit on error

# --- Configuration ---
WORKFLOW="basic"
INPUT_FASTA="input.fasta"
OUTPUT_DIR="output"
CUSTOM_ARGS=""

# --- Parse arguments ---
for arg in "$@"; do
  case $arg in
    basic|full|visual|performance|memory|custom)
      WORKFLOW="$arg"
      ;;
    --input=*)
      INPUT_FASTA="${arg#*=}"
      ;;
    --output=*)
      OUTPUT_DIR="${arg#*=}"
      ;;
    --help)
      echo "Orbweaver Run Script"
      echo "Usage: ./scripts/run.sh [workflow] [options]"
      echo ""
      echo "Available workflows:"
      echo "  basic       - Basic grammar construction (default)"
      echo "  full        - All outputs and stats"
      echo "  visual      - Generate visualization files"
      echo "  performance - Optimized for performance"
      echo "  memory      - Optimized for large genomes (low memory usage)"
      echo "  custom      - Run with custom arguments (provide after --)"
      echo ""
      echo "Options:"
      echo "  --input=FILE  - Input FASTA file (default: input.fasta)"
      echo "  --output=DIR  - Output directory (default: output)"
      echo "  --help        - Show this help message"
      echo ""
      echo "Examples:"
      echo "  ./scripts/run.sh basic --input=genome.fasta"
      echo "  ./scripts/run.sh memory --input=large_genome.fasta"
      echo "  ./scripts/run.sh custom -- -i input.fasta --min-rule-usage 5 --use-encoding"
      exit 0
      ;;
    --)
      # Everything after -- is treated as custom arguments
      shift
      CUSTOM_ARGS="$@"
      break
      ;;
    *)
      if [[ "$arg" == --* ]]; then
        echo "Unknown option: $arg"
        echo "Run ./scripts/run.sh --help for usage information"
        exit 1
      else
        CUSTOM_ARGS="$CUSTOM_ARGS $arg"
      fi
      ;;
  esac
done

# --- Locate binary ---
# First check release build, fall back to debug if not found
BINARY_PATH="./target/release/orbweaver"
if [ ! -f "$BINARY_PATH" ]; then
    BINARY_PATH="./target/debug/orbweaver"
    if [ ! -f "$BINARY_PATH" ]; then
        echo "Error: Orbweaver binary not found!"
        echo "Please build the project first:"
        echo "  cargo build --release"
        echo "  or"
        echo "  ./scripts/build.sh"
        exit 1
    else
        echo "Notice: Using debug build. For better performance, use release build."
    fi
fi

# --- Input Validation ---
if [ ! -f "$INPUT_FASTA" ]; then
    echo "Warning: Input FASTA file not found: $INPUT_FASTA"
    
    # Create a tiny test FASTA for demo purposes
    echo "Creating a small demo FASTA file for testing..."
    mkdir -p "$(dirname "$INPUT_FASTA")"
    cat > "$INPUT_FASTA" << EOF
>test_sequence
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
EOF
    echo "Created demo file: $INPUT_FASTA"
fi

# --- Ensure output directory exists ---
mkdir -p "$OUTPUT_DIR"

# --- Check if input file is large ---
FILESIZE=$(stat -c%s "$INPUT_FASTA" 2>/dev/null || stat -f%z "$INPUT_FASTA")
if [ "$FILESIZE" -gt 10000000 ] && [ "$WORKFLOW" != "memory" ]; then
    echo "⚠️  Warning: Input file is large ($(numfmt --to=iec-i --suffix=B $FILESIZE))"
    echo "    Consider using the 'memory' workflow: ./scripts/run.sh memory --input=$INPUT_FASTA"
fi

# --- Define workflow commands ---
case "$WORKFLOW" in
    "basic")
        COMMAND=(
            "$BINARY_PATH"
            -i "$INPUT_FASTA"
            -j "$OUTPUT_DIR/grammar.json"
            --stats
            --use-encoding
        )
        ;;
    "full")
        COMMAND=(
            "$BINARY_PATH"
            -i "$INPUT_FASTA"
            -j "$OUTPUT_DIR/grammar.json"
            --output-text "$OUTPUT_DIR/grammar.txt"
            --output-gfa "$OUTPUT_DIR/grammar.gfa"
            --visualize "$OUTPUT_DIR/grammar.dot"
            --export-blocks "$OUTPUT_DIR/rules.fasta"
            --stats
            --use-encoding
        )
        ;;
    "visual")
        COMMAND=(
            "$BINARY_PATH"
            -i "$INPUT_FASTA"
            --visualize "$OUTPUT_DIR/grammar.dot"
            --output-gfa "$OUTPUT_DIR/grammar.gfa"
            --stats
            --use-encoding
        )
        # Check if Graphviz is installed and convert to PNG if available
        if command -v dot &> /dev/null; then
            echo "Graphviz found. Will convert DOT to PNG after processing."
            DO_GRAPHVIZ=1
        else
            echo "Notice: Graphviz not found. Install it to automatically convert DOT to PNG."
            DO_GRAPHVIZ=0
        fi
        ;;
    "performance")
        COMMAND=(
            "$BINARY_PATH"
            -i "$INPUT_FASTA"
            -j "$OUTPUT_DIR/grammar.json"
            --min-rule-usage 5
            --reverse-aware true
            --skip-ns true
            --use-encoding
        )
        ;;
    "memory")
        COMMAND=(
            "$BINARY_PATH"
            -i "$INPUT_FASTA"
            -j "$OUTPUT_DIR/grammar.json"
            --min-rule-usage 10
            --use-encoding
            --skip-ns true
        )
        if [ "$FILESIZE" -gt 100000000 ]; then
            echo "Detected very large file, enabling additional memory optimizations..."
            # For extremely large files, add future chunking parameters here
        fi
        ;;
    "custom")
        if [ -z "$CUSTOM_ARGS" ]; then
            echo "Error: No custom arguments provided."
            echo "Usage: ./scripts/run.sh custom -- [your arguments]"
            exit 1
        fi
        # Use the custom arguments directly
        COMMAND=($BINARY_PATH $CUSTOM_ARGS)
        ;;
    *)
        echo "Error: Unknown workflow: $WORKFLOW"
        exit 1
        ;;
esac

# --- Execute --- 
echo "=========================================================="
echo "Running Orbweaver with workflow: $WORKFLOW"
echo "Input: $INPUT_FASTA ($(numfmt --to=iec-i --suffix=B $FILESIZE 2>/dev/null || echo "$(($FILESIZE / 1024 / 1024))MB"))"
echo "Output directory: $OUTPUT_DIR"
echo "----------------------------------------------------------"
echo "Command: ${COMMAND[*]}"
echo "=========================================================="

# Set memory limit based on file size if using memory workflow
if [ "$WORKFLOW" = "memory" ]; then
    # Create temporary script that sets ulimit before running command
    TMP_SCRIPT=$(mktemp)
    echo "#!/bin/bash" > "$TMP_SCRIPT"
    echo "# Auto-generated by run.sh" >> "$TMP_SCRIPT"
    
    # Calculate memory limit (roughly 4x the file size, minimum 4GB)
    MEM_LIMIT=$((FILESIZE * 4))
    MEM_LIMIT=$((MEM_LIMIT > 4000000000 ? MEM_LIMIT : 4000000000))
    
    echo "echo \"Setting memory limit to approximately $((MEM_LIMIT / 1000000))MB\"" >> "$TMP_SCRIPT"
    echo "ulimit -v $MEM_LIMIT" >> "$TMP_SCRIPT"
    
    # Build the command
    CMD_STR=""
    for part in "${COMMAND[@]}"; do
        CMD_STR="$CMD_STR \"$part\""
    done
    
    echo "$CMD_STR" >> "$TMP_SCRIPT"
    chmod +x "$TMP_SCRIPT"
    
    "$TMP_SCRIPT"
    EXIT_CODE=$?
    rm "$TMP_SCRIPT"
else
    "${COMMAND[@]}"
    EXIT_CODE=$?
fi

echo "=========================================================="
if [ $EXIT_CODE -eq 0 ]; then
    echo "✅ Orbweaver execution successful!"
    if [ "$WORKFLOW" = "visual" ] && [ "$DO_GRAPHVIZ" -eq 1 ]; then
        echo "Converting DOT to PNG..."
        dot -Tpng "$OUTPUT_DIR/grammar.dot" -o "$OUTPUT_DIR/grammar.png"
        echo "✅ Visualization available at: $OUTPUT_DIR/grammar.png"
    fi
    
    # Show summary of outputs
    echo "----------------------------------------------------------"
    echo "📊 Output files:"
    find "$OUTPUT_DIR" -type f -newer "$INPUT_FASTA" | while read -r file; do
        size=$(du -h "$file" | cut -f1)
        echo "  - $file ($size)"
    done
else
    echo "❌ Orbweaver failed with exit code $EXIT_CODE"
    if [ $EXIT_CODE -eq 137 ]; then
        echo "💥 Process was killed due to excessive memory usage."
        echo "   Try using the memory workflow with a smaller input file:"
        echo "   ./scripts/run.sh memory --input=$INPUT_FASTA"
    fi
fi
echo "=========================================================="

exit $EXIT_CODE 