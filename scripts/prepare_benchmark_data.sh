#!/bin/bash
# ======================================================
# Prepare Benchmark Data Script
# ======================================================
# Takes an input FASTA file, removes N/n bases, and writes
# a cleaned FASTA file.

set -e

INPUT_FASTA=$1
OUTPUT_FASTA=$2
MAX_SIZE=${3:-0} # Optional: truncate sequence to MAX_SIZE bases (0 = no limit)

if [ -z "$INPUT_FASTA" ] || [ -z "$OUTPUT_FASTA" ]; then
    echo "Usage: $0 <input_fasta> <output_fasta> [max_size]" >&2
    exit 1
fi

if [ ! -f "$INPUT_FASTA" ]; then
    echo "Error: Input FASTA file not found: $INPUT_FASTA" >&2
    exit 1
fi

# Check if output file exists and has content (basic check)
# If max_size is specified, we might need to regenerate even if it exists
if [ "$MAX_SIZE" -eq 0 ] && [ -f "$OUTPUT_FASTA" ] && [ -s "$OUTPUT_FASTA" ]; then
     output_size=$(grep -v '>' "$OUTPUT_FASTA" | tr -d '\n' | wc -c)
     # If output exists and is non-empty (and no size limit), assume it's okay
     if [ "$output_size" -gt 100 ]; then # Heuristic: assume okay if > 100 bases
         echo "Cleaned file $OUTPUT_FASTA already exists. Skipping preparation." >&2
         exit 0
     else
         echo "Cleaned file $OUTPUT_FASTA exists but seems small or empty. Regenerating." >&2
     fi
fi

echo "Preparing benchmark data from $INPUT_FASTA to $OUTPUT_FASTA..." >&2
if [ "$MAX_SIZE" -gt 0 ]; then
    echo "Truncating sequence to first $MAX_SIZE bases (excluding Ns)." >&2
fi

# Create parent directory if it doesn't exist
mkdir -p "$(dirname "$OUTPUT_FASTA")"

# Process the FASTA file:
# 1. Keep the first header line found.
# 2. Filter out subsequent header lines.
# 3. Remove N/n characters (case-insensitive).
# 4. Concatenate sequence lines.
# 5. If MAX_SIZE > 0, truncate the sequence.
# 6. Reformat the sequence with 80 chars per line.

# Get the first header
header=$(grep -m 1 '>' "$INPUT_FASTA")

# Process the sequence
sequence_cleaned=$(grep -v '>' "$INPUT_FASTA" | tr -d 'Nn\n[:space:]') # Remove Ns, newlines, spaces

# Truncate if needed
if [ "$MAX_SIZE" -gt 0 ] && [ ${#sequence_cleaned} -gt $MAX_SIZE ]; then
    echo "Original cleaned length: ${#sequence_cleaned}, Truncating to: $MAX_SIZE" >&2
    sequence_cleaned=${sequence_cleaned:0:$MAX_SIZE}
else
    MAX_SIZE=${#sequence_cleaned} # Use actual length if not truncating or smaller
fi

# Write header to output file
echo "$header" > "$OUTPUT_FASTA"

# Reformat sequence with 80 chars per line
echo "$sequence_cleaned" | fold -w 80 >> "$OUTPUT_FASTA"

# Ensure final newline if sequence is not empty
if [ ${#sequence_cleaned} -gt 0 ]; then
    # Check if the file already ends with a newline
    if [ "$(tail -c1 "$OUTPUT_FASTA" | wc -l)" -eq 0 ]; then
        echo "" >> "$OUTPUT_FASTA"
    fi
fi


echo "Cleaned data written to $OUTPUT_FASTA ($MAX_SIZE bases)" >&2

exit 0 