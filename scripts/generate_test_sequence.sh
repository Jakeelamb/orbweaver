#!/bin/bash
# Generate a test sequence of specified size

if [ -z "$1" ]; then
  echo "Usage: $0 <size_in_bases> [output_file]" >&2
  exit 1
fi

SIZE=$1
OUTPUT_FILE=${2:-} # If $2 is not set, OUTPUT_FILE is empty

echo "Generating a random DNA sequence of $SIZE bases..." >&2

# Determine output: stdout or file
if [ -n "$OUTPUT_FILE" ]; then
    # Create parent directory if it doesn't exist
    mkdir -p "$(dirname "$OUTPUT_FILE")"
    exec > "$OUTPUT_FILE" # Redirect stdout to file
fi

echo ">random_test_sequence_${SIZE}"
chars_generated=0
max_line_len=80
current_line_len=0

while [ $chars_generated -lt $SIZE ]; do
    # Determine how many characters to generate in this iteration
    needed_for_line=$((max_line_len - current_line_len))
    remaining_total=$((SIZE - chars_generated))
    to_generate=$( (($needed_for_line < $remaining_total)) && echo $needed_for_line || echo $remaining_total )
    
    # Generate bases and ensure no newlines from tr/head within the sequence data itself
    tr -dc 'ACGT' < /dev/urandom | head -c $to_generate | tr -d '\n'
    
    chars_generated=$((chars_generated + to_generate))
    current_line_len=$((current_line_len + to_generate))
    
    if [ $current_line_len -ge $max_line_len ] && [ $chars_generated -lt $SIZE ]; then
        echo "" # Add a newline to break the line
        current_line_len=0
    fi
done

# Add a final newline if any characters were generated and the last line had content
# This ensures the FASTA file ends correctly.
if [ $chars_generated -gt 0 ] && [ $current_line_len -gt 0 ]; then
    echo ""
fi

if [ -n "$OUTPUT_FILE" ]; then
    echo "Done! Sequence written to $OUTPUT_FILE" >&2
else
    echo "Done! Sequence generated to stdout." >&2
fi 