#!/bin/bash
# Generate a test sequence of specified size

if [ -z "$1" ]; then
  SIZE=1000000  # Default 1M bases
else
  SIZE=$1
fi

echo "Generating a random DNA sequence of $SIZE bases..."
echo ">test_sequence" > test.fasta
for i in $(seq 1 $((SIZE / 100))); do
  # Generate 100 random bases at a time
  tr -dc 'ACGT' < /dev/urandom | head -c 100 >> test.fasta
  echo "" >> test.fasta
done

echo "Done! Sequence written to test.fasta" 