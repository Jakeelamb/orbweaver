# Memory Optimization in Orbweaver

This document describes the memory optimization strategies implemented in Orbweaver for processing large genomic sequences.

## Overview

Orbweaver implements four major memory efficiency optimizations:

1. Streaming Input Processing
2. 2-Bit DNA Encoding
3. Rule Eviction Strategy
4. Adaptive Chunking

These techniques allow Orbweaver to process very large genomes with a controlled memory footprint.

## 1. Streaming Input Processing

Instead of loading entire FASTA files into memory, sequences are processed in chunks.

### Implementation

- `FastaStream` class handles incremental file reading
- `GrammarBuilder.process_sequence_chunk()` builds grammar incrementally
- Final grammar is produced with `GrammarBuilder.finalize_grammar()`

### Memory Benefits

- Reduces peak memory usage by >90% for large genomes
- Allows processing of sequences larger than available RAM
- Memory usage becomes primarily dependent on grammar complexity, not input size

### Usage

```bash
orbweaver -i large_genome.fasta --streaming
```

## 2. 2-Bit DNA Encoding

DNA sequences are stored using 2 bits per base instead of 8 bits per ASCII character.

### Implementation

- `EncodedBase` struct encodes A=00, C=01, G=10, T=11
- `BitVector` provides efficient storage and operations on bit-packed data
- All internal operations handle encoded bases directly

### Memory Benefits

- 75% reduction in sequence storage requirements
- Improved cache utilization and memory bandwidth
- Better memory locality for sequence operations

### Usage

```bash
# Enabled by default, but can be explicitly set
orbweaver -i input.fasta --use-encoding
```

## 3. Rule Eviction Strategy

A priority queue approach tracks rule usage and can evict least-used rules.

### Implementation

- `GrammarBuilder` maintains usage counts for each rule
- Binary heap (priority queue) organizes rules by usage frequency
- When rule count exceeds threshold, least-used rules are inlined and removed

### Memory Benefits

- Bounds the total memory used by grammar rules
- Prevents memory growth for highly repetitive sequences
- Allows processing of arbitrarily complex grammars with fixed memory

### Eviction Algorithm

1. Track usage count each time a rule is referenced
2. When rule count exceeds threshold:
   - Create priority queue ordered by (usage_count, rule_depth, rule_id)
   - Inline rules with lowest priority at their usage points
   - Remove inlined rules from the grammar
   - Update affected digrams and continue construction

### Usage

```bash
orbweaver -i input.fasta --max-rule-count 5000
```

## 4. Adaptive Chunking

Dynamically adjusts chunk sizes based on sequence complexity and available memory.

### Implementation

- Analyzes sequence entropy to determine optimal chunk size
- Adjusts boundaries to preserve patterns across chunks
- Monitors memory usage during processing
- Automatically scales chunk size based on sequence features

### Memory Benefits

- Optimizes chunk sizes for efficient processing
- Prevents out-of-memory errors on complex sequences
- Smaller chunks for complex regions, larger for simple/repetitive ones
- Balances memory usage with processing efficiency

### Entropy Analysis

- Calculates information entropy on sequence samples
- Higher entropy (more complex/random) = smaller chunks
- Lower entropy (more repetitive) = larger chunks
- Adjusts for available system memory

### Usage

```bash
orbweaver -i input.fasta --adaptive-chunking
```

With memory constraints:

```bash
orbweaver -i input.fasta --adaptive-chunking --max-memory-per-chunk-mb 1000
```

## 5. Multi-Sequence Processing

Process only selected sequences from multi-FASTA files to reduce memory requirements.

### Implementation

- User can specify which sequences to process by index
- Skips unneeded sequences to avoid memory and processing overhead
- Separate grammar construction for each sequence, with optional merging

### Memory Benefits

- Avoids loading unnecessary sequences into memory
- Focuses computation on sequences of interest
- Enables selective processing of very large multi-FASTA files

### Usage

```bash
# Process only sequences 0, 2, and 5
orbweaver -i multi_sequences.fasta --sequence-indices 0,2,5
```

## Combined Optimization

For extremely large or complex genomes, all techniques can be combined:

```bash
orbweaver -i huge_genome.fasta \
  --streaming \
  --adaptive-chunking \
  --max-memory-per-chunk-mb 1000 \
  --max-rule-count 5000 \
  --min-rule-usage 5 \
  --use-encoding
```

## Performance Comparison

| Genome Size | Without Optimization | With All Optimizations | Memory Reduction |
|-------------|---------------------|------------------------|------------------|
| 10 MB       | ~100 MB             | ~40 MB                 | ~60%             |
| 100 MB      | ~1 GB               | ~200 MB                | ~80%             |
| 1 GB        | ~10 GB              | ~500 MB                | ~95%             |
| 10 GB       | Out of memory       | ~1 GB                  | Enables processing |
| 100 GB      | Out of memory       | ~2-5 GB                | Enables processing |

*Note: Actual figures will vary based on sequence complexity and rule parameters.*

## Memory Usage Monitoring

Orbweaver can estimate and report memory usage during processing:

- 2-bit encoding memory savings automatically calculated
- Memory consumption reported during chunked processing
- Memory estimates provided with `--stats` flag

Example memory output:
```
2-Bit Encoding Memory Savings:
  Reduced memory usage by approximately 75.0% (encoding savings)
  Original: 1000000000 bytes, With 2-bit encoding: ~250000000 bytes
```

## Recommendations

| Genome Size | Recommended Flags |
|-------------|------------------|
| < 10 MB     | Default settings |
| 10-100 MB   | `--use-encoding` (default) |
| 100 MB-1 GB | `--max-rule-count 5000` |
| 1-10 GB     | `--streaming --max-rule-count 5000` |
| 10-100 GB   | `--streaming --adaptive-chunking --max-memory-per-chunk-mb 1000` |
| > 100 GB    | `--streaming --adaptive-chunking --max-memory-per-chunk-mb 1000 --max-rule-count 2000 --min-rule-usage 10` |

## Advanced Configuration

For sequences with varying complexity or specialized requirements:

### Highly Repetitive Genomes

With many repeated patterns (e.g., plant genomes with many repeats):
```bash
orbweaver -i repetitive_genome.fasta --max-rule-count 10000 --min-rule-usage 3
```

### Highly Complex Genomes

With few repeated patterns (high entropy):
```bash
orbweaver -i complex_genome.fasta --min-rule-usage 5 --streaming
```

### Ultra-Large Genomes

For processing extremely large genomes (100GB+):
```bash
orbweaver -i ultra_large.fasta --streaming --adaptive-chunking --max-memory-per-chunk-mb 500 --min-rule-usage 10
```

### Low-Memory Systems

For systems with limited RAM (e.g., 8GB or less):
```bash
orbweaver -i genome.fasta --streaming --max-rule-count 1000 --min-rule-usage 5 --max-memory-per-chunk-mb 500
``` 