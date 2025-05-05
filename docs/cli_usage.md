# Orbweaver CLI Usage Guide

This document provides comprehensive information about using the Orbweaver command-line interface.

## Basic Usage

```bash
./orbweaver -i <input_fasta> -j <output_json>
```

## Command-Line Options

### Input and Output

| Option | Description |
|--------|-------------|
| `-i, --input <FILE>` | Input FASTA file path (required) |
| `-j, --output-json <FILE>` | Output JSON file for the grammar |
| `--output-text <FILE>` | Output human-readable text representation |
| `--output-gfa <FILE>` | Output GFA format for visualization in tools like Bandage |
| `--visualize <FILE>` | Generate a DOT file for visualization with Graphviz |
| `--export-blocks <FILE>` | Export grammar rules as sequences in FASTA format |

### Grammar Construction Parameters

| Option | Description | Default |
|--------|-------------|---------|
| `-k, --kmer-size <SIZE>` | K-mer size for processing | 2 |
| `--min-rule-usage <COUNT>` | Minimum usage count for a rule to be kept | 2 |
| `--max-rule-count <COUNT>` | Maximum number of rules allowed (enables rule eviction) | No limit |
| `--reverse-aware <BOOL>` | Perform reverse complement canonicalization | true |

### Processing Options

| Option | Description | Default |
|--------|-------------|---------|
| `--skip-ns <BOOL>` | Skip N bases in FASTA input | true |
| `--chunk-size <SIZE>` | Process genome in chunks of this size | No chunking |
| `--chunk-overlap <SIZE>` | Overlap between chunks when chunking is enabled | 1000 |
| `--use-encoding <BOOL>` | Use 2-bit encoding to reduce memory usage | true |

### Miscellaneous

| Option | Description |
|--------|-------------|
| `--stats` | Print statistics about the generated grammar |
| `-h, --help` | Show help information |
| `-V, --version` | Display version information |

## Feature Details

### Rule Inlining

Rule inlining reduces grammar complexity by recursively expanding rules that are used only once. This optimization happens automatically during grammar construction and:

- Improves representation efficiency
- Removes unnecessary levels of indirection
- Simplifies the resulting grammar

### Rule Eviction

When `--max-rule-count` is specified, Orbweaver will maintain a limit on the number of rules:

- Rules with the lowest usage counts are evicted first
- Evicted rules are inlined at their usage points
- This prevents memory growth in large genomes
- Recommended for genomes > 100MB

Example:
```bash
orbweaver -i large_genome.fasta --max-rule-count 10000
```

### Chunked Grammar Construction

For very large genomes, Orbweaver can process the sequence in chunks:

- Specify `--chunk-size` to enable chunked processing
- Each chunk is processed while maintaining the global grammar
- Chunks overlap by `--chunk-overlap` bases to preserve patterns at boundaries
- Significantly reduces memory usage for large genomes

Example:
```bash
orbweaver -i large_genome.fasta --chunk-size 1000000 --chunk-overlap 10000
```

## Memory Usage Considerations

Orbweaver can process large genomic sequences efficiently using several memory-optimization techniques:

### 2-Bit Encoding

By default, Orbweaver uses 2-bit encoding for DNA sequences (A=00, C=01, G=10, T=11), reducing memory usage by up to 75% compared to storing sequences as ASCII characters. This encoding is controlled by the `--use-encoding` flag.

For extremely large genomes (100MB+), use the memory-optimized workflow:

```bash
./scripts/run.sh memory --input=large_genome.fasta
```

This workflow:
- Uses 2-bit encoding
- Sets conservative rule parameters
- Automatically calculates appropriate memory limits
- Provides additional error handling for out-of-memory situations

### When to Use Memory Optimization

| Genome Size | Recommended Workflow |
|-------------|---------------------|
| < 10 MB | `basic` workflow (default) |
| 10-100 MB | Use `--use-encoding true` flag |
| > 100 MB | `memory` workflow with `--max-rule-count` |
| > 1 GB | `memory` workflow with `--chunk-size` and `--max-rule-count` |

## Grammar Hierarchy Analysis

Orbweaver tracks rule hierarchy depth information:

- Maximum rule depth (height of the grammar tree)
- Average rule depth (complexity measure)
- This information is displayed when using `--stats`

## Examples

### Basic Usage
```bash
orbweaver -i input.fasta -j output.json
```

### Complete Analysis with All Outputs
```bash
orbweaver -i input.fasta -j grammar.json --output-text grammar.txt --output-gfa grammar.gfa --visualize grammar.dot --export-blocks rules.fasta --stats
```

### Memory-Efficient Processing of Large Genomes
```bash
orbweaver -i large_genome.fasta -j grammar.json --use-encoding true --min-rule-usage 5 --stats
```

### Rule Eviction for Memory Optimization
```bash
orbweaver -i large_genome.fasta -j grammar.json --max-rule-count 5000 --stats
```

### Chunked Processing for Very Large Genomes
```bash
orbweaver -i huge_genome.fasta -j grammar.json --chunk-size 5000000 --chunk-overlap 10000
```

### Using the Run Script for Common Workflows
```bash
# Basic usage
./scripts/run.sh basic --input=input.fasta

# Memory-efficient processing
./scripts/run.sh memory --input=large_genome.fasta

# Complete output formats
./scripts/run.sh full --input=input.fasta --output=results
```

## Troubleshooting

### Memory Issues

If you encounter "out of memory" errors:

1. Enable rule eviction with `--max-rule-count`
2. Use chunked processing with `--chunk-size`
3. Ensure you're using the `--use-encoding` flag or the `memory` workflow
4. Increase the `--min-rule-usage` value to reduce the number of rules generated

### Visualization

To generate visualizations from the DOT file:
```bash
# Generate PNG image
dot -Tpng grammar.dot -o grammar.png

# Generate SVG (for better scalability)
dot -Tsvg grammar.dot -o grammar.svg
```

## Performance Tips

1. Always use a release build for production: `cargo build --release`
2. For large genomes, the `--use-encoding` flag significantly reduces memory usage
3. The `--reverse-aware` option can improve compression but slightly increases processing time
4. Increase `--min-rule-usage` to reduce grammar size at the cost of compression ratio
5. Set an appropriate `--max-rule-count` based on your system's memory limitations
6. For processing multiple files, start with smaller genomes to tune parameters

```bash
./target/release/orbweaver --help
```

```text
Orbweaver: A tool for building grammars from genomic data.

Usage: orbweaver [OPTIONS]

Options:
  -i, --input <INPUT>
          Input FASTA file path (.fa, .fasta, .fna)
          [mandatory]

  -j, --output-json <OUTPUT_JSON>
          Output JSON file path for the grammar (optional)

  -k, --kmer-size <KMER_SIZE>
          K-mer size for processing [default: 2]

      --skip-ns
          Skip N bases in FASTA input
          [default: true]

      --chunk-size <CHUNK_SIZE>
          Process genome in chunks of this size (optional)

      --chunk-overlap <CHUNK_OVERLAP>
          Overlap between chunks when chunking is enabled
          [default: 1000]

      --min-rule-usage <MIN_RULE_USAGE>
          Minimum usage count for a rule to be kept (for eviction)
          [default: 2]

      --max-rule-count <MAX_RULE_COUNT>
          Maximum number of rules allowed (triggers rule eviction)

      --reverse-aware
          Perform reverse complement canonicalization
          [default: true]

      --export-blocks <EXPORT_BLOCKS>
          Export grammar rules as sequences in FASTA format (optional)

      --output-gfa <OUTPUT_GFA>
          Output GFAv1 representation of the grammar (optional)

      --output-text <OUTPUT_TEXT>
          Output human-readable text representation of the grammar (optional)

      --visualize <VISUALIZE>
          Generate a .dot file for visualizing the grammar (optional)

      --stats
          Print statistics about the generated grammar

      --use-encoding
          Use 2-bit encoding to reduce memory usage
          [default: true]

  -h, --help
          Print help (see more with '--help')

  -V, --version
          Print version
```

## Argument Details

*   `-i, --input <INPUT>`: **Required**. Path to the input FASTA file. Can contain one or more sequences, but currently only the first sequence is processed.
*   `-j, --output-json <PATH>`: Optional. Path to write the resulting grammar (final sequence and rules) in JSON format.
*   `-k, --kmer-size <INT>`: K-mer size for processing. Default is 2 (digrams).
*   `--skip-ns`: If present, skips 'N' or 'n' bases during FASTA reading. Default is true (skip Ns).
*   `--chunk-size <INT>`: Size of chunks when processing large genomes. When specified, the genome is processed in overlapping segments of this size.
*   `--chunk-overlap <INT>`: Number of bases to overlap between chunks when chunked processing is enabled. Default is 1000.
*   `--min-rule-usage <INT>`: The minimum number of times a digram must appear to be replaced by a rule. Default is 2.
*   `--max-rule-count <INT>`: Maximum number of rules to maintain before evicting least-used rules. Helps control memory usage.
*   `--reverse-aware`: If present, treats a digram and its reverse complement as equivalent when finding frequent pairs for rule creation. Default is true.
*   `--export-blocks <PATH>`: Optional. Path to write a FASTA file where each record represents a rule and its fully expanded base sequence.
*   `--output-gfa <PATH>`: Optional. Path to write the grammar in GFAv1 format.
*   `--output-text <PATH>`: Optional. Path to write the grammar in a human-readable text format, showing the final sequence and rule definitions.
*   `--visualize <PATH>`: Optional. Path to write the grammar structure in DOT format, suitable for visualization with Graphviz.
*   `--stats`: If present, calculates and prints statistics about the grammar (rule count, depth, compression ratio) to standard output after processing.
*   `--use-encoding`: If present, uses 2-bit encoding for DNA sequences to reduce memory usage. Default is true.
*   `-h, --help`: Prints the help message.
*   `-V, --version`: Prints the version of the tool. 