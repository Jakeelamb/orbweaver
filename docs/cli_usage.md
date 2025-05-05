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
| `-k, --kmer-size <SIZE>` | K-mer size for processing | 21 |
| `--min-rule-usage <COUNT>` | Minimum usage count for a rule to be kept | 2 |
| `--max-rule-count <COUNT>` | Maximum number of rules allowed | No limit |
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
| > 100 MB | `memory` workflow |
| > 1 GB | `memory` workflow with custom chunk size |

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

1. Ensure you're using the `--use-encoding` flag or the `memory` workflow
2. Increase the `--min-rule-usage` value to reduce the number of rules generated
3. For extremely large genomes, use chunking with `--chunk-size` parameter

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
          K-mer size for processing [default: 21]
          (Note: The grammar builder currently uses digrams (k=2) internally,
           this k-mer size is not directly used by the core builder yet,
           but kept for potential future k-mer based analysis integration)

      --skip-ns
          Skip N bases in FASTA input
          [default: true]

      --chunk-size <CHUNK_SIZE>
          Process genome in chunks of this size (optional)
          (Note: Chunking is not yet implemented)

      --chunk-overlap <CHUNK_OVERLAP>
          Overlap between chunks when chunking is enabled
          [default: 1000]
          (Note: Chunking is not yet implemented)

      --min-rule-usage <MIN_RULE_USAGE>
          Minimum usage count for a rule to be kept (for eviction)
          [default: 2]

      --max-rule-count <MAX_RULE_COUNT>
          Maximum number of rules allowed (triggers eviction)
          (Note: Eviction is not yet implemented)

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

  -h, --help
          Print help (see more with '--help')

  -V, --version
          Print version
```

## Argument Details

*   `-i, --input <INPUT>`: **Required**. Path to the input FASTA file. Can contain one or more sequences, but currently only the first sequence is processed.
*   `-j, --output-json <PATH>`: Optional. Path to write the resulting grammar (final sequence and rules) in JSON format.
*   `-k, --kmer-size <INT>`: K-mer size. Currently **not used** by the core grammar builder, which operates on digrams. Default is 21.
*   `--skip-ns`: If present, skips 'N' or 'n' bases during FASTA reading. Default is true (skip Ns).
*   `--chunk-size <INT>`: *Not yet implemented*. Intended for processing large genomes in chunks.
*   `--chunk-overlap <INT>`: *Not yet implemented*. Specifies overlap for chunked processing. Default is 1000.
*   `--min-rule-usage <INT>`: The minimum number of times a digram must appear to be replaced by a rule. Default is 2.
*   `--max-rule-count <INT>`: *Not yet implemented*. Intended to trigger rule eviction when the total number of rules exceeds this limit.
*   `--reverse-aware`: If present, treats a digram and its reverse complement as equivalent when finding frequent pairs for rule creation. Default is true.
*   `--export-blocks <PATH>`: Optional. Path to write a FASTA file where each record represents a rule and its fully expanded base sequence.
*   `--output-gfa <PATH>`: Optional. Path to write the grammar in GFAv1 format. The representation is experimental.
*   `--output-text <PATH>`: Optional. Path to write the grammar in a human-readable text format, showing the final sequence and rule definitions.
*   `--visualize <PATH>`: Optional. Path to write the grammar structure in DOT format, suitable for visualization with Graphviz.
*   `--stats`: If present, calculates and prints statistics about the grammar (rule count, depth, compression ratio) to standard output after processing.
*   `-h, --help`: Prints the help message.
*   `-V, --version`: Prints the version of the tool. 