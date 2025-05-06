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
| `--max-rule-count <COUNT>` | Maximum number of rules allowed (enables rule eviction) | No limit |
| `--reverse-aware <BOOL>` | Perform reverse complement canonicalization | true |

### Processing Options

| Option | Description | Default |
|--------|-------------|---------|
| `--skip-ns <BOOL>` | Skip N bases in FASTA input | true |
| `--chunk-size <SIZE>` | Process genome in chunks of this size | No chunking |
| `--chunk-overlap <SIZE>` | Overlap between chunks when chunking is enabled | 1000 |
| `--streaming` | Process FASTA in streaming mode (low memory usage) | false |
| `--adaptive-chunking` | Dynamically adjust chunk sizes based on sequence complexity | false |
| `--sequence-indices <LIST>` | Process only specific sequences by index (0-based, comma-separated) | All sequences |
| `--max-memory-per-chunk-mb <SIZE>` | Maximum memory usage per chunk in MB | System dependent |
| `--use-encoding <BOOL>` | Use 2-bit encoding to reduce memory usage | true |
| `--threads <COUNT>` | Number of threads to use for parallel processing | Logical CPU count |

### Miscellaneous

| Option | Description |
|--------|-------------|
| `--stats` | Print statistics about the generated grammar |
| `--profile` | Enable performance profiling |
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

### Streaming Mode

Streaming mode processes sequences without loading the entire file into memory:

- Enabled with the `--streaming` flag
- Reads and processes the FASTA file incrementally
- Ideal for very large genomes where memory is limited
- Can be combined with chunking and rule eviction

Example:
```bash
orbweaver -i huge_genome.fasta --streaming --max-rule-count 5000
```

### Adaptive Chunking

Adaptive chunking dynamically adjusts chunk sizes based on sequence complexity:

- Enabled with the `--adaptive-chunking` flag
- Analyzes sequence entropy to determine optimal chunk sizes
- Balances processing speed and memory usage
- Works best with `--max-memory-per-chunk-mb` to set memory constraints

Example:
```bash
orbweaver -i complex_genome.fasta --adaptive-chunking --max-memory-per-chunk-mb 1000
```

### Multi-Sequence Processing

Orbweaver can process specific sequences from multi-FASTA files:

- Use `--sequence-indices` with a comma-separated list of indices (0-based)
- Processes only the specified sequences, skipping others
- Useful for selective analysis of multi-genome files

Example:
```bash
orbweaver -i multiple_genomes.fasta --sequence-indices 0,2,5 -j output.json
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
| 100 MB - 1 GB | `memory` workflow with `--max-rule-count` |
| 1 GB - 10 GB | `memory` workflow with `--streaming` and `--max-rule-count` |
| > 10 GB | `memory` workflow with `--streaming`, `--adaptive-chunking`, and `--max-memory-per-chunk-mb` |

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

### Streaming Mode for Memory Efficiency
```bash
orbweaver -i huge_genome.fasta -j grammar.json --streaming --max-rule-count 5000
```

### Adaptive Chunking for Complex Genomes
```bash
orbweaver -i complex_genome.fasta -j grammar.json --adaptive-chunking --max-memory-per-chunk-mb 2000
```

### Processing Specific Sequences from Multi-FASTA
```bash
# Process only the first and third sequences (0-based indexing)
orbweaver -i multiple_genomes.fasta -j output.json --sequence-indices 0,2 --stats
```

### Maximum Memory Efficiency for Extremely Large Genomes
```bash
orbweaver -i ultra_large_genome.fasta -j grammar.json \
  --streaming \
  --adaptive-chunking \
  --max-memory-per-chunk-mb 1000 \
  --max-rule-count 10000 \
  --min-rule-usage 5
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
2. Use streaming mode with `--streaming`
3. Enable adaptive chunking with `--adaptive-chunking`
4. Set a memory limit with `--max-memory-per-chunk-mb`
5. Increase the `--min-rule-usage` value to reduce the number of rules generated

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
7. Consider using `--adaptive-chunking` when processing genomes with varying complexity
8. Use `--sequence-indices` to process only the specific sequences you need

```bash
./target/release/orbweaver --help
```

```text
Orbweaver: A grammar-based approach to genomic sequence analysis.

Usage: orbweaver [OPTIONS] --input <INPUT>

Options:
  -i, --input <INPUT>
          Input FASTA file path (.fa, .fasta, .fna).
          
          Path to the file containing DNA sequences in FASTA format.
          Processes all sequences in multi-sequence files by default.
          [required]

  -j, --output-json <OUTPUT_JSON>
          Output JSON file path for the grammar.
          
          Writes a detailed JSON representation of the grammar including all
          rules and the final compressed sequence.

      --output-text <OUTPUT_TEXT>
          Output human-readable text representation of the grammar.
          
          Writes a simplified text format showing rules and the final sequence
          in a more readable format than JSON.

      --output-gfa <OUTPUT_GFA>
          Output GFAv1 representation of the grammar.
          
          Exports the grammar as a graph in GFA (Graphical Fragment Assembly) format,
          which can be visualized with tools like Bandage.

      --visualize <VISUALIZE>
          Generate a .dot file for visualizing the grammar.
          
          Creates a DOT file for visualization with Graphviz tools like dot, neato, etc.
          Example usage: dot -Tpng grammar.dot -o grammar.png

      --export-blocks <EXPORT_BLOCKS>
          Export grammar rules as sequences in FASTA format.
          
          Writes each grammar rule as a separate FASTA record, where each record
          contains the fully expanded DNA sequence for that rule.

      --stats
          Print statistics about the generated grammar.
          
          Displays metrics about the grammar: rule count, depth, compression ratio, etc.

  -k, --kmer-size <KMER_SIZE>
          K-mer size for processing.
          
          Used in some analysis algorithms. The grammar builder currently uses
          digrams (k=2) internally regardless of this setting.
          [default: 21]

      --min-rule-usage <MIN_RULE_USAGE>
          Minimum usage count for a rule to be kept.
          
          A digram must appear at least this many times to be replaced by a rule.
          Higher values lead to fewer rules with more usage each.
          [default: 2]

      --max-rule-count <MAX_RULE_COUNT>
          Maximum number of rules allowed (triggers eviction).
          
          When the number of rules exceeds this limit, less frequently used rules
          are evicted and inlined. If not specified, no limit is enforced.

      --reverse-aware <REVERSE_AWARE>
          Perform reverse complement canonicalization.
          
          When enabled, a digram and its reverse complement are treated as equivalent.
          For example, "AC" and "GT" (reverse complement) are considered the same pattern.
          [default: true]

      --skip-ns <SKIP_NS>
          Skip N bases in FASTA input.
          
          When enabled, 'N' or 'n' characters in the input are skipped,
          as they typically represent unknown bases.
          [default: true]

      --chunk-size <CHUNK_SIZE>
          Process genome in chunks of this size.
          
          Divides large sequences into smaller chunks for processing.
          Useful for very large genomes that exceed available memory.
          If set to 0, will determine size dynamically.
          [default: 0]

      --chunk-overlap <CHUNK_OVERLAP>
          Overlap between chunks when chunking is enabled.
          
          Specifies the number of bases that overlap between adjacent chunks.
          Helps ensure patterns spanning chunk boundaries are detected.
          Note: Only relevant when chunking is enabled.
          [default: 1000]
      
      --threads <THREADS>
          Number of threads to use for parallel processing.
          
          Sets the number of worker threads for parallelized operations.
          Default: number of logical CPU cores.

      --profile
          Enable performance profiling.
          
          Generates a CPU profile and flame graph of the application execution.
          Results are saved to the 'profile' directory.

      --use-encoding <USE_ENCODING>
          Enable 2-bit encoding to reduce memory usage.
          
          Use 2-bit encoding for DNA sequences (A=00, C=01, G=10, T=11),
          reducing memory usage by up to 75% compared to ASCII storage.
          [default: true]
          
      --streaming
          Enable streaming mode for low memory usage.
          
          Process the FASTA file in a streaming fashion, without loading
          the entire sequence into memory at once. Recommended for large genomes.
          
      --adaptive-chunking
          Enable adaptive chunk sizing.
          
          Dynamically adjust chunk sizes based on sequence complexity
          and available memory. Requires chunked mode.
          
      --sequence-indices <SEQUENCE_INDICES>
          Process only specific sequences by index (0-based).
          
          For multi-sequence FASTA files, process only the sequences
          with the specified indices. Default: process all sequences.
          
      --max-memory-per-chunk-mb <MAX_MEMORY_PER_CHUNK_MB>
          Maximum memory usage per chunk (in MB).
          
          Limits the memory usage of each chunk when adaptive chunking
          is enabled. Helps prevent out-of-memory errors.

  -h, --help
          Print help (see more with '--help')

  -V, --version
          Print version
```

## Argument Details

*   `-i, --input <INPUT>`: **Required**. Path to the input FASTA file. Can contain one or more sequences, but currently only the first sequence is processed.
*   `-j, --output-json <PATH>`: Optional. Path to write the resulting grammar (final sequence and rules) in JSON format.
*   `-k, --kmer-size <INT>`: K-mer size for processing. Default is 21 (digrams).
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
*   `--streaming`: If present, processes the FASTA file in streaming mode (low memory usage). Default is false.
*   `--adaptive-chunking`: If present, enables adaptive chunking (dynamically adjusts chunk sizes based on sequence complexity). Default is false.
*   `--sequence-indices <LIST>`: Comma-separated list of indices of sequences to process (0-based).
*   `--max-memory-per-chunk-mb <SIZE>`: Maximum memory usage per chunk in MB.
*   `--threads <COUNT>`: Number of threads to use for parallel processing.
*   `-h, --help`: Prints the help message.
*   `-V, --version`: Prints the version of the tool. 