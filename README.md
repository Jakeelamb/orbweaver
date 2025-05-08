# Orbweaver: Grammar-Based Genomic Sequence Analysis

[![Rust](https://github.com/your-username/orbweaver/actions/workflows/rust.yml/badge.svg)](https://github.com/your-username/orbweaver/actions/workflows/rust.yml)

Orbweaver is a command-line tool written in Rust for building and analyzing context-free grammars (CFGs) from genomic sequences. It identifies repeating patterns in DNA sequences and replaces them with non-terminal symbols (rules), effectively compressing the sequence into a hierarchical grammar representation.

## Table of Contents

- [Features](#features)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [How It Works](#how-it-works)
- [Usage Guide](#usage-guide)
- [Output Formats](#output-formats)
- [Examples](#examples)
- [Configuration](#configuration)
- [Architecture](#architecture)
- [Known Limitations](#known-limitations)
- [Development](#development)
- [Contributing](#contributing)
- [License](#license)
- [Memory Efficiency Improvements](#memory-efficiency-improvements)

## Features

- **FASTA Processing**: Reads sequences from FASTA files with support for skipping ambiguous bases (N).
- **Grammar Construction**: Builds a grammar by iteratively replacing the most frequent digrams.
- **Reverse Complement Awareness**: Can treat a digram and its reverse complement as equivalent.
- **Multiple Output Formats**:
  - **JSON**: Serialized representation of the grammar structure
  - **Text**: Human-readable format showing rules and final sequence
  - **GFA**: Graph-based representation (Graphical Fragment Assembly)
  - **DOT**: For visualization using Graphviz
  - **FASTA**: Exports expanded sequences for each rule
  - **Tabular Repeat Summary**: TSV file with rule ID, usage, length, and sequence preview (`--output-repeats`)
- **Statistics**: Calculates compression ratio, rule depth, and other metrics
- **Parallelization**: Configurable multi-threading for parallel chunk processing on large genomes.
- **GPU Acceleration**: OpenCL-based acceleration for computationally intensive operations. Enabled by default; use `--no-gpu` for CPU fallback.
- **Advanced Memory Management**:
  - **2-bit Encoding**: Efficient DNA storage using 2 bits per base
  - **Streaming Mode**: Process very large genomes without loading entirely into memory
  - **Rule Eviction**: Prevents memory growth by removing less useful rules (`--max-rule-count`)
- **Multi-Sequence Support**: Process selected sequences from multi-FASTA files using sequence indices.

## Installation

### Prerequisites

- Rust and Cargo (1.70.0+)
- For visualization: Graphviz (optional)
- For GPU acceleration: OpenCL-compatible device and drivers (optional)

### From Source

```bash
# Clone the repository
git clone https://github.com/your-username/orbweaver.git
cd orbweaver

# Build the project
cargo build --release

# The executable will be in target/release/orbweaver
```

## Quick Start

```bash
# Process a FASTA file using GPU (default) and output the grammar as JSON
./target/release/orbweaver --input-files input.fasta -j output.json

# Process multiple FASTA files using CPU and output a repeat summary table
./target/release/orbweaver --input-files genome1.fna,genome2.fna --no-gpu --output-repeats repeats_summary.tsv

# Generate all outputs and print statistics (using GPU by default)
./target/release/orbweaver --input-files genome.fna \
  -j output/grammar.json \
  --output-text output/grammar.txt \
  --output-gfa output/grammar.gfa \
  --visualize output/grammar.dot \
  --export-blocks output/rules.fasta \
  --output-repeats output/repeats_summary.tsv \
  --stats

# Memory-efficient processing of large genomes with streaming (GPU default)
./target/release/orbweaver --input-files large_genome.fasta \
  -j output/grammar.json \
  --streaming \
  --max-rule-count 10000 \
  --min-rule-usage 5

# See all available options
./target/release/orbweaver --help
```

## How It Works

Orbweaver implements a variant of the Sequitur algorithm, which iteratively identifies the most frequent adjacent pairs (digrams) in a sequence and replaces them with non-terminal symbols. These replacements create a hierarchical grammar that compactly represents the original sequence.

The process works as follows:

1. **Initialization**: Convert the input DNA sequence to a series of terminal symbols.
2. **Iteration**:
   - Find the most frequent digram in the current sequence
   - If the digram appears at least `min_rule_usage` times:
     - Create a new rule that expands to this digram
     - Replace all occurrences of the digram with the new rule
   - Repeat until no digram appears frequently enough
3. **Output**: The resulting grammar consists of:
   - A final sequence (now containing both terminals and non-terminals)
   - A set of production rules

When "reverse-aware" mode is enabled, digrams that are reverse complements of each other (e.g., "AC" and "GT") are treated as equivalent, which can lead to more effective compression in biological sequences.

## Usage Guide

### Basic Usage

```bash
orbweaver --input-files <file1.fa[,file2.fa,...]> [options]
```

### Required Arguments

- `-i, --input-files <FILES>`: Comma-separated paths to input FASTA files (.fa, .fasta, .fna).

### Output Options

- `-j, --output-json <FILE>`: Write grammar to JSON file
- `--output-text <FILE>`: Write grammar in human-readable text format
- `--output-gfa <FILE>`: Write grammar in GFA (graph) format
- `--output-repeats <FILE>`: Write tabular summary of repeats (TSV format)
- `--visualize <FILE>`: Generate DOT file for visualization
- `--export-blocks <FILE>`: Export grammar rules as FASTA sequences
- `--stats`: Print statistics about the generated grammar

### Grammar Construction Options

- `-k, --kmer-size <INT>`: K-mer size used in some analysis algorithms (default: 21)
- `--min-rule-usage <INT>`: Minimum digram frequency for rule creation (default: 2)
- `--max-rule-count <INT>`: Maximum number of rules allowed (enables rule eviction)
- `--reverse-aware <BOOL>`: Perform reverse complement canonicalization (default: true)

### Input Processing Options

- `--skip-ns <BOOL>`: Skip 'N' bases in input (default: true)
- `--streaming`: Process the FASTA file in streaming mode for lower memory usage
- `--no-gpu`: Disable GPU acceleration and use CPU fallback (GPU is enabled by default).

For full details, run `orbweaver --help`.

## Output Formats

### JSON Format

Contains the complete grammar structure:
- `final_sequence`: Array of symbols in the compressed sequence
- `rules`: Map of rule definitions

```json
{
  "final_sequence": [
    {"id": 10, "symbol_type": {"NonTerminal": 0}, "strand": "+"},
    {"id": 11, "symbol_type": {"NonTerminal": 1}, "strand": "-"}
  ],
  "rules": {
    "0": {
      "id": 0,
      "symbols": [
        {"id": 0, "symbol_type": {"Terminal": 65}, "strand": "+"},
        {"id": 1, "symbol_type": {"Terminal": 67}, "strand": "+"}
      ],
      "usage_count": 4
    }
  }
}
```

### Text Format

Human-readable representation:

```
== Final Sequence (3 symbols) ==
R0+ R1- R0+

== Rules (2 total) ==
R0 [Usage=2] -> A+ C+
R1 [Usage=1] -> R0+ G-
```

### FASTA Format

Each rule expanded to its full sequence:

```
>Rule_0 [Usage=4]
ACGTACGT
>Rule_1 [Usage=2]
TTGC
```

See [Output Formats Documentation](docs/output_formats.md) for more details.

## Examples

### Basic Grammar Construction

```bash
orbweaver --input-files sample.fasta -j grammar.json --stats
```

### Visualizing the Grammar

```bash
# Generate DOT file
orbweaver --input-files sample.fasta --visualize grammar.dot

# Convert to PNG using Graphviz
dot -Tpng grammar.dot -o grammar.png
```

### Memory-Efficient Processing for Large Genomes

```bash
# Use streaming mode
orbweaver --input-files large_genome.fasta -j grammar.json --streaming

# Limit memory usage via rule count with streaming
orbweaver --input-files huge_genome.fasta -j grammar.json --streaming --max-rule-count 10000
```

### GPU-Accelerated Processing (Default)

```bash
# GPU acceleration is enabled by default for faster digram finding
orbweaver --input-files genome.fasta -j grammar.json

# To disable GPU and use CPU fallback:
orbweaver --input-files genome.fasta -j grammar.json --no-gpu
```

### Processing Multiple Input Files

```bash
# Provide a comma-separated list of FASTA files
orbweaver --input-files file1.fasta,file2.fasta,file3.fasta -j grammar.json
```

## Configuration

(Note: Configuration via environment variables or config files is not currently implemented. Use command-line flags.)

### Environment Variables

Orbweaver can be configured using environment variables:
- `ORBWEAVER_THREADS`: Number of threads to use for parallelization
- `ORBWEAVER_LOG_LEVEL`: Set logging verbosity (debug, info, warn, error)

### Configuration File (Coming Soon)

You can also use a YAML/TOML configuration file to specify complex settings:

```yaml
# orbweaver.yaml
input:
  path: "genome.fasta"
  skip_ns: true
output:
  json: "output/grammar.json"
  text: "output/grammar.txt"
grammar:
  min_rule_usage: 3
  reverse_aware: true
```

## Architecture

Orbweaver follows a modular architecture:

- **fasta**: FASTA file reading and sequence handling
- **grammar**: Core grammar construction logic
- **io**: Input/output operations for different formats
- **analysis**: Statistical analysis of the grammar
- **parallel**: Parallelization strategies for large genomes
- **encode**: Efficient DNA encoding and memory optimization
- **gpu**: OpenCL-based GPU acceleration

See [Architecture Documentation](docs/architecture.md) for more details.

## Known Limitations

- **Performance**: Large genomes may require significant time and memory
- **GFA Output**: The GFA representation is experimental and may not be optimal for all visualization tools
- **GPU Acceleration**: Currently requires OpenCL-compatible hardware and drivers

## Development

### Building from Source

```bash
# Clone repository
git clone https://github.com/your-username/orbweaver.git
cd orbweaver

# Run tests
cargo test

# Build in debug mode
cargo build

# Build in release mode
cargo build --release
```

### Using Helper Scripts

```bash
# Build the project
./scripts/build.sh

# Run the application
./scripts/run.sh

# Run all tests
./scripts/test_all.sh

# Test memory efficiency
./scripts/test_memory_efficiency.sh input.fasta results
```

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

1. Fork the repository
2. Create your feature branch (`git checkout -b feature/amazing-feature`)
3. Commit your changes (`git commit -m 'Add some amazing feature'`)
4. Push to the branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request

Please ensure your code follows the project's coding style and includes appropriate tests.

### Development Roadmap

- [x] Implement chunking for processing large genomes
- [x] Add rule eviction strategy for limiting memory usage
- [x] Add streaming mode for large genomes
- [x] Implement adaptive chunking
- [x] Support multiple sequences in FASTA files
- [x] Add GPU acceleration support
- [ ] Improve GFA output format
- [ ] Add benchmarking tools

## License

This project is licensed under the MIT License - see the LICENSE file for details.

### Memory Efficiency Improvements

This project implements three major memory efficiency optimizations:

1. **Streaming Input Processing**: Instead of loading entire FASTA files into memory, sequences can be processed in chunks using a streaming approach, dramatically reducing memory footprint for large genomes.

2. **2-Bit Encoding Optimization**: DNA sequences are efficiently stored using 2-bit encoding (A=00, C=01, G=10, T=11), reducing memory usage by up to 75% compared to ASCII storage.

3. **Rule Eviction Strategy**: A priority queue-based approach tracks usage of grammar rules and can evict least-used rules when memory constraints are reached, ensuring bounded memory usage even for complex sequences (`--max-rule-count`).

To use these memory optimizations:

```bash
# Use streaming and rule eviction
orbweaver --input-files large_genome.fasta --streaming --max-rule-count 5000

# Just use streaming mode with default settings
orbweaver --input-files large_genome.fasta --streaming

# Just use rule eviction with a limit of 10,000 rules
orbweaver --input-files large_genome.fasta --max-rule-count 10000
```

These improvements allow processing of large genomes (>10GB) with controlled memory usage. 