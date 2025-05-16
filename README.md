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
# Process a FASTA file. Species ID is derived from 'input.fasta' (-> 'input').
# Output is saved in ./input/<timestamp_or_run_id>/ by default.
# Default files generated in that directory: grammar.json, rules.fasta, repeats_summary.tsv, grammar.dot, grammar.graphml
./target/release/orbweaver --input-files input.fasta --stats

# Process multiple FASTA files, deriving species ID from the first file ('genome1').
# Output in ./genome1/<timestamp_or_run_id>/ (or specify --run-id for a custom <run_id>).
./target/release/orbweaver --input-files genome1.fna,genome2.fna --stats --no-gpu

# Generate specific outputs, specifying a run ID.
# Output in ./genome/<my_specific_run>/ 
# Explicitly named files are created in that directory. Default-named files also go there.
./target/release/orbweaver --input-files genome.fna --run-id my_specific_run \
  --output-json my_grammar.json \
  --output-text my_grammar.txt \
  --stats

# Memory-efficient processing of large genomes with streaming
# Output in ./large_genome/<timestamp_or_run_id>/
./target/release/orbweaver --input-files large_genome.fasta \
  --streaming \
  --max-rule-count 10000 \
  --min-rule-usage 5 --stats

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
# Species ID derived from first file in list. Output dir: ./<derived_species_id>/<run_id_or_timestamp>/
./orbweaver --input-files <file1.fa[,file2.fa,...]> [options]
```

### Required Arguments

- `-i, --input-files <FILES>`: Comma-separated paths to input FASTA files (.fa, .fasta, .fna).

### Output Directory and Naming

The output directory structure is automatically managed:
- **Species ID**: Derived from the filename of the *first* input file (e.g., `input.fasta` leads to a species ID `input`).
- **Run ID**: Can be specified with `--run-id <ID>`. If not provided, a timestamp (e.g., `YYYYMMDD_HHMMSS`) is used.
- **Output Path**: All files are saved into `./<derived_species_id>/<run_id_or_timestamp>/` in the current working directory.

### Output Options

Many common output files are generated by default inside the run-specific directory if a full path (which is interpreted as a filename) is not given for them:
- `grammar.json` (see `--output-json`)
- `rules.fasta` (see `--export-blocks`)
- `repeats_summary.tsv` (see `--output-repeats`)
- `grammar.dot` (see `--visualize`)
- `grammar.graphml` (see `--output-graphml`)

If you provide a filename to an output option (e.g., `--output-json my_custom_name.json`), that file will be created with that name inside the run-specific directory.

- `-j, --output-json <FILENAME>`: Write grammar to JSON file.
- `--output-text <FILENAME>`: Write grammar in human-readable text format.
- `--output-gfa <FILENAME>`: Write grammar in GFA (graph) format.
- `--output-repeats <FILENAME>`: Write tabular summary of repeats (TSV format).
- `--visualize <FILENAME>`: Generate DOT file for visualization.
- `--export-blocks <FILENAME>`: Export grammar rules as FASTA sequences.
- `--output-graphml <FILENAME>`: Output GraphML representation.
- `--stats`: Print statistics about the generated grammar.

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
# Species ID derived from sample.fasta, output in ./sample/<timestamp_or_run_id>/
# Includes default outputs like grammar.json, rules.fasta, etc., in that directory.
orbweaver --input-files sample.fasta --stats
```

### Visualizing the Grammar

```bash
# Generate DOT file (e.g., ./sample/<timestamp_or_run_id>/grammar.dot by default)
orbweaver --input-files sample.fasta --visualize

# If you specified a filename: orbweaver --input-files sample.fasta --visualize my_grammar_viz.dot
# (This creates ./sample/<timestamp_or_run_id>/my_grammar_viz.dot)

# Convert to PNG using Graphviz (assuming default path & name)
# You'll need to know the exact <timestamp_or_run_id> used by orbweaver for that run.
dot -Tpng sample/<run_id_here>/grammar.dot -o sample/<run_id_here>/grammar.png
```

### Memory-Efficient Processing for Large Genomes

```bash
# Use streaming mode. Output in ./large_genome/<timestamp_or_run_id>/
orbweaver --input-files large_genome.fasta --streaming --stats

# Limit memory usage via rule count with streaming
orbweaver --input-files huge_genome.fasta --streaming --max-rule-count 10000 --stats
```

### GPU-Accelerated Processing (Default)

```bash
# GPU acceleration is enabled by default. Output in ./genome/<timestamp_or_run_id>/
orbweaver --input-files genome.fasta --stats

# To disable GPU and use CPU fallback:
orbweaver --input-files genome.fasta --no-gpu --stats
```

### Processing Multiple Input Files

```bash
# Provide a comma-separated list of FASTA files.
# Species ID derived from file1.fasta. Output in ./file1/<timestamp_or_run_id>/
orbweaver --input-files file1.fasta,file2.fasta,file3.fasta --stats
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