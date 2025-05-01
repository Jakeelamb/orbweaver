# Orbweaver

Orbweaver is a Rust-based tool for computing assembly indices, constructing substring motif graphs for genome assemblies, and performing comparative genomic analysis.

## Overview

Orbweaver processes genome assemblies at the chromosome level to:

1. Compute the assembly index—a measure of structural complexity based on assembly theory—for individual chromosomes and scaffolds.
2. Construct substring motif graphs (representing Straight-Line Programs - SLPs) per chromosome to represent hierarchical sequence construction.
3. Store results (per-chromosome metrics, graphs, and genome-level summaries) in a Neo4j graph database for comparative analysis across multiple genomes.
4. Generate visualization-ready data for comparative analysis across taxonomic groups.

## Features

- **VBQ Format Support**: Uses the vbinseq (VBQ) format via bqtools for efficient encoding and decoding of genomic sequences.
- **Straight-Line Program Construction**: Implements a greedy, Byte Pair Encoding (BPE)-like algorithm to build hierarchical motif representations of genomic sequences.
- **Graph Representation**: Converts SLPs to directed graphs for visualization and analysis.
- **Parallel Processing**: Processes chromosomes in parallel for improved performance.
- **Neo4j Integration**: Stores genome metadata, chromosome metrics, and motif graphs in a Neo4j graph database.
- **Comparative Analysis**: Generates statistics and data files for comparing genomes across taxonomic groups.

## Installation

### Prerequisites

- Rust 1.56 or later
- Neo4j 4.x or later (if using database features)
- bqtools (for VBQ format handling)

### Building from Source

```bash
git clone https://github.com/yourusername/orbweaver.git
cd orbweaver
cargo build --release
```

The compiled binary will be available at `target/release/orbweaver`.

## Usage

```bash
orbweaver --input genome.fasta --genome-id "GCF_000001405.39" --species-name "Homo sapiens" --related-group "Primates" --neo4j-password "yourpassword"
```

### Command Line Options

```
USAGE:
    orbweaver [OPTIONS] --input <INPUT> --genome-id <GENOME_ID> --species-name <SPECIES_NAME> --related-group <RELATED_GROUP> --neo4j-password <NEO4J_PASSWORD>

OPTIONS:
    -i, --input <INPUT>                      Path to input FASTA file containing chromosomes/scaffolds
    -o, --output <OUTPUT>                    Output directory for storing results [default: output]
        --neo4j-uri <NEO4J_URI>              Neo4j connection URI [default: bolt://localhost:7687]
        --neo4j-user <NEO4J_USER>            Neo4j username [default: neo4j]
        --neo4j-password <NEO4J_PASSWORD>    Neo4j password
        --genome-id <GENOME_ID>              Genome ID (e.g., accession number)
        --species-name <SPECIES_NAME>        Species name
        --related-group <RELATED_GROUP>      Related group (taxonomic group, clade, etc.)
        --threshold-factor <THRESHOLD>        SLP threshold factor (higher means more aggressive compression) [default: 1000000]
    -t, --threads <THREADS>                  Number of threads to use for parallel processing [default: 0]
        --skip-db                            Skip database operations
    -d, --debug                              Enable debug logging
    -h, --help                               Print help information
    -V, --version                            Print version information
```

## Technical Details

### Genome Preprocessing

Orbweaver processes genome assemblies using bqtools to convert FASTA files to VBQ format. It uses the `-p r` policy (random assignment for non-ATCG nucleotides) to handle ambiguous bases.

### SLP Construction

The tool constructs a Straight-Line Program for each chromosome using a greedy algorithm:

1. Initialize with individual nucleotides (A, C, G, T).
2. Iteratively find and replace the most frequent pair of adjacent motifs.
3. When multiple pairs have the same frequency, select the pair whose combined sequence comes first lexicographically.
4. Stop when the highest frequency of any adjacent pair falls below a normalized threshold relative to sequence length.

### Assembly Index Calculation

The assembly index for each chromosome is calculated as:
```
Assembly Index = (Total number of unique motifs generated in the chromosome's SLP) - 4
```

The genome-level assembly index is calculated as a weighted average of chromosome assembly indices, where the weight is the chromosome length.

### Database Schema

Orbweaver uses the following Neo4j schema:

**Nodes**:
- (:Genome {genome_id, species_name, related_group, original_N_percent_avg, total_processed_size_bp, assembly_index_avg})
- (:Chromosome {chromosome_id, processed_length_bp, original_N_percent, assembly_index, graph_depth, avg_motif_length})
- (:Motif {sequence, length})

**Relationships**:
- (:Chromosome)-[:PART_OF]->(:Genome)
- (:Motif)-[:FOUND_IN {frequency}]->(:Chromosome)
- (:Motif)-[:DERIVED_FROM]->(:Motif)

## Output Files

Orbweaver generates several output files:

- `output/vbq/` - Directory containing VBQ files for each chromosome
- `output/temp/` - Temporary files used during processing
- `output/group_statistics.json` - Summary statistics for each taxonomic group
- `output/correlation_analysis.json` - Correlation analysis of genome size vs. assembly index
- `output/chromosome_metrics.csv` - CSV file with chromosome metrics for visualization

## License

[MIT License](LICENSE)

## Citation

If you use Orbweaver in your research, please cite our paper:

```
[Citation information - to be added]
``` 