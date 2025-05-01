# Orbweaver

Orbweaver is a Rust-based tool for computing assembly indices, constructing substring motif graphs for genome assemblies, and performing comparative genomic analysis.

## Overview

Orbweaver processes genome assemblies at the chromosome level to:

1. Compute the assembly index—a measure of structural complexity based on assembly theory—for individual chromosomes and scaffolds.
2. Construct substring motif graphs (representing Straight-Line Programs - SLPs) per chromosome to represent hierarchical sequence construction.
3. Store results (per-chromosome metrics, graphs, and genome-level summaries) in a Neo4j graph database for comparative analysis across multiple genomes.
4. Export CSV files for analysis across taxonomic groups.

## Features

- **VBQ Format Support**: Uses the vbinseq (VBQ) format via bqtools for efficient encoding and decoding of genomic sequences.
- **Straight-Line Program Construction**: Implements a greedy, Byte Pair Encoding (BPE)-like algorithm to build hierarchical motif representations of genomic sequences.
- **FM-Index Integration**: Uses FM-Index for efficient substring queries during SLP construction, optimizing the process for large chromosomes.
- **Graph Representation**: Converts SLPs to directed graphs for visualization and analysis.
- **Parallel Processing**: Processes chromosomes in parallel for improved performance.
- **Metadata Acquisition**: Integrates with NCBI taxonomy to automatically determine taxonomic groups and lineage information.
- **Neo4j Integration**: Stores genome metadata, chromosome metrics, and motif graphs in a Neo4j graph database.
- **Batch Processing**: Efficiently processes batches of genomes using a single command.

## Installation

### Prerequisites

- Rust 1.56 or later
- Neo4j 4.x or later (if using database features)
- bqtools (for VBQ format handling)
- Python 3.6+ with required packages (for metadata acquisition)

### Required Python Packages

For metadata acquisition, the following Python packages are required:
```bash
pip install pandas numpy requests dotenv
```

### Building from Source

```bash
git clone https://github.com/yourusername/orbweaver.git
cd orbweaver
cargo build --release
```

The compiled binary will be available at `target/release/orbweaver`.

## Usage

Orbweaver provides three main commands:

### Process a Single Genome

```bash
orbweaver process --input genome.fasta --genome-id "GCF_000001405.39" --species-name "Homo sapiens" --related-group "Primates" --neo4j-password "yourpassword"
```

### Process a Batch of Genomes

```bash
orbweaver batch --tsv-file genomes.tsv --neo4j-password "yourpassword"
```

The TSV file should contain genome IDs and paths to FASTA files, with optional species names and taxonomic groups:
```
# genome_id   fasta_path
GCF_000001405.39  /path/to/human_genome.fasta
GCF_002880755.1   /path/to/chimp_genome.fasta
```

### Run Analysis Only

```bash
orbweaver analyze --output results --neo4j-password "yourpassword"
```

### Command Line Options

#### Common Options
```
    -o, --output <OUTPUT>                  Output directory for storing results [default: output]
        --neo4j-uri <NEO4J_URI>            Neo4j connection URI [default: bolt://localhost:7687]
        --neo4j-user <NEO4J_USER>          Neo4j username [default: neo4j]
        --neo4j-password <NEO4J_PASSWORD>  Neo4j password
    -t, --threads <THREADS>                Number of threads to use for parallel processing [default: 0]
        --threshold-factor <THRESHOLD>     SLP threshold factor (higher means more aggressive compression) [default: 1000000]
        --use-fm-index                     Use FM-Index for SLP construction
        --metadata-file <FILE>             Path to metadata CSV file [default: metazoan_chromosome_assemblies_with_lineage.csv]
        --metadata-script <SCRIPT>         Path to Python script for metadata generation [default: src/ncbi_scrape.py]
    -d, --debug                            Enable debug logging
    -h, --help                             Print help information
    -V, --version                          Print version information
```

#### Process Command Options
```
    -i, --input <INPUT>                    Path to input FASTA file containing chromosomes/scaffolds
        --genome-id <GENOME_ID>            Genome ID (e.g., accession number)
        --species-name <SPECIES_NAME>      Species name
        --related-group <RELATED_GROUP>    Related group (taxonomic group, clade, etc.)
        --skip-db                          Skip database operations
```

#### Batch Command Options
```
    -i, --tsv-file <TSV_FILE>              Path to TSV file with genome information
```

## Technical Details

### Genome Preprocessing

Orbweaver processes genome assemblies using bqtools to convert FASTA files to VBQ format. It uses the `-p r` policy (random assignment for non-ATCG nucleotides) to handle ambiguous bases.

### SLP Construction

The tool constructs a Straight-Line Program for each chromosome using a greedy algorithm:

1. Initialize with individual nucleotides (A, C, G, T).
2. Iteratively find and replace the most frequent pair of adjacent motifs.
   - Uses FM-Index for efficient substring queries on large sequences.
3. When multiple pairs have the same frequency, select the pair whose combined sequence comes first lexicographically.
4. Stop when the highest frequency of any adjacent pair falls below a normalized threshold relative to sequence length.

### Metadata Acquisition

Orbweaver can automatically fetch and process metadata from NCBI, including:
1. Taxonomic lineage information (kingdom, phylum, class, order, family, genus)
2. Assembly information (release date, assembly level, genome size)
3. Species information

The metadata is used to automatically determine the most appropriate taxonomic grouping for comparative analysis.

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
- (:Motif {sequence, length, creation_iteration})

**Relationships**:
- (:Chromosome)-[:PART_OF]->(:Genome)
- (:Motif)-[:FOUND_IN {frequency}]->(:Chromosome)
- (:Motif)-[:DERIVED_FROM]->(:Motif)

## Output Files

Orbweaver generates several output files:

- `output/vbq/` - Directory containing VBQ files for each chromosome
- `output/temp/` - Temporary files used during processing
- `output/group_statistics.json` - Summary statistics for each taxonomic group
- `output/correlation_analysis.json` - Basic correlation analysis of genome size vs. assembly index
- `output/chromosome_metrics.csv` - CSV file with chromosome metrics
- `output/group_metrics.csv` - Metrics grouped by taxonomic classification

## License

[MIT License](LICENSE)

## Citation

If you use Orbweaver in your research, please cite our paper:

```
[Citation information - to be added]
``` 