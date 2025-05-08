# Orbweaver Architecture

This document outlines the module structure, data flow, algorithms, and design decisions of the Orbweaver grammar builder.

## Table of Contents

1. [Module Structure](#module-structure)
2. [Core Data Types](#core-data-types)
3. [Main Algorithm](#main-algorithm)
4. [Data Flow](#data-flow)
5. [Design Patterns](#design-patterns)
6. [Memory Management](#memory-management)
7. [Parallelization Strategy](#parallelization-strategy)
8. [Error Handling](#error-handling)
9. [Testing Strategy](#testing-strategy)
10. [Future Design Considerations](#future-design-considerations)
11. [GPU Acceleration](#gpu-acceleration)

## Module Structure

Orbweaver is organized into several logical modules, each with a specific responsibility:

```
orbweaver/
├── src/                    # Source code
│   ├── main.rs             # Entry point, CLI handling (clap)
│   ├── lib.rs              # Library exports
│   ├── analysis/           # Analysis algorithms
│   │   ├── stats.rs        # Grammar statistics calculation
│   │   └── mod.rs          # Module definitions
│   ├── encode/             # Encoding and memory optimization
│   │   ├── dna_2bit.rs     # 2-bit encoding for DNA bases
│   │   ├── bitvec.rs       # Bit-level operations and memory optimization
│   │   ├── kmer.rs         # K-mer representation and operations
│   │   └── mod.rs          # Module definitions
│   ├── fasta/              # FASTA handling
│   │   ├── reader.rs       # FASTA file reading (with streaming support)
│   │   ├── encoder.rs      # Sequence encoding utilities
│   │   └── mod.rs          # Module definitions
│   ├── grammar/            # Core grammar components
│   │   ├── symbol.rs       # Symbol representation (Terminal/NonTerminal)
│   │   ├── rule.rs         # Rule data structure
│   │   ├── digram.rs       # Digram representation and operations
│   │   ├── digram_table.rs # Efficient digram tracking
│   │   ├── builder.rs      # Grammar construction algorithm
│   │   ├── engine.rs       # Grammar types and interfaces
│   │   └── mod.rs          # Module definitions
│   ├── io/                 # Input/Output operations
│   │   ├── output_json.rs  # JSON serialization
│   │   ├── output_gfa.rs   # GFA format output
│   │   ├── output_text.rs  # Human-readable text output
│   │   ├── output_fasta.rs # FASTA export of rules
│   │   ├── output_dot.rs   # DOT graph format for visualization
│   │   └── mod.rs          # Module definitions
│   ├── parallel/           # Parallelization strategies
│   │   ├── engine.rs       # Parallel execution engine 
│   │   ├── chunking.rs     # Adaptive chunking algorithm
│   │   └── mod.rs          # Module definitions
│   └── utils/              # Shared utilities
│       ├── export.rs       # Export utilities for grammar
│       ├── hash.rs         # Hashing utilities
│       ├── io.rs           # I/O utilities
│       ├── progress.rs     # Progress tracking
│       ├── stats.rs        # Statistical utilities
│       ├── visualization.rs# Visualization utilities
│       └── mod.rs          # Module definitions
├── tests/                  # Integration tests
│   ├── integration_tests.rs# End-to-end tests
│   ├── memory_opt_tests.rs # Memory optimization tests
│   └── memory_opt_integration.rs # Memory optimization integration tests
```

## Core Data Types

### Symbol

The `Symbol` struct represents a single unit in the grammar. It can be either a terminal (a DNA base) or a non-terminal (a reference to a rule). Each symbol has:

```rust
pub struct Symbol {
    pub id: usize,             // Unique identifier for this instance
    pub symbol_type: SymbolType, // Terminal or NonTerminal
    pub strand: Direction,     // Forward or Reverse strand
}

pub enum SymbolType {
    Terminal(EncodedBase),  // 2-bit encoded DNA base
    NonTerminal(usize),     // Rule ID
}

pub enum Direction {
    Forward,
    Reverse,
}
```

### EncodedBase

The `EncodedBase` struct provides memory-efficient storage of DNA bases:

```rust
pub struct EncodedBase(pub u8); // A=0, C=1, G=2, T=3
```

### Rule

The `Rule` struct represents a grammar production rule:

```rust
pub struct Rule {
    pub id: usize,           // Unique identifier
    pub symbols: Vec<Symbol>, // The expansion sequence (typically 2 symbols initially)
    pub usage_count: usize,  // Number of times this rule appears
    pub positions: Vec<usize>, // Original positions where this rule was found
    pub depth: Option<usize>, // Hierarchical depth of this rule
}
```

### DigramTable

The `DigramTable` is a specialized hash map that efficiently tracks digrams (pairs of adjacent symbols) in the sequence:

```rust
pub type DigramKey = u64; // A 64-bit hash representing a canonical digram

pub struct DigramTable {
    // Maps canonical digram key -> Vec of (position, original digram instance)
    occurrences: HashMap<DigramKey, Vec<(usize, (Symbol, Symbol))>>,
}
```

### Grammar

The `Grammar` struct represents the final output of the grammar construction process:

```rust
pub struct Grammar {
    pub sequence: Vec<Symbol>,   // The compressed sequence
    pub rules: HashMap<usize, Rule>, // Rules discovered during inference
    pub max_depth: usize,       // Maximum rule depth (for hierarchy analysis)
}
```

### GrammarBuilder

The `GrammarBuilder` orchestrates the construction process:

```rust
pub struct GrammarBuilder {
    sequence: Vec<Symbol>,        // The working sequence
    rules: HashMap<usize, Rule>,  // Rules created during processing
    digram_table: DigramTable,    // Tracks digram occurrences
    next_rule_id: usize,          // Counter for assigning rule IDs
    min_rule_usage: usize,        // Configuration settings
    reverse_aware: bool,          // Configuration settings
    max_rule_count: Option<usize>, // Rule count limit for eviction
    rule_depths: HashMap<usize, usize>, // Track rule depths
    metrics: PerformanceMetrics,  // Performance tracking
    stream_mode: bool,            // Streaming mode flag
    chunk_count: usize,           // Number of chunks processed
    total_bases_processed: usize, // Total bases processed
}
```

## Main Algorithm

Orbweaver implements a variant of the Sequitur algorithm with extensions for DNA sequences. The core algorithm can be summarized as:

1. **Initialize**: Convert the input DNA sequence into a list of Terminal symbols.
2. **Iterate**:
   - Scan the sequence to build the digram table, tracking all adjacent symbol pairs
   - Find the most frequent digram (above min_rule_usage threshold)
   - Create a new rule for this digram
   - Replace all occurrences of this digram with the new rule
   - Repeat until no digram occurs frequently enough

### Pseudocode

```
procedure BuildGrammar(sequence, min_rule_usage, reverse_aware):
    symbols = ConvertToSymbols(sequence)
    rule_id = 0
    
    while true:
        digram_table = BuildDigramTable(symbols, reverse_aware)
        most_frequent = FindMostFrequentDigram(digram_table)
        
        if most_frequent.count < min_rule_usage:
            break
            
        new_rule = CreateRule(rule_id, most_frequent.digram)
        rule_id += 1
        
        symbols = ReplaceOccurrences(symbols, most_frequent.occurrences, new_rule)
    
    return (symbols, rules)
```

### Memory-Optimized Extensions

Several extensions to the core algorithm improve memory efficiency:

1. **Streaming Processing**: Process the sequence in chunks rather than all at once
2. **Adaptive Chunking**: Dynamically adjust chunk sizes based on sequence entropy
3. **Rule Eviction**: Remove less useful rules when exceeding memory constraints
4. **2-bit Encoding**: Store DNA bases in 2 bits instead of 8-bit ASCII
5. **Suffix Array Optimization**: Use suffix arrays for initial digram finding in large sequences

### Reverse Complement Handling

When `reverse_aware` is enabled, the digram table treats a digram and its reverse complement as equivalent. For example, "AC" (forward strand) and "GT" (reverse strand) are tracked together as a single canonical form. This is particularly useful for DNA where patterns can appear on either strand.

## Data Flow

```
┌────────────────┐     ┌─────────────────┐     ┌──────────────────┐
│                │     │                 │     │                  │
│  FASTA Input   │────►│  EncodedBase    │────►│  Vec<Symbol>     │
│                │     │  Sequence       │     │                  │
└────────────────┘     └─────────────────┘     └──────────────────┘
                                                         │
                                                         ▼
┌────────────────┐     ┌─────────────────┐     ┌──────────────────┐
│                │     │                 │     │                  │
│  Output        │◄────│  Final Grammar  │◄────│  Iterative Rule  │
│  Formats       │     │  Structure      │     │  Creation        │
│                │     │                 │     │                  │
└────────────────┘     └─────────────────┘     └──────────────────┘
```

### Streaming Mode Flow

```
┌────────────────┐     ┌─────────────────┐     ┌──────────────────┐
│                │     │                 │     │                  │
│  FASTA Input   │────►│  Chunk 1        │────►│  Process Chunk   │
│  Stream        │     │                 │     │                  │
└────────────────┘     └─────────────────┘     └──────────────────┘
       │                                                 │
       │                                                 ▼
       │               ┌─────────────────┐     ┌──────────────────┐
       │               │                 │     │                  │
       └──────────────►│  Chunk 2...N    │────►│  Process Chunk   │
                       │                 │     │                  │
                       └─────────────────┘     └──────────────────┘
                                                         │
                                                         ▼
                                               ┌──────────────────┐
                                               │                  │
                                               │  Finalize        │
                                               │  Grammar         │
                                               │                  │
                                               └──────────────────┘
```

1. **Input Processing** (`src/main.rs` -> `fasta::reader`):
   - Reads FASTA file using streaming or memory-mapped approaches
   - Optionally filters out 'N' bases
   - Processes selected sequences by index (if specified)
   - Returns sequence as `Vec<EncodedBase>`

2. **Symbol Conversion** (`grammar::builder::initialize_sequence`):
   - Converts encoded bases to `Vec<Symbol>` with Terminal symbols
   - Assigns position IDs and strand information

3. **Digram Tracking** (`grammar::builder::rebuild_digram_table` -> `digram_table::add_digram`):
   - Scans the current sequence to find all adjacent pairs
   - Handles reverse complement canonicalization if enabled
   - Creates efficiency map of digrams to their positions
   - Uses suffix arrays for initial optimization on large sequences

4. **Rule Creation Loop** (`grammar::builder::step`):
   - Queries the digram table for most frequent pair
   - Creates a new rule if frequency meets threshold
   - Replaces all occurrences with the non-terminal
   - Rebuilds the digram table and repeats
   - Performs rule eviction if exceeding memory limits

5. **Output Generation** (`io::output_*`):
   - Creates various representations of the final grammar
   - Some formats (FASTA, GFA) require rule expansion

## Design Patterns

### Builder Pattern

The `GrammarBuilder` class implements the Builder pattern, handling the step-by-step construction of the grammar. It abstracts the complexity of the grammar construction process while providing flexibility through configuration parameters.

### Strategy Pattern

The canonicalization strategy in `DigramTable` implements a form of the Strategy pattern. The behavior for determining canonical forms changes based on whether reverse complement awareness is enabled.

### Factory Methods

The `Symbol` struct includes factory methods (`terminal` and `non_terminal`) that encapsulate the creation logic, making it more maintainable and preventing invalid states.

### Observer Pattern

The `ProgressTracker` uses an observer-like pattern to monitor and report progress, with periodic callbacks to display status.

## Memory Management

### Memory Efficiency

- **2-bit Encoding**: DNA bases stored as 2 bits (A=00, C=01, G=10, T=11) instead of 8-bit ASCII
- **Streaming Processing**: Processes files incrementally without loading entire sequences
- **Rule Eviction**: Priority-based eviction of least-used rules to limit memory growth
- **Adaptive Chunking**: Dynamically sizes chunks based on sequence complexity
- **Memory-Mapped I/O**: Efficient file access without loading entire contents
- **Parallel Chunking**: Processes large genomes in overlapping chunks
- **BitVector**: Efficient bit-level storage and operations for encoded sequences

### Memory Optimizations

- The `DigramTable` efficiently tracks digram occurrences without duplicating symbol data
- Symbols are stored with minimal overhead (type, ID, strand only)
- The **replace_occurrences** method works in-place on the sequence
- Entropy analysis adjusts chunk sizes for complex vs. repetitive regions
- Memory usage monitoring and adaptive limits per chunk
- Selective sequence processing to avoid loading unnecessary sequences

## Parallelization Strategy

Orbweaver employs multiple parallelization strategies:

1. **Chunk-Based Parallelism**:
   - Divides large sequences into overlapping chunks
   - Processes chunks in parallel using thread pools
   - Uses rayon's work-stealing thread pool for load balancing
   - Merges chunk results to create a unified grammar

2. **Digram Table Parallelism**:
   - Parallel population of digram table for large sequences
   - Concurrent processing of non-overlapping digrams
   - Parallel canonicalization of digrams

3. **Multi-Sequence Parallelism**:
   - Processes separate sequences in parallel
   - Applies parallelism at both sequence and chunk levels
   - Merges grammars from different sequences

4. **Adaptive Parallelism**:
   - Adjusts thread allocation based on available cores
   - Optimizes parallelism based on sequence characteristics
   - Thread count configurable via CLI

### Implementation Details

The `parallel` module coordinates parallelism:

- `parallel::engine`: Manages parallel grammar construction
- `parallel::chunking`: Handles chunk creation and processing
- Uses channels for efficient communication between worker threads
- Thread synchronization via atomic operations and barriers
- Progress tracking for parallel operations

## Error Handling

Orbweaver uses the `anyhow` crate for flexible error handling:

- Error types propagate upward with context
- Multi-level error context provides helpful debugging information
- The CLI presents user-friendly error messages

Example:

```rust
read_fasta_sequences(&args.input, args.skip_ns)
    .context("Failed to read FASTA file")?;
```

### Error Recovery

- Graceful handling of memory allocation failures
- Recovery from chunk processing errors
- Adaptive reduction of memory usage when constraints are reached

## Testing Strategy

The codebase employs a multi-level testing approach:

1. **Unit Tests**: Within module files using `#[cfg(test)]`
   - Tests for `Symbol`, `Rule`, `DigramTable` creation and manipulation
   - Algorithm validations for smaller inputs
   - Memory optimization verification

2. **Integration Tests**: In the `tests/` directory
   - End-to-end tests that verify the CLI and file outputs
   - Tests for edge cases like empty files, invalid inputs, etc.
   - Memory optimization integration tests

3. **Performance Tests**: Found in `memory_opt_tests.rs`
   - Verifies memory reduction techniques
   - Benchmarks performance of different strategies

4. **Continuous Integration**:
   - GitHub Actions workflow runs tests on each commit
   - Verifies functionality across platforms

## Future Design Considerations

1. **Further Chunking Optimizations**:
   - Improve chunk boundary handling for better pattern detection
   - Dynamic adaptation of chunking strategies based on runtime metrics

2. **Advanced Rule Eviction**:
   - Develop more sophisticated eviction policies based on rule utility
   - Consider compression impact in eviction decisions
   - Explore probabilistic data structures for rule tracking

3. **Rule Optimization**:
   - Identify rules that can be merged or decomposed
   - Consider optimal rule extraction based on compression metrics
   - Implement rule hierarchy pruning

4. **Distributed Processing**:
   - Support for distributed computation across multiple machines
   - Cloud-native deployment options

5. **GPU Acceleration**:
   - Explore GPU-based acceleration for digram finding
   - Optimize suffix array construction for GPUs 

## GPU Acceleration

Orbweaver implements GPU acceleration for computationally intensive operations using OpenCL, significantly improving performance for large genome sequences.

### GPU Module Structure

```
orbweaver/
├── src/
│   ├── gpu/                    # GPU acceleration module
│   │   ├── mod.rs              # Core GPU context handling
│   │   ├── digram.rs           # Digram finding on GPU
│   │   └── suffix_array.rs     # Suffix array construction on GPU
├── orbweaver-kernels/          # OpenCL kernel code
│   ├── src/
│   │   └── lib.rs              # Embedded kernel source code
```

### Implementation Details

#### OpenCL Kernels

1. **Digram Finding Kernel**: 
   - Processes DNA sequences in parallel to identify and count adjacent symbol pairs
   - Uses atomic operations for thread-safe counting
   - Implements reverse complement canonicalization directly on the GPU
   - Handles 2-bit encoded sequence data for memory efficiency

2. **Suffix Array Construction Kernel**:
   - Implements prefix doubling algorithm for efficient suffix array construction
   - Uses a multi-stage approach with host-side sorting coordination
   - Initializes, sorts, and updates rank information in parallel

#### GPU Context Management

The `GpuContext` struct encapsulates an OpenCL environment:

```rust
pub struct GpuContext {
    pub platform: Platform,
    pub device: Device,
    pub context: Context,
    pub queue: Queue,
    pub program: Option<Program>,
}
```

- Automatically handles device selection and fallback to CPU if no GPU is available
- Manages OpenCL program compilation and kernel execution
- Provides memory management for device buffers

#### Integration with Grammar Construction

The GPU-accelerated workflow integrates with the core grammar construction:

1. Sequence data is uploaded to GPU memory once
2. Iterative digram finding is performed on the GPU
3. Results are transferred back to host memory for rule creation
4. The process continues until no more frequent patterns are found
5. Final grammar is constructed from the discovered rules

### Performance Considerations

- **Memory Transfer Overhead**: Data transfer between host and GPU memory can become a bottleneck for smaller sequences
- **Batch Processing**: Operations are batched to minimize PCIe bus transfers
- **Kernel Fusion**: Multiple operations are combined in a single kernel where possible
- **Fallback Mechanisms**: Automatic fallback to CPU implementation if GPU operations fail

### Future Enhancements

- **Streaming Processing**: Support for processing sequences larger than GPU memory
- **Multi-GPU Support**: Distribution of work across multiple GPUs
- **Advanced Algorithms**: Implementation of more specialized algorithms like skew algorithm for suffix array construction
- **Compressed Data Structures**: Keeping data compressed in GPU memory 