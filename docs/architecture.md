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

## Module Structure

Orbweaver is organized into several logical modules, each with a specific responsibility:

```
orbweaver/
├── src/                    # Source code
│   ├── main.rs             # Entry point, CLI handling (clap)
│   ├── lib.rs              # Library exports
│   ├── config.rs           # Configuration utilities (placeholder)
│   ├── analysis/           # Analysis algorithms
│   │   ├── motifs.rs       # Pattern finding (suffix arrays, rolling hash, minimizers)
│   │   ├── stats.rs        # Grammar statistics calculation
│   │   └── mod.rs          # Module definitions
│   ├── fasta/              # FASTA handling
│   │   ├── reader.rs       # FASTA file reading (via bio crate)
│   │   ├── encoder.rs      # 2-bit encoding, reverse complementation
│   │   └── mod.rs          # Module definitions
│   ├── grammar/            # Core grammar components
│   │   ├── symbol.rs       # Symbol representation (Terminal/NonTerminal)
│   │   ├── rule.rs         # Rule data structure
│   │   ├── digram_table.rs # Efficient digram tracking
│   │   ├── builder.rs      # Grammar construction algorithm
│   │   └── mod.rs          # Module definitions
│   ├── io/                 # Input/Output operations
│   │   ├── output_json.rs  # JSON serialization
│   │   ├── output_gfa.rs   # GFA format output
│   │   ├── output_text.rs  # Human-readable text output
│   │   ├── output_fasta.rs # FASTA export of rules
│   │   ├── output_dot.rs   # DOT graph format for visualization
│   │   └── mod.rs          # Module definitions
│   └── utils/              # Shared utilities (placeholder)
├── tests/                  # Integration tests
│   └── integration_tests.rs # End-to-end tests
```

## Core Data Types

### Symbol

The `Symbol` struct represents a single unit in the grammar. It can be either a terminal (a DNA base) or a non-terminal (a reference to a rule). Each symbol has:

```rust
pub struct Symbol {
    pub id: usize,             // Unique identifier for this instance
    pub symbol_type: SymbolType, // Terminal or NonTerminal
    pub strand: char,          // '+' or '-' strand
}

pub enum SymbolType {
    Terminal(u8),      // ASCII value of DNA base (e.g., 65 for 'A')
    NonTerminal(usize), // Rule ID
}
```

### Rule

The `Rule` struct represents a grammar production rule:

```rust
pub struct Rule {
    pub id: usize,           // Unique identifier
    pub symbols: Vec<Symbol>, // The expansion sequence (typically 2 symbols initially)
    pub usage_count: usize,  // Number of times this rule appears
}
```

### DigramTable

The `DigramTable` is a specialized hash map that efficiently tracks digrams (pairs of adjacent symbols) in the sequence:

```rust
pub type DigramKey = ((SymbolType, char), (SymbolType, char));

pub struct DigramTable {
    // Maps canonical digram key -> Vec of (position, original digram instance)
    occurrences: HashMap<DigramKey, Vec<(usize, (Symbol, Symbol))>>,
}
```

### GrammarBuilder

The `GrammarBuilder` orchestrates the construction process:

```rust
pub struct GrammarBuilder {
    sequence: Vec<Symbol>,   // The working sequence
    rules: HashMap<usize, Rule>, // Rules created during processing
    digram_table: DigramTable, // Tracks digram occurrences
    next_rule_id: usize,     // Counter for assigning rule IDs
    min_rule_usage: usize,   // Configuration settings
    reverse_aware: bool,     // Configuration settings
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

### Reverse Complement Handling

When `reverse_aware` is enabled, the digram table treats a digram and its reverse complement as equivalent. For example, "AC" (forward strand) and "GT" (reverse strand) are tracked together as a single canonical form. This is particularly useful for DNA where patterns can appear on either strand.

## Data Flow

```
┌────────────────┐     ┌─────────────────┐     ┌──────────────────┐
│                │     │                 │     │                  │
│  FASTA Input   │────►│  Vec<u8> Bytes  │────►│  Vec<Symbol>     │
│                │     │                 │     │                  │
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

1. **Input Processing** (`src/main.rs` -> `fasta::reader`):
   - Reads FASTA file using the `bio` crate
   - Optionally filters out 'N' bases
   - Returns sequence as `Vec<u8>`

2. **Symbol Conversion** (`grammar::builder::initialize_sequence`):
   - Converts raw bytes to `Vec<Symbol>` with Terminal symbols
   - Assigns position IDs and strand information

3. **Digram Tracking** (`grammar::builder::rebuild_digram_table` -> `digram_table::add_digram`):
   - Scans the current sequence to find all adjacent pairs
   - Handles reverse complement canonicalization if enabled
   - Creates efficiency map of digrams to their positions

4. **Rule Creation Loop** (`grammar::builder::step`):
   - Queries the digram table for most frequent pair
   - Creates a new rule if frequency meets threshold
   - Replaces all occurrences with the non-terminal
   - Rebuilds the digram table and repeats

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

## Memory Management

### Memory Efficiency

- The `DigramTable` efficiently tracks digram occurrences without duplicating symbol data
- Symbols are stored with minimal overhead (type, ID, strand only)
- The **replace_occurrences** method works in-place on the sequence

### Potential Improvements

- Chunking for large genomes (not yet implemented) would process segments sequentially
- Rule eviction could limit memory usage by removing infrequently used rules (not yet implemented)

## Parallelization Strategy

The current implementation has parallelization in the `analysis::motifs` module:

- Motif finding algorithms use Rayon's work-stealing thread pool
- Memory-mapped I/O with parallel chunk processing
- Channel-based communication between worker threads

Future enhancement opportunities:
- Parallelized digram scanning
- Concurrent rule replacements for non-overlapping regions

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

## Testing Strategy

The codebase employs a multi-level testing approach:

1. **Unit Tests**: Within module files using `#[cfg(test)]`
   - Tests for `Symbol`, `Rule`, `DigramTable` creation and manipulation
   - Algorithm validations for smaller inputs

2. **Integration Tests**: In the `tests/` directory
   - End-to-end tests that verify the CLI and file outputs
   - Tests for edge cases like empty files, invalid inputs, etc.

3. **Continuous Integration**:
   - GitHub Actions workflow runs tests on each commit
   - Verifies functionality across platforms

## Future Design Considerations

1. **Chunking Strategy**:
   - Divide large genomes into overlapping chunks
   - Process each chunk separately
   - Merge local grammars into a global grammar

2. **Rule Eviction**:
   - Maintain a rule usage priority queue
   - When rule count exceeds threshold, remove least-used rules
   - Inline removed rules back into the sequence

3. **Rule Optimization**:
   - Identify rules that can be merged or decomposed
   - Consider optimal rule extraction based on compression metrics
   - Implement rule hierarchy pruning

4. **Multiple Sequence Handling**:
   - Concatenate sequences with delimiter symbols
   - Build separate grammars for each sequence
   - Create a meta-grammar that references shared rules 