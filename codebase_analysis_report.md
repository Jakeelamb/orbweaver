# Orbweaver Codebase Analysis Report

This report summarizes the findings from an automated analysis of the Orbweaver codebase, focusing on its structure, core algorithms, performance characteristics, and potential areas for optimization. **Significant refactoring has been undertaken to address key CPU performance bottlenecks, primarily by shifting from full pattern table rebuilds to incremental updates.**

## 1. Project Structure

The Orbweaver project is organized as a Rust workspace with a main binary application and several library crates.

**Workspace Root:**
- `Cargo.toml`, `Cargo.lock`: Project manifest and dependencies.
- `src/`: Main source code for the `orbweaver` binary.
- `orbweaver-utils/`: Workspace member, likely containing command-line utility tools.
  - `src/main.rs`: Indicates it builds one or more binaries.
  - Depends on the main `orbweaver` crate, which is an unusual pattern that might suggest a potential for circular dependencies if not managed carefully.
- `orbweaver-kernels/`: Workspace member, providing OpenCL kernel code as a library.
  - `src/lib.rs`: Defines functions to retrieve kernel source strings.
- `tests/`, `test_data/`: Integration tests and associated data.
- `genomes/`: Likely stores genome data for analysis.
- `scripts/`: Utility scripts.
- `docs/`: Project documentation.
- `build.rs`: Rust build script for pre-compilation tasks.

**Main Application (`src/` directory):**
- `main.rs`: Entry point for the `orbweaver` binary; orchestrates the application flow.
- `lib.rs`: Potentially defines a library component for testing or module re-export.
- `config.rs`: Handles application configuration.
- `args.rs`, `cli.rs`: Define command-line argument parsing and handling (likely using `clap`).
- `run_management.rs`: Code related to managing application runs.
- **Module Directories:**
  - `grammar/`: **Central module** for grammar induction, digram/k-mer processing, and SLP representation. This is a key area for performance.
  - `parallel/`: Code related to parallel processing (e.g., using `rayon`).
  - `utils/`: General utility functions and structures.
  - `io/`: Input/output operations.
  - `fasta/`: Code for handling FASTA formatted files.
  - `gpu/`: Code related to GPU acceleration using OpenCL.
  - `analysis/`: Modules for performing various analyses.
  - `encode/`: Deals with sequence encoding (e.g., 2-bit DNA encoding via `bitnuc`).

## 2. Core Grammar Construction Logic (`src/grammar/`)

This module is at the heart of Orbweaver's functionality.

### 2.1. Fundamental Types

- **`Symbol` (`symbol.rs`):**
  - Represents a terminal (e.g., `EncodedBase` from `dna_2bit`) or a non-terminal (rule ID).
  - Contains an `id` (instance ID in sequence), `symbol_type` (Terminal/NonTerminal), and `strand` (`Direction::Forward`/`Reverse`).
  - **Crucially, `PartialEq`, `Eq`, `Hash`, `Ord` implementations for `Symbol` IGNORE the `id` field**, considering only `symbol_type` and `strand`. This is appropriate for grouping symbols of the same type.
  - `reverse_complement()` method correctly handles terminals and non-terminals (flips strand, complements base for terminals).
  - Efficient due to small size and `Copy` trait.

- **`SymbolKey` (`symbol.rs`):**
  - A struct containing `symbol_type` and `direction`.
  - Explicitly drops the `id` from `Symbol`, serving as the clean key for map lookups.
  - Efficient (`Copy`, `Hash`, `Ord`).

- **`Digram` (`digram.rs`):**
  - Simple struct holding `first: Symbol` and `second: Symbol`.
  - Implements standard traits (`Hash`, `Eq`, `Ord`).
  - Contains helper functions like `count_digrams`, `find_most_frequent_digram` (CPU-based, using `HashMap`), and a notable `find_most_frequent_terminal_digram_suffix_array` (uses `suffix_array` crate for initial raw sequence analysis).

- **`DigramKeyTuple` (`digram_table.rs`):**
  - Struct `(SymbolKey, SymbolKey)`.
  - This is the actual key type used in `DigramTable`'s `DashMap`.
  - Efficient (`Copy`, `Hash`, `Ord`).

- **`Rule` (`rule.rs`):**
  - `id: usize`, `symbols: Vec<Symbol>` (definition of the rule), `usage_count: usize`, `positions: Vec<usize>`, `depth: Option<usize>`.
  - Provides methods for creation and expansion (`expand_recursive`, `expand_for_inlining`) which are important for rule management.

### 2.2. `DigramTable` (`digram_table.rs`)

- **Data Structure:** `occurrences: dashmap::DashMap<DigramKeyTuple, Vec<(usize, DigramSource)>>`.
  - `DigramSource` indicates if the digram is from the original sequence or a derived rule.
  - `DashMap` allows for concurrent population.
- **Key Generation (`canonical_key`):**
  - Takes `(&Symbol, &Symbol)` and `reverse_aware` flag.
  - Creates `DigramKeyTuple` for forward (`s1, s2`) and reverse complemented (`s2_rc, s1_rc`) versions.
  - Canonicalization applies if both symbols are Terminals OR if both are NonTerminals from the *same rule ID*.
  - Returns the lexicographically smaller `DigramKeyTuple`. This logic is sound and relies on efficient `SymbolKey` and `DigramKeyTuple` comparisons.
- **Population:**
  - `add_digram`: Adds a single occurrence.
  - `add_digrams_from_sequence`, `build_parallel`: Uses `rayon::par_windows` to process sequence chunks in parallel and populate the `DashMap`.
- **`find_most_frequent_digram()`:**
  - **Performance Issue:** Iterates over the entire `DashMap` to find the entry with the longest occurrence vector. This is a linear scan and becomes slow as the number of unique digrams (keys in DashMap) increases.
  - **Update:** The implementation has been optimized to clone the `Vec` of occurrences only for the most frequent entry at the end of the scan, rather than for every entry during the iteration. However, the linear scan itself remains.
- **Other Methods:** `remove_occurrence`, `merge`, `clear`.

### 2.3. `KmerTable` (`kmer_table.rs`)

- Generalizes `DigramTable` for k-mers of length `k >= 2`.
- **Data Structure:** `occurrences: std::collections::HashMap<KmerKey, Vec<(usize, Vec<Symbol>)>>`.
  - Uses a standard `HashMap`.
  - `KmerKey = Vec<Symbol>`: The key is a vector of `Symbol`s.
  - **Performance Note:** Using `Vec<Symbol>` as a key is less efficient for hashing and equality checks than fixed-size tuples, especially for larger `k`. Cloning these keys is also more expensive.
- **Canonicalization (`get_canonical_kmer`):**
  - More complex than for digrams.
  - Normalizes k-mers (sets symbol IDs to 0, strands to Forward) before comparison with their reverse complement (if applicable).
  - Involves multiple iterations and cloning, making it costlier per k-mer.
- **Population (`add_kmers_from_sequence`):**
  - Uses `par_windows` but collects results into a temporary `HashMap` before merging into the main `occurrences` map (standard pattern for non-concurrent `HashMap`).
- **`find_most_frequent_kmer()`:**
  - **Performance Issue:** Similar to `DigramTable`, performs a linear scan of the `HashMap`.
  - **Note:** This function already employed an efficient strategy of cloning the occurrence vector only for the most frequent k-mer. The linear scan aspect remains.
- **Memory Usage:** Potentially high due to storing `Vec<Symbol>` for keys and for each occurrence.

### 2.4. `GrammarBuilder` (`builder.rs`) - The Iterative Process

- Orchestrates grammar construction by iteratively finding and replacing frequent patterns.
- **Core Loop (`build_grammar` calling `step()` repeatedly):**
  1.  **`rebuild_pattern_table()`:**
      - **Previously a MAJOR CPU BOTTLENECK:** This function was called in *every single `step()`* and after rule eviction, clearing and repopulating `DigramTable` (or `KmerTable`) from scratch by iterating the entire current `self.sequence`.
      - **Update:** This has been significantly refactored. `rebuild_pattern_table()` is no longer called in every `step()`. Instead, it's called once at the beginning of main grammar construction methods (`build_grammar`, `build_grammar_with_gpu`), at the start of `finalize_grammar`, and when the first chunk is processed in `process_sequence_chunk`.
  2.  **Find Most Frequent Pattern:** Uses `digram_table.find_most_frequent_digram()` (or k-mer equivalent).
  3.  **Rule Creation:** If pattern count >= `min_rule_usage`, a new `Rule` is created.
  4.  **Sequence Replacement:**
      - `replace_digram_occurrences()`: Replaces digrams with the new non-terminal.
          - **Update:** This function now performs incremental updates to the `DigramTable`. It removes affected old digrams around the replacement sites and adds newly formed digrams, avoiding a full table rebuild.
      - `replace_kmer_occurrences()`:
          - **Performance Issue (Original):** Used `Vec::contains()` in a loop for checking positions, which was O(N*M).
          - **Update:** This has been **fixed**. The function now uses a `HashSet` for efficient O(1) average time complexity lookups of positions.
          - **Update:** This function now also performs incremental updates to the `KmerTable`, similar to `replace_digram_occurrences`.
          - Strand for new non-terminal from k-mer defaults to `Forward`, which might be an oversimplification.
  5.  **Rule Eviction/Inlining:**
      - `evict_least_used_rules()`: Sorts rules by usage, calls `inline_rule()` for least used ones.
          - **Update:** No longer triggers `rebuild_pattern_table()`. The `inline_rule` function now handles its own incremental updates.
      - `inline_rule()`: Replaces occurrences of a rule's non-terminal with its definition.
          - **Update:** This function now performs incremental updates to the relevant pattern table (`DigramTable` or `KmerTable`) after inlining a rule, reflecting the changes in the sequence.
      - `inline_single_use_rules()`: Similar, called at the end.

- **Streaming Mode:** `process_sequence_chunk` and `finalize_grammar` attempt to handle large sequences by chunking. The initial call to `rebuild_pattern_table` for the first chunk is now the main full build in this mode, with subsequent operations relying on incremental updates.

## 3. GPU Acceleration (`src/gpu/` & `orbweaver-kernels/`)

Aims to speed up parts of the grammar construction, primarily digram finding.

### 3.1. `GpuContext` (`src/gpu/mod.rs`)

- Manages OpenCL state: `Platform`, `Device`, `Context`, `Queue`, `Program`.
- `new()`: Initializes OpenCL, selects device (GPU preferred), loads and compiles kernels from `orbweaver-kernels`.
- **`Clone` Implementation:** Re-runs `GpuContext::new()`, meaning full OpenCL re-initialization (including kernel compilation) on each clone. This can be very slow if `GpuContext` is cloned frequently (e.g., per iteration or per GPU call).

### 3.2. `GpuSequence` (`src/gpu/digram.rs`)

- Holds sequence data as `Vec<u8>` (2-bit encoded bases) for GPU.
- `buffer: Option<Buffer<u8>>`: OpenCL buffer for the sequence on GPU.
- `upload_to_gpu()`: Transfers data from host `Vec<u8>` to GPU `Buffer<u8>`.
- `Clone` implementation only clones `data` (host `Vec<u8>`), `buffer` becomes `None`.
- `from_symbols()`: Converts `Vec<Symbol>` to `Vec<u8>`, but fails if non-terminals are present (GPU path primarily for raw terminal sequences).

### 3.3. GPU Digram Finding (`src/gpu/digram.rs` & `orbweaver-kernels/`)

- **Orchestration (`GpuSequence::find_most_frequent_digram`):**
  - If GPU context and data buffer exist, calls `find_most_frequent_digram_opencl()`.
  - Falls back to CPU (`find_most_frequent_digram_cpu()`) on failure or if GPU is not ready.
  - The CPU fallback uses `DigramTable::canonical_key` and `custom_hash` to produce a `u64 DigramKey`.

- **OpenCL Kernels (`orbweaver-kernels/src/lib.rs` -> `get_digram_kernel()`):**
  1.  **`compute_digram_hashes`:**
      - Takes `uchar *sequence`.
      - Outputs `ulong *hashes`. Each `hashes[gid]` is a 64-bit canonical key for the digram at `gid` (e.g., `(base1 << 32) | base2`, with reverse complement handling). This kernel appears sound.
  2.  **`count_digrams`:**
      - Takes `ulong *hashes` (from previous kernel).
      - Outputs `uint *counts` (histogram).
      - Uses `index = hash % counts_size` and `atomic_inc(&counts[index])`.
      - **Issue:** Relies on modulo for indexing, which can lead to collisions if `counts_size` is too small or not well-chosen relative to the distribution of 64-bit hash keys. CPU-side collision resolution would be needed for accuracy.
  3.  **`find_digrams` (Enhanced kernel):**
      - Takes `uint *symbols`.
      - Computes a 32-bit `digram_key = (s1 << 16) | s2`.
      - **MAJOR ISSUE:** Applies `digram_key = digram_key % 1000000`. This fixed, small modulo severely limits the key space, causing massive collisions and making it unable to accurately find the true most frequent digram for typical genomic sequences.
      - Stores the collided key in `digram_indices[gid]`.

- **Suffix Array Kernel (`get_suffix_array_kernel()`):**
  - The `sort_suffixes` kernel uses a basic insertion sort (O(N^2)), which is not practical for large sequences. It's likely a placeholder.

### 3.4. GPU Path Integration in `GrammarBuilder` (`builder.rs`)

- `build_grammar_with_gpu()`:
  - **MAJOR ISSUE (Data Synchronization):** The loop finds frequent digrams using `gpu_seq.find_most_frequent_digram()`. It then calls `self.replace_digram_occurrences()`, which modifies the CPU-side `self.sequence: Vec<Symbol>`. However, the `gpu_seq` (holding GPU data) is **not updated** within the loop. Subsequent iterations of `gpu_seq.find_most_frequent_digram()` operate on stale GPU data, leading to incorrect results after the first replacement.

## 4. Utility Modules

- **`src/utils/hash.rs`:**
  - `custom_hash<T: Hash>(value: &T) -> u64`: Generic hashing utility using `twox_hash::XxHash64`. Used by `GpuSequence::find_most_frequent_digram_cpu` to get a u64 key.
- **`src/encode/dna_2bit.rs` (Inferred):**
  - Provides `EncodedBase` (likely a `u8` wrapper for 0-3 A,C,G,T) and `revcomp()` functionality. Critical for memory efficiency and symbol representation.

## 5. Summary of Key Performance Issues & Optimization Opportunities

### Critical CPU Performance Bottlenecks:
1.  **Repeated Full `rebuild_pattern_table()`:** In `GrammarBuilder::step()` and after rule evictions.
    - **Status: Largely Addressed.**
    - **Fix Applied:** The `rebuild_pattern_table()` function is no longer called repeatedly within the main loop or after every eviction. It is now called once at the beginning of major processing phases. Core sequence modification functions (`replace_digram_occurrences`, `replace_kmer_occurrences`, `inline_rule`) have been refactored to perform incremental updates to the `DigramTable` and `KmerTable`.
2.  **Linear Scan for Most Frequent Pattern:** In `DigramTable::find_most_frequent_digram()` and `KmerTable` equivalent.
    - **Status: Partially Addressed (for `DigramTable` cloning).**
    - **Fix Applied (Partial):** `DigramTable::find_most_frequent_digram` was optimized to only clone the occurrence vector for the winning candidate. `KmerTable` already did this.
    - **Remaining Suggestion:** Maintain an auxiliary data structure (e.g., a priority queue or a sorted list updated incrementally) for O(1) or O(log N) retrieval of the most frequent item, avoiding the full scan.
3.  **`Vec::contains()` in `replace_kmer_occurrences()`:** Inefficient O(N*M) lookup.
    - **Status: Addressed.**
    - **Fix Applied:** Replaced with `HashSet::contains()` for efficient lookups.

### Major GPU Path Issues (Correctness & Performance):
1.  **Stale GPU Data:** `GpuSequence` data is not updated after CPU-side sequence replacements in `GrammarBuilder`'s GPU loop.
    - **Fix:** Implement data synchronization (e.g., re-upload sequence, or GPU-side replacements).
2.  **Faulty `find_digrams` Kernel:** The fixed modulo (`% 1000000`) in the `find_digrams` OpenCL kernel causes excessive collisions, making it inaccurate.
    - **Fix:** Prioritize the two-pass `compute_digram_hashes` + `count_digrams` approach with careful handling of its (potentially smaller) collision issues, or a more robust single-pass counting kernel.
3.  **`GpuContext::clone()` Overhead:** Re-initializes OpenCL, potentially slow if cloned often.
    - **Fix:** Refactor `GpuContext` cloning (e.g., `Arc<Inner>`) or usage patterns.
4.  **Basic Suffix Array Kernel:** Not viable for large sequences.
    - **Fix:** Implement or integrate a high-performance parallel suffix array algorithm.

### Other Potential Optimizations:
-   **`KmerTable` Keying/Storage:** `Vec<Symbol>` as key and for storing occurrences is memory/CPU intensive for large `k`. Consider specialized k-mer hashing or representation.
-   **Memory Management:** Profile for excessive allocations or cloning, especially in main loops and data structure manipulations.

This report should serve as a good reference for understanding the current state of the codebase and for planning future development and optimization efforts. 