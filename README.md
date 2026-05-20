# Orbweaver

Orbweaver is a Rust research tool for building compact, grammar-derived motif graphs from genomic FASTA files.

The narrowed goal is not general compression, not a genome browser, and not a GPU demo. The goal is:

> Build stable grammar motifs from assemblies, export bounded graph artifacts, and compare those motif graphs across species with honest memory and runtime measurements.

## Current Focus

The production research path is:

```bash
cargo run --release -- \
  --input-files genome.fna \
  --mmap \
  --hierarchical-merge \
  --auto-tune-lcg \
  --streaming-output \
  --stats
```

This path uses memory-mapped FASTA input, locally consistent grammar parsing, chunk overlap handling, hierarchical merge, and streaming JSON output.

GPU support exists, but is experimental and opt-in with `--gpu`. It is not the baseline until an end-to-end benchmark proves it improves the focused path.

## Outputs

Each run writes a run directory containing:

- `grammar.json`: compressed final sequence plus grammar rules
- `rules.fasta`: expanded sequence for each grammar rule
- `grammar.dot`: Graphviz grammar graph
- `grammar.graphml`: graph exchange format
- `grammar.gfa`: GFA representation
- `motif_table.tsv`: rule-level motif table for comparison
- `grammar_stats.txt`: compression and assembly-index statistics when `--stats` is set
- `run_metadata.json`: run arguments, status, and checkpoints

Graph output is for bounded motif inspection. Do not try to render whole-genome terminal-heavy graphs and call that visualization. Compare stable rule/motif fingerprints first; visualize selected neighborhoods later.

`motif_table.tsv` is the main comparison surface. It contains `rule_id`, canonical sequence `fingerprint`, expanded `length`, `usage_count`, `depth`, `assembly_index`, and `sequence_preview`.

## Modes

Supported research path:

- `--mmap --hierarchical-merge --auto-tune-lcg --streaming-output`

Compatibility paths:

- default in-memory mode for tiny fixtures and tests
- `--streaming` legacy chunk streaming
- `--chunk-size` / `--adaptive-chunking` legacy parallel Sequitur path

The compatibility paths should not receive new features unless they are needed to preserve tests while migrating behavior into the LCG mmap path.

## Development Baseline

```bash
cargo check
cargo test --lib
cargo test --test integration_tests
cargo test --test gpu_digram_tests
```

Use the benchmark harness for real memory and artifact checks:

```bash
cargo build --release --bins
target/release/orbweaver-bench \
  --input genome.fna \
  --output-dir benchmark_results \
  --run-id bacteria_baseline
```

The harness runs the focused mmap LCG path, records wall time, kernel-reported peak RSS, output artifact sizes, motif row count, and GraphML node/edge counts into `benchmark_summary.json` and `benchmark_summary.tsv`.

Large-sequence tests are a benchmark gate, not a correctness oracle. They should assert process success, output existence, rule counts, runtime, and memory/RSS through this harness.

## Near-Term Engineering Rules

- One serious code path: mmap + LCG + hierarchical merge.
- CPU first. GPU only when measured end-to-end.
- Delete fake outputs instead of shipping placeholders.
- Compare motifs by content fingerprints, not by raw node IDs.
- Keep memory claims tied to benchmark artifacts.
- Prefer bacterial and 10-50 MB slices until the baseline is stable.

## Useful Flags

```bash
# Focused path, CPU baseline
orbweaver -i genome.fna --mmap --hierarchical-merge --auto-tune-lcg --streaming-output --stats

# Smaller chunks for lower peak memory
orbweaver -i genome.fna --mmap --mmap-chunk-size 1000000 --hierarchical-merge --streaming-output

# Experimental GPU path
orbweaver -i genome.fna --mmap --hierarchical-merge --gpu
```
