# Orbweaver CLI

This document reflects the narrowed research path. Older in-memory, streaming, and chunked modes still exist for compatibility, but new work should target the mmap LCG path.

## Focused Command

```bash
orbweaver \
  --input-files genome.fna \
  --mmap \
  --hierarchical-merge \
  --auto-tune-lcg \
  --streaming-output \
  --stats
```

## Required Input

| Option | Meaning |
| --- | --- |
| `-i, --input-files <FILES>` | One or more FASTA files, comma-separated. |

## Focused Path Options

| Option | Meaning |
| --- | --- |
| `--mmap` | Memory-map FASTA input and process it in chunks. |
| `--mmap-chunk-size <BASES>` | Chunk size for mmap processing. Default: `10000000`. |
| `--chunk-overlap <BASES>` | Boundary overlap between chunks. Default: `1000`. |
| `--hierarchical-merge` | Merge chunk grammars as a tree to reduce peak grammar fan-in. |
| `--auto-tune-lcg` | Tune LCG phrase parameters from a sequence sample. |
| `--streaming-output` | Write JSON incrementally. |
| `--streaming-output-flush-interval <N>` | Rule flush interval for streaming JSON. Default: `1000`. |
| `--stats` | Write grammar statistics and assembly-index output. |

## Output Options

| Option | Meaning |
| --- | --- |
| `--run-id <ID>` | Run directory name. Defaults to timestamp. |
| `-j, --output-json <FILE>` | JSON filename inside the run directory. Defaults to `grammar.json`. |
| `--output-text <FILE>` | Optional human-readable grammar text. |
| `--output-gfa <FILE>` | GFA filename. Defaults to `grammar.gfa`. |

Orbweaver always writes `grammar.json`, `rules.fasta`, `grammar.dot`, `grammar.graphml`, `grammar.gfa`, `motif_table.tsv`, and `run_metadata.json` into the run directory unless a specific output filename overrides the default for that format.

## Experimental GPU

GPU support is opt-in:

```bash
orbweaver -i genome.fna --mmap --hierarchical-merge --gpu
```

Do not treat GPU as the baseline until an end-to-end benchmark shows better wall time or memory behavior for the focused path.

## Compatibility Modes

These modes remain for tests and old scripts:

| Option | Status |
| --- | --- |
| default in-memory mode | Tiny fixtures only. |
| `--streaming` | Legacy chunk streaming. |
| `--chunk-size` / `--adaptive-chunking` | Legacy parallel Sequitur path. |
| `--use-lcg` without `--mmap` | Transitional path. Prefer `--mmap`. |

## Practical Examples

```bash
# CPU baseline for real work
orbweaver -i genome.fna --mmap --hierarchical-merge --auto-tune-lcg --streaming-output --stats

# Lower memory pressure
orbweaver -i genome.fna --mmap --mmap-chunk-size 1000000 --hierarchical-merge --streaming-output

# Tiny fixture sanity check
orbweaver -i tiny.fasta --stats
```

## Benchmark Harness

Build both binaries:

```bash
cargo build --release --bins
```

Run the focused path with measured artifacts:

```bash
target/release/orbweaver-bench \
  --input genome.fna \
  --output-dir benchmark_results \
  --run-id bacteria_baseline
```

The harness invokes `orbweaver` with `--mmap --hierarchical-merge --auto-tune-lcg --streaming-output --stats`. It writes `benchmark_summary.json` and `benchmark_summary.tsv` into the run directory with wall time, peak RSS, grammar stats, motif row count, and graph size.
