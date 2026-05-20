# Orbweaver Goal

## Distinct Goal

Turn Orbweaver into a credible Rust tool for cross-species grammar motif comparison.

Success means Orbweaver can take one or more FASTA assemblies, produce stable locally consistent grammar motifs, export bounded graph artifacts, and report enough runtime/memory/statistical evidence to compare motif structure across species.

## Non-Goals

- General-purpose genome assembly
- General-purpose compression benchmark
- Whole-genome visual rendering
- Database integration
- Network data acquisition
- GPU-first architecture
- Parallel implementations of the same grammar semantics

## Primary Path

The only path that should receive new feature work is:

```bash
orbweaver \
  --input-files <assembly.fna> \
  --mmap \
  --hierarchical-merge \
  --auto-tune-lcg \
  --streaming-output \
  --stats
```

Legacy in-memory, streaming, and chunked modes may remain while tests depend on them, but they are compatibility paths.

## First Usable Artifact

For each assembly:

- `grammar.json`
- `rules.fasta`
- `grammar.graphml`
- `grammar_stats.txt`
- a motif table with rule id, content fingerprint, length, usage count, depth, assembly index, and expanded sequence preview

The motif table is now the first-class comparison artifact. Keep it stable before adding new algorithms.

## Benchmark Gate

Every serious change should preserve or improve:

- process success on a tiny fixture
- process success on a 10 MB repetitive fixture
- process success on a real 10-50 MB FASTA slice
- peak RSS
- wall time
- output rule count
- graph node/edge count
- compression ratio

No benchmark claim is accepted without a saved command and artifact.

## Cleanup Order

1. Remove fake outputs and dead modules.
2. Make CPU LCG mmap the documented baseline.
3. Keep the motif table stable and covered by tests.
4. Add a small benchmark harness that records wall time and RSS.
5. Fix assembly-index calculation for streaming/mmap output.
6. Add cross-species comparison over motif fingerprints.
7. Reconsider GPU only after profiling shows a real hot path.
