# Memory Strategy

The memory strategy is now intentionally narrow: make the mmap LCG path work, measure it, then optimize only the measured bottlenecks.

## Baseline Path

```bash
orbweaver \
  --input-files genome.fna \
  --mmap \
  --hierarchical-merge \
  --auto-tune-lcg \
  --streaming-output \
  --stats
```

## Mechanisms

- Memory-mapped FASTA input avoids loading the full file into owned memory.
- 2-bit base representation reduces sequence working-set size.
- Locally consistent grammar parsing creates content-stable phrases.
- Chunk overlap gives boundary-spanning motifs a chance to survive.
- Hierarchical merge avoids holding every chunk grammar in one flat merge set.
- Streaming JSON reduces output memory spikes.

## Compatibility Mechanisms

These exist but are not the research baseline:

- legacy `--streaming`
- legacy `--chunk-size`
- legacy `--adaptive-chunking`
- `--max-rule-count` rule eviction
- opt-in `--gpu`

## Measurement Rules

Do not claim memory wins without saving:

- exact command
- input file or fixture name
- wall time
- peak RSS
- rule count
- compressed sequence length
- output artifact paths

The next benchmark harness should record wall time and RSS for a tiny fixture, a 10 MB repetitive fixture, and one real 10-50 MB FASTA slice.

