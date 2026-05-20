# Orbweaver Output Formats

Orbweaver writes run artifacts into the current working directory under the selected `--run-id` or a timestamp.

## `grammar.json`

JSON serialization of the compressed sequence and grammar rules.

Core fields:

- `final_sequence` / `sequence`: symbols remaining after grammar construction
- `rules`: rule id to rule definition
- `max_depth`: maximum known grammar depth

Symbols are terminals or non-terminals. Rules contain their expansion symbols, usage count, optional depth, and optional assembly index.

## `rules.fasta`

Expanded sequence for every grammar rule.

Use this for direct motif inspection and for building the planned motif table.

## `grammar.dot`

Graphviz DOT representation of the grammar.

This is for bounded motif neighborhoods, not whole-genome terminal-heavy rendering.

## `grammar.graphml`

GraphML representation of the grammar for graph tools.

This is the preferred graph exchange format for downstream motif graph comparison.

## `grammar.gfa`

GFA representation of rule segments and symbol adjacencies.

## `grammar_stats.txt`

Written when `--stats` is set. Includes rule count, compressed sequence length, total symbols in rule definitions, depth, and assembly-index output when available.

## `motif_table.tsv`

Rule-level motif table for cross-species comparison:

```tsv
rule_id	fingerprint	length	usage_count	depth	assembly_index	sequence_preview
```

`fingerprint` is a canonical content hash over the expanded motif sequence and its reverse complement. This is the first field to join on across species.

Do not reintroduce placeholder repeat output. This table must be derived from real grammar rules.
