# Orbweaver Output Formats

Orbweaver can output the generated grammar in several formats using command-line flags.

## 1. JSON (`--output-json <PATH>`)

This format provides a direct serialization of the internal grammar representation using `serde_json`.

```json
{
  "final_sequence": [
    {
      "id": 10,
      "symbol_type": {
        "NonTerminal": 0
      },
      "strand": "+"
    },
    {
      "id": 11,
      "symbol_type": {
        "NonTerminal": 0
      },
      "strand": "+"
    }
    // ... more symbols in the final compressed sequence
  ],
  "rules": {
    "0": {
      "id": 0,
      "symbols": [
        {
          "id": 0,
          "symbol_type": {
            "Terminal": 65 // ASCII for 'A'
          },
          "strand": "+"
        },
        {
          "id": 1,
          "symbol_type": {
            "Terminal": 67 // ASCII for 'C'
          },
          "strand": "+"
        }
      ],
      "usage_count": 4
    },
    "1": {
       // ... definition for Rule 1 ...
    }
    // ... more rules 
  }
}
```

*   `final_sequence`: An array of `Symbol` objects representing the compressed sequence after all rule replacements.
*   `rules`: A map where keys are rule IDs (as strings) and values are `Rule` objects.
*   `Symbol`: Contains the instance `id` (its original position or ID context), `symbol_type` (`Terminal` with ASCII value or `NonTerminal` with rule ID), and `strand` (`+` or `-`).
*   `Rule`: Contains the rule `id`, the `symbols` it expands to, and its `usage_count`.

## 2. GFAv1 (`--output-gfa <PATH>`)

Graphical Fragment Assembly format, version 1. This output is **experimental** and aims to represent the grammar structure.

*   **H Line**: Standard GFA header (`H\tVN:Z:1.0`).
*   **S Lines**: Each rule is represented as a segment (`S`).
    *   `<sid>`: The rule ID (e.g., `0`, `1`).
    *   `<sequence>`: The fully reconstructed DNA sequence of the rule.
    *   `LN:i:<len>`: Optional tag indicating the length of the reconstructed sequence.
    *   Example: `S\t0\tACGT\tLN:i:4`
*   **L Lines**: Links represent the adjacency of symbols within a rule definition (`Rule R -> S1 S2`).
    *   `<sid1>`: Segment ID of the first symbol (e.g., `R0`, `T_A+`).
    *   `<orient1>`: Strand of the first symbol (`+` or `-`).
    *   `<sid2>`: Segment ID of the second symbol.
    *   `<orient2>`: Strand of the second symbol.
    *   `<overlap>`: `0M` (zero overlap, indicating adjacency).
    *   Example: `L\tR0\t+\tR1\t-\t0M` (Link from rule R0 to rule R1 with opposite strand)
    *   Example: `L\tR0\t+\tT_C+\t+\t0M` (Link from rule R0 to terminal C+)

**Note:** The interpretation of grammar rules in GFA can vary. This representation focuses on showing the components of each rule. Visualizers may or may not display this intuitively as a grammar hierarchy.

## 3. Text (`--output-text <PATH>`)

A human-readable text format describing the final sequence and rules.

```text
== Final Sequence (3 symbols) ==
R0+ R1- R0+

== Rules (2 total) ==
R0 [Usage=2] -> A+ C+
R1 [Usage=1] -> R0+ G-
```

*   **Final Sequence**: Shows the sequence of symbols (Terminals like `A+`, NonTerminals like `R0-`) remaining after compression. Line wrapping may occur.
*   **Rules**: Lists each rule definition.
    *   `R<id>`: The rule identifier.
    *   `[Usage=<count>]`: The number of times this rule was used (replaced digrams).
    *   `-> Symbol1 Symbol2 ...`: The sequence of symbols the rule expands to.

## 4. FASTA (`--export-blocks <PATH>`)

Exports the fully expanded DNA sequence for each rule in standard FASTA format.

```fasta
>Rule_0 [Usage=4]
ACGTACGT
>Rule_1 [Usage=2]
TTGC
>Rule_2_ERROR [Usage=1]
Error reconstructing sequence: Cycle detected involving rule 2
```

*   Each record corresponds to one rule.
*   The header line includes `Rule_<id>` and the `Usage=<count>`.
*   The sequence is the fully expanded DNA sequence for that rule.
*   If an error occurs during reconstruction (e.g., a cycle is detected), an error record is written instead. 