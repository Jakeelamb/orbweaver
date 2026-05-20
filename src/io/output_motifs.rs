use crate::grammar::engine::Grammar;
use crate::utils::export::expand_rule_to_string;
use crate::utils::hash::hash_sequence;
use anyhow::{Context, Result};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

const PREVIEW_LEN: usize = 80;

#[derive(Debug, Clone, Eq, PartialEq)]
struct MotifRow {
    rule_id: usize,
    fingerprint: u64,
    length: usize,
    usage_count: usize,
    depth: usize,
    assembly_index: Option<usize>,
    sequence_preview: String,
}

pub fn write_motif_table(path: &Path, grammar: &Grammar) -> Result<()> {
    let mut rows = Vec::with_capacity(grammar.rules.len());
    for rule in grammar.rules.values() {
        let expanded = expand_rule_to_string(rule, &grammar.rules, crate::grammar::symbol::Direction::Forward);
        let fingerprint = canonical_sequence_fingerprint(expanded.as_bytes());
        let preview = expanded.chars().take(PREVIEW_LEN).collect();

        rows.push(MotifRow {
            rule_id: rule.id,
            fingerprint,
            length: expanded.len(),
            usage_count: rule.usage_count,
            depth: rule.depth.unwrap_or(0),
            assembly_index: rule.assembly_index,
            sequence_preview: preview,
        });
    }

    rows.sort_by(|a, b| {
        b.usage_count
            .cmp(&a.usage_count)
            .then_with(|| b.length.cmp(&a.length))
            .then_with(|| a.rule_id.cmp(&b.rule_id))
    });

    let file = File::create(path)
        .with_context(|| format!("failed to create motif table: {}", path.display()))?;
    let mut writer = BufWriter::new(file);

    writeln!(
        writer,
        "rule_id\tfingerprint\tlength\tusage_count\tdepth\tassembly_index\tsequence_preview"
    )?;

    for row in rows {
        writeln!(
            writer,
            "{}\t{:016x}\t{}\t{}\t{}\t{}\t{}",
            row.rule_id,
            row.fingerprint,
            row.length,
            row.usage_count,
            row.depth,
            row.assembly_index
                .map(|v| v.to_string())
                .unwrap_or_else(|| "NA".to_string()),
            row.sequence_preview
        )?;
    }

    Ok(())
}

fn canonical_sequence_fingerprint(sequence: &[u8]) -> u64 {
    let forward = hash_sequence(sequence);
    let reverse_complement = reverse_complement_ascii(sequence);
    let reverse = hash_sequence(&reverse_complement);
    forward.min(reverse)
}

fn reverse_complement_ascii(sequence: &[u8]) -> Vec<u8> {
    sequence
        .iter()
        .rev()
        .map(|base| match base {
            b'A' | b'a' => b'T',
            b'C' | b'c' => b'G',
            b'G' | b'g' => b'C',
            b'T' | b't' => b'A',
            _ => b'N',
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::analysis::assembly_index::calculate_rule_assembly_indices;
    use crate::encode::dna_2bit::EncodedBase;
    use crate::grammar::rule::Rule;
    use crate::grammar::symbol::{Direction, Symbol};
    use std::collections::HashMap;

    #[test]
    fn canonical_fingerprint_matches_reverse_complement() {
        let forward = canonical_sequence_fingerprint(b"ACGTACG");
        let reverse = canonical_sequence_fingerprint(b"CGTACGT");
        assert_eq!(forward, reverse);
    }

    #[test]
    fn motif_table_contains_real_rows() {
        let mut rules = HashMap::new();
        rules.insert(
            0,
            Rule {
                id: 0,
                symbols: vec![
                    Symbol::terminal(0, EncodedBase(0), Direction::Forward, None, None),
                    Symbol::terminal(1, EncodedBase(1), Direction::Forward, None, None),
                ],
                usage_count: 3,
                positions: vec![],
                depth: Some(1),
                assembly_index: None,
            },
        );
        let mut grammar = Grammar {
            sequence: vec![Symbol::non_terminal(0, 0, Direction::Forward)],
            rules,
            max_depth: 1,
            origins: HashMap::new(),
        };
        calculate_rule_assembly_indices(&mut grammar.rules).unwrap();

        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("motif_table.tsv");
        write_motif_table(&path, &grammar).unwrap();
        let contents = std::fs::read_to_string(path).unwrap();

        assert!(contents.contains("rule_id\tfingerprint\tlength"));
        assert!(contents.contains("0\t"));
        assert!(contents.contains("\t2\t3\t1\t1\tAC"));
    }
}
