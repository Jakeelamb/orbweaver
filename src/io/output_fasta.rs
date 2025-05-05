use crate::grammar::builder::GrammarBuilder;
use crate::grammar::rule::Rule;
use crate::grammar::symbol::SymbolType;
use anyhow::{Context, Result};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

/// Recursively reconstructs the base sequence of a rule, handling potential cycles.
/// Needs access to the full rule map.
fn reconstruct_rule_sequence_fasta(
    rule_id: usize,
    all_rules: &HashMap<usize, Rule>,
    visited_stack: &mut std::collections::HashSet<usize>,
) -> Result<String> {
    if !visited_stack.insert(rule_id) {
        return Err(anyhow::anyhow!(
            "Cycle detected during FASTA sequence reconstruction involving rule {}",
            rule_id
        ));
    }

    let rule = all_rules.get(&rule_id).ok_or_else(|| {
        anyhow::anyhow!(
            "Rule ID {} not found during FASTA sequence reconstruction",
            rule_id
        )
    })?;

    let mut sequence = String::new();
    for &symbol in &rule.symbols {
        match symbol.symbol_type {
            SymbolType::Terminal(base) => {
                sequence.push(base as char);
            }
            SymbolType::NonTerminal(sub_rule_id) => {
                let sub_sequence = reconstruct_rule_sequence_fasta(sub_rule_id, all_rules, visited_stack)?;
                sequence.push_str(&sub_sequence);
            }
        }
    }

    visited_stack.remove(&rule_id); // Remove after successfully processing children
    Ok(sequence)
}

/// Writes the rules (blocks) of the grammar to a FASTA file.
/// Each rule's fully expanded sequence is written as a FASTA record.
///
/// Args:
///     grammar_builder: The GrammarBuilder instance after build_grammar() has run.
///     output_path: The path to the output FASTA file.
pub fn export_rules_as_fasta(grammar_builder: &GrammarBuilder, output_path: &Path) -> Result<()> {
    println!("Exporting grammar rules as FASTA: {}", output_path.display());

    let (_final_sequence, rules) = grammar_builder.get_grammar();

    let file = File::create(output_path)
        .with_context(|| format!("Failed to create FASTA export file: {}", output_path.display()))?;
    let mut writer = BufWriter::new(file);

    // Sort rules by ID for consistent output
    let mut sorted_rule_ids: Vec<_> = rules.keys().copied().collect();
    sorted_rule_ids.sort_unstable();

    for rule_id in sorted_rule_ids {
        let rule = &rules[&rule_id];
        let mut visited = std::collections::HashSet::new(); // For cycle detection
        match reconstruct_rule_sequence_fasta(rule.id, rules, &mut visited) {
            Ok(sequence) => {
                writeln!(writer, ">Rule_{} [Usage={}]", rule.id, rule.usage_count)
                    .with_context(|| format!("Failed to write FASTA header for rule {}", rule.id))?;
                
                // Write sequence, wrapping lines
                let max_line_len = 70;
                for chunk in sequence.as_bytes().chunks(max_line_len) {
                    writeln!(writer, "{}", std::str::from_utf8(chunk).unwrap_or("ERROR"))
                        .with_context(|| format!("Failed to write FASTA sequence for rule {}", rule.id))?;
                }
            }
            Err(e) => {
                eprintln!(
                    "Warning: Skipping FASTA export for rule {} due to error: {}",
                    rule.id,
                    e
                );
                // Optionally write an error marker in the FASTA file
                writeln!(writer, ">Rule_{}_ERROR [Usage={}]", rule.id, rule.usage_count)
                    .context("Failed to write FASTA error header")?;
                writeln!(writer, "Error reconstructing sequence: {}", e)
                    .context("Failed to write FASTA error message")?;
            }
        }
    }

    println!("Successfully exported rules to FASTA.");
    Ok(())
} 