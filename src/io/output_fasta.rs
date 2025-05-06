use crate::grammar::engine::Grammar;
use crate::grammar::symbol::Direction;
use crate::utils::export::expand_rule_to_string;
use anyhow::{Context, Result};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

/// Exports all rules as FASTA sequences.
/// Each rule becomes a separate FASTA entry.
pub fn write_grammar_fasta(path: &Path, grammar: &Grammar) -> Result<()> {
    let file = File::create(path)
        .context(format!("Failed to create FASTA file: {}", path.display()))?;
    let mut writer = BufWriter::new(file);
    
    // Get the sequence and rules from the grammar
    let sequence = &grammar.sequence;
    let rules = &grammar.rules;
    
    // Export each rule as a separate FASTA entry
    for (rule_id, rule) in rules {
        // Generate the header line
        let header = format!(">Rule_{} usage={} depth={}\n", 
            rule_id, 
            rule.usage_count,
            rule.depth.unwrap_or(0));
        writer.write_all(header.as_bytes())?;
        
        // Expand the rule to its literal DNA representation
        let expanded = expand_rule_to_string(rule, rules, Direction::Forward);
        
        // Write the sequence in FASTA format (80 chars per line)
        for chunk in expanded.as_bytes().chunks(80) {
            writer.write_all(chunk)?;
            writer.write_all(b"\n")?;
        }
    }
    
    println!("Successfully wrote {} rules to FASTA file", rules.len());
    Ok(())
} 