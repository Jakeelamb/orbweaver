use anyhow::{Context, Result};
use crate::grammar::engine::Grammar;
use crate::grammar::rule::Rule;
use crate::grammar::symbol::{Symbol, SymbolType, Direction};
use crate::utils::io::OutputFormat;
use crate::utils::visualization::{grammar_to_dot, DotOptions, grammar_to_gfa, grammar_to_graphml_direct};
use std::collections::HashMap;
use std::path::Path;
use std::fs::File;
use std::io::{BufWriter, Write};
use serde::Serialize;
use crate::encode::dna_2bit::EncodedBase;

/// Export grammar to a file based on format
pub fn export_grammar(
    grammar: &Grammar,
    path: &Path,
    format: OutputFormat,
) -> Result<()> {
    let file = File::create(path)
        .with_context(|| format!("Failed to create output file: {}", path.display()))?;
    let mut writer = BufWriter::new(file);
    
    match format {
        OutputFormat::Json => export_grammar_json(grammar, &mut writer),
        OutputFormat::Text => export_grammar_text(grammar, &mut writer),
        OutputFormat::Gfa => export_grammar_gfa(grammar, &mut writer),
        OutputFormat::Dot => export_grammar_dot(grammar, &mut writer),
        OutputFormat::Fasta => export_grammar_fasta(grammar, &mut writer),
    }
}

/// Export grammar to JSON format
pub fn export_grammar_json<W: Write>(grammar: &Grammar, writer: &mut W) -> Result<()> {
    let json = grammar_to_json(&grammar.sequence, &grammar.rules)?;
    writer.write_all(json.as_bytes())?;
    Ok(())
}

/// Export grammar to text format
pub fn export_grammar_text<W: Write>(grammar: &Grammar, writer: &mut W) -> Result<()> {
    // Output the final sequence
    writeln!(writer, "== Final Sequence ({} symbols) ==", grammar.sequence.len())?;
    
    let sequence_str: Vec<String> = grammar.sequence.iter()
        .map(|sym| format_symbol(sym, grammar))
        .collect();
    
    writeln!(writer, "{}", sequence_str.join(" "))?;
    writeln!(writer)?;
    
    // Output the rules
    writeln!(writer, "== Rules ({} total) ==", grammar.rules.len())?;
    
    // Sort rules by ID for consistent output
    let mut rule_ids: Vec<usize> = grammar.rules.keys().cloned().collect();
    rule_ids.sort();
    
    for rule_id in rule_ids {
        let rule = &grammar.rules[&rule_id];
        
        let symbols_str: Vec<String> = rule.symbols.iter()
            .map(|sym| format_symbol(sym, grammar))
            .collect();
        
        writeln!(
            writer,
            "R{} [Usage={}] -> {}",
            rule.id,
            rule.usage_count,
            symbols_str.join(" ")
        )?;
    }
    
    Ok(())
}

/// Export grammar to DOT format
pub fn export_grammar_dot<W: Write>(grammar: &Grammar, writer: &mut W) -> Result<()> {
    let options = DotOptions::default();
    let dot = grammar_to_dot(grammar, &options)?;
    writer.write_all(dot.as_bytes())?;
    Ok(())
}

/// Export grammar to GFA format
pub fn export_grammar_gfa<W: Write>(grammar: &Grammar, writer: &mut W) -> Result<()> {
    let gfa = grammar_to_gfa(grammar)?;
    writer.write_all(gfa.as_bytes())?;
    Ok(())
}

/// Export grammar rules as FASTA sequences
pub fn export_grammar_fasta<W: Write>(grammar: &Grammar, writer: &mut W) -> Result<()> {
    // Sort rules by ID for consistent output
    let mut rule_ids: Vec<usize> = grammar.rules.keys().cloned().collect();
    rule_ids.sort();
    
    for rule_id in rule_ids {
        let rule = &grammar.rules[&rule_id];
        
        // Expand rule to get the full sequence
        let sequence = expand_rule_to_string(rule, &grammar.rules, Direction::Forward);
        
        // Write FASTA entry
        writeln!(
            writer,
            ">Rule_{} [Usage={}] [Symbols={}]",
            rule.id,
            rule.usage_count,
            rule.symbols.len()
        )?;
        
        // Write sequence with line wrapping (80 chars per line)
        for i in 0..sequence.len() {
            write!(writer, "{}", sequence.chars().nth(i).unwrap())?;
            if (i + 1) % 80 == 0 {
                writeln!(writer)?;
            }
        }
        
        // Add final newline if not already added
        if sequence.len() % 80 != 0 {
            writeln!(writer)?;
        }
    }
    
    // Add the final sequence
    writeln!(writer, ">Final_Sequence [Length={}]", grammar.sequence.len())?;
    let final_sequence = expand_sequence_to_string(&grammar.sequence, &grammar.rules);
    
    // Write sequence with line wrapping (80 chars per line)
    for i in 0..final_sequence.len() {
        write!(writer, "{}", final_sequence.chars().nth(i).unwrap())?;
        if (i + 1) % 80 == 0 {
            writeln!(writer)?;
        }
    }
    
    // Add final newline if not already added
    if final_sequence.len() % 80 != 0 {
        writeln!(writer)?;
    }
    
    Ok(())
}

/// Format a symbol for display
pub fn format_symbol(symbol: &Symbol, _grammar: &Grammar) -> String {
    match symbol.symbol_type {
        SymbolType::Terminal(base) => {
            // Map terminal ID to a DNA base
            let base_char = match base.0 {
                0 => 'A',
                1 => 'C',
                2 => 'G',
                3 => 'T',
                _ => 'N',
            };
            
            let strand = symbol.strand.to_char();
            
            format!("{}{}", base_char, strand)
        },
        SymbolType::NonTerminal(rule_id) => {
            let strand = symbol.strand.to_char();
            
            format!("R{}{}", rule_id, strand)
        }
    }
}

/// Expands a rule to its DNA string representation
pub fn expand_rule_to_string(rule: &Rule, grammar: &HashMap<usize, Rule>, parent_strand: Direction) -> String {
    let mut result = String::new();
    
    for symbol in &rule.symbols {
        // Combine parent and symbol directions
        let effective_strand = if parent_strand == Direction::Reverse {
            // Flip the symbol's strand if parent is reversed
            symbol.strand.flip()
        } else {
            // Keep the symbol's strand if parent is forward
            symbol.strand
        };
        
        match &symbol.symbol_type {
            SymbolType::Terminal(base) => {
                let base_char = match base.0 {
                    0 => 'A',
                    1 => 'C',
                    2 => 'G',
                    3 => 'T',
                    _ => 'N',
                };
                
                if effective_strand == Direction::Reverse {
                    // Use lowercase for reverse complement
                    result.push(base_char.to_lowercase().next().unwrap());
                } else {
                    result.push(base_char);
                }
            },
            SymbolType::NonTerminal(rule_id) => {
                if let Some(nested_rule) = grammar.get(rule_id) {
                    let expanded = expand_rule_to_string(nested_rule, grammar, effective_strand);
                    result.push_str(&expanded);
                } else {
                    // Rule not found, use a placeholder
                    result.push_str(&format!("[R{}?]", rule_id));
                }
            }
        }
    }
    
    result
}

/// Expands a compressed sequence to its full representation using the grammar rules
pub fn expand_sequence_to_string(sequence: &[Symbol], rules: &HashMap<usize, Rule>) -> String {
    let mut result = String::new();
    
    for symbol in sequence {
        match &symbol.symbol_type {
            SymbolType::Terminal(base) => {
                let base_char = match base.0 {
                    0 => 'A',
                    1 => 'C',
                    2 => 'G',
                    3 => 'T',
                    _ => 'N',
                };
                
                if symbol.strand == Direction::Reverse {
                    // Use lowercase for reverse complement
                    result.push(base_char.to_lowercase().next().unwrap());
                } else {
                    result.push(base_char);
                }
            },
            SymbolType::NonTerminal(rule_id) => {
                if let Some(rule) = rules.get(rule_id) {
                    let expanded = expand_rule_to_string(rule, rules, symbol.strand);
                    result.push_str(&expanded);
                } else {
                    // Rule not found, use a placeholder
                    result.push_str(&format!("[R{}?]", rule_id));
                }
            }
        }
    }
    
    result
}

/// Calculate compression statistics for a grammar
pub fn calculate_compression_stats(grammar: &Grammar) -> CompressionStats {
    // Calculate original size (expanded sequence)
    let expanded_size = expand_sequence_to_string(&grammar.sequence, &grammar.rules).len();
    
    // Calculate compressed size
    let rule_sizes: usize = grammar.rules.values()
        .map(|rule| rule.symbols.len())
        .sum();
    
    let compressed_size = grammar.sequence.len() + rule_sizes;
    
    // Calculate compression ratio
    let compression_ratio = if compressed_size > 0 {
        expanded_size as f64 / compressed_size as f64
    } else {
        0.0
    };
    
    // Find the most used rule
    let most_used_rule = grammar.rules.values()
        .max_by_key(|rule| rule.usage_count)
        .map(|rule| (rule.id, rule.usage_count));
    
    // Calculate average rule usage
    let avg_rule_usage = if !grammar.rules.is_empty() {
        grammar.rules.values()
            .map(|rule| rule.usage_count as f64)
            .sum::<f64>() / grammar.rules.len() as f64
    } else {
        0.0
    };
    
    CompressionStats {
        original_size: expanded_size,
        compressed_size,
        compression_ratio,
        rule_count: grammar.rules.len(),
        max_depth: grammar.max_depth,
        most_used_rule,
        avg_rule_usage,
    }
}

/// Compression statistics
#[derive(Debug)]
pub struct CompressionStats {
    /// Original (expanded) size in bases
    pub original_size: usize,
    /// Compressed size (grammar representation)
    pub compressed_size: usize,
    /// Compression ratio (original / compressed)
    pub compression_ratio: f64,
    /// Number of rules
    pub rule_count: usize,
    /// Maximum recursion depth
    pub max_depth: usize,
    /// Most used rule (id, count)
    pub most_used_rule: Option<(usize, usize)>,
    /// Average rule usage
    pub avg_rule_usage: f64,
}

/// Print compression statistics
pub fn print_compression_stats(stats: &CompressionStats) -> Result<()> {
    println!("Compression Statistics:");
    println!("  Original Size: {} bases", stats.original_size);
    println!("  Compressed Size: {} symbols", stats.compressed_size);
    println!("  Compression Ratio: {:.2}x", stats.compression_ratio);
    println!("  Rule Count: {}", stats.rule_count);
    println!("  Maximum Depth: {}", stats.max_depth);
    
    if let Some((rule_id, usage)) = stats.most_used_rule {
        println!("  Most Used Rule: R{} (used {} times)", rule_id, usage);
    }
    
    println!("  Average Rule Usage: {:.2}", stats.avg_rule_usage);
    
    Ok(())
}

/// Converts a grammar to a JSON representation for export
pub fn grammar_to_json(
    sequence: &[Symbol],
    rules: &HashMap<usize, Rule>,
) -> Result<String> {
    #[derive(Serialize)]
    struct GrammarJson {
        #[serde(rename = "final_sequence")]
        sequence: Vec<SymbolJson>,
        rules: Vec<RuleJson>,
    }

    #[derive(Serialize)]
    struct SymbolJson {
        id: usize,
        #[serde(rename = "type")]
        type_: String,
        value: String,
        strand: char,
    }

    #[derive(Serialize)]
    struct RuleJson {
        id: usize,
        symbols: Vec<SymbolJson>,
        usage_count: usize,
        depth: Option<usize>,
    }

    // Helper function to convert Symbol to SymbolJson
    fn symbol_to_json(symbol: &Symbol) -> SymbolJson {
        match symbol.symbol_type {
            SymbolType::Terminal(base) => SymbolJson {
                id: symbol.id,
                type_: "terminal".to_string(),
                value: base.to_char().to_string(),
                strand: symbol.strand.to_char(),
            },
            SymbolType::NonTerminal(rule_id) => SymbolJson {
                id: symbol.id,
                type_: "non-terminal".to_string(),
                value: rule_id.to_string(),
                strand: symbol.strand.to_char(),
            },
        }
    }

    // Convert sequence
    let sequence_json: Vec<SymbolJson> = sequence.iter().map(symbol_to_json).collect();

    // Convert rules
    let mut rules_json: Vec<RuleJson> = rules
        .iter()
        .map(|(&id, rule)| RuleJson {
            id,
            symbols: rule.symbols.iter().map(symbol_to_json).collect(),
            usage_count: rule.usage_count,
            depth: rule.depth,
        })
        .collect();

    // Sort rules by ID for consistent output
    rules_json.sort_by_key(|rule| rule.id);

    // Create final structure
    let grammar_json = GrammarJson {
        sequence: sequence_json,
        rules: rules_json,
    };

    // Serialize to JSON
    serde_json::to_string_pretty(&grammar_json).context("Failed to serialize grammar to JSON")
}

/// Writes a grammar to a file in JSON format
pub fn write_grammar_json<P: AsRef<Path>>(
    path: P,
    sequence: &[Symbol],
    rules: &HashMap<usize, Rule>,
) -> Result<()> {
    let json = grammar_to_json(sequence, rules)?;
    let file = File::create(path).context("Failed to create JSON output file")?;
    let mut writer = BufWriter::new(file);
    writer.write_all(json.as_bytes())?;
    Ok(())
}

/// Convert a Symbol to its plain text representation
fn symbol_to_string(symbol: &Symbol) -> String {
    match symbol.symbol_type {
        SymbolType::Terminal(base) => format!("{}{}", base.to_char(), symbol.strand.to_char()),
        SymbolType::NonTerminal(rule_id) => format!("R{}{}", rule_id, symbol.strand.to_char()),
    }
}

/// Convert a rule to its plain text representation
fn rule_to_string(rule: &Rule) -> String {
    let mut result = String::new();
    
    for symbol in &rule.symbols {
        result.push_str(&symbol_to_string(symbol));
    }
    
    result
}

/// Writes a grammar in a simplified text format suitable for analysis
pub fn write_grammar_text<P: AsRef<Path>>(
    path: P,
    sequence: &[Symbol],
    rules: &HashMap<usize, Rule>,
) -> Result<()> {
    let file = File::create(path).context("Failed to create text output file")?;
    let mut writer = BufWriter::new(file);
    
    // Write the sequence
    writeln!(&mut writer, "== Sequence ==")?;
    for symbol in sequence {
        write!(&mut writer, "{}", symbol_to_string(symbol))?;
    }
    writeln!(&mut writer)?;
    writeln!(&mut writer)?;
    
    // Write the rules sorted by ID
    writeln!(&mut writer, "== Rules ==")?;
    let mut rule_ids: Vec<usize> = rules.keys().cloned().collect();
    rule_ids.sort();
    
    for &rule_id in &rule_ids {
        if let Some(rule) = rules.get(&rule_id) {
            writeln!(
                &mut writer,
                "Rule {}: {} (usage: {})",
                rule_id,
                rule_to_string(rule),
                rule.usage_count
            )?;
        }
    }
    
    Ok(())
}

/// Extracts the underlying DNA sequence from a grammar
pub fn extract_dna_sequence(sequence: &[Symbol], rules: &HashMap<usize, Rule>) -> Vec<(EncodedBase, char)> {
    fn expand_symbol(
        symbol: &Symbol,
        rules: &HashMap<usize, Rule>,
        result: &mut Vec<(EncodedBase, char)>,
        depth: usize,
    ) {
        // Prevent infinite recursion
        if depth > 1000 {
            return;
        }
        
        match &symbol.symbol_type {
            SymbolType::Terminal(base) => {
                // For terminals, just add the base and strand
                result.push((*base, symbol.strand.to_char()));
            }
            SymbolType::NonTerminal(rule_id) => {
                // For non-terminals, expand the rule
                if let Some(rule) = rules.get(rule_id) {
                    for sym in &rule.symbols {
                        // Create a new symbol with the strand combined with the parent strand
                        let effective_strand = if symbol.strand == Direction::Reverse {
                            // Flip the strand of the child when parent is '-'
                            sym.strand.flip()
                        } else {
                            // Keep the strand of the child when parent is '+'
                            sym.strand
                        };
                        
                        let mut effective_symbol = *sym; // Dereference sym to create a mutable copy
                        effective_symbol.strand = effective_strand;
                        
                        // Recursively expand this symbol
                        expand_symbol(&effective_symbol, rules, result, depth + 1);
                    }
                }
            }
        }
    }
    
    let mut result = Vec::new();
    
    for symbol in sequence {
        expand_symbol(symbol, rules, &mut result, 0);
    }
    
    result
}

/// Writes the expanded DNA sequence to a file in FASTA format
pub fn write_dna_fasta<P: AsRef<Path>>(
    path: P,
    sequence: &[Symbol],
    rules: &HashMap<usize, Rule>,
    header: &str,
) -> Result<()> {
    let dna_sequence = extract_dna_sequence(sequence, rules);
    
    let file = File::create(path).context("Failed to create FASTA output file")?;
    let mut writer = BufWriter::new(file);
    
    // Write the FASTA header
    writeln!(&mut writer, ">{}", header)?;
    
    // Write the sequence with line wrapping
    let line_width = 80;
    let mut line_count = 0;
    
    for (base, _strand) in &dna_sequence {
        // For FASTA, we don't include strand information directly
        write!(&mut writer, "{}", base.to_char())?;
        
        line_count += 1;
        if line_count >= line_width {
            writeln!(&mut writer)?;
            line_count = 0;
        }
    }
    
    // Ensure the file ends with a newline
    if line_count > 0 {
        writeln!(&mut writer)?;
    }
    
    Ok(())
}

/// Writes a tabular summary of identified repeats (rules) to a file.
pub fn write_repeat_summary(path: &Path, grammar: &Grammar) -> Result<()> {
    let file = File::create(path)
        .with_context(|| format!("Failed to create repeat summary file: {}", path.display()))?;
    let mut writer = BufWriter::new(file);

    // Write header
    writeln!(writer, "RuleID\tUsageCount\tExpandedLength\tSequence(Preview)")?;

    // Identify repeats (e.g., rules with usage >= 2)
    let mut repeats: Vec<_> = grammar.rules.values()
        .filter(|rule| rule.usage_count >= 2) // Simple definition of repeat for now
        .collect();

    // Sort repeats by usage count (descending)
    repeats.sort_by(|a, b| b.usage_count.cmp(&a.usage_count));

    // Write data for each repeat
    for rule in repeats {
        let expanded_sequence = expand_rule_to_string(rule, &grammar.rules, Direction::Forward);
        let preview_len = 50; // Show first 50 bases as preview
        let sequence_preview = if expanded_sequence.len() > preview_len {
            format!("{}...", &expanded_sequence[..preview_len])
        } else {
            expanded_sequence.clone()
        };

        writeln!(writer, "{}\t{}\t{}\t{}",
            rule.id,
            rule.usage_count,
            expanded_sequence.len(),
            sequence_preview
        )?;
    }

    println!("Successfully wrote repeat summary.");
    Ok(())
}

/// Export grammar to GraphML format
pub fn write_grammar_graphml(path: &Path, grammar: &Grammar, dot_options: &DotOptions) -> Result<()> {
    // Use the direct GraphML generation function
    let graphml_string = grammar_to_graphml_direct(grammar, dot_options)
        .map_err(|e| anyhow::anyhow!("Failed to generate GraphML content: {}", e))?;

    let file = File::create(path)
        .map_err(|e| anyhow::anyhow!("Failed to create GraphML file '{}': {}", path.display(), e))?;
    let mut writer = BufWriter::new(file);
    writer.write_all(graphml_string.as_bytes())
        .map_err(|e| anyhow::anyhow!("Failed to write to GraphML file '{}': {}", path.display(), e))?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::grammar::rule::Rule;
    use crate::grammar::symbol::{Symbol, SymbolType, Direction};
    use crate::encode::dna_2bit::EncodedBase;
    use std::collections::HashMap;

    fn create_test_grammar() -> Grammar {
        let mut rules = HashMap::new();
        let rule0 = Rule {
            id: 0,
            symbols: vec![
                Symbol::terminal(0, EncodedBase(0), Direction::Forward),
                Symbol::terminal(1, EncodedBase(1), Direction::Forward),
            ],
            usage_count: 4,
            positions: vec![0, 2, 4, 6],
            depth: Some(1),
        };
        let rule1 = Rule {
            id: 1,
            symbols: vec![
                Symbol::non_terminal(2, 0, Direction::Forward),
                Symbol::terminal(3, EncodedBase(2), Direction::Forward),
            ],
            usage_count: 2,
            positions: vec![1, 5],
            depth: Some(2),
        };
        rules.insert(0, rule0);
        rules.insert(1, rule1);
        Grammar {
            sequence: vec![
                Symbol::non_terminal(10, 0, Direction::Forward),
                Symbol::non_terminal(11, 1, Direction::Forward),
            ],
            rules,
            max_depth: 2,
            origins: HashMap::new(),
        }
    }

    #[test]
    fn test_expand_rule_to_string() {
        let grammar = create_test_grammar();
        let rule = grammar.rules.get(&1).unwrap();
        let expanded = expand_rule_to_string(rule, &grammar.rules, Direction::Forward);
        // Rule 1 expands to Rule 0 + G, Rule 0 is AC
        assert_eq!(expanded, "ACG");
    }

    #[test]
    fn test_expand_sequence_to_string() {
        let grammar = create_test_grammar();
        let expanded = expand_sequence_to_string(&grammar.sequence, &grammar.rules);
        // Sequence is [R0, R1], R0 expands to AC, R1 expands to ACG
        assert_eq!(expanded, "ACACG");
    }

    #[test]
    fn test_calculate_compression_stats() {
        let grammar = create_test_grammar();
        let stats = calculate_compression_stats(&grammar);
        assert_eq!(stats.original_size, 5); // AC + ACG = 5
        // Compressed size = final_sequence.len + sum(rule_symbol_lens)
        // final_sequence = [R0, R1] -> len 2
        // Rule 0 symbols = [A, C] -> len 2
        // Rule 1 symbols = [R0, G] -> len 2
        // Compressed = 2 + 2 + 2 = 6
        assert_eq!(stats.compressed_size, 6); // Corrected assertion from 2 to 6
        assert_eq!(stats.rule_count, 2);
        assert_eq!(stats.max_depth, 2);
    }

    #[test]
    fn test_export_grammar_text() {
        let grammar = create_test_grammar();
        let mut buf = Vec::new();
        export_grammar_text(&grammar, &mut buf).unwrap();
        let text = String::from_utf8(buf).unwrap();
        assert!(text.contains("== Final Sequence"));
        assert!(text.contains("R0"));
        assert!(text.contains("R1"));
    }
} 