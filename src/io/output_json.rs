use crate::grammar::symbol::{Symbol, SymbolType};
use crate::grammar::rule::Rule;
use crate::grammar::engine::Grammar;
use anyhow::{Context, Result};
use std::collections::HashMap;
use std::path::Path;
use std::fs::File;
use std::io::{BufWriter, Write};
use serde_json::{json, Value};

/// Writes the grammar to a JSON file.
/// Updated signature to accept sequence and rules directly
pub fn write_grammar_json(path: &Path, sequence: &[Symbol], rules: &HashMap<usize, Rule>) -> Result<()> {
    println!("Writing grammar to JSON: {}", path.display());
    
    let file = File::create(path)
        .context(format!("Failed to create file: {}", path.display()))?;
    let mut writer = BufWriter::new(file);
    
    // Convert the grammar to JSON
    let json_value = grammar_to_json(sequence, rules);
    
    // Serialize to pretty-printed JSON
    let json_str = serde_json::to_string_pretty(&json_value)
        .context("Failed to serialize grammar to JSON")?;
    
    // Write to file
    writer.write_all(json_str.as_bytes())
        .context("Failed to write to file")?;
    
    println!("Successfully wrote grammar to JSON.");
    Ok(())
}

/// Overload that accepts a Grammar struct directly
pub fn write_grammar_json_from_grammar(path: &Path, grammar: &Grammar) -> Result<()> {
    write_grammar_json(path, &grammar.sequence, &grammar.rules)
}

/// Convert the grammar (sequence and rules) to a JSON value.
fn grammar_to_json(sequence: &[Symbol], rules: &HashMap<usize, Rule>) -> Value {
    // Create a JSON object for the grammar
    let mut json_grammar = json!({
        "final_sequence": serialize_sequence(sequence),
        "rules": {},
    });

    // Add each rule to the JSON object
    if let Some(rules_obj) = json_grammar["rules"].as_object_mut() {
        for (rule_id, rule) in rules {
            let rule_json = json!({
                "symbols": serialize_sequence(&rule.symbols),
                "usage_count": rule.usage_count,
                "depth": rule.depth.unwrap_or(0)
            });
            rules_obj.insert(rule_id.to_string(), rule_json);
        }
    }

    // Extract and clone the sequence array to avoid borrow issues
    let sequence_array = json_grammar["final_sequence"]
        .as_array()
        .cloned()
        .unwrap_or_default();

    // Validate that all rule IDs in the sequence exist in the rules map
    // If they don't exist, add placeholder empty rules to avoid broken references
    if let Some(rules_obj) = json_grammar["rules"].as_object_mut() {
        for symbol_value in sequence_array {
            if let Some(symbol_map) = symbol_value.as_object() {
                if symbol_map.get("type").and_then(Value::as_str) == Some("non_terminal") {
                    if let Some(rule_id_val) = symbol_map.get("rule_id") {
                        if let Some(rule_id_num) = rule_id_val.as_u64() {
                            let rule_id_str = rule_id_num.to_string();
                            if !rules_obj.contains_key(&rule_id_str) {
                                // Add a placeholder rule with empty symbols
                                rules_obj.insert(rule_id_str.clone(), json!({
                                    "symbols": [],
                                    "usage_count": 1,
                                    "depth": 0
                                }));
                                eprintln!("Warning: Added placeholder for missing rule ID {} referenced in sequence", rule_id_str);
                            }
                        }
                    }
                }
            }
        }
    }

    json_grammar
}

/// Serialize a sequence of symbols to a JSON array.
fn serialize_sequence(sequence: &[Symbol]) -> Value {
    let mut json_sequence = Vec::new();

    for symbol in sequence {
        let json_symbol = match &symbol.symbol_type {
            SymbolType::Terminal(base) => {
                let base_char = match base.0 {
                    0 => 'A',
                    1 => 'C',
                    2 => 'G',
                    3 => 'T',
                    _ => 'N',
                };

                let strand_char = symbol.strand.to_char();

                json!({
                    "type": "terminal",
                    "value": format!("{}{}", base_char, strand_char)
                })
            },
            SymbolType::NonTerminal(rule_id) => {
                let strand_char = symbol.strand.to_char();

                json!({
                    "type": "non_terminal",
                    "rule_id": rule_id,
                    "strand": format!("{}", strand_char)
                })
            }
        };

        json_sequence.push(json_symbol);
    }

    Value::Array(json_sequence)
}

/// Streaming JSON writer that writes rules incrementally to reduce memory usage.
/// This is useful for very large grammars where building the full JSON structure
/// in memory would be prohibitive.
pub struct StreamingJsonWriter<W: Write> {
    writer: BufWriter<W>,
    rules_written: usize,
    flush_interval: usize,
}

impl<W: Write> StreamingJsonWriter<W> {
    /// Create a new streaming JSON writer with the given flush interval.
    /// The writer will flush to disk every `flush_interval` rules.
    pub fn new(writer: W, flush_interval: usize) -> Self {
        Self {
            writer: BufWriter::new(writer),
            rules_written: 0,
            flush_interval,
        }
    }

    /// Begin the JSON output, writing the opening structure.
    pub fn begin(&mut self, sequence: &[Symbol]) -> Result<()> {
        // Write opening brace and final_sequence
        writeln!(self.writer, "{{")?;

        // Write the final sequence
        let sequence_json = serialize_sequence(sequence);
        let sequence_str = serde_json::to_string(&sequence_json)?;
        writeln!(self.writer, "  \"final_sequence\": {},", sequence_str)?;

        // Begin the rules object
        writeln!(self.writer, "  \"rules\": {{")?;

        Ok(())
    }

    /// Write a single rule to the output.
    pub fn write_rule(&mut self, rule_id: usize, rule: &Rule) -> Result<()> {
        let rule_json = json!({
            "symbols": serialize_sequence(&rule.symbols),
            "usage_count": rule.usage_count,
            "depth": rule.depth.unwrap_or(0)
        });

        let rule_str = serde_json::to_string(&rule_json)?;

        // Add comma before rule if not the first
        if self.rules_written > 0 {
            writeln!(self.writer, ",")?;
        }
        write!(self.writer, "    \"{}\": {}", rule_id, rule_str)?;

        self.rules_written += 1;

        // Flush periodically
        if self.rules_written % self.flush_interval == 0 {
            self.writer.flush()?;
        }

        Ok(())
    }

    /// End the JSON output, writing the closing structure.
    pub fn end(mut self) -> Result<()> {
        // Close rules object
        writeln!(self.writer)?;
        writeln!(self.writer, "  }}")?;

        // Close main object
        writeln!(self.writer, "}}")?;

        self.writer.flush()?;
        Ok(())
    }
}

/// Write grammar to JSON using streaming mode.
/// This writes rules incrementally and flushes periodically to reduce memory usage.
pub fn write_grammar_json_streaming(
    path: &Path,
    sequence: &[Symbol],
    rules: &HashMap<usize, Rule>,
    flush_interval: usize,
) -> Result<()> {
    println!("Writing grammar to JSON (streaming): {}", path.display());

    let file = File::create(path)
        .context(format!("Failed to create file: {}", path.display()))?;

    let mut writer = StreamingJsonWriter::new(file, flush_interval);

    // Write opening and sequence
    writer.begin(sequence)?;

    // Write rules in sorted order for deterministic output
    let mut rule_ids: Vec<_> = rules.keys().collect();
    rule_ids.sort();

    for rule_id in rule_ids {
        if let Some(rule) = rules.get(rule_id) {
            writer.write_rule(*rule_id, rule)?;
        }
    }

    // Close the JSON
    writer.end()?;

    println!("Successfully wrote grammar to JSON (streaming).");
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::grammar::symbol::Direction;
    use crate::encode::dna_2bit::EncodedBase;
    use tempfile::tempdir;

    #[test]
    fn test_streaming_json_writer() {
        let dir = tempdir().unwrap();
        let path = dir.path().join("test_streaming.json");

        // Create a simple sequence
        let sequence = vec![
            Symbol::terminal(0, EncodedBase(0), Direction::Forward, None, None), // A
            Symbol::non_terminal(1, 1, Direction::Forward),
        ];

        // Create a simple rule
        let mut rules = HashMap::new();
        rules.insert(1, Rule {
            id: 1,
            symbols: vec![
                Symbol::terminal(0, EncodedBase(1), Direction::Forward, None, None), // C
                Symbol::terminal(1, EncodedBase(2), Direction::Forward, None, None), // G
            ],
            usage_count: 5,
            depth: Some(1),
            positions: vec![],
            assembly_index: None,
        });

        // Write using streaming
        write_grammar_json_streaming(&path, &sequence, &rules, 100).unwrap();

        // Read and verify it's valid JSON
        let content = std::fs::read_to_string(&path).unwrap();
        let parsed: serde_json::Value = serde_json::from_str(&content)
            .expect("Streaming output should be valid JSON");

        assert!(parsed.get("final_sequence").is_some());
        assert!(parsed.get("rules").is_some());
        assert!(parsed["rules"]["1"].is_object());
    }

    #[test]
    fn test_streaming_vs_regular_output_equivalent() {
        let dir = tempdir().unwrap();
        let path_streaming = dir.path().join("streaming.json");
        let path_regular = dir.path().join("regular.json");

        // Create test data
        let sequence = vec![
            Symbol::terminal(0, EncodedBase(0), Direction::Forward, None, None),
            Symbol::non_terminal(1, 1, Direction::Forward),
            Symbol::non_terminal(2, 2, Direction::Reverse),
        ];

        let mut rules = HashMap::new();
        rules.insert(1, Rule {
            id: 1,
            symbols: vec![
                Symbol::terminal(0, EncodedBase(1), Direction::Forward, None, None),
            ],
            usage_count: 3,
            depth: Some(1),
            positions: vec![],
            assembly_index: None,
        });
        rules.insert(2, Rule {
            id: 2,
            symbols: vec![
                Symbol::terminal(0, EncodedBase(2), Direction::Forward, None, None),
            ],
            usage_count: 2,
            depth: Some(1),
            positions: vec![],
            assembly_index: None,
        });

        // Write both ways
        write_grammar_json_streaming(&path_streaming, &sequence, &rules, 100).unwrap();
        write_grammar_json(&path_regular, &sequence, &rules).unwrap();

        // Parse both
        let streaming_content = std::fs::read_to_string(&path_streaming).unwrap();
        let regular_content = std::fs::read_to_string(&path_regular).unwrap();

        let streaming_json: serde_json::Value = serde_json::from_str(&streaming_content).unwrap();
        let regular_json: serde_json::Value = serde_json::from_str(&regular_content).unwrap();

        // Compare key fields
        assert_eq!(
            streaming_json["final_sequence"],
            regular_json["final_sequence"],
            "Sequences should match"
        );

        // Both should have the same rules
        assert!(streaming_json["rules"]["1"].is_object());
        assert!(streaming_json["rules"]["2"].is_object());
        assert!(regular_json["rules"]["1"].is_object());
        assert!(regular_json["rules"]["2"].is_object());
    }
} 