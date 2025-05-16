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
    let rules_obj = json_grammar["rules"].as_object_mut().unwrap();
    for (rule_id, rule) in rules {
        let rule_json = json!({
            "symbols": serialize_sequence(&rule.symbols),
            "usage_count": rule.usage_count,
            "depth": rule.depth.unwrap_or(0)
        });
        rules_obj.insert(rule_id.to_string(), rule_json);
    }
    
    // Extract and clone the sequence array to avoid borrow issues
    let sequence_array = json_grammar["final_sequence"].as_array().unwrap().clone();
    
    // Validate that all rule IDs in the sequence exist in the rules map
    // If they don't exist, add placeholder empty rules to avoid broken references
    let rules_obj = json_grammar["rules"].as_object_mut().unwrap();
    for symbol_value in sequence_array { // symbol_value is a Value from json_grammar["final_sequence"]
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