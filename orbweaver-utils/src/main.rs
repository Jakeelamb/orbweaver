use anyhow::{Context, Result};
use clap::{Parser, Subcommand};
use orbweaver::encode::EncodedBase;
use orbweaver::grammar::engine::Grammar;
use orbweaver::grammar::rule::Rule;
use orbweaver::grammar::symbol::{Direction, Symbol, SymbolType};
use serde::Deserialize;
use serde_json::Value;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufReader, BufWriter, Write};
use std::path::PathBuf;

#[derive(Parser, Debug)]
#[clap(name = "orbweaver-utils", version = "0.1.0", about = "Utility tools for Orbweaver grammars")]
struct Cli {
    #[clap(subcommand)]
    command: Commands,
}

#[derive(Subcommand, Debug)]
enum Commands {
    /// Extracts the top N longest rules from a grammar JSON file and outputs them as a GFA file.
    ExtractTopGfa {
        /// Path to the input grammar.json file from an orbweaver run.
        #[clap(short, long, value_parser)]
        input_json: PathBuf,

        /// Path for the output GFA file containing the top N rules.
        #[clap(short, long, value_parser)]
        output_gfa: PathBuf,

        /// The number of longest rules to extract (by expanded sequence length).
        #[clap(short, long, value_parser, default_value_t = 5)]
        top_n: usize,
    },
}

#[derive(Deserialize, Debug)]
struct JsonSymbolUtil {
    #[serde(rename = "type")]
    symbol_kind: String,
    value: Option<String>,   // For terminals, e.g., "A+"; for non-terminals, this will hold the rule_id as string
    strand: Option<String>,
}

#[derive(Deserialize, Debug)]
struct JsonRuleUtil {
    id: usize, // ADDED: Expect rule ID directly in the rule object
    symbols: Vec<JsonSymbolUtil>,
    usage_count: usize,
    depth: Option<usize>, // Already optional
}

#[derive(Deserialize, Debug)]
struct JsonGrammarUtil {
    final_sequence: Vec<JsonSymbolUtil>,
    rules: Vec<JsonRuleUtil>, // MODIFIED: Was Vec<HashMap<String, JsonRuleUtil>>
    max_depth: Option<usize>,             // Not in current JSON, will be None
}

fn convert_json_symbol_to_orb_symbol(json_symbol: &JsonSymbolUtil, symbol_instance_id: usize) -> Result<Symbol> {
    let symbol_kind_cleaned = json_symbol.symbol_kind.trim().to_lowercase();
    let sk_str = symbol_kind_cleaned.as_str();

    let orb_symbol_type;
    let orb_direction;

    if sk_str == "terminal" {
        let base_str = json_symbol.value.as_ref()
            .context("Terminal symbol missing 'value' field (expected base character)")?;
        if base_str.len() != 1 {
            anyhow::bail!("Terminal 'value' (base) string '{}' must be a single character.", base_str);
        }
        let base_char = base_str.chars().next().unwrap();

        let strand_str = json_symbol.strand.as_ref()
            .context("Terminal symbol missing 'strand' field")?;
        if strand_str.len() != 1 {
            anyhow::bail!("Terminal 'strand' string '{}' must be a single character.", strand_str);
        }
        let strand_char = strand_str.chars().next().unwrap();

        let base = EncodedBase::from_base(base_char as u8)
            .ok_or_else(|| anyhow::anyhow!("Invalid base character '{}' from terminal value", base_char))?;
        orb_symbol_type = SymbolType::Terminal(base);
        orb_direction = Direction::from_char(strand_char);

    } else if sk_str == "non-terminal" {
        let val_str = json_symbol.value.as_ref()
            .context("Non-terminal symbol missing 'value' field (expected to contain rule_id)")?;
        let rule_id = val_str.parse::<usize>()
            .with_context(|| format!("Failed to parse rule_id from non-terminal value: '{}'", val_str))?;

        let strand_str = json_symbol.strand.as_ref().context("Non-terminal symbol missing 'strand'")?;
        if strand_str.is_empty() {
            anyhow::bail!("Non-terminal strand string is empty for rule_id {}", rule_id);
        }
        let strand_char = strand_str.chars().next().unwrap();
        orb_symbol_type = SymbolType::NonTerminal(rule_id);
        orb_direction = Direction::from_char(strand_char);
    } else {
        anyhow::bail!("Unknown symbol type: '{}' (cleaned from '{}')", sk_str, json_symbol.symbol_kind);
    }

    Ok(Symbol {
        id: symbol_instance_id, // Assign a unique instance ID
        symbol_type: orb_symbol_type,
        strand: orb_direction,
    })
}

fn convert_json_grammar_to_orb_grammar(json_grammar_util: JsonGrammarUtil) -> Result<Grammar> {
    let mut grammar_rules = HashMap::new();
    let mut symbol_instance_counter = 0; // Counter for unique symbol instance IDs

    for json_rule in json_grammar_util.rules.iter() {
        let rule_id = json_rule.id;

        let mut rule_symbols = Vec::new();
        for json_sym in &json_rule.symbols {
            let orb_sym = convert_json_symbol_to_orb_symbol(json_sym, symbol_instance_counter)?;
            rule_symbols.push(orb_sym);
            symbol_instance_counter += 1;
        }

        grammar_rules.insert(rule_id, Rule {
            id: rule_id,
            symbols: rule_symbols,
            usage_count: json_rule.usage_count,
            positions: Vec::new(), // Positions are not in JSON, so initialize as empty
            depth: json_rule.depth,
        });
    }

    let mut grammar_sequence = Vec::new();
    for json_sym in json_grammar_util.final_sequence {
        let orb_sym = convert_json_symbol_to_orb_symbol(&json_sym, symbol_instance_counter)?;
        grammar_sequence.push(orb_sym);
        symbol_instance_counter += 1;
    }
    
    let mut max_depth_calc = 0;
    for rule in grammar_rules.values() {
        if let Some(d) = rule.depth {
            if d > max_depth_calc {
                max_depth_calc = d;
            }
        }
    }


    Ok(Grammar {
        sequence: grammar_sequence,
        rules: grammar_rules,
        max_depth: json_grammar_util.max_depth.unwrap_or(max_depth_calc),
        origins: HashMap::new(), // Origins are not in JSON
    })
}

// (Removed get_segment_id_and_orientation as it's from the main lib's output_gfa.rs)
// This utility binary will reconstruct its own GFA segments/links.

fn write_gfa_for_grammar(grammar: &Grammar, output_path: &PathBuf, top_n: usize) -> Result<()> {
    let mut rule_lengths: Vec<(usize, usize, String)> = Vec::new(); // rule_id, length, expanded_sequence

    for (id, rule) in &grammar.rules {
        let expanded_seq = orbweaver::utils::export::expand_rule_to_string(rule, &grammar.rules, Direction::Forward);
        rule_lengths.push((*id, expanded_seq.len(), expanded_seq));
    }

    rule_lengths.sort_by(|a, b| b.1.cmp(&a.1)); // Sort by length descending

    // Collect info for top N rules (ID and their pre-expanded forward sequence for their S-lines)
    let top_n_rules_s_lines: Vec<(String, String)> = rule_lengths.iter().take(top_n)
        .map(|(id, _len, seq_str)| (format!("R{}", id), seq_str.clone()))
        .collect();
    
    // Get just the IDs of the top_n rules for iterating their structure
    let top_n_rule_ids: Vec<usize> = rule_lengths.iter().take(top_n).map(|(id, _, _)| *id).collect();


    if top_n_rule_ids.is_empty() && !grammar.rules.is_empty() {
        println!("Warning: Top N rules list is empty (N={}), but grammar has {} rules. GFA will be minimal.", top_n, grammar.rules.len());
    } else if grammar.rules.is_empty() {
        println!("No rules found in the grammar. GFA file will be effectively empty (header only).");
    }

    let file = File::create(output_path)
        .with_context(|| format!("Failed to create GFA output file: {}", output_path.display()))?;
    let mut writer = BufWriter::new(file);

    writeln!(writer, "H\tVN:Z:1.0").context("Failed to write GFA header")?;

    // Keep track of all S-lines written (seg_id -> sequence)
    let mut segments_written: HashMap<String, String> = HashMap::new();

    // 1. Write S-Lines for the top N rules themselves
    for (seg_id, expanded_seq_str) in &top_n_rules_s_lines {
        if !expanded_seq_str.is_empty() {
            writeln!(writer, "S\t{}\t{}\tLN:i:{}", seg_id, expanded_seq_str, expanded_seq_str.len())
                .with_context(|| format!("Failed to write S-line for rule {}", seg_id))?;
            segments_written.insert(seg_id.clone(), expanded_seq_str.clone());
        } else {
            eprintln!("Warning: Top-N rule {} has an empty expanded sequence. Skipping its S-line.", seg_id);
        }
    }
    
    // 2. Write L-Lines for the structure of top N rules, and S-Lines for their components if not already written
    for rule_id_from_top_n in &top_n_rule_ids {
        if let Some(rule_obj) = grammar.rules.get(rule_id_from_top_n) {
            if rule_obj.symbols.len() >= 2 {
                for i in 0..(rule_obj.symbols.len() - 1) {
                    let sym1 = &rule_obj.symbols[i];
                    let sym2 = &rule_obj.symbols[i + 1];

                    // Process Symbol 1 to get its GFA identity and sequence for S-line
                    let (s1_id_str, s1_orient_char, s1_seq_for_s_line) = 
                        symbol_to_gfa_segment_parts(sym1, &grammar.rules, &segments_written)?;
                    
                    // Write S-line for sym1 if not already written
                    if !segments_written.contains_key(&s1_id_str) {
                        if s1_seq_for_s_line.is_empty() {
                             eprintln!("Warning: Segment {} (component of R{}) has empty sequence. Skipping its S-line.", s1_id_str, rule_id_from_top_n);
                        } else {
                            writeln!(writer, "S\t{}\t{}\tLN:i:{}", s1_id_str, s1_seq_for_s_line, s1_seq_for_s_line.len())
                                .with_context(|| format!("Failed to write S-line for segment {}", s1_id_str))?;
                            segments_written.insert(s1_id_str.clone(), s1_seq_for_s_line);
                        }
                    }

                    // Process Symbol 2
                    let (s2_id_str, s2_orient_char, s2_seq_for_s_line) = 
                        symbol_to_gfa_segment_parts(sym2, &grammar.rules, &segments_written)?;

                    // Write S-line for sym2 if not already written
                     if !segments_written.contains_key(&s2_id_str) {
                         if s2_seq_for_s_line.is_empty() {
                             eprintln!("Warning: Segment {} (component of R{}) has empty sequence. Skipping its S-line.", s2_id_str, rule_id_from_top_n);
                        } else {
                            writeln!(writer, "S\t{}\t{}\tLN:i:{}", s2_id_str, s2_seq_for_s_line, s2_seq_for_s_line.len())
                                .with_context(|| format!("Failed to write S-line for segment {}", s2_id_str))?;
                            segments_written.insert(s2_id_str.clone(), s2_seq_for_s_line);
                        }
                    }

                    // Write L-line, but only if both segments are valid (have been written or were already present)
                    if segments_written.contains_key(&s1_id_str) && segments_written.contains_key(&s2_id_str) {
                        writeln!(writer, "L\t{}\t{}\t{}\t{}\t0M", s1_id_str, s1_orient_char, s2_id_str, s2_orient_char)
                            .with_context(|| format!("Failed to write L-line between {} and {} for rule R{}", s1_id_str, s2_id_str, rule_id_from_top_n))?;
                    } else {
                        eprintln!("Warning: Skipping L-line between {} and {} due to missing S-line for one or both segments.", s1_id_str, s2_id_str);
                    }
                }
            }
        }
    }

    println!("Successfully wrote GFA for top {} rules (and their components) to {}", top_n_rule_ids.len(), output_path.display());
    Ok(())
}

// Helper to get GFA segment parts and manage writing S-lines for terminals
fn symbol_to_gfa_segment_parts<'a>(
    symbol: &'a Symbol,
    all_rules: &'a HashMap<usize, Rule>,
    // Read-only map of segments whose S-lines are already written (or planned from top-N)
    // and their canonical forward sequences.
    segments_already_defined: &HashMap<String, String> 
) -> Result<(String, char, String)> { // Returns: seg_id, orientation_in_link, sequence_for_s_line (always forward)
    match symbol.symbol_type {
        SymbolType::Terminal(base) => {
            let base_char = base.to_char();
            // Terminals like 'A', 'C', 'G', 'T' will be unique segments T_A, T_C etc.
            let seg_id = format!("T_{}", base_char); 
            let orientation = symbol.strand.to_char();
            let s_line_sequence = base_char.to_string(); // S-line sequence is just the base itself
            Ok((seg_id, orientation, s_line_sequence))
        }
        SymbolType::NonTerminal(rule_id) => {
            let seg_id = format!("R{}", rule_id);
            let orientation = symbol.strand.to_char();
            
            // Determine the sequence for this non-terminal's S-line.
            // It's always the forward expansion of the rule.
            let s_line_sequence = if let Some(existing_seq) = segments_already_defined.get(&seg_id) {
                // This rule's S-line was already written (e.g., it was a top-N rule) or its sequence determined.
                existing_seq.clone()
            } else if let Some(rule) = all_rules.get(&rule_id) {
                // This rule is a component, not a top-N rule itself. Expand it for its S-line.
                orbweaver::utils::export::expand_rule_to_string(rule, all_rules, Direction::Forward)
            } else {
                // This should ideally not happen if grammar is consistent.
                eprintln!("Warning: Rule R{} referenced in GFA structure not found in grammar. Using empty sequence for its S-line.", rule_id);
                String::new() 
            };
            Ok((seg_id, orientation, s_line_sequence))
        }
    }
}


fn main() -> Result<()> {
    let cli = Cli::parse();

    match cli.command {
        Commands::ExtractTopGfa { input_json, output_gfa, top_n } => {
            println!(
                "Extracting top {} rules from {} to {}",
                top_n, input_json.display(), output_gfa.display()
            );

            let file = File::open(&input_json)
                .with_context(|| format!("Failed to open input JSON file: {}", input_json.display()))?;
            let reader = BufReader::new(file);

            let generic_json_value: Value = serde_json::from_reader(reader)
                .with_context(|| format!("Failed to parse input JSON into a generic Value: {}. Check if the file is valid JSON.", input_json.display()))?;
            
            println!("Successfully parsed grammar.json into a generic serde_json::Value.");

            match &generic_json_value {
                Value::Object(map) => {
                    println!("DEBUG: Top-level JSON is an Object with keys: {:?}", map.keys().collect::<Vec<_>>());
                    if let Some(rules_value) = map.get("rules") {
                        print!("DEBUG: 'rules' field is of type: ");
                        match rules_value {
                            Value::Object(rules_map) => {
                                println!("Object with keys (first 5): {:?}", rules_map.keys().take(5).collect::<Vec<_>>());
                            }
                            Value::Array(rules_arr) => {
                                print!("Array with length: {}. ", rules_arr.len());
                                if let Some(first_elem) = rules_arr.first() {
                                    let first_elem_type_str = if first_elem.is_object() { "Object" } 
                                                            else if first_elem.is_array() { "Array" } 
                                                            else if first_elem.is_string() { "String" } 
                                                            else if first_elem.is_number() { "Number" } 
                                                            else if first_elem.is_boolean() { "Boolean" } 
                                                            else if first_elem.is_null() { "Null" } 
                                                            else { "Unknown" };
                                    println!("First element type: {}", first_elem_type_str);
                                } else {
                                    println!("Array is empty.");
                                }
                            }
                            Value::String(_) => println!("String"),
                            Value::Number(_) => println!("Number"),
                            Value::Bool(_) => println!("Boolean"),
                            Value::Null => println!("Null"),
                        }
                    } else {
                        println!("DEBUG: 'rules' field not found in top-level Object.");
                    }
                }
                Value::Array(arr) => {
                    println!("DEBUG: Top-level JSON is an Array with length: {}", arr.len());
                }
                Value::String(s) => {
                    println!("DEBUG: Top-level JSON is a String (length {}): '{}...'", s.len(), s.chars().take(50).collect::<String>());
                }
                Value::Number(n) => {
                    println!("DEBUG: Top-level JSON is a Number: {}", n);
                }
                Value::Bool(b) => {
                    println!("DEBUG: Top-level JSON is a Bool: {}", b);
                }
                Value::Null => {
                    println!("DEBUG: Top-level JSON is Null");
                }
            }

            let json_grammar_util: JsonGrammarUtil = serde_json::from_value(generic_json_value)
                .with_context(|| format!("Failed to deserialize generic JSON Value into JsonGrammarUtil struct for file: {}. Check struct definitions against JSON content.", input_json.display()))?;
            
            println!("Successfully deserialized generic Value into JsonGrammarUtil.");

            let orb_grammar = convert_json_grammar_to_orb_grammar(json_grammar_util)
                .context("Failed to convert JSON values to Orbweaver Grammar")?;
            
            println!("Successfully converted to Orbweaver Grammar struct. Total rules: {}", orb_grammar.rules.len());

            write_gfa_for_grammar(&orb_grammar, &output_gfa, top_n)
                .context("Failed to write GFA output for top N rules")?;
        }
    }

    Ok(())
} 