
// Placeholder for stats module

use anyhow::{Context, Result};
use crate::grammar::engine::Grammar;
use crate::grammar::rule::Rule;
use crate::grammar::symbol::{Symbol, SymbolType};
use std::collections::HashMap;
use std::fs::File;
use std::io::Write;
// use crate::encode::dna_2bit::EncodedBase; // Required for EncodedBase.to_char() - Unused
use crate::analysis::assembly_index::calculate_chromosome_assembly_index;
use std::fs::OpenOptions; // For appending to file
use std::path::PathBuf;

// Helper function to calculate uncompressed length of a sequence of symbols
fn get_sequence_uncompressed_length(
    sequence: &[Symbol],
    rules: &HashMap<usize, Rule>,
    memo: &mut HashMap<usize, usize>, // Memoization for rule_id -> uncompressed_length
) -> usize {
    let mut total_length = 0;
    for symbol in sequence {
        match symbol.symbol_type {
            SymbolType::Terminal(_) => {
                total_length += 1; // Assuming each EncodedBase is 1 bp
            }
            SymbolType::NonTerminal(rule_id) => {
                if let Some(cached_len) = memo.get(&rule_id) {
                    total_length += *cached_len;
                } else if let Some(rule) = rules.get(&rule_id) {
                    // To prevent infinite recursion on malformed/cyclic grammars,
                    // this assumes an acyclic grammar (max_depth calculation helps ensure this).
                    let rule_uncompressed_len =
                        get_sequence_uncompressed_length(&rule.symbols, rules, memo);
                    memo.insert(rule_id, rule_uncompressed_len);
                    total_length += rule_uncompressed_len;
                } else {
                    // This case indicates an invalid grammar or an issue in rule mapping.
                    // eprintln!("Warning: Rule ID {} not found during uncompressed length calculation.", rule_id);
                }
            }
        }
    }
    total_length
}

pub fn calculate_and_print_stats(
    grammar: &Grammar,
    output_path: &PathBuf,
    print_to_console: bool,
    // _verbose: bool, // Replaced by chromosome_details
    chromosome_details: &[(String, usize)], // List of (chromosome_name, original_base_length)
) -> Result<()> {
    let num_rules = grammar.rules.len();
    let compressed_sequence_len_symbols = grammar.sequence.len();

    let mut memo_uncompressed_len = HashMap::new();
    let compressed_sequence_len_bp = get_sequence_uncompressed_length(
        &grammar.sequence,
        &grammar.rules,
        &mut memo_uncompressed_len,
    );

    let max_depth = grammar.max_depth;
    let total_symbols_in_rules = grammar
        .rules
        .values()
        .map(|r| r.symbols.len())
        .sum::<usize>();
    
    let avg_symbols_per_rule = if num_rules > 0 {
        total_symbols_in_rules as f64 / num_rules as f64
    } else {
        0.0
    };

    let rules_with_assembly_index = grammar
        .rules
        .values()
        .filter(|r| r.assembly_index.is_some())
        .count();

    // Note: True compression ratio would require original_sequence_length_bp,
    // which is not passed to this function currently. 
    // This could be added as a new parameter if needed.

    let stats_string = format!(
        "--- Grammar Statistics ---\n\
        Number of Rules: {}\n\
        Compressed Sequence Length (Symbols): {}\n\
        Compressed Sequence Length (Base Pairs): {}\n\
        Maximum Rule Depth: {}\n\
        Total Symbols in Rule Definitions: {}\n\
        Average Symbols per Rule: {:.2}\n\
        Rules with Assembly Index calculated: {}\n\
        Number of Chromosomes/Sequences: {}\n\
        --- End of Statistics ---",
        num_rules,
        compressed_sequence_len_symbols,
        compressed_sequence_len_bp,
        max_depth,
        total_symbols_in_rules,
        avg_symbols_per_rule,
        rules_with_assembly_index,
        chromosome_details.len()
    );

    if print_to_console {
        println!("\n{}\n", stats_string);
    }

    let mut file = File::create(output_path)
        .with_context(|| format!("Failed to create stats file at {:?}", output_path))?;
    writeln!(file, "{}", stats_string)
        .with_context(|| format!("Failed to write stats to file {:?}", output_path))?;

    // --- Add detailed rule assembly index output --- 
    let details_filename = output_path.file_name().map_or_else(
        || "rule_assembly_details.txt".into(),
        |name| {
            let mut new_name = name.to_os_string();
            new_name.push(".rule_details.txt"); // Append to original stats filename
            PathBuf::from(new_name)
        }
    );
    let rule_details_path = output_path.with_file_name(details_filename);

    let mut details_file = File::create(&rule_details_path)
        .with_context(|| format!("Failed to create rule assembly details file at {:?}", rule_details_path))?;

    writeln!(details_file, "--- Rule Assembly Index Details ---")
        .with_context(|| format!("Failed to write header to rule details file {:?}", rule_details_path))?;

    let mut sorted_rules: Vec<&Rule> = grammar.rules.values().collect();
    sorted_rules.sort_by_key(|r| r.id);

    for rule in sorted_rules {
        let symbol_strings: Vec<String> = rule.symbols.iter().map(|symbol| {
            match symbol.symbol_type {
                SymbolType::Terminal(encoded_base) => format!("T({})", encoded_base.to_char()),
                SymbolType::NonTerminal(id) => format!("R{}", id),
            }
        }).collect();
        
        let rule_detail_string = format!(
            "Rule {}: [{}] -> AI: {}",
            rule.id,
            symbol_strings.join(", "),
            rule.assembly_index.map_or_else(|| "None".to_string(), |ai| ai.to_string())
        );
        writeln!(details_file, "{}", rule_detail_string)
            .with_context(|| format!("Failed to write rule detail for rule {} to file {:?}", rule.id, rule_details_path))?;
    }

    writeln!(details_file, "--- End of Rule Assembly Index Details ---")
        .with_context(|| format!("Failed to write footer to rule details file {:?}", rule_details_path))?;

    // --- Detailed Assembly Index Calculations (Per-Chromosome and Whole Genome) ---
    let assembly_index_filename = "assembly_index.txt";
    let assembly_index_path = output_path
        .parent()
        .ok_or_else(|| anyhow::anyhow!("Failed to get parent directory for output path: {:?}", output_path))?
        .join(assembly_index_filename);

    let mut assembly_index_file = File::create(&assembly_index_path)
        .with_context(|| format!("Failed to create assembly index file at {:?}", assembly_index_path))?;

    writeln!(assembly_index_file, "Number of Chromosomes/Sequences: {}", chromosome_details.len())
        .with_context(|| format!("Failed to write chromosome count to assembly index file {:?}", assembly_index_path))?;

    for (chr_idx, (chr_name, original_base_len)) in chromosome_details.iter().enumerate() {
        // Filter symbols belonging to the current chromosome based on source_grammar_id.
        // This assumes merge_grammars (if used) sets symbol.source_grammar_id to the index of the grammar
        // in the list of grammars it merged, corresponding to chr_idx here.
        // Filter symbols belonging to the current chromosome based on source_grammar_id.
        // TEMP: Commenting out filter for now as source_grammar_id is not yet on Symbol.
        // This means all symbols will be processed for each chromosome detail section currently.
        let chromosome_symbols_refs: Vec<&Symbol> = grammar.sequence.iter()
            // .filter(|s| (**s).source_grammar_id == Some(chr_idx)) // Corrected dereference, but field still missing
            .collect();

        if chromosome_symbols_refs.is_empty() {
            if print_to_console {
                println!("Chromosome: {} (Original Base Length: {}) - No symbols found in final grammar with source_grammar_id {}. Skipping detailed AI.", 
                    chr_name, original_base_len, chr_idx);
            }
            writeln!(assembly_index_file, "\n--- Chromosome: {} (Original Base Length: {}, Symbols in Final Grammar: 0) ---", chr_name, original_base_len)
                 .with_context(|| format!("Failed to write header for empty chr {} to assembly index file", chr_name))?;
            writeln!(assembly_index_file, "(No symbols uniquely attributed to this chromosome in the final sequence via source_grammar_id for detailed AI calculation)")
                .with_context(|| format!("Failed to write note for empty chr {} to assembly index file", chr_name))?;
            writeln!(assembly_index_file, "--- End of Chromosome: {} ---", chr_name)
                .with_context(|| format!("Failed to write footer for empty chr {} to assembly index file", chr_name))?;
            continue;
        }

        let chromosome_symbols: Vec<Symbol> = chromosome_symbols_refs.into_iter().cloned().collect();

        match calculate_chromosome_assembly_index(&chromosome_symbols, &grammar.rules) {
            Ok(chr_ai) => {
                if print_to_console {
                    println!("Chromosome: {} (Original Base Length: {}) - Assembly Index: {}", chr_name, original_base_len, chr_ai);
                }

                writeln!(assembly_index_file, "\n--- Chromosome: {} (Original Base Length: {}, Symbols in Final Grammar: {}, Assembly Index: {}) ---", 
                         chr_name, original_base_len, chromosome_symbols.len(), chr_ai)
                    .with_context(|| format!("Failed to write header for chr {} to assembly index file", chr_name))?;

                for symbol in &chromosome_symbols {
                    let symbol_repr: String;
                    let symbol_ai_contrib: usize;
                    match symbol.symbol_type {
                        SymbolType::Terminal(eb) => {
                            symbol_repr = format!("T({})", eb.to_char()); // TEMP: (Original Pos: {}) removed, field not on Symbol
                            symbol_ai_contrib = 0;
                        }
                        SymbolType::NonTerminal(id) => {
                            symbol_repr = format!("R{}", id); // TEMP: (Original Pos: {}) removed, field not on Symbol
                            symbol_ai_contrib = grammar.rules.get(&id)
                                .and_then(|rule| rule.assembly_index)
                                .unwrap_or(0);
                        }
                    }
                    let detail_line = format!("Symbol: {} (AI Contribution: {})", symbol_repr, symbol_ai_contrib);
                    writeln!(assembly_index_file, "{}", detail_line)
                        .with_context(|| format!("Failed to write detail line for chr {} to assembly index file", chr_name))?;
                }
                writeln!(assembly_index_file, "--- End of Chromosome: {} ---", chr_name)
                    .with_context(|| format!("Failed to write footer for chr {} to assembly index file", chr_name))?;
            }
            Err(e) => {
                if print_to_console {
                    println!("Warning: Failed to calculate Assembly Index for chromosome {}: {}. Statistics might be incomplete.", chr_name, e);
                }
                writeln!(assembly_index_file, "\n--- Chromosome: {} (Original Base Length: {}) ---", chr_name, original_base_len)
                    .with_context(|| format!("Failed to write error header for chr {} to assembly index file", chr_name))?;
                writeln!(assembly_index_file, "Warning: Failed to calculate Assembly Index: {}.", e)
                    .with_context(|| format!("Failed to write error for chr {} to assembly index file", chr_name))?;
                writeln!(assembly_index_file, "--- End of Chromosome: {} ---", chr_name)
                    .with_context(|| format!("Failed to write error footer for chr {} to assembly index file", chr_name))?;
            }
        }
    }

    // --- Whole Genome Summary ---    
    writeln!(assembly_index_file, "\n--- Whole Genome Summary ---")
        .with_context(|| format!("Failed to write WGA summary header to assembly index file {:?}", assembly_index_path))?;

    match calculate_chromosome_assembly_index(&grammar.sequence, &grammar.rules) {
        Ok(genome_ai) => {


            let total_line = format!("Total Whole Genome Ensemble Assembly Index: {}", genome_ai);
            
            writeln!(assembly_index_file, "{}", total_line)
                .with_context(|| format!("Failed to write total WGA to assembly index file {:?}", assembly_index_path))?;
            writeln!(assembly_index_file, "--- End of Whole Genome Summary ---")
                .with_context(|| format!("Failed to write WGA footer to assembly index file {:?}", assembly_index_path))?; 

            if print_to_console {
                println!("Detailed per-chromosome and Whole Genome Ensemble Assembly Index calculations saved to: {}", assembly_index_path.display());
                println!("{}", total_line);
            }

            let mut stats_file_appender = OpenOptions::new()
                .append(true)
                .open(output_path)
                .with_context(|| format!("Failed to open stats file {:?} for appending WGAI total", output_path))?;
            
            writeln!(stats_file_appender, "\n{}", total_line)
                .with_context(|| format!("Failed to write WGAI total to stats file {:?}", output_path))?;
        }
        Err(e) => {
            eprintln!("Warning: Failed to calculate Whole Genome Ensemble Assembly Index: {}. Statistics might be incomplete.", e);
            if let Ok(mut file) = OpenOptions::new().append(true).open(output_path) {
                let _ = writeln!(file, "\nWarning: Failed to calculate Whole Genome Ensemble Assembly Index: {}.", e);
            }
        }
    }

    Ok(())
} 