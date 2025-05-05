use crate::grammar::builder::GrammarBuilder;
use crate::grammar::rule::Rule;
use crate::grammar::symbol::SymbolType;
use anyhow::Result;
use std::collections::{HashMap, HashSet};

/// Calculates the depth of a rule (longest path of non-terminal expansions).
/// Handles cycles by returning an error or a max depth indicator.
fn get_rule_depth(
    rule_id: usize,
    all_rules: &HashMap<usize, Rule>,
    depth_cache: &mut HashMap<usize, Option<usize>>,
    visited_stack: &mut HashSet<usize>,
) -> Result<usize> {
    // Check cache first
    if let Some(cached_depth) = depth_cache.get(&rule_id) {
        return match cached_depth {
            Some(d) => Ok(*d),
            None => Err(anyhow::anyhow!("Cycle detected involving rule {}", rule_id)),
        };
    }

    // Check for cycle
    if !visited_stack.insert(rule_id) {
        depth_cache.insert(rule_id, None); // Mark as cyclic in cache
        return Err(anyhow::anyhow!("Cycle detected involving rule {}", rule_id));
    }

    let rule = all_rules.get(&rule_id).ok_or_else(|| {
        anyhow::anyhow!("Rule ID {} not found during depth calculation", rule_id)
    })?;

    let mut max_sub_depth = 0;
    for symbol in &rule.symbols {
        if let SymbolType::NonTerminal(sub_rule_id) = symbol.symbol_type {
            match get_rule_depth(sub_rule_id, all_rules, depth_cache, visited_stack) {
                Ok(depth) => max_sub_depth = max_sub_depth.max(depth),
                Err(e) => {
                    // Propagate cycle error
                    visited_stack.remove(&rule_id);
                    depth_cache.insert(rule_id, None);
                    return Err(e);
                }
            }
        }
    }

    visited_stack.remove(&rule_id);
    let current_depth = max_sub_depth + 1;
    depth_cache.insert(rule_id, Some(current_depth)); // Cache result
    Ok(current_depth)
}

/// Calculates the length of a rule definition (number of symbols in its expansion).
fn get_rule_definition_length(rule_id: usize, all_rules: &HashMap<usize, Rule>) -> usize {
    all_rules.get(&rule_id).map_or(0, |r| r.symbols.len())
}

/// Calculates and prints statistics about the generated grammar.
pub fn calculate_and_print_stats(
    grammar_builder: &GrammarBuilder,
    initial_sequence_len: usize,
) -> Result<()> {
    println!("Calculating Grammar Statistics...");

    let (final_sequence, rules) = grammar_builder.get_grammar();
    let num_rules = rules.len();
    let final_seq_len = final_sequence.len();

    // Calculate Max Rule Depth
    let mut max_depth = 0;
    let mut depth_cache: HashMap<usize, Option<usize>> = HashMap::new();
    for rule_id in rules.keys() {
        let mut visited_stack = HashSet::new();
        match get_rule_depth(*rule_id, rules, &mut depth_cache, &mut visited_stack) {
            Ok(depth) => max_depth = max_depth.max(depth),
            Err(e) => {
                // Report cycle and continue, max_depth will be based on non-cyclic paths
                eprintln!("Warning calculating depth: {}", e);
            }
        }
    }

    // Calculate Compression Ratio Component: Sum of Rule Definition Lengths
    let total_rule_def_len: usize = rules.keys().map(|&id| get_rule_definition_length(id, rules)).sum();
    
    // Size of compressed representation = final sequence + rule definitions
    let compressed_size = final_seq_len + total_rule_def_len;
    let compression_ratio = if initial_sequence_len > 0 {
        compressed_size as f64 / initial_sequence_len as f64
    } else {
        0.0 // Avoid division by zero
    };

    // --- Print Stats --- 
    println!("----------------------------------------");
    println!("Grammar Statistics:");
    println!("  Initial Sequence Length: {}", initial_sequence_len);
    println!("  Final Sequence Length:   {}", final_seq_len);
    println!("  Number of Rules:         {}", num_rules);
    println!("  Max Rule Depth:          {}", max_depth);
    println!("  Sum Rule Def Lengths:  {}", total_rule_def_len);
    println!("  Compressed Size:       {}", compressed_size);
    println!("  Compression Ratio:     {}", compression_ratio);
    println!("----------------------------------------");

    Ok(())
} 