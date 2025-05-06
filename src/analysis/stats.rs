use crate::grammar::symbol::SymbolType;
use crate::grammar::rule::Rule;
use crate::grammar::engine::Grammar;
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

/// Calculate and print statistics about the generated grammar
pub fn calculate_and_print_stats(grammar: &Grammar) -> Result<()> {
    let (sequence, rules) = (&grammar.sequence, &grammar.rules);
    
    // Initialize counters
    let mut total_terminals = 0;
    let mut total_non_terminals = 0;
    let mut rule_usage_stats = HashMap::new();
    let mut rule_depths = HashMap::new();
    let mut rule_lengths = HashMap::new();
    
    // Count terminal and non-terminal symbols in the final sequence
    for symbol in sequence {
        match symbol.symbol_type {
            SymbolType::Terminal(_) => total_terminals += 1,
            SymbolType::NonTerminal(_) => total_non_terminals += 1,
        }
    }
    
    // Analyze rules
    let mut total_rule_size = 0;
    for (rule_id, rule) in rules {
        // Count usage
        rule_usage_stats.insert(*rule_id, rule.usage_count);
        
        // Count rule length (number of symbols)
        let rule_len = rule.symbols.len();
        rule_lengths.insert(*rule_id, rule_len);
        total_rule_size += rule_len;
        
        // Get rule depth
        rule_depths.insert(*rule_id, rule.depth.unwrap_or(0));
    }
    
    // Calculate some derived statistics
    let total_rules = rules.len();
    let avg_rule_length = if total_rules > 0 {
        total_rule_size as f64 / total_rules as f64
    } else {
        0.0
    };
    
    let avg_rule_usage = if total_rules > 0 {
        rule_usage_stats.values().sum::<usize>() as f64 / total_rules as f64
    } else {
        0.0
    };
    
    let max_depth = grammar.max_depth;
    
    // Print statistics
    println!("\n--- Grammar Statistics ---");
    println!("Total number of rules: {}", total_rules);
    println!("Final sequence statistics:");
    println!("  Terminal symbols: {}", total_terminals);
    println!("  Non-terminal symbols (rule references): {}", total_non_terminals);
    println!("  Total length: {}", total_terminals + total_non_terminals);
    println!("Rule statistics:");
    println!("  Average rule length: {:.2} symbols", avg_rule_length);
    println!("  Average rule usage: {:.2} times", avg_rule_usage);
    println!("  Maximum rule depth: {}", max_depth);
    
    // Find most frequently used rules
    if !rule_usage_stats.is_empty() {
        let mut usage_vec: Vec<(&usize, &usize)> = rule_usage_stats.iter().collect();
        usage_vec.sort_by(|a, b| b.1.cmp(a.1));
        
        println!("\nTop 5 most frequently used rules:");
        for (i, (rule_id, usage)) in usage_vec.iter().take(5).enumerate() {
            println!("  {}. Rule {} - used {} times, depth {}, length {} symbols",
                i + 1, rule_id, usage, 
                rule_depths.get(rule_id).unwrap_or(&0),
                rule_lengths.get(rule_id).unwrap_or(&0));
        }
    }
    
    Ok(())
} 