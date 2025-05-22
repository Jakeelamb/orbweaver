use crate::grammar::rule::Rule;
use crate::grammar::symbol::{Symbol, SymbolType};
// use crate::encode::dna_2bit::EncodedBase; // This is not used in this function
use std::collections::{HashMap, HashSet, VecDeque};
use anyhow::{Result, bail};

/// Calculates and caches the assembly index for each rule in the grammar.
/// The assembly index is stored in Rule::assembly_index.
pub fn calculate_rule_assembly_indices(rules: &mut HashMap<usize, Rule>) -> Result<()> {
    // Ensure anyhow is imported if not already: use anyhow::{Result, anyhow, bail};
    // Track dependencies: rule_id -> set of child rule_ids
    let mut dependencies: HashMap<usize, HashSet<usize>> = HashMap::new();
    // Track reverse dependencies: rule_id -> set of parent rule_ids
    let mut reverse_dependencies: HashMap<usize, HashSet<usize>> = HashMap::new();
    // Track in-degree for topological sort
    let mut in_degree: HashMap<usize, usize> = HashMap::new();

    // Initialize dependencies and in-degree
    for (rule_id, rule) in rules.iter() {
        let mut deps = HashSet::new();
        for symbol in &rule.symbols {
            if let SymbolType::NonTerminal(child_id) = symbol.symbol_type {
                deps.insert(child_id);
                reverse_dependencies.entry(child_id).or_default().insert(*rule_id);
            }
        }
        dependencies.insert(*rule_id, deps.clone());
        in_degree.insert(*rule_id, deps.len());
    }

    // Topological sort: queue of rules with in-degree 0 (no dependencies)
    let mut queue: VecDeque<usize> = in_degree
        .iter()
        .filter_map(|(&rule_id, &deg)| if deg == 0 { Some(rule_id) } else { None })
        .collect();

    // Assembly index cache
    let mut ai_cache: HashMap<usize, usize> = HashMap::new();

    while let Some(rule_id) = queue.pop_front() {
        let rule = rules.get_mut(&rule_id)
            .ok_or_else(|| anyhow::anyhow!("Rule ID {} not found in rules map during AI calculation", rule_id))?;
        // Compute assembly index for this rule
        let mut sum_child_ai = 0;
        for symbol in &rule.symbols {
            if let SymbolType::NonTerminal(child_id) = symbol.symbol_type {
                let child_ai_val = ai_cache.get(&child_id).copied()
                    .ok_or_else(|| anyhow::anyhow!(
                        "Assembly index for child rule {} not found when calculating for parent rule {}",
                        child_id, rule_id
                    ))?;
                sum_child_ai += child_ai_val;
            }
            // Terminals contribute 0 to sum_child_ai as per plan: ai(Terminal) = 0
        }

        let num_symbols_in_rule = rule.symbols.len();
        let rule_ai = if num_symbols_in_rule > 0 {
            // ai(Rk) = (sum ai(Si)) + (m-1)
            sum_child_ai + (num_symbols_in_rule - 1)
        } else {
            0 // An empty rule has an AI of 0; (m-1) would be problematic.
        };
        rule.assembly_index = Some(rule_ai);
        ai_cache.insert(rule_id, rule_ai);

        // Update parents' in-degree
        if let Some(parents) = reverse_dependencies.get(&rule_id) {
            for &parent_id in parents {
                let entry = in_degree.get_mut(&parent_id)
                    .ok_or_else(|| anyhow::anyhow!("Parent rule ID {} not found in in_degree map", parent_id))?;
                *entry -= 1;
                if *entry == 0 {
                    queue.push_back(parent_id);
                }
            }
        }
    }

    // Check for cycles (rules with nonzero in-degree)
    let remaining: Vec<_> = in_degree.iter().filter(|(_, &deg)| deg > 0).map(|(&id, _)| id).collect();
    if !remaining.is_empty() {
        bail!("Cycle detected in rule dependencies: {:?}", remaining);
    }
    Ok(())
}

/// Calculates the assembly index for the final sequence using cached rule assembly indices.
pub fn calculate_chromosome_assembly_index(sequence: &[Symbol], rules: &HashMap<usize, Rule>) -> Result<usize> {
    let mut total_chromosome_ai = 0;
    for symbol in sequence {
        match symbol.symbol_type {
            SymbolType::Terminal(_) => {
                // Terminals contribute 0 to chromosome AI as per plan: ai(Terminal) = 0
            }
            SymbolType::NonTerminal(rule_id) => {
                let rule = rules.get(&rule_id)
                    .ok_or_else(|| anyhow::anyhow!("Rule ID {} referenced in sequence not found in rules map", rule_id))?;
                let rule_ai = rule.assembly_index
                    .ok_or_else(|| anyhow::anyhow!("Assembly index for rule {} must be computed before chromosome AI calculation", rule_id))?;
                total_chromosome_ai += rule_ai;
            }
        }
    }
    Ok(total_chromosome_ai)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::grammar::rule::Rule;
    use crate::grammar::symbol::{Symbol, SymbolType, Direction};
    use crate::encode::dna_2bit::EncodedBase;
    use std::collections::HashMap;

    // Helper function to create a simple terminal symbol
    fn term(id: usize, val: u8) -> Symbol {
        Symbol::terminal(id, EncodedBase(val), Direction::Forward, None, None)
    }

    // Helper function to create a non-terminal symbol
    fn non_term(id: usize, rule_id: usize) -> Symbol {
        Symbol::non_terminal(id, rule_id, Direction::Forward)
    }
    
    // Helper to create a rule for testing
    fn create_rule(id: usize, symbols: Vec<Symbol>, assembly_index: Option<usize>) -> Rule {
        Rule {
            id,
            symbols,
            usage_count: 1, // Dummy value
            positions: vec![], // Dummy value
            depth: None, // Dummy value
            assembly_index,
        }
    }


    #[test]
    fn test_calculate_rule_assembly_indices_simple() {
        let mut grammar = Grammar {
            sequence: vec![],
            rules: HashMap::new(),
            max_depth: 0,
            origins: HashMap::new(),
        };
        // R0 -> A T
        grammar.rules.insert(0, create_rule(0, vec![term(0,0), term(1,3)], None));
        // R1 -> C R0
        grammar.rules.insert(1, create_rule(1, vec![term(2,1), non_term(3,0)], None));
        // R2 -> R1 G
        grammar.rules.insert(2, create_rule(2, vec![non_term(4,1), term(5,2)], None));

        calculate_rule_assembly_indices(&mut grammar.rules).unwrap();

        assert_eq!(grammar.rules.get(&0).unwrap().assembly_index, Some(1)); // 1 + 0 + 0
        assert_eq!(grammar.rules.get(&1).unwrap().assembly_index, Some(1 + 1)); // 1 + 0 + ai(R0)
        assert_eq!(grammar.rules.get(&2).unwrap().assembly_index, Some(1 + (1+1)));// 1 + ai(R1) + 0
    }
    
    #[test]
    fn test_calculate_rule_assembly_indices_multiple_dependencies() {
        let mut grammar = Grammar {
            sequence: vec![],
            rules: HashMap::new(),
            max_depth: 0,
            origins: HashMap::new(),
        };
        // R0 -> A T
        grammar.rules.insert(0, create_rule(0, vec![term(0,0), term(1,3)], None));
        // R1 -> C G
        grammar.rules.insert(1, create_rule(1, vec![term(2,1), term(3,2)], None));
        // R2 -> R0 R1
        grammar.rules.insert(2, create_rule(2, vec![non_term(4,0), non_term(5,1)], None));

        calculate_rule_assembly_indices(&mut grammar.rules).unwrap();
        assert_eq!(grammar.rules.get(&0).unwrap().assembly_index, Some(1));
        assert_eq!(grammar.rules.get(&1).unwrap().assembly_index, Some(1));
        assert_eq!(grammar.rules.get(&2).unwrap().assembly_index, Some(1 + 1 + 1)); // 1 + ai(R0) + ai(R1)
    }

    #[test]
    fn test_calculate_rule_assembly_indices_no_rules() {
        let mut grammar = Grammar {
            sequence: vec![],
            rules: HashMap::new(),
            max_depth: 0,
            origins: HashMap::new(),
        };
        calculate_rule_assembly_indices(&mut grammar.rules).unwrap();
        assert!(grammar.rules.is_empty());
    }
    
    #[test]
    fn test_calculate_rule_assembly_indices_cyclic_dependency_warning() {
        // This test primarily checks if it handles cycles gracefully (doesn't panic, prints warning)
        // The AI values might be None for rules involved in a cycle.
        let mut grammar = Grammar {
            sequence: vec![],
            rules: HashMap::new(),
            max_depth: 0,
            origins: HashMap::new(),
        };
        // R0 -> R1
        grammar.rules.insert(0, create_rule(0, vec![non_term(0,1)], None));
        // R1 -> R0
        grammar.rules.insert(1, create_rule(1, vec![non_term(1,0)], None));

        calculate_rule_assembly_indices(&mut grammar.rules).unwrap();
        // Expectation is that indices will be None and warnings printed (captured via stderr or logs if set up)
        assert!(grammar.rules.get(&0).unwrap().assembly_index.is_none());
        assert!(grammar.rules.get(&1).unwrap().assembly_index.is_none());
    }

    #[test]
    fn test_calculate_chromosome_assembly_index_simple() {
        let mut rules = HashMap::new();
        rules.insert(0, create_rule(0, vec![term(0,0), term(1,1)], Some(1))); // R0 -> A C, ai=1
        rules.insert(1, create_rule(1, vec![term(2,2), term(3,3)], Some(1))); // R1 -> G T, ai=1
        rules.insert(2, create_rule(2, vec![non_term(4,0), non_term(5,1)], Some(3))); // R2 -> R0 R1, ai=3 (1+ai(R0)+ai(R1))

        let grammar = Grammar {
            sequence: vec![non_term(10, 2), term(11,0)], // R2 A
            rules,
            max_depth: 0,
            origins: HashMap::new(),
        };
        
        let chromosome_ai = calculate_chromosome_assembly_index(&grammar.sequence, &grammar.rules).unwrap();
        assert_eq!(chromosome_ai, 3); // ai(R2) + ai(A) = 3 + 0
    }

    #[test]
    fn test_calculate_chromosome_assembly_index_all_terminals() {
        let grammar = Grammar {
            sequence: vec![term(0,0), term(1,1), term(2,2)], // A C G
            rules: HashMap::new(),
            max_depth: 0,
            origins: HashMap::new(),
        };
        let chromosome_ai = calculate_chromosome_assembly_index(&grammar.sequence, &grammar.rules).unwrap();
        assert_eq!(chromosome_ai, 0); // 0+0+0
    }

    #[test]
    fn test_calculate_chromosome_assembly_index_empty_sequence() {
        let grammar = Grammar {
            sequence: vec![],
            rules: HashMap::new(),
            max_depth: 0,
            origins: HashMap::new(),
        };
        let chromosome_ai = calculate_chromosome_assembly_index(&grammar.sequence, &grammar.rules).unwrap();
        assert_eq!(chromosome_ai, 0);
    }
}