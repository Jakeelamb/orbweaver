use crate::grammar::symbol::{Symbol, SymbolType};
// use std::fmt; // Remove unused import
use serde::{Serialize, Deserialize}; // Import both Serialize and Deserialize
use std::collections::HashMap;

/// Represents a grammar rule containing a sequence of symbols
/// 
/// Each rule has a unique ID, a list of symbols (terminals and non-terminals),
/// and metadata such as a usage counter and original positions in the sequence.
#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct Rule {
    /// Unique identifier for this rule
    pub id: usize,
    
    /// The list of symbols (terminals and non-terminals) in this rule
    pub symbols: Vec<Symbol>,
    
    /// How many times this rule is used in the grammar
    pub usage_count: usize,
    
    /// The original positions where this rule was found (for reporting)
    pub positions: Vec<usize>,
}

impl Rule {
    /// Creates a new rule resulting from replacing a digram.
    pub fn new(id: usize, symbol1: Symbol, symbol2: Symbol) -> Self {
        Rule {
            id,
            symbols: vec![symbol1, symbol2],
            usage_count: 0, // Initial usage count is 0, incremented when used
            positions: Vec::new(), // Initialize empty positions
        }
    }
    
    /// Recursively expands a rule, substituting all non-terminal symbols with their corresponding rules.
    /// This is used for rule inlining when optimizing the grammar.
    /// 
    /// Args:
    ///    rule_map: A mapping from rule IDs to Rules.
    ///    max_depth: Maximum recursion depth to prevent infinite expansion (due to cycles)
    ///
    /// Returns:
    ///    A new Vec<Symbol> with all references expanded to their full form.
    pub fn expand_recursive(&self, rule_map: &HashMap<usize, Rule>, max_depth: usize) -> Vec<Symbol> {
        if max_depth == 0 {
            // Safety check against infinite recursion
            return self.symbols.clone();
        }
        
        let mut result = Vec::new();
        
        for symbol in &self.symbols {
            match symbol.symbol_type {
                SymbolType::NonTerminal(rule_id) => {
                    if let Some(referenced_rule) = rule_map.get(&rule_id) {
                        // Recursively expand the referenced rule
                        let expanded = referenced_rule.expand_recursive(rule_map, max_depth - 1);
                        
                        // Apply proper strand propagation
                        for exp_symbol in expanded {
                            let mut new_symbol = exp_symbol;
                            if symbol.strand == '-' {
                                // Flip the strand when expanding a rule referenced with a negative strand
                                new_symbol.strand = if new_symbol.strand == '+' { '-' } else { '+' };
                            }
                            result.push(new_symbol);
                        }
                    } else {
                        // If rule not found, keep the original symbol
                        result.push(*symbol);
                    }
                },
                SymbolType::Terminal(_) => {
                    // Terminal symbols are added as-is
                    result.push(*symbol);
                }
            }
        }
        
        result
    }
    
    /// Expand a rule that is used only once, for inlining purposes.
    /// Unlike expand_recursive, this is used when replacing a rule in the grammar.
    ///
    /// Args:
    ///    rule_map: Mapping from rule IDs to Rules.
    ///
    /// Returns:
    ///    The expanded symbols with correct strand information.
    pub fn expand_for_inlining(&self, rule_map: &HashMap<usize, Rule>) -> Vec<Symbol> {
        // Set a reasonable recursion limit
        self.expand_recursive(rule_map, 100)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::grammar::symbol::Symbol;

    #[test]
    fn test_rule_creation() {
        let s1 = Symbol::terminal(0, b'A', '+');
        let s2 = Symbol::non_terminal(1, 100, '-');
        let rule = Rule::new(101, s1, s2); // Rule 101 replaces digram (s1, s2)

        assert_eq!(rule.id, 101);
        assert_eq!(rule.usage_count, 0);
        assert_eq!(rule.symbols.len(), 2);
        assert_eq!(rule.symbols[0], s1);
        assert_eq!(rule.symbols[1], s2);
    }
} 