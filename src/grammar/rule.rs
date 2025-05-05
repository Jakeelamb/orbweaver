use crate::grammar::symbol::Symbol;
// use std::fmt; // Remove unused import
use serde::Serialize; // Import Serialize

/// Represents a production rule in the grammar.
/// Typically, rule R -> S1 S2, where S1 and S2 can be Terminal or NonTerminal.
#[derive(Debug, Serialize)]
pub struct Rule {
    pub id: usize,           // Unique identifier for this rule (corresponds to NonTerminal Symbol IDs)
    pub symbols: Vec<Symbol>, // The sequence this rule produces (initially 2 symbols)
    pub usage_count: usize,  // How many times this rule is used in the current sequence/other rules
                             // Strand usage tracking might be added here or managed by the GrammarBuilder
}

impl Rule {
    /// Creates a new rule resulting from replacing a digram.
    pub fn new(id: usize, symbol1: Symbol, symbol2: Symbol) -> Self {
        Rule {
            id,
            symbols: vec![symbol1, symbol2],
            usage_count: 0, // Initial usage count is 0, incremented when used
        }
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