use crate::encode::dna_2bit::EncodedBase;
use crate::grammar::rule::Rule;
use crate::grammar::symbol::{Direction, Symbol};
use anyhow::Result;
use std::collections::{HashMap};
use rayon::prelude::*;
use serde::{Serialize, Deserialize};

/// A Grammar is the result of the Sequitur algorithm,
/// containing a sequence of symbols and a set of rules
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Grammar {
    /// The compressed sequence
    pub sequence: Vec<Symbol>,
    
    /// The rules discovered during grammar inference
    pub rules: HashMap<usize, Rule>,
    
    /// Maximum rule depth (for hierarchy analysis)
    pub max_depth: usize,
}

/// Sequitur algorithm implementation for inferring hierarchical structure in sequences
pub struct Sequitur {
    /// Minimum number of occurrences required to create a rule
    min_rule_usage: usize,
    
    /// Whether to consider reverse complements
    reverse_aware: bool,
    
    /// Counter for assigning rule IDs
    next_rule_id: usize,
}

impl Sequitur {
    /// Create a new Sequitur instance
    pub fn new(min_rule_usage: usize, reverse_aware: bool) -> Self {
        Self {
            min_rule_usage,
            reverse_aware,
            next_rule_id: 0,
        }
    }
    
    /// Build a grammar from a DNA sequence
    pub fn build_grammar(&mut self, sequence: &[EncodedBase]) -> Result<Grammar> {
        // This is just a stub implementation
        // The real implementation would use the digram table and build a grammar
        
        // For now, just return a simple grammar with the original sequence
        let mut symbols = Vec::with_capacity(sequence.len());
        
        for (i, &base) in sequence.iter().enumerate() {
            symbols.push(Symbol::terminal(i, base, Direction::Forward));
        }
        
        Ok(Grammar {
            sequence: symbols,
            rules: HashMap::new(),
            max_depth: 0,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_basic_sequitur() -> Result<()> {
        let seq = vec![
            EncodedBase(0), // A
            EncodedBase(1), // C
            EncodedBase(0), // A
            EncodedBase(1), // C
        ];
        
        let mut sequitur = Sequitur::new(2, false);
        let grammar = sequitur.build_grammar(&seq)?;
        
        // The stub implementation doesn't create rules, so verify the sequence
        assert_eq!(grammar.sequence.len(), seq.len());
        assert!(grammar.rules.is_empty());
        
        Ok(())
    }
} 