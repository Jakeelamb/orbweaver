use std::fmt;
use serde::{Serialize, Deserialize};

/// Represents the type of a symbol in the grammar.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum SymbolType {
    Terminal(u8), // Represents a single base (e.g., A, C, G, T as u8)
    NonTerminal(usize), // Represents a rule ID
}

/// Represents a symbol instance in the grammar sequence or rules.
#[derive(Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct Symbol {
    pub id: usize, // Unique identifier for this *instance* of the symbol (e.g., its position in the original sequence or a rule expansion)
    pub symbol_type: SymbolType,
    pub strand: char, // '+' or '-'
}

impl Symbol {
    /// Creates a new terminal symbol instance.
    pub fn terminal(id: usize, base: u8, strand: char) -> Self {
        // Optionally validate base and strand here
        Symbol {
            id,
            symbol_type: SymbolType::Terminal(base),
            strand,
        }
    }

    /// Creates a new non-terminal symbol instance.
    pub fn non_terminal(id: usize, rule_id: usize, strand: char) -> Self {
        // Optionally validate strand here
        Symbol {
            id,
            symbol_type: SymbolType::NonTerminal(rule_id),
            strand,
        }
    }
}

// Custom Debug implementation for cleaner output
impl fmt::Debug for Symbol {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self.symbol_type {
            SymbolType::Terminal(base) => write!(
                f,
                "Sym(id:{}, T:{}, {})",
                self.id, base as char, self.strand
            ),
            SymbolType::NonTerminal(rule_id) => {
                write!(f, "Sym(id:{}, R:{}, {})", self.id, rule_id, self.strand)
            }
        }
    }
}


#[cfg(test)]
mod tests {
    use super::*; // Imports items from the parent module (symbol.rs)
    use std::collections::HashSet;

    #[test]
    fn test_symbol_creation() {
        let term = Symbol::terminal(0, b'A', '+');
        assert_eq!(term.id, 0);
        assert_eq!(term.strand, '+');
        assert_eq!(term.symbol_type, SymbolType::Terminal(b'A'));

        let nonterm = Symbol::non_terminal(1, 100, '-');
        assert_eq!(nonterm.id, 1);
        assert_eq!(nonterm.strand, '-');
        assert_eq!(nonterm.symbol_type, SymbolType::NonTerminal(100));
    }

    #[test]
    fn test_symbol_equality() {
        let term1a = Symbol::terminal(0, b'A', '+');
        let term1b = Symbol::terminal(0, b'A', '+');
        let term2 = Symbol::terminal(1, b'A', '+'); // Different ID
        let term3 = Symbol::terminal(0, b'C', '+'); // Different base
        let term4 = Symbol::terminal(0, b'A', '-'); // Different strand

        assert_eq!(term1a, term1b);
        assert_ne!(term1a, term2);
        assert_ne!(term1a, term3);
        assert_ne!(term1a, term4);

        let nonterm1a = Symbol::non_terminal(5, 100, '-');
        let nonterm1b = Symbol::non_terminal(5, 100, '-');
        let nonterm2 = Symbol::non_terminal(6, 100, '-'); // Different ID
        let nonterm3 = Symbol::non_terminal(5, 101, '-'); // Different rule ID
        let nonterm4 = Symbol::non_terminal(5, 100, '+'); // Different strand

        assert_eq!(nonterm1a, nonterm1b);
        assert_ne!(nonterm1a, nonterm2);
        assert_ne!(nonterm1a, nonterm3);
        assert_ne!(nonterm1a, nonterm4);
    }

    #[test]
    fn test_symbol_hashing() {
        let mut set = HashSet::new();
        let term1a = Symbol::terminal(0, b'A', '+');
        let term1b = Symbol::terminal(0, b'A', '+');
        let term2 = Symbol::terminal(1, b'A', '+');
        let nonterm1 = Symbol::non_terminal(5, 100, '-');
        let nonterm2 = Symbol::non_terminal(5, 100, '-');

        assert!(set.insert(term1a));
        assert!(!set.insert(term1b)); // Should not insert duplicate
        assert!(set.insert(term2));
        assert!(set.insert(nonterm1));
        assert!(!set.insert(nonterm2)); // Should not insert duplicate

        assert_eq!(set.len(), 3);
    }
} 