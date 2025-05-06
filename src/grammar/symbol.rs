use std::fmt;
use std::cmp::Ordering;
use serde::{Serialize, Deserialize};
use crate::encode::dna_2bit::EncodedBase;

/// Represents the strand direction in DNA sequences.
#[derive(Debug, Clone, Copy, Eq, Hash, Serialize, Deserialize, PartialEq, PartialOrd)]
pub enum Direction {
    Forward,
    Reverse,
}

impl Direction {
    /// Convert a strand character ('+' or '-') to a Direction.
    pub fn from_char(c: char) -> Self {
        match c {
            '-' => Direction::Reverse,
            _ => Direction::Forward, // Default to forward for anything but '-'
        }
    }
    
    /// Convert Direction to a strand character ('+' or '-').
    pub fn to_char(&self) -> char {
        match self {
            Direction::Forward => '+',
            Direction::Reverse => '-',
        }
    }
    
    /// Get the opposite direction.
    pub fn flip(&self) -> Self {
        match self {
            Direction::Forward => Direction::Reverse,
            Direction::Reverse => Direction::Forward,
        }
    }
}

// Implement Display for Direction
impl fmt::Display for Direction {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.to_char())
    }
}

/// Type of symbol in a grammar.
#[derive(Debug, Clone, Copy, Eq, Hash, Serialize, Deserialize, PartialEq, PartialOrd)]
pub enum SymbolType {
    Terminal(EncodedBase), // Represents a single encoded base
    NonTerminal(usize), // Represents a rule ID
}

/// Represents a symbol in a grammar, either a terminal (DNA base) or
/// a non-terminal (reference to a rule).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct Symbol {
    pub id: usize, // Unique identifier for this *instance* of the symbol (e.g., its position in the original sequence or a rule expansion)
    pub symbol_type: SymbolType,
    pub strand: Direction, // Direction (Forward or Reverse)
}

impl PartialOrd for Symbol {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Symbol {
    fn cmp(&self, other: &Self) -> Ordering {
        // First compare by symbol type
        match (&self.symbol_type, &other.symbol_type) {
            (SymbolType::Terminal(a), SymbolType::Terminal(b)) => {
                // For terminals, compare by base value
                match a.0.cmp(&b.0) {
                    Ordering::Equal => {},
                    ord => return ord,
                }
            },
            (SymbolType::NonTerminal(a), SymbolType::NonTerminal(b)) => {
                // For non-terminals, compare by rule ID
                match a.cmp(b) {
                    Ordering::Equal => {},
                    ord => return ord,
                }
            },
            (SymbolType::Terminal(_), SymbolType::NonTerminal(_)) => {
                // Terminals come before non-terminals
                return Ordering::Less;
            },
            (SymbolType::NonTerminal(_), SymbolType::Terminal(_)) => {
                // Non-terminals come after terminals
                return Ordering::Greater;
            },
        }
        
        // If symbol types are equal, compare by strand (Forward < Reverse)
        let strand_ord = match (&self.strand, &other.strand) {
            (Direction::Forward, Direction::Reverse) => Ordering::Less,
            (Direction::Reverse, Direction::Forward) => Ordering::Greater,
            _ => Ordering::Equal,
        };
        
        if strand_ord != Ordering::Equal {
            return strand_ord;
        }
        
        // Finally, compare by ID
        self.id.cmp(&other.id)
    }
}

impl Symbol {
    /// Create a new terminal symbol.
    pub fn terminal(id: usize, base: EncodedBase, strand: Direction) -> Self {
        Self {
            id,
            symbol_type: SymbolType::Terminal(base),
            strand,
        }
    }
    
    /// Create a new non-terminal symbol.
    pub fn non_terminal(id: usize, rule_id: usize, strand: Direction) -> Self {
        Self {
            id,
            symbol_type: SymbolType::NonTerminal(rule_id),
            strand,
        }
    }
    
    /// Get the reverse complement of this symbol.
    pub fn reverse_complement(&self) -> Self {
        match &self.symbol_type {
            SymbolType::Terminal(base) => {
                // For terminals, flip the base and strand
                Symbol {
                    id: self.id,
                    symbol_type: SymbolType::Terminal(base.revcomp()),
                    strand: self.strand.flip(),
                }
            },
            SymbolType::NonTerminal(rule_id) => {
                // For non-terminals, just flip the strand
                Symbol {
                    id: self.id,
                    symbol_type: SymbolType::NonTerminal(*rule_id),
                    strand: self.strand.flip(),
                }
            }
        }
    }
}

impl fmt::Display for Symbol {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match &self.symbol_type {
            SymbolType::Terminal(base) => {
                write!(f, "{}{}", base.to_char(), self.strand.to_char())
            },
            SymbolType::NonTerminal(rule_id) => {
                write!(f, "R{}{}", rule_id, self.strand.to_char())
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::encode::dna_2bit::EncodedBase;
    
    #[test]
    fn test_symbol_creation() {
        // Test terminal symbol
        let a_plus = Symbol::terminal(0, EncodedBase(0), Direction::Forward);
        assert_eq!(a_plus.id, 0);
        assert!(matches!(a_plus.symbol_type, SymbolType::Terminal(EncodedBase(0))));
        assert_eq!(a_plus.strand, Direction::Forward);
        
        // Test non-terminal symbol
        let rule_minus = Symbol::non_terminal(5, 10, Direction::Reverse);
        assert_eq!(rule_minus.id, 5);
        assert!(matches!(rule_minus.symbol_type, SymbolType::NonTerminal(10)));
        assert_eq!(rule_minus.strand, Direction::Reverse);
    }
    
    #[test]
    fn test_symbol_equality() {
        let sym1 = Symbol::terminal(1, EncodedBase(0), Direction::Forward);
        let sym2 = Symbol::terminal(1, EncodedBase(0), Direction::Forward);
        let sym3 = Symbol::terminal(2, EncodedBase(0), Direction::Forward);
        let sym4 = Symbol::terminal(1, EncodedBase(1), Direction::Forward);
        let sym5 = Symbol::terminal(1, EncodedBase(0), Direction::Reverse);
        
        assert_eq!(sym1, sym2);
        assert_ne!(sym1, sym3); // Different ID
        assert_ne!(sym1, sym4); // Different base
        assert_ne!(sym1, sym5); // Different strand
        
        let nt1 = Symbol::non_terminal(1, 5, Direction::Forward);
        let nt2 = Symbol::non_terminal(1, 5, Direction::Forward);
        let nt3 = Symbol::non_terminal(2, 5, Direction::Forward);
        let nt4 = Symbol::non_terminal(1, 6, Direction::Forward);
        let nt5 = Symbol::non_terminal(1, 5, Direction::Reverse);
        
        assert_eq!(nt1, nt2);
        assert_ne!(nt1, nt3); // Different ID
        assert_ne!(nt1, nt4); // Different rule
        assert_ne!(nt1, nt5); // Different strand
    }
    
    #[test]
    fn test_symbol_hashing() {
        use std::collections::HashMap;
        
        let mut map = HashMap::new();
        
        let sym1 = Symbol::terminal(1, EncodedBase(0), Direction::Forward);
        let sym2 = Symbol::terminal(1, EncodedBase(0), Direction::Forward);
        
        map.insert(sym1, "symbol1");
        assert_eq!(map.get(&sym2), Some(&"symbol1"));
        
        // Different symbols should have different hash values
        let sym3 = Symbol::terminal(2, EncodedBase(0), Direction::Forward);
        assert!(map.get(&sym3).is_none());
    }
    
    #[test]
    fn test_reverse_complement() {
        // A+ -> T-
        let sym_a_plus = Symbol::terminal(1, EncodedBase(0), Direction::Forward);
        let rev_comp = sym_a_plus.reverse_complement();
        
        assert_eq!(rev_comp.id, 1); // ID stays the same
        assert!(matches!(rev_comp.symbol_type, SymbolType::Terminal(EncodedBase(3)))); // A -> T
        assert_eq!(rev_comp.strand, Direction::Reverse); // + -> -
        
        // For non-terminal
        let nt = Symbol::non_terminal(2, 5, Direction::Forward);
        let nt_rev = nt.reverse_complement();
        
        assert_eq!(nt_rev.id, 2);
        assert!(matches!(nt_rev.symbol_type, SymbolType::NonTerminal(5)));
        assert_eq!(nt_rev.strand, Direction::Reverse);
    }
    
    #[test]
    fn test_symbol_comparison() {
        // Terminal < NonTerminal
        let term = Symbol::terminal(1, EncodedBase(0), Direction::Forward);
        let nonterm = Symbol::non_terminal(1, 5, Direction::Forward);
        assert!(term < nonterm);
        
        // Compare terminals by base
        let a = Symbol::terminal(1, EncodedBase(0), Direction::Forward);
        let c = Symbol::terminal(1, EncodedBase(1), Direction::Forward);
        assert!(a < c);
        
        // Compare non-terminals by rule ID
        let r1 = Symbol::non_terminal(1, 1, Direction::Forward);
        let r2 = Symbol::non_terminal(1, 2, Direction::Forward);
        assert!(r1 < r2);
        
        // Compare by strand (Forward < Reverse)
        let a_fwd = Symbol::terminal(1, EncodedBase(0), Direction::Forward);
        let a_rev = Symbol::terminal(1, EncodedBase(0), Direction::Reverse);
        assert!(a_fwd < a_rev);
        
        // Compare by ID if everything else is equal
        let id1 = Symbol::terminal(1, EncodedBase(0), Direction::Forward);
        let id2 = Symbol::terminal(2, EncodedBase(0), Direction::Forward);
        assert!(id1 < id2);
    }
} 