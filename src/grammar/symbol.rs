use std::fmt;
use std::cmp::Ordering;
use serde::{Serialize, Deserialize};
use crate::encode::dna_2bit::EncodedBase;
use std::hash::{Hash, Hasher};

/// Represents the strand direction in DNA sequences.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Ord, PartialOrd, Serialize, Deserialize)]
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
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Ord, PartialOrd, Serialize, Deserialize)]
pub enum SymbolType {
    Terminal(EncodedBase), // Represents a single encoded base
    NonTerminal(usize), // Represents a rule ID
}

/// Represents a symbol in a grammar, either a terminal (DNA base) or
/// a non-terminal (reference to a rule).
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct Symbol {
    pub id: usize, // Instance ID of the symbol in a sequence
    pub symbol_type: SymbolType,
    pub strand: Direction,
}

impl Symbol {
    /// Create a new terminal symbol.
    pub fn terminal(id: usize, base: EncodedBase, strand: Direction) -> Self {
        Symbol::new(id, SymbolType::Terminal(base), strand)
    }
    
    /// Create a new non-terminal symbol.
    pub fn non_terminal(id: usize, rule_id: usize, strand: Direction) -> Self {
        Symbol::new(id, SymbolType::NonTerminal(rule_id), strand)
    }
    
    /// Get the reverse complement of this symbol.
    pub fn reverse_complement(&self) -> Self {
        match self.symbol_type {
            SymbolType::Terminal(base) => Symbol {
                id: self.id, // Keep the same instance ID for now, or use a convention?
                symbol_type: SymbolType::Terminal(base.revcomp()),
                strand: self.strand.flip(),
            },
            SymbolType::NonTerminal(rule_id) => Symbol {
                id: self.id,
                symbol_type: SymbolType::NonTerminal(rule_id), // Rule ID doesn't change
                strand: self.strand.flip(),
            },
        }
    }

    pub fn new(id: usize, symbol_type: SymbolType, strand: Direction) -> Self {
        Symbol { id, symbol_type, strand }
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

// For canonical keys, we only care about type and strand, not instance ID.
#[derive(Debug, Clone, Copy, Eq, PartialEq, Hash, Ord, PartialOrd, Serialize, Deserialize)]
pub struct SymbolKey {
    pub symbol_type: SymbolType, // Type (Terminal base or NonTerminal rule_id)
    pub direction: Direction,    // Strand
}

impl SymbolKey {
    pub fn new(symbol: &Symbol) -> Self {
        SymbolKey {
            symbol_type: symbol.symbol_type, // This already ignores id
            direction: symbol.strand,
        }
    }
}

// Implementations for Symbol that ignore 'id' for comparison, ordering, and hashing.
impl PartialEq for Symbol {
    fn eq(&self, other: &Self) -> bool {
        self.symbol_type == other.symbol_type && self.strand == other.strand
    }
}
impl Eq for Symbol {}

impl PartialOrd for Symbol {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Symbol {
    fn cmp(&self, other: &Self) -> Ordering {
        (self.symbol_type, self.strand).cmp(&(other.symbol_type, other.strand))
    }
}

impl Hash for Symbol {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.symbol_type.hash(state);
        self.strand.hash(state);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashMap;

    fn create_terminal(id: usize, base_val: u8, strand: Direction) -> Symbol {
        Symbol::terminal(id, EncodedBase(base_val), strand)
    }

    fn create_non_terminal(id: usize, rule_id: usize, strand: Direction) -> Symbol {
        Symbol::non_terminal(id, rule_id, strand)
    }

    #[test]
    fn test_symbol_creation() {
        let term_fwd = create_terminal(0, 0, Direction::Forward);
        assert_eq!(term_fwd.id, 0);
        assert_eq!(term_fwd.symbol_type, SymbolType::Terminal(EncodedBase(0)));
        assert_eq!(term_fwd.strand, Direction::Forward);

        let nt_rev = create_non_terminal(1, 10, Direction::Reverse);
        assert_eq!(nt_rev.id, 1);
        assert_eq!(nt_rev.symbol_type, SymbolType::NonTerminal(10));
        assert_eq!(nt_rev.strand, Direction::Reverse);
    }

    #[test]
    fn test_symbol_equality() {
        let sym1 = create_terminal(1, 0, Direction::Forward); // A+ id 1
        let sym2 = create_terminal(2, 0, Direction::Forward); // A+ id 2
        let sym3 = create_terminal(1, 1, Direction::Forward); // C+ id 1
        let sym4 = create_terminal(1, 0, Direction::Reverse); // A- id 1

        // Symbols with different IDs but same type and strand ARE now equal
        assert_eq!(sym1, sym2, "Symbols with different IDs but same type/strand should be equal");

        assert_ne!(sym1, sym3, "Symbols with different bases should not be equal");
        assert_ne!(sym1, sym4, "Symbols with different strands should not be equal");

        let nt1 = create_non_terminal(5, 100, Direction::Forward);
        let nt2 = create_non_terminal(6, 100, Direction::Forward);
        let nt3 = create_non_terminal(5, 101, Direction::Forward);

        // Non-terminals with different IDs but same rule_id and strand ARE now equal
        assert_eq!(nt1, nt2, "Non-terminals with different IDs but same rule/strand should be equal");
        assert_ne!(nt1, nt3, "Non-terminals with different rule_ids should not be equal");
    }

    #[test]
    fn test_reverse_complement() {
        let term_a_fwd = create_terminal(0, 0, Direction::Forward); // A+
        let term_t_rev = create_terminal(0, 3, Direction::Reverse); // T-
        assert_eq!(term_a_fwd.reverse_complement(), term_t_rev);

        // Check if ID is preserved by reverse_complement (it should be, for now)
        let term_g_fwd_id5 = create_terminal(5, 2, Direction::Forward); // G+ id 5
        let term_c_rev_id5 = create_terminal(5, 1, Direction::Reverse); // C- id 5
        assert_eq!(term_g_fwd_id5.reverse_complement(), term_c_rev_id5);
        assert_eq!(term_g_fwd_id5.reverse_complement().id, 5);

        let nt_fwd = create_non_terminal(1, 10, Direction::Forward);
        let nt_rev = create_non_terminal(1, 10, Direction::Reverse);
        assert_eq!(nt_fwd.reverse_complement(), nt_rev);
        assert_eq!(nt_fwd.reverse_complement().id, 1);
    }

    #[test]
    fn test_symbol_hashing() {
        let sym1 = create_terminal(1, 0, Direction::Forward); // A+ id 1
        let sym2 = create_terminal(2, 0, Direction::Forward); // A+ id 2 (same content as sym1 for hashing)
        let sym3 = create_terminal(3, 1, Direction::Forward); // C+ id 3

        let mut map = HashMap::new();
        map.insert(sym1, "value1");

        // sym2 should map to the same entry as sym1 because their hash is the same
        assert_eq!(map.get(&sym2), Some(&"value1"));
        // sym3 is different and should not be found / map to value1
        assert!(map.get(&sym3).is_none() || map.get(&sym3) != Some(&"value1"));

        // Test with non-terminals
        let nt1 = create_non_terminal(10, 100, Direction::Forward);
        let nt2 = create_non_terminal(20, 100, Direction::Forward); // Same content as nt1 for hashing
        map.clear();
        map.insert(nt1, "value_nt1");
        assert_eq!(map.get(&nt2), Some(&"value_nt1"));
    }

    #[test]
    fn test_symbol_comparison() {
        // Order: Terminal < NonTerminal
        // For Terminals: by base, then by strand
        // For NonTerminals: by rule_id, then by strand
        // ID is ignored.

        let t_a_f_1 = create_terminal(1, 0, Direction::Forward); // A+
        let t_a_f_2 = create_terminal(2, 0, Direction::Forward); // A+ (same as t_a_f_1 for comparison)
        let t_c_f_1 = create_terminal(1, 1, Direction::Forward); // C+
        let t_a_r_1 = create_terminal(1, 0, Direction::Reverse); // A-

        let nt_10_f_1 = create_non_terminal(1, 10, Direction::Forward); // R10+
        let nt_10_f_2 = create_non_terminal(2, 10, Direction::Forward); // R10+ (same as nt_10_f_1)
        let nt_20_f_1 = create_non_terminal(1, 20, Direction::Forward); // R20+
        let nt_10_r_1 = create_non_terminal(1, 10, Direction::Reverse); // R10-

        // Equality (ID ignored)
        assert_eq!(t_a_f_1.cmp(&t_a_f_2), Ordering::Equal);
        assert_eq!(nt_10_f_1.cmp(&nt_10_f_2), Ordering::Equal);

        // Terminals vs Terminals
        assert_eq!(t_a_f_1.cmp(&t_c_f_1), Ordering::Less, "A+ < C+");
        assert_eq!(t_c_f_1.cmp(&t_a_f_1), Ordering::Greater, "C+ > A+");
        assert_eq!(t_a_r_1.cmp(&t_a_f_1), Ordering::Greater, "A- > A+ (Reverse > Forward)");
        assert_eq!(t_a_f_1.cmp(&t_a_r_1), Ordering::Less, "A+ < A-");

        // NonTerminals vs NonTerminals
        assert_eq!(nt_10_f_1.cmp(&nt_20_f_1), Ordering::Less, "R10+ < R20+");
        assert_eq!(nt_20_f_1.cmp(&nt_10_f_1), Ordering::Greater, "R20+ > R10+");
        assert_eq!(nt_10_r_1.cmp(&nt_10_f_1), Ordering::Greater, "R10- > R10+ (Reverse > Forward)");
        assert_eq!(nt_10_f_1.cmp(&nt_10_r_1), Ordering::Less, "R10+ < R10-");

        // Terminals vs NonTerminals
        assert_eq!(t_a_f_1.cmp(&nt_10_f_1), Ordering::Less, "Terminal < NonTerminal");
        assert_eq!(nt_10_f_1.cmp(&t_a_f_1), Ordering::Greater, "NonTerminal > Terminal");
    }
} 