// Placeholder for hashing utilities

use crate::grammar::symbol::Symbol;
use twox_hash::XxHash64;
use std::hash::{Hash, Hasher};

/// Hash a vector of symbols using xxHash algorithm for fast hashing
pub fn hash_symbols(symbols: &[Symbol]) -> u64 {
    let mut hasher = XxHash64::default();
    for symbol in symbols {
        symbol.hash(&mut hasher);
    }
    hasher.finish()
}

/// Hash a DNA sequence represented as bytes
pub fn hash_sequence(sequence: &[u8]) -> u64 {
    let mut hasher = XxHash64::default();
    sequence.hash(&mut hasher);
    hasher.finish()
}

/// Hash a sequence of symbols using XXH3
/// This is a fast non-cryptographic hash function suitable for rule hashing
pub fn hash_symbols_xxh3(symbols: &[Symbol]) -> u64 {
    let mut hasher = twox_hash::xxh3::Hash64::default();
    symbols.hash(&mut hasher);
    hasher.finish()
}

/// Hash a sequence of symbols considering reverse complement equivalence
/// Returns the canonical hash (minimum of forward and reverse complement)
pub fn canonical_hash_symbols(symbols: &[Symbol]) -> u64 {
    // Hash the forward sequence
    let forward_hash = hash_symbols_xxh3(symbols);
    
    // Create the reverse complement sequence
    let revcomp: Vec<Symbol> = symbols.iter().rev().map(|s| s.reverse_complement()).collect();
    
    // Hash the reverse complement
    let reverse_hash = hash_symbols_xxh3(&revcomp);
    
    // Return the minimum hash as canonical
    std::cmp::min(forward_hash, reverse_hash)
}

/// A trait for types that can be canonicalized (have a canonical representation)
pub trait Canonicalizable {
    /// Returns the canonical representation of self
    fn canonical(&self) -> Self;
}

/// Helper for creating a canonical representation of a symbol pair (digram)
/// considering reverse complement
pub fn canonical_digram(sym1: Symbol, sym2: Symbol) -> (Symbol, Symbol) {
    // Original digram
    let forward = (sym1, sym2);
    
    // Reverse complement digram (in reverse order)
    let reverse = (sym2.reverse_complement(), sym1.reverse_complement());
    
    // Use lexicographic comparison to determine canonical form
    if forward < reverse {
        forward
    } else {
        reverse
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::encode::dna_2bit::EncodedBase;
    use crate::grammar::symbol::{Symbol, Direction};
    
    #[test]
    fn test_hash_symbols() {
        let symbols = vec![
            Symbol::terminal(0, EncodedBase(0), Direction::Forward),
            Symbol::terminal(1, EncodedBase(1), Direction::Forward),
            Symbol::non_terminal(2, 1, Direction::Forward),
        ];
        
        let hash1 = hash_symbols(&symbols);
        let hash2 = hash_symbols(&symbols);
        assert_eq!(hash1, hash2, "Same input should produce same hash");
        
        let different_symbols = vec![
            Symbol::terminal(0, EncodedBase(2), Direction::Forward),
            Symbol::terminal(1, EncodedBase(3), Direction::Forward),
        ];
        
        let hash3 = hash_symbols(&different_symbols);
        assert_ne!(hash1, hash3, "Different input should produce different hash");
    }
    
    #[test]
    fn test_hash_sequence() {
        let seq1 = b"ACGTACGT";
        let seq2 = b"ACGTACGT";
        let seq3 = b"ACGTTCGT"; // Different
        
        let hash1 = hash_sequence(seq1);
        let hash2 = hash_sequence(seq2);
        let hash3 = hash_sequence(seq3);
        
        assert_eq!(hash1, hash2, "Same input should produce same hash");
        assert_ne!(hash1, hash3, "Different input should produce different hash");
    }
    
    #[test]
    fn test_canonical_hash_symbols() {
        let base_a = EncodedBase::from_base(b'A').unwrap();
        let base_t = EncodedBase::from_base(b'T').unwrap();
        
        // A+ T+ should be equivalent to A- T- (reverse complement)
        let sym_a_fwd = Symbol::terminal(0, base_a, Direction::Forward);
        let sym_t_fwd = Symbol::terminal(1, base_t, Direction::Forward);
        let sym_a_rev = sym_a_fwd.reverse_complement();
        let sym_t_rev = sym_t_fwd.reverse_complement();
        
        let symbols1 = vec![sym_a_fwd, sym_t_fwd];
        let symbols2 = vec![sym_t_rev, sym_a_rev];
        
        // These should produce the same canonical hash
        assert_eq!(canonical_hash_symbols(&symbols1), canonical_hash_symbols(&symbols2));
    }
    
    #[test]
    fn test_canonical_digram() {
        let base_a = EncodedBase::from_base(b'A').unwrap();
        let base_t = EncodedBase::from_base(b'T').unwrap();
        let base_c = EncodedBase::from_base(b'C').unwrap();
        let base_g = EncodedBase::from_base(b'G').unwrap();
        
        let sym_a_fwd = Symbol::terminal(0, base_a, Direction::Forward);
        let sym_t_fwd = Symbol::terminal(1, base_t, Direction::Forward);
        let sym_c_fwd = Symbol::terminal(2, base_c, Direction::Forward);
        let sym_g_fwd = Symbol::terminal(3, base_g, Direction::Forward);
        
        // Test AT/TA case
        let digram1 = (sym_a_fwd, sym_t_fwd);
        let canonical1 = canonical_digram(sym_a_fwd, sym_t_fwd);
        let revcomp1 = (sym_t_fwd.reverse_complement(), sym_a_fwd.reverse_complement());
        
        // The canonical form should be one of the two
        assert!(canonical1 == digram1 || canonical1 == revcomp1);
        
        // Test CG/GC case
        let digram2 = (sym_c_fwd, sym_g_fwd);
        let canonical2 = canonical_digram(sym_c_fwd, sym_g_fwd);
        let revcomp2 = (sym_g_fwd.reverse_complement(), sym_c_fwd.reverse_complement());
        
        // The canonical form should be one of the two
        assert!(canonical2 == digram2 || canonical2 == revcomp2);
    }
} 