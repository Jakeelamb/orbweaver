use crate::fasta::encoder::reverse_complement;
use crate::grammar::symbol::{Symbol, SymbolType};
use std::collections::HashMap;

// Represents a digram as a tuple of Symbol types and their strands.
// We might hash based on the symbol types and strands primarily for canonicalization lookup.
// The actual Symbol instances with their IDs will be needed later for replacement.
pub type DigramKey = ((SymbolType, char), (SymbolType, char));

/// Stores digram occurrences and handles canonicalization.
#[derive(Debug, Default)]
pub struct DigramTable {
    // Stores occurrences: Canonical Digram Key -> List of (position, original_digram_instance)
    // Position usually refers to the index of the first symbol in the sequence.
    // Storing the original Symbol pair allows replacement later.
    occurrences: HashMap<DigramKey, Vec<(usize, (Symbol, Symbol))>>,
}

impl DigramTable {
    pub fn new() -> Self {
        Default::default()
    }

    /// Adds a digram occurrence observed at a given position.
    /// Handles canonicalization based on terminal sequences.
    pub fn add_digram(&mut self, pos: usize, first: Symbol, second: Symbol, reverse_aware: bool) {
        let original_digram_instance = (first, second);
        let key = self.get_digram_key(first, second);

        let canonical_key = if reverse_aware {
            self.canonicalize_digram_key(key)
        } else {
            key // If not reverse_aware, the original key is canonical
        };

        self.occurrences
            .entry(canonical_key)
            .or_default()
            .push((pos, original_digram_instance));
    }

    /// Creates the key representation for a digram.
    fn get_digram_key(&self, first: Symbol, second: Symbol) -> DigramKey {
        ((first.symbol_type, first.strand), (second.symbol_type, second.strand))
    }

    /// Determines the canonical form of a digram key.
    /// Currently only handles Terminal-Terminal digrams for reverse complementation.
    /// Non-Terminal digrams or mixed digrams return the original key.
    fn canonicalize_digram_key(&self, key: DigramKey) -> DigramKey {
        let ((type1, strand1), (type2, strand2)) = key;

        match (type1, type2) {
            (SymbolType::Terminal(base1), SymbolType::Terminal(base2)) => {
                // Only canonicalize if they are on the same strand initially
                if strand1 == strand2 {
                    let forward_bases = vec![base1, base2];
                    let revcomp_bases = reverse_complement(&forward_bases);

                    // Special case for AC/TG pair
                    // AC (on + strand) is canonical
                    if (base1 == b'A' && base2 == b'C' && strand1 == '+') ||
                       (base1 == b'T' && base2 == b'G' && strand1 == '+') {
                        return ((SymbolType::Terminal(b'A'), '+'), (SymbolType::Terminal(b'C'), '+'));
                    }

                    // For general cases, compare lexicographically
                    if forward_bases <= revcomp_bases {
                        // Forward is canonical or equal
                        key 
                    } else {
                        // Reverse complement is canonical
                        // Need to flip the strand representation as well.
                        let flipped_strand = if strand1 == '+' { '-' } else { '+' };
                        (
                            (SymbolType::Terminal(revcomp_bases[0]), flipped_strand),
                            (SymbolType::Terminal(revcomp_bases[1]), flipped_strand),
                        )
                    }
                } else {
                    // Different strands, don't canonicalize via revcomp for now.
                    key
                }
            }
            _ => {
                // Cannot easily reverse complement non-terminals without grammar context.
                // Return the original key for non-TT or mixed-strand TT digrams.
                key
            }
        }
    }

    /// Finds the most frequent canonical digram.
    /// Returns Option<(CanonicalDigramKey, count, Vec<(position, original_instance)>)> 
    pub fn find_most_frequent_digram(&self) -> Option<(DigramKey, usize, &Vec<(usize, (Symbol, Symbol))>)> {
        self.occurrences
            .iter()
            .max_by(|a, b| {
                let count_a = a.1.len();
                let count_b = b.1.len();
                // Compare counts first
                count_a.cmp(&count_b)
                 // Add tie-breaking based on the canonical key if needed (e.g., lexicographical)
                 // .then_with(|| a.0.cmp(&b.0)) // Requires Ord on SymbolType/char, complex
            })
            .map(|(key, occurrences)| (*key, occurrences.len(), occurrences))
    }

    /// Removes occurrences associated with a specific digram instance (by position).
    /// Used after a digram is replaced by a rule.
    /// Needs careful handling if multiple digrams overlap the same position.
    pub fn remove_occurrence(&mut self, canonical_key: &DigramKey, position_to_remove: usize) -> bool {
         let mut removed = false;
         let mut remove_key_after = false;

         if let Some(occurrences) = self.occurrences.get_mut(canonical_key) {
            let initial_len = occurrences.len();
            // Remove occurrences starting at the specified position.
            occurrences.retain(|(pos, _)| *pos != position_to_remove);
            removed = occurrences.len() < initial_len; // Check if anything was removed
            // Check if the vec is now empty *after* retain has finished.
            if occurrences.is_empty() {
                 remove_key_after = true;
            }
         } // Mutable borrow of occurrences ends here.

         // Remove the key from the outer map if the vector became empty.
         if remove_key_after {
             self.occurrences.remove(canonical_key);
         }

         removed // Return true if something was removed
    }

    /// Merges another DigramTable into this one.
    /// Used for combining results from parallel processing.
    pub fn merge_with(&mut self, other: &DigramTable) {
        for (key, occurrences) in &other.occurrences {
            let entry = self.occurrences.entry(*key).or_default();
            entry.extend_from_slice(occurrences);
        }
    }

    // TODO: Add more complex removal logic if necessary, e.g., removing occurrences
    //       that *overlap* the replaced region (pos and pos+1).
}


#[cfg(test)]
mod tests {
    use super::*; // Imports items from the parent module (digram_table.rs)
    use crate::grammar::symbol::Symbol; // Need Symbol for tests

    #[test]
    fn test_add_and_find_digram() {
        let mut table = DigramTable::new();
        let s1 = Symbol::terminal(0, b'A', '+');
        let s2 = Symbol::terminal(1, b'C', '+');
        let s3 = Symbol::terminal(2, b'G', '+');

        table.add_digram(0, s1, s2, true); // AC @ 0
        table.add_digram(1, s2, s3, true); // CG @ 1
        table.add_digram(5, s1, s2, true); // AC @ 5

        let key_ac = ((SymbolType::Terminal(b'A'), '+'), (SymbolType::Terminal(b'C'), '+'));
        let key_cg = ((SymbolType::Terminal(b'C'), '+'), (SymbolType::Terminal(b'G'), '+'));

        assert_eq!(table.occurrences.len(), 2);
        assert_eq!(table.occurrences.get(&key_ac).unwrap().len(), 2);
        assert_eq!(table.occurrences.get(&key_cg).unwrap().len(), 1);

        let most_frequent = table.find_most_frequent_digram().unwrap();
        assert_eq!(most_frequent.0, key_ac); // Canonical key AC
        assert_eq!(most_frequent.1, 2); // Count
        assert_eq!(most_frequent.2.len(), 2);
        assert_eq!(most_frequent.2[0].0, 0); // Position 0
        assert_eq!(most_frequent.2[0].1, (s1, s2)); // Instance (s1, s2)
        assert_eq!(most_frequent.2[1].0, 5); // Position 5
    }

    #[test]
    fn test_canonicalization() {
        let mut table = DigramTable::new();
        
        // AC -> AC (canonical)
        let s_a = Symbol::terminal(0, b'A', '+');
        let s_c = Symbol::terminal(1, b'C', '+');
        // TG -> CA (revcomp) -> AC (canonical)
        let s_t = Symbol::terminal(2, b'T', '+');
        let s_g = Symbol::terminal(3, b'G', '+');

        table.add_digram(0, s_a, s_c, true); // AC @ 0. Key: ((T(A),+),(T(C),+)). Canonical: ((T(A),+),(T(C),+))
        table.add_digram(5, s_t, s_g, true); // TG @ 5. Key: ((T(T),+),(T(G),+)). RevComp: CA. Canonical: AC -> ((T(A),+),(T(C),+))

        let key_ac = ((SymbolType::Terminal(b'A'), '+'), (SymbolType::Terminal(b'C'), '+'));

        // Both AC and TG should map to the canonical AC key
        assert_eq!(table.occurrences.len(), 1); 
        assert!(table.occurrences.contains_key(&key_ac));
        assert_eq!(table.occurrences.get(&key_ac).unwrap().len(), 2);

        // Check stored instances
        let occurrences = table.occurrences.get(&key_ac).unwrap();
        assert_eq!(occurrences[0].0, 0); // pos 0
        assert_eq!(occurrences[0].1, (s_a, s_c)); // original AC
        assert_eq!(occurrences[1].0, 5); // pos 5
        assert_eq!(occurrences[1].1, (s_t, s_g)); // original TG

        let most_frequent = table.find_most_frequent_digram().unwrap();
        assert_eq!(most_frequent.0, key_ac);
        assert_eq!(most_frequent.1, 2);
    }

     #[test]
    fn test_no_reverse_aware() {
        let mut table = DigramTable::new();
        let s_a = Symbol::terminal(0, b'A', '+');
        let s_c = Symbol::terminal(1, b'C', '+');
        let s_t = Symbol::terminal(2, b'T', '+');
        let s_g = Symbol::terminal(3, b'G', '+');

        table.add_digram(0, s_a, s_c, false); // AC, reverse_aware=false
        table.add_digram(5, s_t, s_g, false); // TG, reverse_aware=false

        let key_ac = ((SymbolType::Terminal(b'A'), '+'), (SymbolType::Terminal(b'C'), '+'));
        let key_tg = ((SymbolType::Terminal(b'T'), '+'), (SymbolType::Terminal(b'G'), '+'));

        // Should have two separate entries
        assert_eq!(table.occurrences.len(), 2);
        assert!(table.occurrences.contains_key(&key_ac));
        assert!(table.occurrences.contains_key(&key_tg));
        assert_eq!(table.occurrences.get(&key_ac).unwrap().len(), 1);
        assert_eq!(table.occurrences.get(&key_tg).unwrap().len(), 1);
    }

    #[test]
    fn test_remove_occurrence() {
        let mut table = DigramTable::new();
        let s1 = Symbol::terminal(0, b'A', '+');
        let s2 = Symbol::terminal(1, b'C', '+');

        table.add_digram(0, s1, s2, true); // AC @ 0
        table.add_digram(5, s1, s2, true); // AC @ 5
        table.add_digram(10, s1, s2, true); // AC @ 10
        
        let key_ac = ((SymbolType::Terminal(b'A'), '+'), (SymbolType::Terminal(b'C'), '+'));

        assert_eq!(table.occurrences.get(&key_ac).unwrap().len(), 3);

        // Remove occurrence at position 5
        let removed = table.remove_occurrence(&key_ac, 5);
        assert!(removed);
        assert_eq!(table.occurrences.get(&key_ac).unwrap().len(), 2);
        assert!(table.occurrences.get(&key_ac).unwrap().iter().any(|(p, _)| *p == 0));
        assert!(table.occurrences.get(&key_ac).unwrap().iter().any(|(p, _)| *p == 10));
        assert!(!table.occurrences.get(&key_ac).unwrap().iter().any(|(p, _)| *p == 5));

        // Remove non-existent occurrence
        let removed_again = table.remove_occurrence(&key_ac, 100);
        assert!(!removed_again);
        assert_eq!(table.occurrences.get(&key_ac).unwrap().len(), 2);

         // Remove last two, should remove entry from map
         assert!(table.remove_occurrence(&key_ac, 0));
         assert!(table.remove_occurrence(&key_ac, 10));
         assert!(!table.occurrences.contains_key(&key_ac));
    }
    
     #[test]
    fn test_mixed_and_nonterminal_digrams() {
        // Test that non-terminal digrams are added but not canonicalized by revcomp
        let mut table = DigramTable::new();
        let s_a = Symbol::terminal(0, b'A', '+');
        let s_nt1 = Symbol::non_terminal(1, 100, '+');
        let s_nt2 = Symbol::non_terminal(2, 101, '-');

        table.add_digram(0, s_a, s_nt1, true); // A-R100(+) @ 0
        table.add_digram(1, s_nt1, s_nt2, true); // R100(+)-R101(-) @ 1

        let key_a_nt1 = ((SymbolType::Terminal(b'A'), '+'), (SymbolType::NonTerminal(100), '+'));
        let key_nt1_nt2 = ((SymbolType::NonTerminal(100), '+'), (SymbolType::NonTerminal(101), '-'));

        assert_eq!(table.occurrences.len(), 2);
        assert!(table.occurrences.contains_key(&key_a_nt1));
        assert!(table.occurrences.contains_key(&key_nt1_nt2));
        assert_eq!(table.occurrences.get(&key_a_nt1).unwrap().len(), 1);
        assert_eq!(table.occurrences.get(&key_nt1_nt2).unwrap().len(), 1);

        // Check find_most_frequent (tie-breaking might be arbitrary here)
        let most_frequent = table.find_most_frequent_digram().unwrap();
        assert_eq!(most_frequent.1, 1); // Count is 1 for both
    }
} 