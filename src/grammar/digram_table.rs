use crate::grammar::symbol::{Symbol, SymbolType, Direction};
use crate::encode::dna_2bit::EncodedBase;
use dashmap;
use rayon::prelude::*;

/// The key used for indexing digrams in the digram table
/// Format: ((first_symbol_type, first_strand), (second_symbol_type, second_strand))
pub type DigramKey = ((SymbolType, Direction), (SymbolType, Direction));

/// Stores digram occurrences and handles canonicalization.
#[derive(Debug, Default)]
pub struct DigramTable {
    // Stores occurrences: Canonical Digram Key -> List of (position, original_digram_instance)
    // Position usually refers to the index of the first symbol in the sequence.
    // Storing the original Symbol pair allows replacement later.
    occurrences: dashmap::DashMap<DigramKey, Vec<(usize, (Symbol, Symbol))>>,
}

impl DigramTable {
    pub fn new() -> Self {
        Self {
            occurrences: dashmap::DashMap::new(),
        }
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

    /// Add digrams in parallel from a sequence
    pub fn add_digrams_from_sequence(&mut self, sequence: &[Symbol], reverse_aware: bool) {
        if sequence.len() < 2 {
            return;
        }

        if sequence.len() < 10000 {
            // For small sequences, process sequentially to avoid parallel overhead
            for i in 0..(sequence.len() - 1) {
                let sym1 = sequence[i];
                let sym2 = sequence[i + 1];
                self.add_digram(i, sym1, sym2, reverse_aware);
            }
        } else {
            // For large sequences, use parallel processing with rayon
            sequence.par_windows(2)
                .enumerate()
                .for_each(|(i, window)| {
                    let sym1 = window[0];
                    let sym2 = window[1];
                    let key = self.get_digram_key(sym1, sym2);
                    let canonical_key = if reverse_aware {
                        self.canonicalize_digram_key(key)
                    } else {
                        key
                    };
                    
                    // Thread-safe concurrent update using DashMap
                    self.occurrences
                        .entry(canonical_key)
                        .or_default()
                        .push((i, (sym1, sym2)));
                });
        }
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
                    let revcomp_bases_enc: Vec<EncodedBase> = forward_bases.iter().rev().map(|b| b.revcomp()).collect();

                    // Special case for AC/TG pair
                    // AC (on + strand) is canonical
                    if (base1 == EncodedBase(0) && base2 == EncodedBase(1) && strand1 == Direction::Forward) ||
                       (base1 == EncodedBase(3) && base2 == EncodedBase(2) && strand1 == Direction::Forward) {
                        return ((SymbolType::Terminal(EncodedBase(0)), Direction::Forward), (SymbolType::Terminal(EncodedBase(1)), Direction::Forward));
                    }

                    // For general cases, compare lexicographically using EncodedBase
                    if forward_bases <= revcomp_bases_enc {
                        // Forward is canonical or equal
                        key 
                    } else {
                        // Reverse complement is canonical
                        let flipped_strand = strand1.flip();
                        (
                            (SymbolType::Terminal(revcomp_bases_enc[0]), flipped_strand),
                            (SymbolType::Terminal(revcomp_bases_enc[1]), flipped_strand),
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
    pub fn find_most_frequent_digram(&self) -> Option<(DigramKey, usize, Vec<(usize, (Symbol, Symbol))>)> {
        let mut max_key = None;
        let mut max_count = 0;
        let mut max_occurrences = Vec::new();

        // Iterate over entries and find the one with the most occurrences
        for entry in self.occurrences.iter() {
            let key = *entry.key();
            let occurrences = entry.value();
            let count = occurrences.len();

            if count > max_count {
                max_key = Some(key);
                max_count = count;
                max_occurrences = occurrences.clone();
            }
        }

        max_key.map(|key| (key, max_count, max_occurrences))
    }

    /// Removes occurrences associated with a specific digram instance (by position).
    /// Used after a digram is replaced by a rule.
    /// Needs careful handling if multiple digrams overlap the same position.
    pub fn remove_occurrence(&mut self, canonical_key: &DigramKey, position_to_remove: usize) -> bool {
        if let Some(mut entry) = self.occurrences.get_mut(canonical_key) {
            let occurrences = entry.value_mut();
            let initial_len = occurrences.len();
            
            // Remove occurrences starting at the specified position.
            occurrences.retain(|(pos, _)| *pos != position_to_remove);
            let removed = occurrences.len() < initial_len; // Check if anything was removed
            
            // If the vector became empty, remove the entry
            if occurrences.is_empty() {
                drop(entry); // Drop the mutable reference first
                self.occurrences.remove(canonical_key);
            }
            
            removed
        } else {
            false
        }
    }

    /// Merges another DigramTable into this one.
    /// Used for combining results from parallel processing.
    pub fn merge_with(&mut self, other: &DigramTable) {
        for entry in other.occurrences.iter() {
            let key = *entry.key();
            let other_occurrences = entry.value().clone();
            
            // Get or create the entry in this table
            if let Some(mut entry) = self.occurrences.get_mut(&key) {
                // Append the occurrences
                entry.value_mut().extend_from_slice(&other_occurrences);
            } else {
                // Insert a new entry
                self.occurrences.insert(key, other_occurrences);
            }
        }
    }

    /// Get the number of distinct digrams in the table
    pub fn len(&self) -> usize {
        self.occurrences.len()
    }

    /// Check if the table is empty
    pub fn is_empty(&self) -> bool {
        self.occurrences.is_empty()
    }

    /// Clear all entries from the table
    pub fn clear(&mut self) {
        self.occurrences.clear();
    }
}

#[cfg(test)]
mod tests {
    use super::*; // Imports items from the parent module (digram_table.rs)
    use crate::grammar::symbol::{Symbol, Direction, SymbolType}; // Ensure Direction is imported
    use crate::encode::dna_2bit::EncodedBase; // Import EncodedBase

    #[test]
    fn test_add_and_find_digram() {
        let mut table = DigramTable::new();
        let s1 = Symbol::terminal(0, EncodedBase(0), Direction::Forward); // A+
        let s2 = Symbol::terminal(1, EncodedBase(1), Direction::Forward); // C+
        let s3 = Symbol::terminal(2, EncodedBase(2), Direction::Forward); // G+

        table.add_digram(0, s1, s2, true); // AC @ 0
        table.add_digram(1, s2, s3, true); // CG @ 1
        table.add_digram(5, s1, s2, true); // AC @ 5

        let key_ac = ((SymbolType::Terminal(EncodedBase(0)), Direction::Forward), (SymbolType::Terminal(EncodedBase(1)), Direction::Forward));
        let key_cg = ((SymbolType::Terminal(EncodedBase(1)), Direction::Forward), (SymbolType::Terminal(EncodedBase(2)), Direction::Forward));

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
        let s_a = Symbol::terminal(0, EncodedBase(0), Direction::Forward); // A+
        let s_c = Symbol::terminal(1, EncodedBase(1), Direction::Forward); // C+
        // TG -> CA (revcomp) -> AC (canonical)
        let s_t = Symbol::terminal(2, EncodedBase(3), Direction::Forward); // T+
        let s_g = Symbol::terminal(3, EncodedBase(2), Direction::Forward); // G+

        table.add_digram(0, s_a, s_c, true); // AC @ 0.
        table.add_digram(5, s_t, s_g, true); // TG @ 5.

        let key_ac = ((SymbolType::Terminal(EncodedBase(0)), Direction::Forward), (SymbolType::Terminal(EncodedBase(1)), Direction::Forward));

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
        let s_a = Symbol::terminal(0, EncodedBase(0), Direction::Forward); // A+
        let s_c = Symbol::terminal(1, EncodedBase(1), Direction::Forward); // C+
        let s_t = Symbol::terminal(2, EncodedBase(3), Direction::Forward); // T+
        let s_g = Symbol::terminal(3, EncodedBase(2), Direction::Forward); // G+

        table.add_digram(0, s_a, s_c, false); // AC, reverse_aware=false
        table.add_digram(5, s_t, s_g, false); // TG, reverse_aware=false

        let key_ac = ((SymbolType::Terminal(EncodedBase(0)), Direction::Forward), (SymbolType::Terminal(EncodedBase(1)), Direction::Forward));
        let key_tg = ((SymbolType::Terminal(EncodedBase(3)), Direction::Forward), (SymbolType::Terminal(EncodedBase(2)), Direction::Forward));

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
        let s1 = Symbol::terminal(0, EncodedBase(0), Direction::Forward); // A+
        let s2 = Symbol::terminal(1, EncodedBase(1), Direction::Forward); // C+

        table.add_digram(0, s1, s2, true); // AC @ 0
        table.add_digram(5, s1, s2, true); // AC @ 5
        table.add_digram(10, s1, s2, true); // AC @ 10
        
        let key_ac = ((SymbolType::Terminal(EncodedBase(0)), Direction::Forward), (SymbolType::Terminal(EncodedBase(1)), Direction::Forward));

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
        let s_a = Symbol::terminal(0, EncodedBase(0), Direction::Forward); // A+
        let s_nt1 = Symbol::non_terminal(1, 100, Direction::Forward); // R100+
        let s_nt2 = Symbol::non_terminal(2, 101, Direction::Reverse); // R101-

        table.add_digram(0, s_a, s_nt1, true); // A-R100(+) @ 0
        table.add_digram(1, s_nt1, s_nt2, true); // R100(+)-R101(-) @ 1

        let key_a_nt1 = ((SymbolType::Terminal(EncodedBase(0)), Direction::Forward), (SymbolType::NonTerminal(100), Direction::Forward));
        let key_nt1_nt2 = ((SymbolType::NonTerminal(100), Direction::Forward), (SymbolType::NonTerminal(101), Direction::Reverse));

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