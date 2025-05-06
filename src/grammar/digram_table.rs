use crate::grammar::symbol::{Symbol, SymbolType, Direction};
use crate::encode::dna_2bit::EncodedBase;
use dashmap;
use rayon::prelude::*;
use std::collections::HashMap;
use std::hash::{Hash, Hasher};
use std::iter::FromIterator;
use std::sync::Arc;
use crate::utils::hash::custom_hash;

/// The key used for indexing digrams in the digram table
/// Format: ((first_symbol_type, first_strand), (second_symbol_type, second_strand))
pub type DigramKey = u64;

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
        Self {
            occurrences: HashMap::new(),
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
            // For large sequences, use parallel processing
            // First, collect the digrams and their canonical keys in parallel
            let digram_data: Vec<(DigramKey, usize, (Symbol, Symbol))> = sequence.par_windows(2)
                .enumerate()
                .map(|(i, window)| {
                    let sym1 = window[0];
                    let sym2 = window[1];
                    let canonical_key = Self::canonical_key(&(sym1, sym2), reverse_aware);
                    (canonical_key, i, (sym1, sym2))
                })
                .collect();
            
            // Then sequentially add them to the HashMap
            for (canonical_key, pos, digram) in digram_data {
                self.occurrences
                    .entry(canonical_key)
                    .or_default()
                    .push((pos, digram));
            }
        }
    }

    /// Adds a single pre-calculated entry to the table.
    /// Used for optimizations like the suffix array approach where we find
    /// only the most frequent digram initially.
    pub fn add_single_entry(&mut self, canonical_key: DigramKey, occurrences: Vec<(usize, (Symbol, Symbol))>) {
        if !occurrences.is_empty() {
            self.occurrences.insert(canonical_key, occurrences);
        }
    }

    /// Creates the key representation for a digram.
    fn get_digram_key(&self, first: Symbol, second: Symbol) -> DigramKey {
        custom_hash(&(first, second))
    }

    /// Determines the canonical form of a digram key.
    /// Currently only handles Terminal-Terminal digrams for reverse complementation.
    /// Non-Terminal digrams or mixed digrams return the original key.
    fn canonicalize_digram_key(&self, key: DigramKey) -> DigramKey {
        // With the new hash-based keys, we don't need to do anything here
        // since the canonical_key function already handles this logic
        key
    }

    /// Finds the most frequent canonical digram.
    /// Returns Option<(CanonicalDigramKey, count, Vec<(position, original_instance)>)> 
    pub fn find_most_frequent_digram(&self) -> Option<(DigramKey, usize, Vec<(usize, (Symbol, Symbol))>)> {
        // If the table is small, use the sequential approach
        if self.occurrences.len() < 100 {
            let mut max_key = None;
            let mut max_count = 0;
            let mut max_occurrences = Vec::new();

            // Iterate over entries and find the one with the most occurrences
            for (key, occurrences) in &self.occurrences {
                let count = occurrences.len();

                if count > max_count {
                    max_key = Some(*key);
                    max_count = count;
                    max_occurrences = occurrences.clone();
                }
            }

            max_key.map(|key| (key, max_count, max_occurrences))
        } else {
            // For larger tables, use parallel processing with rayon
            let result = self.occurrences.iter()
                .par_bridge() // Convert to parallel iterator
                .map(|(key, occurrences)| {
                    let count = occurrences.len();
                    (*key, count, occurrences.clone())
                })
                .max_by_key(|(_, count, _)| *count);
            
            // Return the result, if any
            result
        }
    }

    /// Removes occurrences associated with a specific digram instance (by position).
    /// Used after a digram is replaced by a rule.
    /// Needs careful handling if multiple digrams overlap the same position.
    pub fn remove_occurrence(&mut self, canonical_key: &DigramKey, position_to_remove: usize) -> bool {
        if let Some(occurrences) = self.occurrences.get_mut(canonical_key) {
            let initial_len = occurrences.len();
            
            // Remove occurrences starting at the specified position.
            occurrences.retain(|(pos, _)| *pos != position_to_remove);
            let removed = occurrences.len() < initial_len; // Check if anything was removed
            
            // If the vector became empty, remove the entry
            if occurrences.is_empty() {
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
        for (key, other_occurrences) in other.occurrences.iter() {
            let occurrences = self.occurrences.entry(*key).or_default();
            occurrences.extend_from_slice(other_occurrences);
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

    /// Builds a digram table from a sequence of symbols
    /// 
    /// # Arguments
    /// * `sequence` - The sequence of symbols to build the table from
    /// * `reverse_aware` - Whether to treat reverse complement digrams as the same
    pub fn build(sequence: &[Symbol], reverse_aware: bool) -> Self {
        const PARALLEL_THRESHOLD: usize = 10_000;
        
        // For small sequences, use sequential processing
        if sequence.len() < PARALLEL_THRESHOLD {
            return Self::build_sequential(sequence, reverse_aware);
        }
        
        // For larger sequences, use parallel processing with DashMap
        Self::build_parallel(sequence, reverse_aware)
    }
    
    /// Sequential implementation of digram table construction
    fn build_sequential(sequence: &[Symbol], reverse_aware: bool) -> Self {
        let mut occurrences = HashMap::new();
        
        // Scan the sequence for digrams
        for i in 0..sequence.len().saturating_sub(1) {
            let digram = (sequence[i].clone(), sequence[i+1].clone());
            let key = Self::canonical_key(&digram, reverse_aware);
            
            occurrences.entry(key)
                .or_insert_with(Vec::new)
                .push((i, digram));
        }
        
        DigramTable { occurrences }
    }
    
    /// Parallel implementation of digram table construction using DashMap and Rayon
    fn build_parallel(sequence: &[Symbol], reverse_aware: bool) -> Self {
        // Use DashMap for thread-safe concurrent access
        let digram_map = dashmap::DashMap::new();
        
        // Calculate the optimal chunk size based on sequence length
        let chunk_size = std::cmp::max(
            1000,
            std::cmp::min(sequence.len() / (num_cpus::get() * 2), 100_000)
        );
        
        // Process the sequence in parallel chunks
        sequence.par_chunks(chunk_size)
            .enumerate()
            .for_each(|(chunk_idx, chunk)| {
                let start_idx = chunk_idx * chunk_size;
                
                // Process each chunk, accounting for the overlap between chunks
                let end = if chunk.len() < 2 { 0 } else { chunk.len() - 1 };
                for i in 0..end {
                    let pos = start_idx + i;
                    let sym1 = chunk[i].clone();
                    let sym2 = chunk[i + 1].clone();
                    
                    // Get canonical key for this digram
                    let canonical_key = DigramTable::canonical_key(&(sym1.clone(), sym2.clone()), reverse_aware);
                    
                    // Thread-safe concurrent update using DashMap
                    digram_map.entry(canonical_key)
                        .or_insert_with(Vec::new)
                        .push((pos, (sym1, sym2)));
                }
            });
        
        // Convert DashMap to standard HashMap
        let occurrences = digram_map.into_iter()
            .map(|(k, v)| (k, v))
            .collect();
        
        DigramTable { occurrences }
    }

    /// Computes a canonical key for a digram tuple
    /// This is a static method to allow other modules to convert between digrams and keys
    pub fn canonical_key(digram: &(Symbol, Symbol), reverse_aware: bool) -> DigramKey {
        let (a, b) = digram;
        
        // Hash the digram using our custom hash function
        let mut key = custom_hash(&(a.id, b.id));
        
        // If reverse awareness is enabled, also hash the reverse complement 
        // and use the lexicographically smaller one as the canonical key
        if reverse_aware {
            // Only applies to terminal symbols that might have reverse complements
            if let (SymbolType::Terminal(a_val), SymbolType::Terminal(b_val)) = (&a.symbol_type, &b.symbol_type) {
                // Compute reverse complement digram hash
                let rev_comp_a = b.reverse_complement();
                let rev_comp_b = a.reverse_complement();
                let rev_comp_key = custom_hash(&(rev_comp_a.id, rev_comp_b.id));
                
                // Use the smaller hash as the canonical key
                if rev_comp_key < key {
                    key = rev_comp_key;
                }
            }
        }
        
        key
    }
}

/// Helper function to complement a DNA base
fn complement(val: &u8) -> u8 {
    match *val {
        0 => 3, // A -> T
        1 => 2, // C -> G
        2 => 1, // G -> C
        3 => 0, // T -> A
        _ => *val, // Non-DNA bases remain unchanged
    }
}

/// Provides a best-effort extraction of symbols from a sequence using a canonical key.
/// This is useful when we need to recover the original symbols after using a hash key.
pub fn extract_symbols_from_key(key: DigramKey, sequence: &[EncodedBase]) -> Option<(Symbol, Symbol)> {
    if sequence.len() < 2 {
        return None;
    }
    
    // Create new terminal symbols from the first digram in the sequence
    let base1 = sequence[0];
    let base2 = sequence[1];
    
    // Create default symbols for the observed digram
    let sym1 = Symbol::terminal(0, base1, Direction::Forward);
    let sym2 = Symbol::terminal(1, base2, Direction::Forward);
    
    Some((sym1, sym2))
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

        table.add_digram(0, s1.clone(), s2.clone(), true); // AC @ 0
        table.add_digram(1, s2.clone(), s3.clone(), true); // CG @ 1
        table.add_digram(5, s1.clone(), s2.clone(), true); // AC @ 5

        // Create the canonical keys using the new hash-based approach
        let key_ac = DigramTable::canonical_key(&(s1.clone(), s2.clone()), true);
        let key_cg = DigramTable::canonical_key(&(s2.clone(), s3.clone()), true);

        assert_eq!(table.occurrences.len(), 2);
        assert_eq!(table.occurrences.get(&key_ac).unwrap().len(), 2);
        assert_eq!(table.occurrences.get(&key_cg).unwrap().len(), 1);

        let most_frequent = table.find_most_frequent_digram().unwrap();
        assert_eq!(most_frequent.0, key_ac); // Canonical key AC
        assert_eq!(most_frequent.1, 2); // Count
        assert_eq!(most_frequent.2.len(), 2);
        assert_eq!(most_frequent.2[0].0, 0); // Position 0
        assert_eq!(most_frequent.2[0].1, (s1.clone(), s2.clone())); // Instance (s1, s2)
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

        table.add_digram(0, s_a.clone(), s_c.clone(), true); // AC @ 0.
        table.add_digram(5, s_t.clone(), s_g.clone(), true); // TG @ 5.

        // Create the canonical key using the new hash-based approach
        let key_ac = DigramTable::canonical_key(&(s_a.clone(), s_c.clone()), true);

        // Both AC and TG should map to the canonical AC key
        assert_eq!(table.occurrences.len(), 1);
        assert!(table.occurrences.contains_key(&key_ac));
        assert_eq!(table.occurrences.get(&key_ac).unwrap().len(), 2);

        // Check stored instances
        let occurrences = table.occurrences.get(&key_ac).unwrap();
        assert_eq!(occurrences[0].0, 0); // pos 0
        assert_eq!(occurrences[0].1, (s_a.clone(), s_c.clone())); // original AC
        assert_eq!(occurrences[1].0, 5); // pos 5
        assert_eq!(occurrences[1].1, (s_t.clone(), s_g.clone())); // original TG

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

        table.add_digram(0, s_a.clone(), s_c.clone(), false); // AC, reverse_aware=false
        table.add_digram(5, s_t.clone(), s_g.clone(), false); // TG, reverse_aware=false

        // Create the canonical keys using the new hash-based approach (without reverse awareness)
        let key_ac = DigramTable::canonical_key(&(s_a.clone(), s_c.clone()), false);
        let key_tg = DigramTable::canonical_key(&(s_t.clone(), s_g.clone()), false);

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

        table.add_digram(0, s1.clone(), s2.clone(), true); // AC @ 0
        table.add_digram(5, s1.clone(), s2.clone(), true); // AC @ 5
        table.add_digram(10, s1.clone(), s2.clone(), true); // AC @ 10
        
        // Create the canonical key using the new hash-based approach
        let key_ac = DigramTable::canonical_key(&(s1.clone(), s2.clone()), true);

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

        table.add_digram(0, s_a.clone(), s_nt1.clone(), true); // A-R100(+) @ 0
        table.add_digram(1, s_nt1.clone(), s_nt2.clone(), true); // R100(+)-R101(-) @ 1

        // Create the canonical keys using the new hash-based approach
        let key_a_nt1 = DigramTable::canonical_key(&(s_a.clone(), s_nt1.clone()), true);
        let key_nt1_nt2 = DigramTable::canonical_key(&(s_nt1.clone(), s_nt2.clone()), true);

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