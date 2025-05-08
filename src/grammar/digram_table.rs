use crate::grammar::symbol::{Symbol, SymbolType, Direction, SymbolKey};
use crate::encode::dna_2bit::EncodedBase;
use dashmap;
use rayon::prelude::*;
use std::collections::HashMap;
use std::hash::{Hash, Hasher};
use crate::utils::hash::custom_hash;
use std::cmp::Ordering;

/// The key used for indexing digrams in the digram table
/// For terminals: ((base1, strand1), (base2, strand2))
/// For non-terminals: ((rule_id1, strand1), (rule_id2, strand2))
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Ord, PartialOrd)]
pub struct DigramKeyTuple(pub SymbolKey, pub SymbolKey);

impl DigramKeyTuple {
    pub fn new(s1: &Symbol, s2: &Symbol) -> Self {
        Self(SymbolKey::new(s1), SymbolKey::new(s2))
    }
}

pub type DigramKey = u64;

/// Definition of DigramSource
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum DigramSource {
    Original,
    Rule(usize), // Rule ID from which this digram was derived
}

/// Helper function to compare symbols based on type only (Terminal < NonTerminal, then base/rule_id)
/// Ignores strand and instance ID.
fn cmp_symbols_type_only(s1: &Symbol, s2: &Symbol) -> Ordering {
    match (&s1.symbol_type, &s2.symbol_type) {
        (SymbolType::Terminal(b1), SymbolType::Terminal(b2)) => b1.0.cmp(&b2.0),
        (SymbolType::NonTerminal(r1), SymbolType::NonTerminal(r2)) => r1.cmp(r2),
        (SymbolType::Terminal(_), SymbolType::NonTerminal(_)) => Ordering::Less, // Terminals < NonTerminals
        (SymbolType::NonTerminal(_), SymbolType::Terminal(_)) => Ordering::Greater,
    }
}

/// Stores digram occurrences and handles canonicalization.
#[derive(Debug, Default)]
pub struct DigramTable {
    // Stores occurrences: Canonical Digram Key -> List of (position, original_digram_instance)
    // Position usually refers to the index of the first symbol in the sequence.
    // Storing the original Symbol pair allows replacement later.
    occurrences: dashmap::DashMap<DigramKeyTuple, Vec<(usize, DigramSource)>>,
    pub(crate) rule_references: dashmap::DashMap<usize, Vec<DigramKeyTuple>>,
}

impl DigramTable {
    pub fn new() -> Self {
        Self {
            occurrences: dashmap::DashMap::new(),
            rule_references: dashmap::DashMap::new(),
        }
    }

    /// Adds a digram occurrence observed at a given position.
    /// Handles canonicalization based on terminal sequences.
    pub fn add_digram(&self, pos: usize, first: Symbol, second: Symbol, reverse_aware: bool) {
        let key = Self::canonical_key((&first, &second), reverse_aware);
        self.occurrences.entry(key).or_default().push((pos, DigramSource::Original));
    }

    /// Add digrams in parallel from a sequence
    pub fn add_digrams_from_sequence(&self, sequence: &[Symbol], reverse_aware: bool) {
        if sequence.len() < 2 {
            return;
        }

        if sequence.len() < 1000 { // Threshold for switching to parallel, adjust as needed
            // Sequential for small sequences
            for i in 0..sequence.len() - 1 {
                self.add_digram(i, sequence[i].clone(), sequence[i+1].clone(), reverse_aware);
            }
        } else {
            // Parallel for larger sequences (simplified from original complex rayon code)
            let digram_data: Vec<(DigramKeyTuple, usize, DigramSource)> = sequence // Changed DigramKey to DigramKeyTuple
                .par_windows(2)
                .enumerate()
                .map(|(i, window)| {
                    let sym1 = window[0].clone(); // Clone symbols for ownership
                    let sym2 = window[1].clone();
                    let key = Self::canonical_key((&sym1, &sym2), reverse_aware); // Pass tuple of refs
                    (key, i, DigramSource::Original)
                })
                .collect();

            for (key, pos, source) in digram_data {
                self.occurrences.entry(key).or_default().push((pos, source));
            }
        }
    }

    /// Adds a single pre-calculated entry to the table.
    /// Used for optimizations like the suffix array approach where we find
    /// only the most frequent digram initially.
    pub fn add_single_entry(&mut self, canonical_key: DigramKeyTuple, new_occurrences: Vec<(usize, DigramSource)>) {
        if !new_occurrences.is_empty() {
            self.occurrences.entry(canonical_key).or_default().extend(new_occurrences);
        }
    }

    /// Creates the key representation for a digram.
    fn get_digram_key(&self, first: Symbol, second: Symbol) -> DigramKey {
        let key_tuple = DigramKeyTuple::new(&first, &second);
        custom_hash(&key_tuple)
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
    pub fn find_most_frequent_digram(&self) -> Option<(DigramKeyTuple, usize, Vec<(usize, DigramSource)>)> {
        self.occurrences
            .iter()
            .map(|entry| { // entry is RefMulti
                let key = *entry.key();
                let value = entry.value();
                (key, value.len(), value.clone()) // Clone occurrences for return
            })
            .max_by_key(|&(_, count, _)| count)
    }

    /// Removes occurrences associated with a specific digram instance (by position).
    /// Used after a digram is replaced by a rule.
    /// Needs careful handling if multiple digrams overlap the same position.
    pub fn remove_occurrence(&mut self, canonical_key: &DigramKeyTuple, position_to_remove: usize, source_to_remove: DigramSource) -> bool {
        let mut removed = false;
        if let Some(mut occurrences_list) = self.occurrences.get_mut(canonical_key) {
            let initial_len = occurrences_list.len();
            occurrences_list.retain(|(pos, src)| !(*pos == position_to_remove && *src == source_to_remove));
            removed = occurrences_list.len() < initial_len;
            if occurrences_list.is_empty() {
                // If the vec is empty, remove the key from the map
                drop(occurrences_list); // Release the mutable borrow before removing
                self.occurrences.remove(canonical_key);
            }
        }
        removed
    }

    /// Merges another DigramTable into this one.
    /// Used for combining results from parallel processing.
    pub fn merge(&mut self, other: Self) {
        for entry in other.occurrences.iter() { // Correct way to iterate DashMap
            let key = *entry.key();
            let other_occurrences_list = entry.value();
            self.occurrences.entry(key).or_default().extend_from_slice(other_occurrences_list);
        }
        for entry in other.rule_references.iter() {
            let rule_id = *entry.key();
            let rule_keys = entry.value();
            self.rule_references.entry(rule_id).or_default().extend_from_slice(rule_keys);
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
        let table = Self::new();
        for i in 0..sequence.len().saturating_sub(1) {
            let s1 = &sequence[i]; // Get references
            let s2 = &sequence[i+1]; // Get references
            let key = Self::canonical_key((s1, s2), reverse_aware); // Pass tuple of refs
            table.occurrences.entry(key).or_default().push((i, DigramSource::Original));
        }
        table
    }
    
    /// Parallel implementation of digram table construction using DashMap and Rayon
    fn build_parallel(sequence: &[Symbol], reverse_aware: bool) -> Self {
        let occurrences = dashmap::DashMap::new();
        let rule_references = dashmap::DashMap::new(); // Initialize rule_references

        let num_chunks = rayon::current_num_threads();
        // Ensure chunk_size is at least 1, and at least 2 for digram processing if sequence is long enough
        let chunk_size = if sequence.len() < 2 { 
            sequence.len().max(1) 
        } else { 
            ((sequence.len() as f64 / num_chunks as f64).ceil() as usize).max(2) // Corrected casting and max call
        };

        if sequence.len() < 2 || chunk_size == 0 { 
            return Self { occurrences, rule_references };
        }
        
        let digram_data_chunks: Vec<Vec<(DigramKeyTuple, usize, DigramSource)>> = sequence
            .par_chunks(chunk_size)
            .enumerate()
            .map(|(chunk_idx, chunk)| {
                let mut local_digrams = Vec::new();
                if chunk.len() >= 2 { 
                    for i in 0..chunk.len() - 1 {
                        let s1 = &chunk[i]; // Get references
                        let s2 = &chunk[i+1]; // Get references
                        let key = Self::canonical_key((s1, s2), reverse_aware); // Pass tuple of refs
                        let original_pos = chunk_idx * chunk_size + i; 
                        local_digrams.push((key, original_pos, DigramSource::Original));
                    }
                }
                local_digrams
            })
            .collect();

        for chunk_data in digram_data_chunks {
            for (key, pos, source) in chunk_data {
                occurrences.entry(key).or_default().push((pos, source));
            }
        }
        Self { occurrences, rule_references }
    }

    /// Generates a canonical key for a digram, considering reverse complement if reverse_aware is true.
    /// This key is used for consistent lookup in the occurrence map.
    pub fn canonical_key(symbols: (&Symbol, &Symbol), reverse_aware: bool) -> DigramKeyTuple {
        let (s1, s2) = symbols;
        let key_tuple_fwd = DigramKeyTuple::new(s1, s2); // Key for S1 S2

        if reverse_aware {
            let can_revcomp = match (&s1.symbol_type, &s2.symbol_type) {
                (SymbolType::Terminal(_), SymbolType::Terminal(_)) => true,
                (SymbolType::NonTerminal(r1), SymbolType::NonTerminal(r2)) if r1 == r2 => true,
                _ => false,
            };

            if can_revcomp {
                 let s1_rc = s1.reverse_complement();
                 let s2_rc = s2.reverse_complement();
                 // Key tuple representing the sequence S2_rc S1_rc
                 let key_tuple_representing_rev_seq = DigramKeyTuple::new(&s2_rc, &s1_rc);

                 // The canonical key is DEFINED as the smaller of these two representations.
                 std::cmp::min(key_tuple_fwd, key_tuple_representing_rev_seq)
            } else {
                 key_tuple_fwd
            }
        } else {
            key_tuple_fwd
        }
    }

    /// Populates the digram table from a sequence of symbols
    pub fn populate(&mut self, sequence: &[Symbol], reverse_aware: bool) {
        self.add_digrams_from_sequence(sequence, reverse_aware);
    }

    /// Get the top N digrams by occurrence count
    pub fn get_top_digrams(&self, n: usize) -> Vec<(DigramKeyTuple, Vec<(usize, DigramSource)>)> {
        let mut all_digrams: Vec<(DigramKeyTuple, Vec<(usize, DigramSource)>)> = 
            self.occurrences.iter()
                .map(|entry| (*entry.key(), entry.value().clone())) // Correct iteration and cloning
                .collect();
            
        // Sort by count in descending order
        all_digrams.sort_by(|a, b| b.1.len().cmp(&a.1.len()));
        
        // Take the top N (or fewer if there aren't enough)
        all_digrams.truncate(n);
        
        all_digrams
    }

    pub fn to_hashmap(&self) -> HashMap<DigramKeyTuple, Vec<(usize, DigramSource)>> {
        self.occurrences.iter()
            .map(|entry| (*entry.key(), entry.value().clone())) // Correct way to map DashMap to HashMap
            .collect()
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
        let s1b = Symbol::terminal(10, EncodedBase(0), Direction::Forward); // A+ (different id)

        table.add_digram(0, s1.clone(), s2.clone(), true); // AC @ 0
        table.add_digram(1, s2.clone(), s3.clone(), true); // CG @ 1
        table.add_digram(5, s1b.clone(), s2.clone(), true); // AC @ 5 (should match first AC)

        let key_ac = DigramTable::canonical_key((&s1, &s2), true);
        let key_cg = DigramTable::canonical_key((&s2, &s3), true);

        assert_eq!(table.occurrences.len(), 2);
        assert!(table.occurrences.get(&key_ac).is_some());
        assert_eq!(table.occurrences.get(&key_ac).unwrap().len(), 2);
        assert!(table.occurrences.get(&key_cg).is_some());
        assert_eq!(table.occurrences.get(&key_cg).unwrap().len(), 1);
    }

    #[test]
    fn test_different_ids_same_bases() {
        let mut table = DigramTable::new();
        let s1 = Symbol::terminal(0, EncodedBase(0), Direction::Forward); // A+
        let s2 = Symbol::terminal(1, EncodedBase(1), Direction::Forward); // C+
        let s1b = Symbol::terminal(100, EncodedBase(0), Direction::Forward); // A+ (different id)
        let s2b = Symbol::terminal(101, EncodedBase(1), Direction::Forward); // C+ (different id)

        table.add_digram(0, s1.clone(), s2.clone(), false);
        table.add_digram(1, s1b.clone(), s2b.clone(), false);

        let key1 = DigramTable::canonical_key((&s1, &s2), false);
        let key2 = DigramTable::canonical_key((&s1b, &s2b), false);
        assert_eq!(key1, key2, "Digrams with same bases but different ids should have the same key");
        assert_eq!(table.occurrences.len(), 1);
        assert_eq!(table.occurrences.get(&key1).unwrap().len(), 2);
    }

    #[test]
    fn test_canonicalization() {
        let mut table = DigramTable::new();
        let s_a_fwd = Symbol::terminal(0, EncodedBase(0), Direction::Forward); // A+ id0
        let s_c_fwd = Symbol::terminal(1, EncodedBase(1), Direction::Forward); // C+ id1

        // True reverse complement of (A+ id0, C+ id1) is (G- id1, T- id0)
        // s_c_fwd.reverse_complement() would be G-(id=1)
        // s_a_fwd.reverse_complement() would be T-(id=0)
        let s_g_rev = s_c_fwd.reverse_complement(); 
        let s_t_rev = s_a_fwd.reverse_complement();

        table.add_digram(0, s_a_fwd.clone(), s_c_fwd.clone(), true);      // Add (A+ id0, C+ id1)
        table.add_digram(5, s_g_rev.clone(), s_t_rev.clone(), true);      // Add (G- id1, T- id0)

        assert_eq!(table.occurrences.len(), 1, 
            "Table should contain only one entry for a digram and its reverse complement when reverse_aware is true.");

        // The canonical key for (A+, C+) should be DKT(SK(A,F), SK(C,F))
        let expected_canonical_key = DigramKeyTuple::new(&s_a_fwd, &s_c_fwd);
        
        // We need to ensure this is indeed the smaller one if canonical_key was called on it
        let actual_expected_key_for_ac = DigramTable::canonical_key((&s_a_fwd, &s_c_fwd), true);
        assert_eq!(actual_expected_key_for_ac, expected_canonical_key, "Self-check failed: expected canonical key for A+C+ is not as assumed.");


        assert!(table.occurrences.contains_key(&expected_canonical_key),
            "Occurrences map should contain the canonical key {:?}", expected_canonical_key);
        
        let occurrences = table.occurrences.get(&expected_canonical_key).unwrap();
        assert_eq!(occurrences.len(), 2, "Should have two occurrences for the canonical digram.");
        assert!(occurrences.iter().any(|(pos, _src)| *pos == 0), "Position 0 missing.");
        assert!(occurrences.iter().any(|(pos, _src)| *pos == 5), "Position 5 missing.");
    }

    #[test]
    fn test_strand_and_nonterminal() {
        let mut table = DigramTable::new();
        let s_a_fwd = Symbol::terminal(0, EncodedBase(0), Direction::Forward); // A+
        let s_a_rev = Symbol::terminal(1, EncodedBase(0), Direction::Reverse); // A-
        let s_nt1 = Symbol::non_terminal(2, 100, Direction::Forward); // R100+
        let s_nt2 = Symbol::non_terminal(3, 100, Direction::Reverse); // R100-

        // Strand should matter
        table.add_digram(0, s_a_fwd.clone(), s_a_fwd.clone(), false);
        table.add_digram(1, s_a_fwd.clone(), s_a_rev.clone(), false);
        let key_fwd = DigramTable::canonical_key((&s_a_fwd, &s_a_fwd), false);
        let key_rev = DigramTable::canonical_key((&s_a_fwd, &s_a_rev), false);
        assert_ne!(key_fwd, key_rev);
        assert_eq!(table.occurrences.len(), 2);

        // Non-terminals with same rule id but different strand should be different
        table.add_digram(2, s_nt1.clone(), s_nt2.clone(), false);
        let key_nt = DigramTable::canonical_key((&s_nt1, &s_nt2), false);
        assert!(table.occurrences.contains_key(&key_nt));
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
        let key_ac = DigramTable::canonical_key((&s_a, &s_c), false);
        let key_tg = DigramTable::canonical_key((&s_t, &s_g), false);

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
        
        let key_ac = DigramTable::canonical_key((&s1, &s2), true);

        assert!(table.occurrences.get(&key_ac).is_some());
        assert_eq!(table.occurrences.get(&key_ac).unwrap().len(), 3);

        let removed = table.remove_occurrence(&key_ac, 5, DigramSource::Original);
        assert!(removed);
        assert!(table.occurrences.get(&key_ac).is_some());
        assert_eq!(table.occurrences.get(&key_ac).unwrap().len(), 2);
        assert!(table.occurrences.get(&key_ac).unwrap().iter().any(|(p, _)| *p == 0));
        assert!(!table.occurrences.get(&key_ac).unwrap().iter().any(|(p, _)| *p == 5));

        let removed_again = table.remove_occurrence(&key_ac, 100, DigramSource::Original);
        assert!(!removed_again);
        assert!(table.occurrences.get(&key_ac).is_some());
        assert_eq!(table.occurrences.get(&key_ac).unwrap().len(), 2);

        assert!(table.remove_occurrence(&key_ac, 0, DigramSource::Original));
        assert!(table.occurrences.get(&key_ac).is_some());
        assert_eq!(table.occurrences.get(&key_ac).unwrap().len(), 1);
        assert!(table.remove_occurrence(&key_ac, 10, DigramSource::Original));
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
        let key_a_nt1 = DigramTable::canonical_key((&s_a, &s_nt1), true);
        let key_nt1_nt2 = DigramTable::canonical_key((&s_nt1, &s_nt2), true);

        assert_eq!(table.occurrences.len(), 2);
        assert!(table.occurrences.contains_key(&key_a_nt1));
        assert!(table.occurrences.contains_key(&key_nt1_nt2));
        assert_eq!(table.occurrences.get(&key_a_nt1).unwrap().len(), 1);
        assert_eq!(table.occurrences.get(&key_nt1_nt2).unwrap().len(), 1);

        // Check find_most_frequent (tie-breaking might be arbitrary here)
        let most_frequent = table.find_most_frequent_digram().unwrap();
        assert_eq!(most_frequent.1, 1); // Count is 1 for both
    }

    #[test]
    fn test_new_empty() {
        let table = DigramTable::new();
        assert_eq!(table.occurrences.len(), 0);
        assert!(table.find_most_frequent_digram().is_none()); // Now it should work
    }
} 