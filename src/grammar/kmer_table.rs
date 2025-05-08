//! K-mer table implementation for tracking frequent k-mers in sequences.
//!
//! This module generalizes the DigramTable to support arbitrary k-mer sizes (k ≥ 2),
//! allowing for detection of longer repeating patterns in genomic sequences.

use crate::grammar::symbol::{Symbol, SymbolType, Direction};
use std::collections::HashMap;
use anyhow::Result;
use rayon::prelude::*;
 // Import Ordering

// KmerKey is now the Vec<Symbol> itself representing the canonical k-mer
pub type KmerKey = Vec<Symbol>;

/// Stores k-mer occurrences and handles canonicalization.
#[derive(Debug)]
pub struct KmerTable {
    // Stores occurrences: Canonical K-mer Key (Vec<Symbol>) -> List of (position, original_kmer_instance)
    occurrences: HashMap<KmerKey, Vec<(usize, Vec<Symbol>)>>,
    // The k in k-mer, i.e., the length of patterns we're looking for
    k: usize,
}

impl KmerTable {
    /// Creates a new KmerTable with the specified k-mer size.
    pub fn new(k: usize) -> Self {
        assert!(k >= 2, "k-mer size must be at least 2");
        Self {
            occurrences: HashMap::new(),
            k,
        }
    }

    /// Gets the k value (pattern length) for this table
    pub fn k(&self) -> usize {
        self.k
    }

    /// Gets the canonical representation of a k-mer (Vec<Symbol>).
    /// The returned Vec<Symbol> has all symbol IDs normalized to 0 and strands to Forward
    /// after determining the lexicographically smaller of the k-mer and its reverse complement sequence.
    pub fn get_canonical_kmer(&self, symbols: &[Symbol], reverse_aware: bool) -> KmerKey {
        assert_eq!(symbols.len(), self.k, "Symbol slice must be of length k");

        // 1. Normalize original symbols to have ID 0 and Forward strand for comparison purposes.
        let kmer_fwd_norm: Vec<Symbol> = symbols.iter()
            .map(|s| Symbol::new(0, s.symbol_type, Direction::Forward))
            .collect();

        if !reverse_aware {
            return kmer_fwd_norm; // Return normalized original if not reverse_aware
        }

        // Check if all symbols are terminal (DNA bases) or all are the same non-terminal.
        // This is a simplified check; more complex rules might be needed for mixed/multi-NT kmers.
        let can_revcomp = if symbols.is_empty() { false } else {
            let first_type = &symbols[0].symbol_type;
            symbols.iter().all(|s| {
                match (&s.symbol_type, first_type) {
                    (SymbolType::Terminal(_), SymbolType::Terminal(_)) => true,
                    (SymbolType::NonTerminal(r1), SymbolType::NonTerminal(r2)) if r1 == r2 => true,
                    _ => false,
                }
            })
        };

        if can_revcomp {
            // 2. Create the reverse complement sequence of the original k-mer.
            let rev_comp_seq: Vec<Symbol> = symbols.iter().rev()
                .map(|s| s.reverse_complement()) // This preserves original IDs and flips strands
                .collect();

            // 3. Normalize this reverse complement sequence to ID 0 and Forward strand.
            let kmer_rc_seq_fwd_norm: Vec<Symbol> = rev_comp_seq.iter()
                .map(|s| Symbol::new(0, s.symbol_type, Direction::Forward))
                .collect();
            
            // 4. Compare the two normalized, forward-stranded Vec<Symbol>s and return the smaller.
            // Comparison uses Symbol's Ord, which will compare by type (since strand is now always Fwd).
            if kmer_fwd_norm <= kmer_rc_seq_fwd_norm {
                kmer_fwd_norm
            } else {
                kmer_rc_seq_fwd_norm
            }
        } else {
            // If not all terminals or not same non-terminal, just return the forward-normalized original.
            kmer_fwd_norm
        }
    }


    /// Adds a k-mer occurrence observed at a given position.
    pub fn add_kmer(&mut self, pos: usize, original_symbols: Vec<Symbol>, reverse_aware: bool) {
        assert_eq!(original_symbols.len(), self.k, "Symbol vector must be of length k");

        let canonical_key = self.get_canonical_kmer(&original_symbols, reverse_aware);

        self.occurrences
            .entry(canonical_key) // The key is the canonical Vec<Symbol>
            .or_default()
            .push((pos, original_symbols)); // Store original symbols with the position
    }

    /// Add k-mers from a sequence
    pub fn add_kmers_from_sequence(&mut self, sequence: &[Symbol], reverse_aware: bool) {
        if sequence.len() < self.k {
            return;
        }

        const PARALLEL_THRESHOLD: usize = 10_000;

        if sequence.len() < PARALLEL_THRESHOLD {
            // For small sequences, process sequentially
            for i in 0..=sequence.len() - self.k {
                let kmer_slice = &sequence[i..(i + self.k)];
                let canonical_key = self.get_canonical_kmer(kmer_slice, reverse_aware);
                self.occurrences
                    .entry(canonical_key)
                    .or_default()
                    .push((i, kmer_slice.to_vec())); // Store original slice
            }
        } else {
            // For large sequences, use parallel processing
            let kmer_data: Vec<(KmerKey, usize, Vec<Symbol>)> = sequence
                .par_windows(self.k)
                .enumerate()
                .map(|(i, window)| {
                    let canonical_key = self.get_canonical_kmer(window, reverse_aware);
                    (canonical_key, i, window.to_vec()) // Store original window
                })
                .collect();

            // Add collected k-mers to the table
            // Grouping by key first might be more efficient for locking HashMap less often
            let mut grouped_data: HashMap<KmerKey, Vec<(usize, Vec<Symbol>)>> = HashMap::new();
            for (canonical_key, pos, kmer) in kmer_data {
                grouped_data.entry(canonical_key).or_default().push((pos, kmer));
            }

            for (key, data_vec) in grouped_data {
                 self.occurrences.entry(key).or_default().extend(data_vec);
            }
        }
    }

    /// Builds a k-mer table from a sequence of symbols
    pub fn build(sequence: &[Symbol], k: usize, reverse_aware: bool) -> Self {
        let mut table = Self::new(k);
        table.add_kmers_from_sequence(sequence, reverse_aware);
        table
    }

    /// Finds the most frequent k-mer in the table
    /// Returns Option<(Canonical KmerKey (Vec<Symbol>), count, Vec<(position, original_kmer)>)>
    pub fn find_most_frequent_kmer(&self, min_count: usize) -> Option<(KmerKey, usize, Vec<(usize, Vec<Symbol>)>)> {
        self.occurrences
            .iter()
            .filter(|(_, occurrences)| occurrences.len() >= min_count)
            .max_by_key(|(_, occurrences)| occurrences.len())
            // Clone the key (Vec<Symbol>) and the occurrences Vec
            .map(|(key, occurrences)| (key.clone(), occurrences.len(), occurrences.clone()))
    }


    /// Alias for find_most_frequent_kmer to maintain compatibility
    pub fn find_most_frequent(&self, min_count: usize) -> Option<(KmerKey, usize, Vec<(usize, Vec<Symbol>)>)> {
        self.find_most_frequent_kmer(min_count)
    }

    /// Remove a specific k-mer occurrence identified by its canonical key and starting position
    pub fn remove_occurrence(&mut self, canonical_key: &KmerKey, position_to_remove: usize) -> bool {
        if let Some(occurrences) = self.occurrences.get_mut(canonical_key) {
            let initial_len = occurrences.len();

            // Remove occurrences starting at the specified position
            occurrences.retain(|(pos, _)| *pos != position_to_remove);
            let removed = occurrences.len() < initial_len;

            // If the vector became empty, remove the entry
            if occurrences.is_empty() {
                self.occurrences.remove(canonical_key);
            }

            removed
        } else {
            false
        }
    }


    /// Merges another KmerTable into this one
    pub fn merge_with(&mut self, other: &KmerTable) {
        assert_eq!(self.k, other.k, "Cannot merge KmerTables with different k values");

        for (key, other_occurrences) in other.occurrences.iter() {
            // Clone the key for insertion if it doesn't exist
            let occurrences = self.occurrences.entry(key.clone()).or_default();
            // Clone the occurrences from the other table to extend
            occurrences.extend(other_occurrences.iter().cloned());
        }
    }


    /// Get the number of distinct k-mers in the table
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

    /// Populates the table with k-mers from a sequence.
    pub fn populate(&mut self, sequence: &[Symbol], reverse_aware: bool) -> Result<()> {
        self.add_kmers_from_sequence(sequence, reverse_aware);
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::grammar::symbol::{Symbol, Direction, SymbolType};
    use crate::encode::dna_2bit::EncodedBase;

    // Helper to create terminal symbols
    fn term(id: usize, base_val: u8, dir: Direction) -> Symbol {
        Symbol::terminal(id, EncodedBase(base_val), dir)
    }

    #[test]
    fn test_kmer_table_basic() {
        let s1 = term(0, 0, Direction::Forward); // A+
        let s2 = term(1, 1, Direction::Forward); // C+
        let s3 = term(2, 2, Direction::Forward); // G+
        let s4 = term(3, 3, Direction::Forward); // T+
        let s5 = term(4, 0, Direction::Forward); // A+
        let sequence = vec![s1.clone(), s2.clone(), s3.clone(), s4.clone(), s5.clone()]; // ACGTA

        let mut table = KmerTable::new(3);
        table.populate(&sequence, true).unwrap(); // Populate with k=3, reverse aware

        // Kmers: ACG (pos 0), CGT (pos 1), GTA (pos 2)
        // Reverse Complements: CGT (rc of ACG), ACG (rc of CGT), TAC (rc of GTA)

        // Calculate the canonical key for the ACG/CGT pair
        let kmer_acg = vec![s1.clone(), s2.clone(), s3.clone()]; // A+C+G+
        let canonical_key_acg_cgt = table.get_canonical_kmer(&kmer_acg, true);

        // Calculate the canonical key for GTA
        let kmer_gta = vec![s3.clone(), s4.clone(), s5.clone()]; // G+T+A+
        let canonical_key_gta = table.get_canonical_kmer(&kmer_gta, true);

        // The canonical representation of ACG/CGT might be ACG (A+C+G+) or its reverse complement (T-G-C-)
        // Construct the reverse complement of kmer_acg
        let kmer_acg_rc_seq = vec![
            s3.reverse_complement(), // G+ -> C-
            s2.reverse_complement(), // C+ -> G-
            s1.reverse_complement()  // A+ -> T-
        ];
        // The canonical key should be one of these two forms.
        assert!(
            canonical_key_acg_cgt == kmer_acg || canonical_key_acg_cgt == kmer_acg_rc_seq,
            "Canonical key for ACG was {:?}, expected {:?} or {:?}",
            canonical_key_acg_cgt, kmer_acg, kmer_acg_rc_seq
        );

        // Assertions about the table state
        assert_eq!(table.len(), 2, "Expected 2 distinct canonical k-mers (ACG/CGT group and GTA group)");

        // Check occurrences for the ACG/CGT group using its canonical key
        assert!(table.occurrences.contains_key(&canonical_key_acg_cgt));
        let occurrences_acg_cgt = table.occurrences.get(&canonical_key_acg_cgt).unwrap();
        assert_eq!(occurrences_acg_cgt.len(), 2, "Expected 2 occurrences for the ACG/CGT canonical group");
        assert!(occurrences_acg_cgt.iter().any(|(pos, _)| *pos == 0), "Missing occurrence at pos 0 for ACG/CGT group");
        assert!(occurrences_acg_cgt.iter().any(|(pos, _)| *pos == 1), "Missing occurrence at pos 1 for ACG/CGT group");

        // Check occurrences for the GTA group using its canonical key
        assert!(table.occurrences.contains_key(&canonical_key_gta));
        let occurrences_gta = table.occurrences.get(&canonical_key_gta).unwrap();
        assert_eq!(occurrences_gta.len(), 1, "Expected 1 occurrence for GTA");
        assert_eq!(occurrences_gta[0].0, 2, "Incorrect position for GTA occurrence");

        // Test find_most_frequent
        let (most_freq_key, count, _positions) = table.find_most_frequent(1).unwrap(); // Min count 1
        assert_eq!(most_freq_key, canonical_key_acg_cgt, "Most frequent should be ACG/CGT group key");
        assert_eq!(count, 2, "Most frequent count should be 2");
    }


    #[test]
    fn test_canonicalization_kmer() {
        let table = KmerTable::new(3);
        let s1 = Symbol::terminal(0, EncodedBase(0), Direction::Forward); // A+ id0
        let s2 = Symbol::terminal(1, EncodedBase(1), Direction::Forward); // C+ id1
        let s3 = Symbol::terminal(2, EncodedBase(2), Direction::Forward); // G+ id2

        // kmer_ctg_rev is G.rc, C.rc, A.rc (from s3,s2,s1)
        let kmer_ctg_rev = vec![
            s3.reverse_complement(), // C (id2) Rev
            s2.reverse_complement(), // G (id1) Rev
            s1.reverse_complement(), // T (id0) Rev
        ];
        // println!("DEBUG: kmer_ctg_rev = {:?}", kmer_ctg_rev);

        let canonical_key_ctg_rev = table.get_canonical_kmer(&kmer_ctg_rev, true);
        // println!("DEBUG: canonical_key_ctg_rev = {:?}", canonical_key_ctg_rev);

        // The expected canonical form is the normalized forward k-mer (ACG) because A < C
        let expected_canonical_form: Vec<Symbol> = vec![
            Symbol::new(0, s1.symbol_type, Direction::Forward), // A0F
            Symbol::new(0, s2.symbol_type, Direction::Forward), // C0F
            Symbol::new(0, s3.symbol_type, Direction::Forward)  // G0F
        ];
        // println!("DEBUG: expected_canonical_form = {:?}", expected_canonical_form);

        assert_eq!(canonical_key_ctg_rev, expected_canonical_form, "Canonical key for CTG_rev should be the normalized ACG form");

        // Test ACG part
        let kmer_acg = vec![s1.clone(), s2.clone(), s3.clone()];
        let canonical_key_acg = table.get_canonical_kmer(&kmer_acg, true);
        assert_eq!(canonical_key_acg, expected_canonical_form, "Canonical key for ACG should be the normalized ACG form");
        assert_eq!(canonical_key_acg, canonical_key_ctg_rev, "Canonical keys for kmer and its revcomp sequence should be identical");
    }

}