use crate::grammar::symbol::{Symbol, SymbolType, SymbolKey};
use dashmap::DashMap;
use rayon::prelude::*;
use std::collections::{HashMap, BinaryHeap};
use std::hash::{Hash, Hasher};
use std::cmp::Reverse;
use fxhash::FxBuildHasher;

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
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, Ord, PartialOrd)]
pub enum DigramSource {
    Original,
    Rule(usize), // Rule ID from which this digram was derived
}

/// Stores digram occurrences and handles canonicalization.
#[derive(Debug)]
pub struct DigramTable {
    occurrences_shards: Vec<DashMap<DigramKeyTuple, Vec<(usize, DigramSource)>, FxBuildHasher>>,
    num_shards: usize,
    pub(crate) rule_references: DashMap<usize, Vec<DigramKeyTuple>, FxBuildHasher>,
}

impl DigramTable {
    pub fn new() -> Self {
        let num_shards = rayon::current_num_threads().max(1);
        let mut occurrences_shards = Vec::with_capacity(num_shards);
        for _ in 0..num_shards {
            occurrences_shards.push(DashMap::with_hasher(FxBuildHasher::default()));
        }
        Self {
            occurrences_shards,
            num_shards,
            rule_references: DashMap::with_hasher(FxBuildHasher::default()),
        }
    }

    fn get_shard_index(&self, key: &DigramKeyTuple) -> usize {
        let mut hasher = std::collections::hash_map::DefaultHasher::new();
        key.hash(&mut hasher);
        (hasher.finish() as usize) % self.num_shards
    }

    /// Adds a digram occurrence observed at a given position.
    /// Handles canonicalization based on terminal sequences.
    pub fn add_digram(&self, pos: usize, first: Symbol, second: Symbol, reverse_aware: bool) {
        let key = Self::canonical_key((&first, &second), reverse_aware);
        let shard_index = self.get_shard_index(&key);
        self.occurrences_shards[shard_index].entry(key).or_default().push((pos, DigramSource::Original));
    }

    /// Add digrams in parallel from a sequence
    pub fn add_digrams_from_sequence(&self, sequence: &[Symbol], reverse_aware: bool) {
        if sequence.len() < 2 {
            return;
        }

        if sequence.len() < 1000 {
            for i in 0..sequence.len() - 1 {
                self.add_digram(i, sequence[i].clone(), sequence[i+1].clone(), reverse_aware);
            }
        } else {
            sequence
                .par_windows(2)
                .enumerate()
                .for_each(|(i, window)| {
                    let sym1 = window[0].clone();
                    let sym2 = window[1].clone();
                    let key = Self::canonical_key((&sym1, &sym2), reverse_aware);
                    let shard_index = self.get_shard_index(&key);
                    self.occurrences_shards[shard_index]
                        .entry(key)
                        .or_default()
                        .push((i, DigramSource::Original));
                });
        }
    }

    /// Adds a single pre-calculated entry to the table.
    /// Used for optimizations like the suffix array approach where we find
    /// only the most frequent digram initially.
    pub fn add_single_entry(&mut self, canonical_key: DigramKeyTuple, new_occurrences: Vec<(usize, DigramSource)>) {
        if !new_occurrences.is_empty() {
            let shard_index = self.get_shard_index(&canonical_key);
            self.occurrences_shards[shard_index].entry(canonical_key).or_default().extend(new_occurrences);
        }
    }

    /// Finds the most frequent canonical digram.
    /// Returns Option<(CanonicalDigramKey, count, Vec<(position, original_instance)>)>
    pub fn find_most_frequent_digram(&self) -> Option<(DigramKeyTuple, usize, Vec<(usize, DigramSource)>)> {
        self.occurrences_shards
            .par_iter()
            .map(|shard| {
                shard
                    .iter()
                    .max_by_key(|entry| entry.value().len())
                    .map(|entry| (*entry.key(), entry.value().len(), entry.value().clone()))
            })
            .filter_map(|opt_entry| opt_entry)
            .reduce_with(|entry1, entry2| {
                if entry1.1 >= entry2.1 {
                    entry1
                } else {
                    entry2
                }
            })
    }

    /// Removes occurrences associated with a specific digram instance (by position).
    /// Used after a digram is replaced by a rule.
    /// Needs careful handling if multiple digrams overlap the same position.
    pub fn remove_occurrence(&mut self, canonical_key: &DigramKeyTuple, position_to_remove: usize, source_to_remove: DigramSource) -> bool {
        let shard_index = self.get_shard_index(canonical_key);
        let mut removed = false;
        if let Some(mut occurrences_list) = self.occurrences_shards[shard_index].get_mut(canonical_key) {
            let initial_len = occurrences_list.len();
            occurrences_list.retain(|(pos, src)| !(*pos == position_to_remove && *src == source_to_remove));
            removed = occurrences_list.len() < initial_len;
            if occurrences_list.is_empty() {
                drop(occurrences_list);
                self.occurrences_shards[shard_index].remove(canonical_key);
            }
        }
        removed
    }

    /// Merges another DigramTable into this one.
    /// Used for combining results from parallel processing.
    pub fn merge(&mut self, other: Self) {
        if self.num_shards == other.num_shards {
            for shard_idx in 0..self.num_shards {
                for entry in other.occurrences_shards[shard_idx].iter() {
                    let key = *entry.key();
                    let other_occurrences_list = entry.value();
                    self.occurrences_shards[shard_idx].entry(key).or_default().extend_from_slice(other_occurrences_list);
                }
            }
        } else {
            for other_shard in other.occurrences_shards {
                for entry in other_shard.into_iter() {
                    let (key, occurrences) = entry;
                    let self_shard_idx = self.get_shard_index(&key);
                    self.occurrences_shards[self_shard_idx].entry(key).or_default().extend(occurrences);
                }
            }
        }

        for entry in other.rule_references.iter() {
            let rule_id = *entry.key();
            let rule_keys = entry.value();
            self.rule_references.entry(rule_id).or_default().extend_from_slice(rule_keys);
        }
    }

    /// Get the number of distinct digrams in the table
    pub fn len(&self) -> usize {
        self.occurrences_shards.par_iter().map(|shard| shard.len()).sum()
    }

    /// Check if the table is empty
    pub fn is_empty(&self) -> bool {
        self.occurrences_shards.par_iter().all(|shard| shard.is_empty())
    }

    /// Clear all entries from the table
    pub fn clear(&mut self) {
        for shard in self.occurrences_shards.iter_mut() {
            shard.clear();
        }
        self.rule_references.clear();
    }

    /// Builds a digram table from a sequence of symbols
    /// 
    /// # Arguments
    /// * `sequence` - The sequence of symbols to build the table from
    /// * `reverse_aware` - Whether to treat reverse complement digrams as the same
    pub fn build(sequence: &[Symbol], reverse_aware: bool) -> Self {
        let table = Self::new();
        table.populate_internal(sequence, reverse_aware);
        table
    }
    
    /// Populates the digram table from a sequence of symbols
    pub fn populate(&mut self, sequence: &[Symbol], reverse_aware: bool) {
        self.add_digrams_from_sequence(sequence, reverse_aware);
    }

    /// Populates the digram table from a sequence of symbols (internal version)
    fn populate_internal(&self, sequence: &[Symbol], reverse_aware: bool) {
        self.add_digrams_from_sequence(sequence, reverse_aware);
    }

    /// Generates a canonical key for a digram, considering reverse complement if reverse_aware is true.
    /// This key is used for consistent lookup in the occurrence map.
    pub fn canonical_key(symbols: (&Symbol, &Symbol), reverse_aware: bool) -> DigramKeyTuple {
        let (s1, s2) = symbols;
        let key_tuple_fwd = DigramKeyTuple::new(s1, s2);

        if reverse_aware {
            let can_revcomp = match (&s1.symbol_type, &s2.symbol_type) {
                (SymbolType::Terminal(_), SymbolType::Terminal(_)) => true,
                (SymbolType::NonTerminal(r1), SymbolType::NonTerminal(r2)) if r1 == r2 => true,
                _ => false,
            };

            if can_revcomp {
                 let s1_rc = s1.reverse_complement();
                 let s2_rc = s2.reverse_complement();
                 let key_tuple_representing_rev_seq = DigramKeyTuple::new(&s2_rc, &s1_rc);
                 std::cmp::min(key_tuple_fwd, key_tuple_representing_rev_seq)
            } else {
                 key_tuple_fwd
            }
        } else {
            key_tuple_fwd
        }
    }

    /// Get the top N digrams by occurrence count
    pub fn get_top_digrams(&self, n: usize) -> Vec<(DigramKeyTuple, Vec<(usize, DigramSource)>)> {
        if n == 0 {
            return Vec::new();
        }

        // Use a min-heap for each shard, then merge the top N candidates from all shards.

        // Step 1: Get top N (at most) from each shard in parallel.
        let top_n_per_shard: Vec<Vec<(usize, DigramKeyTuple, Vec<(usize, DigramSource)>)>> = self.occurrences_shards
            .par_iter()
            .map(|shard| {
                let mut shard_min_heap = BinaryHeap::with_capacity(n + 1);
                for entry in shard.iter() {
                    let key = *entry.key();
                    let count = entry.value().len();
                    if count == 0 { continue; }

                    if shard_min_heap.len() < n {
                        shard_min_heap.push(Reverse((count, key, entry.value().clone())));
                    } else if count > shard_min_heap.peek().unwrap().0.0 {
                        shard_min_heap.pop();
                        shard_min_heap.push(Reverse((count, key, entry.value().clone())));
                    }
                }
                shard_min_heap.into_sorted_vec().into_iter().map(|rev_entry| rev_entry.0).collect()
            })
            .collect();

        // Step 2: Combine all candidates and find the global top N.
        let mut all_candidates = Vec::new();
        for shard_top_n in top_n_per_shard {
            for (count, key, occurrences_vec) in shard_top_n {
                all_candidates.push((key, count, occurrences_vec));
            }
        }

        // Sort all candidates by count descending
        all_candidates.sort_by(|a, b| b.1.cmp(&a.1));

        // Take top N and format output
        all_candidates.into_iter().take(n).map(|(key, _count, occurrences_vec)| (key, occurrences_vec)).collect()
    }

    pub fn to_hashmap(&self) -> HashMap<DigramKeyTuple, Vec<(usize, DigramSource)>> {
        let mut combined_map = HashMap::new();
        for shard in &self.occurrences_shards {
            for entry in shard.iter() {
                combined_map.insert(*entry.key(), entry.value().clone());
            }
        }
        combined_map
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::grammar::symbol::{Symbol, Direction, SymbolType};
    use crate::encode::dna_2bit::EncodedBase;

    #[test]
    fn test_add_and_find_digram_sharded() {
        let mut table = new_test_table(4);
        let s1 = Symbol::terminal(0, EncodedBase(0), Direction::Forward);
        let s2 = Symbol::terminal(1, EncodedBase(1), Direction::Forward);
        let s3 = Symbol::terminal(2, EncodedBase(2), Direction::Forward);

        table.add_digram(0, s1.clone(), s2.clone(), true);
        table.add_digram(1, s2.clone(), s3.clone(), true);
        table.add_digram(5, s1.clone(), s2.clone(), true);

        let key_ac = DigramTable::canonical_key((&s1, &s2), true);
        let key_cg = DigramTable::canonical_key((&s2, &s3), true);

        assert_eq!(table.len(), 2);
        assert!(table.occurrences_shards[0].contains_key(&key_ac));
        assert_eq!(table.occurrences_shards[0].get(&key_ac).unwrap().len(), 2);
        assert!(table.occurrences_shards[0].contains_key(&key_cg));
        assert_eq!(table.occurrences_shards[0].get(&key_cg).unwrap().len(), 1);
    }

    #[test]
    fn test_get_top_digrams_sharded() {
        let mut table = new_test_table(2);
        let s_a = Symbol::terminal(0, EncodedBase(0), Direction::Forward);
        let s_c = Symbol::terminal(1, EncodedBase(1), Direction::Forward);
        let s_g = Symbol::terminal(2, EncodedBase(2), Direction::Forward);
        let s_t = Symbol::terminal(3, EncodedBase(3), Direction::Forward);

        table.add_digram(0, s_a.clone(), s_c.clone(), true);
        table.add_digram(1, s_c.clone(), s_g.clone(), true);
        table.add_digram(5, s_g.clone(), s_t.clone(), true);

        let top_2 = table.get_top_digrams(2);
        assert_eq!(top_2.len(), 2);

        let key_ac = DigramTable::canonical_key((&s_a, &s_c), true);
        let key_cg = DigramTable::canonical_key((&s_c, &s_g), true);

        assert_eq!(top_2[0].0, key_ac);
        assert_eq!(top_2[0].1.len(), 2);
        assert_eq!(top_2[1].0, key_cg);
        assert_eq!(top_2[1].1.len(), 1);

        let top_1 = table.get_top_digrams(1);
        assert_eq!(top_1.len(), 1);
        assert_eq!(top_1[0].0, key_ac);

        let top_5 = table.get_top_digrams(5);
        assert_eq!(top_5.len(), 2);
    }

    #[test]
    fn test_remove_occurrence_sharded() {
        let mut table = new_test_table(4);
        let s1 = Symbol::terminal(0, EncodedBase(0), Direction::Forward);
        let s2 = Symbol::terminal(1, EncodedBase(1), Direction::Forward);

        table.add_digram(0, s1.clone(), s2.clone(), true);
        table.add_digram(5, s1.clone(), s2.clone(), true);
        table.add_digram(10, s1.clone(), s2.clone(), true);
        
        let key_ac = DigramTable::canonical_key((&s1, &s2), true);
        let shard_index_ac = table.get_shard_index(&key_ac);

        assert_eq!(table.occurrences_shards[shard_index_ac].get(&key_ac).unwrap().len(), 3);

        let removed = table.remove_occurrence(&key_ac, 5, DigramSource::Original);
        assert!(removed);
        assert_eq!(table.occurrences_shards[shard_index_ac].get(&key_ac).unwrap().len(), 2);
        assert!(table.occurrences_shards[shard_index_ac].get(&key_ac).unwrap().iter().any(|(p, _)| *p == 0));
        assert!(!table.occurrences_shards[shard_index_ac].get(&key_ac).unwrap().iter().any(|(p, _)| *p == 5));

        assert!(table.remove_occurrence(&key_ac, 0, DigramSource::Original));
        assert_eq!(table.occurrences_shards[shard_index_ac].get(&key_ac).unwrap().len(), 1);
        assert!(table.remove_occurrence(&key_ac, 10, DigramSource::Original));
        assert!(!table.occurrences_shards[shard_index_ac].contains_key(&key_ac));
    }

    fn new_test_table(num_shards: usize) -> DigramTable {
        let mut occurrences_shards = Vec::with_capacity(num_shards);
        for _ in 0..num_shards {
            occurrences_shards.push(DashMap::with_hasher(FxBuildHasher::default()));
        }
        DigramTable {
            occurrences_shards,
            num_shards,
            rule_references: DashMap::with_hasher(FxBuildHasher::default()),
        }
    }
} 