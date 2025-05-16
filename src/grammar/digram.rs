use crate::grammar::symbol::{Symbol, Direction, SymbolType};
use std::collections::HashMap;
use std::hash::{Hash, Hasher};
use std::fmt;
use crate::encode::dna_2bit::EncodedBase;
use suffix_array::SuffixArray;
use crate::grammar::digram_table::DigramKeyTuple;
use std::collections::HashSet;

/// A digram is a pair of adjacent symbols in a sequence
#[derive(Debug, Clone, Eq)]
pub struct Digram {
    /// The first symbol
    pub first: Symbol,
    /// The second symbol
    pub second: Symbol,
}

impl Digram {
    /// Create a new digram from two symbols
    pub fn new(first: Symbol, second: Symbol) -> Self {
        Self { first, second }
    }

    /// Get the reverse complement of this digram
    /// 
    /// In DNA, the reverse complement of a sequence is obtained by 
    /// reversing the order and taking the complement of each base:
    /// A->T, C->G, G->C, T->A
    pub fn reverse_complement(&self) -> Self {
        Self::new(
            self.second.reverse_complement(),
            self.first.reverse_complement(),
        )
    }

    /// Check if this digram is identical to its reverse complement
    pub fn is_palindromic(&self) -> bool {
        let revcomp = self.reverse_complement();
        self == &revcomp
    }

    /// Get the canonical form of this digram (the lexicographically smaller of self and reverse_complement)
    pub fn canonical(&self, reverse_aware: bool) -> Self {
        if !reverse_aware {
            return self.clone();
        }
        
        let revcomp = self.reverse_complement();
        if self <= &revcomp {
            self.clone()
        } else {
            revcomp
        }
    }

    /// Create a string representation of the digram
    pub fn to_string(&self) -> String {
        format!("{}{}", self.first, self.second)
    }
}

impl PartialEq for Digram {
    fn eq(&self, other: &Self) -> bool {
        self.first == other.first && self.second == other.second
    }
}

impl Hash for Digram {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.first.hash(state);
        self.second.hash(state);
    }
}

impl fmt::Display for Digram {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}{}", self.first, self.second)
    }
}

impl PartialOrd for Digram {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Digram {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        match self.first.cmp(&other.first) {
            std::cmp::Ordering::Equal => self.second.cmp(&other.second),
            ordering => ordering,
        }
    }
}

/// Counts the occurrences of each digram in a sequence of symbols
pub fn count_digrams(sequence: &[Symbol], reverse_aware: bool) -> HashMap<DigramKeyTuple, usize> {
    let mut counts = HashMap::new();
    if sequence.len() < 2 {
        return counts;
    }
    for i in 0..sequence.len() - 1 {
        let s1 = &sequence[i];
        let s2 = &sequence[i + 1];
        let key = crate::grammar::digram_table::DigramTable::canonical_key((s1, s2), reverse_aware);
        *counts.entry(key).or_insert(0) += 1;
    }
    counts
}

/// Find digram with the highest count
pub fn find_most_frequent_digram(sequence: &[Symbol], reverse_aware: bool) -> Option<(Digram, usize, Vec<usize>)> {
    if sequence.len() < 2 {
        return None;
    }
    let mut counts: HashMap<DigramKeyTuple, (usize, Vec<usize>)> = HashMap::new();
    let mut original_digrams: HashMap<DigramKeyTuple, Digram> = HashMap::new();

    for i in 0..sequence.len() - 1 {
        let s1 = &sequence[i];
        let s2 = &sequence[i + 1];
        let digram = Digram::new(s1.clone(), s2.clone());
        let key = crate::grammar::digram_table::DigramTable::canonical_key((s1, s2), reverse_aware);
        let entry = counts.entry(key).or_insert((0, Vec::new()));
        entry.0 += 1;
        entry.1.push(i);
        original_digrams.entry(key).or_insert(digram);
    }

    counts.into_iter()
        .max_by_key(|(_, (count, _))| *count)
        .map(|(key, (count, positions))| {
            let original_digram = original_digrams.get(&key).expect("Original digram not found for key").clone();
            (original_digram, count, positions)
        })
}

/// Find all digram occurrences in a sequence
pub fn find_digram_positions(sequence: &[Symbol], target_digram: &Digram, reverse_aware: bool) -> Vec<usize> {
    let target_canonical_key = crate::grammar::digram_table::DigramTable::canonical_key((&target_digram.first, &target_digram.second), reverse_aware);
    let mut positions = Vec::new();
    for i in 0..sequence.len().saturating_sub(1) {
        let current_key = crate::grammar::digram_table::DigramTable::canonical_key((&sequence[i], &sequence[i+1]), reverse_aware);
        if current_key == target_canonical_key {
            positions.push(i);
        }
    }
    positions
}

/// Finds the most frequent terminal digram using a suffix array.
/// Returns Option<(CanonicalDigramKey, count, Vec<position>)>
pub fn find_most_frequent_terminal_digram_suffix_array(
    sequence: &[EncodedBase],
    min_count: usize,
    reverse_aware: bool,
) -> Option<(DigramKeyTuple, usize, Vec<usize>)> {
    if sequence.len() < 2 || min_count < 1 {
        return None;
    }

    // Convert EncodedBase sequence to u8 sequence for SuffixArray
    let seq_u8: Vec<u8> = sequence.iter().map(|eb| eb.0).collect();

    // Build the suffix array
    let sa = SuffixArray::new(&seq_u8);

    let mut digram_results: HashMap<(u8, u8), Vec<usize>> = HashMap::new();

    // Iterate through all 16 possible digrams (00, 01, ..., 11)
    for b1_val in 0..4u8 {
        for b2_val in 0..4u8 {
            let digram_bytes = [b1_val, b2_val];
            let positions = sa.search_all(&digram_bytes);
            if positions.len() >= min_count {
                // SuffixArray returns indices typically as u32 or usize, convert to usize
                let positions_usize: Vec<usize> = positions.iter().map(|&idx| idx as usize).collect();
                digram_results.insert((b1_val, b2_val), positions_usize);
            }
        }
    }

    if digram_results.is_empty() {
        return None;
    }

    // Find the most frequent digram considering canonicalization
    let mut best_digram_key: Option<DigramKeyTuple> = None;
    let mut max_freq = 0;
    let mut best_positions_vec: Vec<usize> = Vec::new();

    if reverse_aware {
        // Using a canonical map with DigramKey
        let mut canonical_freqs: HashMap<DigramKeyTuple, (usize, HashSet<usize>)> = HashMap::new();

        for (digram_tuple, positions) in digram_results {
            let count = positions.len();
            let base1 = EncodedBase(digram_tuple.0);
            let base2 = EncodedBase(digram_tuple.1);
            
            // Create temporary symbols for canonical key calculation
            let sym1 = Symbol::terminal(0, base1, Direction::Forward);
            let sym2 = Symbol::terminal(1, base2, Direction::Forward);
            
            // Use DigramTable's canonical_key function to get hash-based key
            let key = crate::grammar::digram_table::DigramTable::canonical_key((&sym1, &sym2), true);

            let entry = canonical_freqs.entry(key).or_insert((0, HashSet::new()));
            entry.0 += count; // Add frequency
            entry.1.extend(positions); // Merge positions
        }

        // Find max among canonical frequencies
        if let Some((key, (count, positions_set))) = canonical_freqs.into_iter().max_by_key(|(_, (c, _))| *c) {
            if count >= min_count { // Re-check min_count after merging
                best_digram_key = Some(key);
                max_freq = count;
                best_positions_vec = positions_set.into_iter().collect();
            }
        }
    } else {
        // Not reverse aware, just find the max frequency
        if let Some((digram_tuple, positions)) = digram_results.into_iter().max_by_key(|(_, pos_vec)| pos_vec.len()) {
            let count = positions.len();
            let base1 = EncodedBase(digram_tuple.0);
            let base2 = EncodedBase(digram_tuple.1);
            
            // Create temporary symbols for canonical key calculation
            let sym1 = Symbol::terminal(0, base1, Direction::Forward);
            let sym2 = Symbol::terminal(1, base2, Direction::Forward);
            
            // Use DigramTable's canonical_key function to get hash-based key
            let key = crate::grammar::digram_table::DigramTable::canonical_key((&sym1, &sym2), false);
            
            best_digram_key = Some(key);
            max_freq = count;
            best_positions_vec = positions;
        }
    }

    if let Some(key) = best_digram_key {
        best_positions_vec.sort_unstable(); // Sort positions for consistency
        Some((key, max_freq, best_positions_vec))
    } else {
        None
    }
}

pub fn count_digrams_with_positions(sequence: &[Symbol], reverse_aware: bool) 
    -> HashMap<DigramKeyTuple, (usize, HashSet<usize>)> {
    let mut canonical_freqs: HashMap<DigramKeyTuple, (usize, HashSet<usize>)> = HashMap::new();

    for i in 0..sequence.len().saturating_sub(1) {
        let s1 = &sequence[i];
        let s2 = &sequence[i+1];
        let key = crate::grammar::digram_table::DigramTable::canonical_key((s1, s2), reverse_aware);
        
        let entry = canonical_freqs.entry(key).or_insert((0, HashSet::new()));
        entry.0 += 1;
        entry.1.insert(i);
    }
    canonical_freqs
}

pub fn find_most_frequent_digram_with_min_count(
    sequence: &[Symbol],
    min_count: usize,
    reverse_aware: bool,
) -> Option<(Digram, usize, Vec<usize>)> {
    if sequence.len() < 2 {
        return None;
    }

    let mut counts: HashMap<DigramKeyTuple, (usize, Vec<usize>)> = HashMap::new();
    let mut original_digrams: HashMap<DigramKeyTuple, Digram> = HashMap::new();
    let mut best_digram_key: Option<DigramKeyTuple> = None;
    let mut max_count = 0;

    for i in 0..sequence.len() - 1 {
        let s1 = &sequence[i];
        let s2 = &sequence[i + 1];
        let current_digram = Digram::new(s1.clone(), s2.clone());
        let key = crate::grammar::digram_table::DigramTable::canonical_key((s1, s2), reverse_aware);

        let entry = counts.entry(key).or_insert((0, Vec::new()));
        entry.0 += 1;
        entry.1.push(i);
        original_digrams.entry(key).or_insert(current_digram);

        if entry.0 > max_count {
            max_count = entry.0;
            best_digram_key = Some(key);
        }
    }

    if max_count >= min_count {
        best_digram_key.map(|key| {
            let (count, positions) = counts.get(&key).unwrap().clone();
            let original_digram = original_digrams.get(&key).unwrap().clone();
            (original_digram, count, positions)
        })
    } else {
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::encode::dna_2bit::EncodedBase;
    use crate::grammar::symbol::{Direction, Symbol};
    use crate::grammar::digram_table::{DigramKeyTuple, DigramTable};

    fn encode_seq(seq_bytes: &[u8]) -> Vec<EncodedBase> {
        seq_bytes.iter().filter_map(|&b| EncodedBase::from_base(b)).collect()
    }

    fn create_terminal(id: u32, base: u8, strand: Direction) -> Symbol {
        Symbol::terminal(id as usize, EncodedBase(base), strand)
    }

    fn create_nonterminal(id: u32, rule_id: u32, strand: Direction) -> Symbol {
        Symbol::non_terminal(id as usize, rule_id as usize, strand)
    }

    #[test]
    fn test_digram_new() {
        let a = create_terminal(0, 0, Direction::Forward);
        let c = create_terminal(1, 1, Direction::Forward);
        let digram = Digram::new(a, c);
        assert_eq!(digram.first, a);
        assert_eq!(digram.second, c);
    }

    #[test]
    fn test_reverse_complement() {
        let a = create_terminal(0, 0, Direction::Forward); // A+
        let c = create_terminal(1, 1, Direction::Forward); // C+
        let digram = Digram::new(a, c);
        let rev = digram.reverse_complement();
        // Reverse complement of A+C+ is G-T-
        match rev.first.symbol_type {
            SymbolType::Terminal(base) => assert_eq!(base.0, 2), // Should be G (2)
            _ => panic!("Expected terminal symbol"),
        }
        match rev.second.symbol_type {
            SymbolType::Terminal(base) => assert_eq!(base.0, 3), // Should be T (3)
            _ => panic!("Expected terminal symbol"),
        }
    }

    #[test]
    fn test_is_palindromic() {
        // Regular digram A+C+ is not palindromic
        let sym1 = create_terminal(0, 0, Direction::Forward); // A+
        let sym2 = create_terminal(1, 1, Direction::Forward); // C+
        let digram1 = Digram::new(sym1, sym2);             // A+C+
        assert!(!digram1.is_palindromic());

        // Palindromic in DNA: A+T+ has revcomp A-T-. With strict equality (incl. strand), it's not palindromic.
        let sym_a = create_terminal(0, 0, Direction::Forward); // A+
        let sym_t = create_terminal(1, 3, Direction::Forward); // T+
        let digram2 = Digram::new(sym_a, sym_t);             // A+T+
        assert!(!digram2.is_palindromic()); // Updated assertion
    }

    #[test]
    fn test_canonical() {
        let a = create_terminal(0, 0, Direction::Forward);
        let c = create_terminal(1, 1, Direction::Forward);
        let digram = Digram::new(a, c);
        let rev = digram.reverse_complement();
        let can = digram.canonical(true);
        let rev_can = rev.canonical(true);
        assert_eq!(can, rev_can);
    }

    #[test]
    fn test_count_digrams() {
        let s_a_f = create_terminal(0, 0, Direction::Forward); 
        let s_c_f = create_terminal(1, 1, Direction::Forward); 
        let s_g_f = create_terminal(2, 2, Direction::Forward); 
        let s_t_f = create_terminal(3, 3, Direction::Forward); 

        let seq = [s_a_f.clone(), s_c_f.clone(), s_g_f.clone(), s_t_f.clone(), s_a_f.clone()];
        let counts = count_digrams(&seq, true);

        let d1 = Digram::new(s_a_f.clone(), s_c_f.clone());
        let d1_canonical = d1.canonical(true); 
        let key_d1_canonical = DigramTable::canonical_key((&d1_canonical.first, &d1_canonical.second), true);

        let d2 = Digram::new(s_g_f.clone(), s_t_f.clone());
        let d2_canonical = d2.canonical(true); 
        let key_d2_canonical = DigramTable::canonical_key((&d2_canonical.first, &d2_canonical.second), true);

        assert_ne!(key_d1_canonical, key_d2_canonical);

        assert_eq!(counts.get(&key_d1_canonical), Some(&1));
        assert_eq!(counts.get(&key_d2_canonical), Some(&1));
        
        let d_gt = Digram::new(s_g_f.clone(), s_t_f.clone());
        let key_d_gt_canonical = DigramTable::canonical_key((&d_gt.canonical(true).first, &d_gt.canonical(true).second), true);
        assert_eq!(counts.get(&key_d_gt_canonical), Some(&1));

        let d_ta = Digram::new(s_t_f.clone(), s_a_f.clone());
        let d_ta_canonical = d_ta.canonical(true);
        let key_d_ta_canonical = DigramTable::canonical_key((&d_ta_canonical.first, &d_ta_canonical.second), true);
        assert_eq!(counts.get(&key_d_ta_canonical), Some(&1));
    }

    #[test]
    fn test_find_most_frequent_digram() {
        let seq = [
            create_terminal(0, 0, Direction::Forward), // A
            create_terminal(1, 1, Direction::Forward), // C
            create_terminal(2, 0, Direction::Forward), // A
            create_terminal(3, 1, Direction::Forward), // C
            create_terminal(4, 2, Direction::Forward), // G
        ];
        let counts = count_digrams(&seq, true);
        let (digram, count, _positions) = find_most_frequent_digram(&seq, true).unwrap();
        assert_eq!(count, 2);
        assert!(digram == Digram::new(seq[0], seq[1]).canonical(true));
    }

    #[test]
    fn test_find_digram_positions() {
        let seq = [
            create_terminal(0, 0, Direction::Forward), // A
            create_terminal(1, 1, Direction::Forward), // C
            create_terminal(2, 0, Direction::Forward), // A
            create_terminal(3, 1, Direction::Forward), // C
            create_terminal(4, 2, Direction::Forward), // G
        ];
        let digram = Digram::new(seq[0], seq[1]).canonical(true);
        let positions = find_digram_positions(&seq, &digram, true);
        assert_eq!(positions, vec![0, 2]);
    }

    #[test]
    fn test_find_most_frequent_terminal_digram_suffix_array() {
        let seq = encode_seq(b"ACGTAC");
        let result = find_most_frequent_terminal_digram_suffix_array(&seq, 2, true);
        assert!(result.is_some());
        let (key, count, positions) = result.unwrap();
        assert_eq!(count, 2);
        assert_eq!(positions.len(), 2);
    }

    #[test]
    fn test_find_most_frequent_complex() {
        let seq = [
            create_terminal(0,0, Direction::Forward), // A+
            create_terminal(1,1, Direction::Forward), // C+
            create_terminal(2,2, Direction::Forward), // G+
            create_terminal(3,1, Direction::Forward), // C+
            create_terminal(4,0, Direction::Forward), // A+
            create_terminal(5,1, Direction::Forward), // C+
        ]; // A+C+G+C+A+C+
           // AC:0, GC:2, CA:3, AC:4 -> AC is most frequent (0,4)

        // find_most_frequent_digram returns Option<(Digram, count, Vec<positions>)> now
        let result = find_most_frequent_digram_with_min_count(&seq, 1, true);
        assert!(result.is_some());
        let (digram, count, positions) = result.unwrap(); // Unpack all three

        let ac_digram = Digram::new(create_terminal(0,0,Direction::Forward), create_terminal(1,1,Direction::Forward));
        let ac_canonical = ac_digram.canonical(true);

        assert_eq!(digram, ac_canonical);
        assert_eq!(count, 2);
        assert_eq!(positions, vec![0, 4]);
    }
} 