use crate::grammar::symbol::Symbol;
use std::collections::HashMap;
use std::hash::{Hash, Hasher};
use std::fmt;
use crate::grammar::symbol::{Direction, SymbolType};
use crate::encode::dna_2bit::EncodedBase;
use suffix_array::SuffixArray;
use crate::grammar::digram_table::DigramKey;
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
pub fn count_digrams(sequence: &[Symbol], reverse_aware: bool) -> HashMap<Digram, usize> {
    let mut counts = HashMap::new();
    
    for i in 0..sequence.len().saturating_sub(1) {
        let digram = Digram::new(sequence[i], sequence[i + 1]);
        let canonical = digram.canonical(reverse_aware);
        
        *counts.entry(canonical).or_insert(0) += 1;
    }
    
    counts
}

/// Find digram with the highest count
pub fn find_most_frequent_digram(
    counts: &HashMap<Digram, usize>,
    min_count: usize,
) -> Option<(Digram, usize)> {
    counts
        .iter()
        .filter(|(_, &count)| count >= min_count)
        .max_by_key(|(_, &count)| count)
        .map(|(digram, &count)| (digram.clone(), count))
}

/// Find all digram occurrences in a sequence
pub fn find_digram_positions(sequence: &[Symbol], target: &Digram, reverse_aware: bool) -> Vec<usize> {
    let mut positions = Vec::new();
    
    for i in 0..sequence.len().saturating_sub(1) {
        let digram = Digram::new(sequence[i], sequence[i + 1]);
        
        if reverse_aware {
            let canonical = digram.canonical(true);
            let target_canonical = target.canonical(true);
            
            if canonical == target_canonical {
                positions.push(i);
            }
        } else {
            if digram == *target {
                positions.push(i);
            }
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
) -> Option<(DigramKey, usize, Vec<usize>)> {
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
    let mut best_digram_key: Option<DigramKey> = None;
    let mut max_freq = 0;
    let mut best_positions_vec: Vec<usize> = Vec::new();

    if reverse_aware {
        let mut canonical_freqs: HashMap<DigramKey, (usize, HashSet<usize>)> = HashMap::new();

        for (digram_tuple, positions) in digram_results {
            let count = positions.len();
            let base1 = EncodedBase(digram_tuple.0);
            let base2 = EncodedBase(digram_tuple.1);
            let key = ((SymbolType::Terminal(base1), Direction::Forward), (SymbolType::Terminal(base2), Direction::Forward));

            // Determine canonical key
            let rev_base1 = base2.revcomp();
            let rev_base2 = base1.revcomp();
            let rev_key_type_only = ((SymbolType::Terminal(rev_base1), Direction::Forward), (SymbolType::Terminal(rev_base2), Direction::Forward));
            let canonical_key = if key.0 <= rev_key_type_only.0 { key } else { rev_key_type_only };

            let entry = canonical_freqs.entry(canonical_key).or_insert((0, HashSet::new()));
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
            let key = ((SymbolType::Terminal(base1), Direction::Forward), (SymbolType::Terminal(base2), Direction::Forward));
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::encode::dna_2bit::EncodedBase;
    use crate::grammar::symbol::{Symbol, Direction, SymbolType};
    use std::collections::HashSet;

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
        let sym1 = create_terminal(0, 0, Direction::Forward);
        let sym2 = create_terminal(1, 1, Direction::Forward);
        let digram = Digram::new(sym1, sym2);
        
        assert_eq!(digram.first.id, 0);
        assert_eq!(digram.second.id, 1);
    }
    
    #[test]
    fn test_reverse_complement() {
        let sym1 = create_terminal(0, 0, Direction::Forward); // A+
        let sym2 = create_terminal(1, 1, Direction::Forward); // C+
        let digram = Digram::new(sym1, sym2);              // A+C+
        
        let rev_comp = digram.reverse_complement();        // G-T-
        
        // The second symbol (C+) becomes the first in reverse complement, and is complemented/reversed
        assert!(matches!(rev_comp.first.symbol_type, SymbolType::Terminal(EncodedBase(2)))); // G
        assert_eq!(rev_comp.first.strand, Direction::Reverse);
        
        // The first symbol (A+) becomes the second in reverse complement, and is complemented/reversed
        assert!(matches!(rev_comp.second.symbol_type, SymbolType::Terminal(EncodedBase(3)))); // T
        assert_eq!(rev_comp.second.strand, Direction::Reverse);
    }
    
    #[test]
    fn test_is_palindromic() {
        // Regular digram A+C+ is not palindromic
        let sym1 = create_terminal(0, 0, Direction::Forward); // A+
        let sym2 = create_terminal(1, 1, Direction::Forward); // C+
        let digram1 = Digram::new(sym1, sym2);             // A+C+
        assert!(!digram1.is_palindromic());
        
        // Palindromic in DNA: A+T+ is reverse complement of itself
        let sym1 = create_terminal(0, 0, Direction::Forward); // A+
        let sym2 = create_terminal(1, 3, Direction::Forward); // T+
        let digram2 = Digram::new(sym1, sym2);             // A+T+
        // This is A+T+ vs T-A- which are different Symbol objects
        assert!(!digram2.is_palindromic()); 
        
        // For non-terminal symbols we can create true palindromes
        let sym1 = create_nonterminal(1, 1, Direction::Forward);  // R1+
        let sym2 = create_nonterminal(2, 1, Direction::Reverse);  // R1-
        let digram3 = Digram::new(sym1, sym2);                // R1+R1-
        
        // The reverse complement of R1+R1- is:
        // - second (R1-) -> first and flip: R1+
        // - first (R1+) -> second and flip: R1-
        // So we get R1+R1- again
        assert!(digram3.is_palindromic());
    }
    
    #[test]
    fn test_canonical() {
        let sym1 = create_terminal(0, 0, Direction::Forward); // A+
        let sym2 = create_terminal(1, 1, Direction::Forward); // C+
        let digram = Digram::new(sym1, sym2);              // A+C+
        
        // When reverse_aware is false, canonical returns the digram unchanged
        let canonical1 = digram.canonical(false);
        assert!(matches!(canonical1.first.symbol_type, SymbolType::Terminal(EncodedBase(0))));
        assert!(matches!(canonical1.second.symbol_type, SymbolType::Terminal(EncodedBase(1))));
        
        // When reverse_aware is true, canonical returns the lexicographically smaller
        // of the digram and its reverse complement
        // A+C+ vs G-T-
        let canonical2 = digram.canonical(true);
        // A+ < G-, so should return original
        assert!(matches!(canonical2.first.symbol_type, SymbolType::Terminal(EncodedBase(0))));
        assert!(matches!(canonical2.second.symbol_type, SymbolType::Terminal(EncodedBase(1))));
        
        // Test where reverse complement is smaller
        let sym1 = create_terminal(3, 3, Direction::Forward); // T+
        let sym2 = create_terminal(4, 2, Direction::Forward); // G+
        let digram = Digram::new(sym1, sym2);              // T+G+
        
        // Reverse complement is C-A- which is lexicographically smaller
        let canonical = digram.canonical(true);
        assert!(matches!(canonical.first.symbol_type, SymbolType::Terminal(EncodedBase(1)))); // C
        assert_eq!(canonical.first.strand, Direction::Reverse);
        assert!(matches!(canonical.second.symbol_type, SymbolType::Terminal(EncodedBase(0)))); // A
        assert_eq!(canonical.second.strand, Direction::Reverse);
    }
    
    #[test]
    fn test_count_digrams() {
        let sequence = vec![
            create_terminal(0, 0, Direction::Forward), // A+
            create_terminal(1, 1, Direction::Forward), // C+
            create_terminal(2, 2, Direction::Forward), // G+
            create_terminal(3, 3, Direction::Forward), // T+
            create_terminal(4, 0, Direction::Forward), // A+
        ];
        
        // Without reverse awareness
        let counts1 = count_digrams(&sequence, false);
        assert_eq!(counts1.len(), 4);
        
        // With reverse awareness, some digrams may be canonicalized to their reverse complement
        let counts2 = count_digrams(&sequence, true);
        
        // For DNA these canonical pairs would be:
        // A+C+ <-> G-T- (canonical is A+C+)
        // C+G+ <-> C-G- (canonical is C+G+)
        // G+T+ <-> A-C- (canonical is A-C-)
        // T+A+ <-> T-A- (canonical is T+A+)
        
        assert_eq!(counts2.len(), 4); // Still 4 unique digrams
    }
    
    #[test]
    fn test_find_most_frequent_digram() {
        let mut counts = HashMap::new();
        
        let digram1 = Digram::new(
            create_terminal(0, 0, Direction::Forward),
            create_terminal(1, 1, Direction::Forward),
        ); // A+C+
        
        let digram2 = Digram::new(
            create_terminal(2, 1, Direction::Forward),
            create_terminal(3, 2, Direction::Forward),
        ); // C+G+
        
        let digram3 = Digram::new(
            create_terminal(4, 2, Direction::Forward),
            create_terminal(5, 3, Direction::Forward),
        ); // G+T+
        
        counts.insert(digram1.clone(), 5);
        counts.insert(digram2.clone(), 10);
        counts.insert(digram3.clone(), 2);
        
        // Find most frequent with min_count = 1
        let result1 = find_most_frequent_digram(&counts, 1);
        assert!(result1.is_some());
        let (digram, count) = result1.unwrap();
        assert_eq!(digram, digram2);
        assert_eq!(count, 10);
        
        // With min_count = 11, no digram is found
        let result2 = find_most_frequent_digram(&counts, 11);
        assert!(result2.is_none());
    }
    
    #[test]
    fn test_find_digram_positions() {
        let sequence = vec![
            create_terminal(0, 0, Direction::Forward), // A+
            create_terminal(1, 1, Direction::Forward), // C+
            create_terminal(2, 2, Direction::Forward), // G+
            create_terminal(3, 0, Direction::Forward), // A+
            create_terminal(4, 1, Direction::Forward), // C+
        ];
        
        let target = Digram::new(
            create_terminal(0, 0, Direction::Forward),
            create_terminal(1, 1, Direction::Forward),
        ); // A+C+
        
        // Without reverse awareness
        let positions1 = find_digram_positions(&sequence, &target, false);
        assert_eq!(positions1, vec![0, 3]);
        
        // Reverse complement of A+C+ is G-T-
        let target2 = Digram::new(
            create_terminal(2, 2, Direction::Forward),
            create_terminal(3, 3, Direction::Forward),
        ); // G+T+
        
        // With reverse awareness, G+T+ and C-A- are treated as the same
        let positions2 = find_digram_positions(&sequence, &target2, true);
        assert!(positions2.is_empty()); // None in this sequence
    }

    #[test]
    fn test_find_most_frequent_terminal_digram_suffix_array() {
        // Sequence: AC AC GT AC GT (AC appears 3 times, GT appears 2 times)
        let seq = encode_seq(b"ACACGTACGT");

        // Test 1: Not reverse aware, min_count = 2
        let result1 = find_most_frequent_terminal_digram_suffix_array(&seq, 2, false);
        assert!(result1.is_some());
        let (key1, count1, pos1) = result1.unwrap();
        // Expect AC (0,1) on Forward strand
        assert_eq!(key1, ((SymbolType::Terminal(EncodedBase(0)), Direction::Forward), (SymbolType::Terminal(EncodedBase(1)), Direction::Forward)));
        assert_eq!(count1, 3); // AC occurs at 0, 2, 6
        assert_eq!(pos1, vec![0, 2, 6]);

        // Test 2: Reverse aware, min_count = 2
        // AC (+) -> canonical AC(+)
        // GT (+) -> revcomp CA(-), canonical AC(+)
        // Total count for canonical AC(+) should be 3 + 2 = 5
        let result2 = find_most_frequent_terminal_digram_suffix_array(&seq, 2, true);
         assert!(result2.is_some());
         let (key2, count2, pos2) = result2.unwrap();
         // Canonical key should still be AC(+)
         assert_eq!(key2, ((SymbolType::Terminal(EncodedBase(0)), Direction::Forward), (SymbolType::Terminal(EncodedBase(1)), Direction::Forward)));
         assert_eq!(count2, 5, "Combined frequency of AC and GT"); 
         // Positions should include those for AC (0, 2, 6) and GT (4, 8)
         let expected_pos2: HashSet<usize> = [0, 2, 4, 6, 8].iter().cloned().collect();
         let actual_pos2: HashSet<usize> = pos2.into_iter().collect();
         assert_eq!(actual_pos2, expected_pos2);


        // Test 3: min_count = 6 (nothing should be found)
        let result3 = find_most_frequent_terminal_digram_suffix_array(&seq, 6, true);
        assert!(result3.is_none());
    }
} 