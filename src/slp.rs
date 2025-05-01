use std::collections::HashMap;
use std::fs::File;
use std::io::BufReader;
use bio::io::fasta;
use std::path::PathBuf;
use anyhow::{Context, Result};

// Represents the SLP rules: non-terminal ID -> (left_symbol, right_symbol)
// Symbols < 4 are terminals (A=0, C=1, G=2, T=3)
// Symbols >= 4 are non-terminals
pub type SLPRules = HashMap<u32, (u32, u32)>;

// Encodes a nucleotide sequence into a vector of u32 symbols.
// A=0, C=1, G=2, T=3. Non-ACGT characters are ignored/skipped.
fn encode_sequence(sequence: &str) -> Vec<u32> {
    sequence
        .chars()
        .filter_map(|c| match c {
            'A' | 'a' => Some(0),
            'C' | 'c' => Some(1),
            'G' | 'g' => Some(2),
            'T' | 't' => Some(3),
            _ => None, // Skip non-ACGT characters
        })
        .collect()
}

// Finds the most frequent adjacent pair of symbols in the sequence.
// Returns Option<((symbol1, symbol2), frequency)>
// Tie-breaking: Chooses the pair (a, b) with the smallest 'a', then smallest 'b'.
fn find_most_frequent_pair(sequence: &[u32]) -> Option<((u32, u32), usize)> {
    if sequence.len() < 2 {
        return None;
    }

    let mut pair_counts: HashMap<(u32, u32), usize> = HashMap::new();
    for window in sequence.windows(2) {
        let pair = (window[0], window[1]);
        *pair_counts.entry(pair).or_insert(0) += 1;
    }

    // Find the pair with the maximum frequency, applying tie-breaking
    pair_counts.into_iter().max_by(|a, b| {
        // Compare frequencies first
        a.1.cmp(&b.1)
            // Then compare the first element of the pair (smaller first)
            .then_with(|| b.0 .0.cmp(&a.0 .0))
            // Then compare the second element of the pair (smaller first)
            .then_with(|| b.0 .1.cmp(&a.0 .1))
    })
}

// Replaces occurrences of a specific pair with a new symbol.
// Performs a single pass from left to right, non-overlapping replacement.
fn replace_pair_in_sequence(
    sequence: &[u32],
    pair_to_replace: (u32, u32),
    new_symbol: u32,
) -> Vec<u32> {
    let mut new_sequence = Vec::with_capacity(sequence.len()); // Approximate capacity
    let mut i = 0;
    while i < sequence.len() {
        if i + 1 < sequence.len() && sequence[i] == pair_to_replace.0 && sequence[i + 1] == pair_to_replace.1 {
            new_sequence.push(new_symbol);
            i += 2; // Skip the next element as it's part of the pair
        } else {
            new_sequence.push(sequence[i]);
            i += 1;
        }
    }
    new_sequence
}


/// Builds a Straight-Line Program (SLP) for a given sequence using a greedy BPE-like approach.
///
/// Args:
///     sequence_str: The input sequence (e.g., a chromosome).
///     min_frequency: The minimum frequency threshold for a pair to be replaced.
///                    The process stops when the most frequent pair count is less than this.
///
/// Returns:
///     A tuple containing the SLP rules and the final compressed sequence representation.
pub fn build_slp_for_sequence(
    sequence_str: &str,
    min_frequency: usize,
) -> (SLPRules, Vec<u32>) {
    let mut current_sequence = encode_sequence(sequence_str);
    let mut rules: SLPRules = HashMap::new();
    let mut next_non_terminal_id = 4u32; // Start non-terminals after A,C,G,T

    loop {
        if current_sequence.len() < 2 {
             break; // Cannot form pairs if sequence is too short
        }

        match find_most_frequent_pair(&current_sequence) {
            Some((pair, frequency)) => {
                // Apply stopping condition
                if frequency < min_frequency {
                    break;
                }

                // Add new rule
                rules.insert(next_non_terminal_id, pair);

                // Replace pair in sequence
                current_sequence = replace_pair_in_sequence(
                    &current_sequence,
                    pair,
                    next_non_terminal_id,
                );

                // Increment for the next rule
                next_non_terminal_id += 1;
            }
            None => {
                // No pairs found (e.g., sequence length < 2)
                break;
            }
        }
    }

    (rules, current_sequence)
}

/// Reads a FASTA file, builds SLP for the first sequence, and prints results.
///
/// Args:
///     fasta_path: Path to the input FASTA file.
///     n_factor: Normalization factor for the stopping condition (e.g., 1,000,000.0).
pub fn run_slp_build(fasta_path: PathBuf, n_factor: f64) -> Result<()> {
    println!("Running SLP construction on FASTA file: {}", fasta_path.display());

    let file = File::open(&fasta_path)
        .with_context(|| format!("Failed to open FASTA file: {}", fasta_path.display()))?;
    let reader = fasta::Reader::new(BufReader::new(file));

    let mut record_count = 0;

    for result in reader.records() {
        let record = result.context("Failed to read FASTA record")?;
        record_count += 1;

        if record_count > 1 {
             println!("Warning: Found more than one record in {}. Using the first record only.", fasta_path.display());
             break; 
        }

        let chromosome_name = record.id().to_string();
        let chromosome_sequence = match std::str::from_utf8(record.seq()) {
            Ok(s) => s.to_string(),
            Err(e) => {
                 anyhow::bail!("Error converting sequence to UTF-8 for record {}: {}", record.id(), e);
            }
        };

        println!("Processing sequence ID: {}", chromosome_name);
        println!("Sequence length: {}", chromosome_sequence.len());

         if chromosome_sequence.is_empty() {
             println!("Sequence for {} is empty. Skipping SLP construction.", chromosome_name);
             continue;
         }

        // Calculate Dynamic Minimum Frequency
        let chr_len = chromosome_sequence.len() as f64;
        let calculated_threshold = (chr_len / n_factor).ceil() as usize;
        let min_freq_dynamic = std::cmp::max(2, calculated_threshold); // Ensure min_freq is at least 2

        println!("Calculated minimum frequency threshold (N={}): {}", n_factor, min_freq_dynamic);

        // Run SLP Construction (using function from the same module)
        let (rules, final_sequence) = build_slp_for_sequence(&chromosome_sequence, min_freq_dynamic);

        println!("\nSLP Construction Complete for {}:", chromosome_name);
        println!("  Number of rules generated: {}", rules.len());
        println!("  Length of final compressed sequence: {}", final_sequence.len());
    }

     if record_count == 0 {
         anyhow::bail!("No records found in FASTA file: {}", fasta_path.display());
     }

    Ok(())
}

// Placeholder for main logic that iterates over chromosomes
// pub fn build_slp_for_genome(genome_data: /* VCF-BQT Genome Type */) -> HashMap<String, (SLPRules, Vec<u32>)> {
//     let mut all_slps = HashMap::new();
//     // Logic to extract sequences per chromosome/scaffold from genome_data
//     // for (chr_name, sequence) in extract_sequences(genome_data) {
//     //     let min_freq = 2; // Or calculate based on user's formula if clarified
//     //     let slp_result = build_slp_for_sequence(&sequence, min_freq);
//     //     all_slps.insert(chr_name, slp_result);
//     // }
//     all_slps
// }


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_encode_sequence() {
        assert_eq!(encode_sequence("ACGT"), vec![0, 1, 2, 3]);
        assert_eq!(encode_sequence("acgtN"), vec![0, 1, 2, 3]);
        assert_eq!(encode_sequence("AGAG"), vec![0, 2, 0, 2]);
        assert_eq!(encode_sequence(""), Vec::<u32>::new());
        assert_eq!(encode_sequence("XYZ"), Vec::<u32>::new());
    }

    #[test]
    fn test_find_most_frequent_pair_simple() {
        let seq = encode_sequence("ABABAB"); // Pairs: AB(3), BA(2)
        assert_eq!(find_most_frequent_pair(&seq), Some(((0, 1), 3)));
    }

    #[test]
    fn test_find_most_frequent_pair_tie_breaking() {
        // AC(2), CA(1), AG(1), GA(1) -> AC wins
        let seq1 = encode_sequence("ACAGA");
        assert_eq!(find_most_frequent_pair(&seq1), Some(((0, 1), 2)));

         // AA(2), AC(2) -> AA wins (0,0) vs (0,1) - smaller second element
         let seq2 = encode_sequence("AACAAC");
         assert_eq!(find_most_frequent_pair(&seq2), Some(((0, 0), 2)));

        // GA(1), AC(1), CT(1) -> AC wins (0,1) vs (2,0) vs (1,3) - smallest first element
        let seq3 = encode_sequence("GACT");
        assert_eq!(find_most_frequent_pair(&seq3), Some(((0, 1), 1)));
    }

    #[test]
    fn test_find_most_frequent_pair_empty_short() {
        assert_eq!(find_most_frequent_pair(&[]), None);
        assert_eq!(find_most_frequent_pair(&[0]), None);
    }

    #[test]
    fn test_replace_pair_in_sequence() {
        // Replace AC (0,1) with 4
        let seq = encode_sequence("ACGTACGT"); // [0,1,2,3,0,1,2,3]
        let expected = vec![4, 2, 3, 4, 2, 3];
        assert_eq!(replace_pair_in_sequence(&seq, (0, 1), 4), expected);

        // Replace AA (0,0) with 4
        let seq2 = encode_sequence("AAACAAAA"); // [0,0,0,1,0,0,0,0]
        let expected2 = vec![4, 0, 1, 4, 4];
        assert_eq!(replace_pair_in_sequence(&seq2, (0, 0), 4), expected2);

        // Replace GG (2,2) with 4 - no occurrences
        let seq3 = encode_sequence("ACGTACGT"); // [0,1,2,3,0,1,2,3]
        let expected3 = vec![0, 1, 2, 3, 0, 1, 2, 3];
        assert_eq!(replace_pair_in_sequence(&seq3, (2, 2), 4), expected3);
    }

    #[test]
    fn test_build_slp_basic() {
        // Sequence: ACACAC -> (AC -> 4) -> 444
        // Rules: {4: (0, 1)}
        let (rules, final_seq) = build_slp_for_sequence("ACACAC", 2);
        assert_eq!(final_seq, vec![4, 4, 4]);
        assert_eq!(rules.len(), 1);
        assert_eq!(rules.get(&4), Some(&(0, 1)));
    }

     #[test]
    fn test_build_slp_multi_step() {
        // Sequence: AAAA -> (AA -> 4) -> 44 -> (44 -> 5) -> 5
        // Rules: {4: (0, 0), 5: (4, 4)}
        let (rules, final_seq) = build_slp_for_sequence("AAAA", 2);
        assert_eq!(final_seq, vec![5]);
        assert_eq!(rules.len(), 2);
        assert_eq!(rules.get(&4), Some(&(0, 0)));
        assert_eq!(rules.get(&5), Some(&(4, 4)));
    }

    #[test]
    fn test_build_slp_no_replacement() {
        // Sequence: ACGT - all pairs frequency 1
        let (rules, final_seq) = build_slp_for_sequence("ACGT", 2);
        assert_eq!(final_seq, vec![0, 1, 2, 3]); // Remains unchanged
        assert!(rules.is_empty());
    }

    #[test]
    fn test_build_slp_complex_tie_break() {
        // Sequence: AGAGAG -> AG (3), GA (2) -> Replace AG with 4 -> 4G4G -> GG (1), G4 (1), 4G (1) -> Stop
        // Rules: {4: (0, 2)}
        let (rules, final_seq) = build_slp_for_sequence("AGAGAG", 2);
        // Step 1: pairs = {(0,2): 3, (2,0): 2}. Replace (0,2) with 4. Seq = [4, 2, 4, 2]
        // Step 2: pairs = {(4,2): 2, (2,4): 1}. Replace (4,2) with 5. Seq = [5, 5]
        // Step 3: pairs = {(5,5): 1}. Freq < 2. Stop.
        // Expected: final_seq = [5, 5], rules = {4: (0,2), 5: (4,2)}
        assert_eq!(final_seq, vec![5, 5]);
        assert_eq!(rules.len(), 2);
        assert_eq!(rules.get(&4), Some(&(0, 2))); // AG
        assert_eq!(rules.get(&5), Some(&(4, 2))); // 4G -> (AG)G
    }

     #[test]
    fn test_build_slp_long_sequence() {
        // Sequence: GTGTGTACACAC -> GT (3), TG (2), AC (3), CA (2)
        // Iter 1: Tie between GT (2,3) and AC (0,1). AC wins tie-break. Replace AC with 4.
        // Seq: GTGTGT 4 4 4 -> [2,3,2,3,2,3, 4, 4, 4]
        // Iter 2: Pairs: GT (3), TG (2), 44 (2). GT wins (2,3). Replace GT with 5.
        // Seq: 5 T 5 T 5 T 4 4 4 -> [5, 3, 5, 3, 5, 3, 4, 4, 4]
        // Iter 3: Pairs: 5T (3), T5 (2), 44 (2). 5T wins (5,3). Replace 5T with 6.
        // Seq: 6 6 6 4 4 4 -> [6, 6, 6, 4, 4, 4]
        // Iter 4: Pairs: 66 (2), 44 (2). 44 wins (4,4) vs (6,6). Replace 44 with 7.
        // Seq: 6 6 6 7 4 -> [6, 6, 6, 7, 4]
        // Iter 5: Pairs: 66 (2). Replace 66 with 8.
        // Seq: 8 6 7 4 -> [8, 6, 7, 4]
        // Iter 6: Pairs: (8,6)(1), (6,7)(1), (7,4)(1). Freq < 2. Stop.
        // Expected Rules: {4:(0,1), 5:(2,3), 6:(5,3), 7:(4,4), 8:(6,6)}
        let (rules, final_seq) = build_slp_for_sequence("GTGTGTACACAC", 2);

        assert_eq!(final_seq, vec![8, 6, 7, 4]);
        assert_eq!(rules.len(), 5);
        assert_eq!(rules.get(&4), Some(&(0, 1))); // AC
        assert_eq!(rules.get(&5), Some(&(2, 3))); // GT
        assert_eq!(rules.get(&6), Some(&(5, 3))); // 5T -> (GT)T
        assert_eq!(rules.get(&7), Some(&(4, 4))); // 44 -> (AC)(AC)
        assert_eq!(rules.get(&8), Some(&(6, 6))); // 66 -> (5T)(5T) -> ((GT)T)((GT)T)

    }
}
