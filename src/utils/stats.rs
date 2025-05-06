use crate::encode::dna_2bit::EncodedBase;
use std::collections::HashMap;
use std::fmt;

/// Statistics for a DNA sequence
#[derive(Debug, Clone)]
pub struct SequenceStats {
    /// Total length of the sequence
    pub length: usize,
    /// Count of each base (A, C, G, T)
    pub base_counts: HashMap<char, usize>,
    /// GC content as a percentage
    pub gc_percentage: f64,
    /// Count of ambiguous bases (N)
    pub ambiguous_count: usize,
    /// Longest homopolymer run
    pub max_homopolymer: usize,
    /// Most common kmer
    pub most_common_kmer: Option<(String, usize)>,
}

impl SequenceStats {
    /// Calculate sequence statistics from a vector of encoded bases
    pub fn from_encoded_bases(bases: &[EncodedBase], k: Option<usize>) -> Self {
        // Initialize base counts
        let mut base_counts = HashMap::new();
        base_counts.insert('A', 0);
        base_counts.insert('C', 0);
        base_counts.insert('G', 0);
        base_counts.insert('T', 0);

        // Calculate base frequencies
        for base in bases {
            let ch = base.to_char();
            *base_counts.entry(ch).or_insert(0) += 1;
        }

        // Calculate GC percentage
        let gc_count = base_counts.get(&'G').unwrap_or(&0) + base_counts.get(&'C').unwrap_or(&0);
        let gc_percentage = if bases.is_empty() {
            0.0
        } else {
            (gc_count as f64 / bases.len() as f64) * 100.0
        };

        // Find the longest homopolymer
        let mut max_homopolymer = 0;
        let mut current_base = None;
        let mut current_run = 0;

        for base in bases {
            let ch = base.to_char();
            if Some(ch) == current_base {
                current_run += 1;
                if current_run > max_homopolymer {
                    max_homopolymer = current_run;
                }
            } else {
                current_base = Some(ch);
                current_run = 1;
            }
        }

        // Calculate k-mer frequencies if requested
        let most_common_kmer = if let Some(k_size) = k {
            if k_size > 0 && k_size <= bases.len() {
                let mut kmer_counts = HashMap::new();
                
                for i in 0..=bases.len().saturating_sub(k_size) {
                    if i + k_size <= bases.len() {
                        let kmer_string: String = bases[i..i+k_size]
                            .iter()
                            .map(|b| b.to_char())
                            .collect();
                        
                        *kmer_counts.entry(kmer_string).or_insert(0) += 1;
                    }
                }
                
                kmer_counts.into_iter()
                    .max_by_key(|(_, count)| *count)
                    .map(|(kmer, count)| (kmer, count))
            } else {
                None
            }
        } else {
            None
        };

        Self {
            length: bases.len(),
            base_counts,
            gc_percentage,
            ambiguous_count: 0, // Encoded bases don't include ambiguous bases
            max_homopolymer,
            most_common_kmer,
        }
    }

    /// Calculate sequence statistics from a string
    pub fn from_string(sequence: &str, k: Option<usize>) -> Self {
        let mut base_counts = HashMap::new();
        base_counts.insert('A', 0);
        base_counts.insert('C', 0);
        base_counts.insert('G', 0);
        base_counts.insert('T', 0);

        let mut ambiguous_count = 0;

        // Calculate base frequencies
        for ch in sequence.chars().map(|c| c.to_ascii_uppercase()) {
            match ch {
                'A' | 'C' | 'G' | 'T' => {
                    *base_counts.entry(ch).or_insert(0) += 1;
                },
                _ => {
                    ambiguous_count += 1;
                }
            }
        }

        // Calculate GC percentage
        let gc_count = base_counts.get(&'G').unwrap_or(&0) + base_counts.get(&'C').unwrap_or(&0);
        let base_count = base_counts.values().sum::<usize>();
        let gc_percentage = if base_count == 0 {
            0.0
        } else {
            (gc_count as f64 / base_count as f64) * 100.0
        };

        // Find the longest homopolymer
        let mut max_homopolymer = 0;
        let mut current_base = None;
        let mut current_run = 0;

        for ch in sequence.chars().map(|c| c.to_ascii_uppercase()) {
            if Some(ch) == current_base {
                current_run += 1;
                if current_run > max_homopolymer {
                    max_homopolymer = current_run;
                }
            } else {
                current_base = Some(ch);
                current_run = 1;
            }
        }

        // Calculate k-mer frequencies if requested
        let most_common_kmer = if let Some(k_size) = k {
            if k_size > 0 && k_size <= sequence.len() {
                let mut kmer_counts = HashMap::new();
                
                for i in 0..=sequence.len().saturating_sub(k_size) {
                    if i + k_size <= sequence.len() {
                        let kmer_string = &sequence[i..i+k_size].to_uppercase();
                        if !kmer_string.chars().any(|c| c != 'A' && c != 'C' && c != 'G' && c != 'T') {
                            *kmer_counts.entry(kmer_string.to_string()).or_insert(0) += 1;
                        }
                    }
                }
                
                kmer_counts.into_iter()
                    .max_by_key(|(_, count)| *count)
                    .map(|(kmer, count)| (kmer, count))
            } else {
                None
            }
        } else {
            None
        };

        Self {
            length: sequence.len(),
            base_counts,
            gc_percentage,
            ambiguous_count,
            max_homopolymer,
            most_common_kmer,
        }
    }
}

impl fmt::Display for SequenceStats {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "Sequence Length: {}", self.length)?;
        writeln!(f, "Base Counts:")?;
        writeln!(f, "  A: {}", self.base_counts.get(&'A').unwrap_or(&0))?;
        writeln!(f, "  C: {}", self.base_counts.get(&'C').unwrap_or(&0))?;
        writeln!(f, "  G: {}", self.base_counts.get(&'G').unwrap_or(&0))?;
        writeln!(f, "  T: {}", self.base_counts.get(&'T').unwrap_or(&0))?;
        writeln!(f, "GC Content: {:.2}%", self.gc_percentage)?;
        writeln!(f, "Ambiguous Bases: {}", self.ambiguous_count)?;
        writeln!(f, "Longest Homopolymer: {}", self.max_homopolymer)?;
        
        if let Some((kmer, count)) = &self.most_common_kmer {
            writeln!(f, "Most Common Kmer: {} (count: {})", kmer, count)?;
        }
        
        Ok(())
    }
}

/// Compute metrics for a grammar
#[derive(Debug, Clone)]
pub struct GrammarStats {
    /// Number of rules in the grammar
    pub rule_count: usize,
    /// Final sequence length (after grammar construction)
    pub final_sequence_length: usize,
    /// Original sequence length (before grammar construction)
    pub original_sequence_length: usize,
    /// Compression ratio (original / final)
    pub compression_ratio: f64,
    /// Maximum rule depth (recursion level)
    pub max_rule_depth: usize,
    /// Average rule usage
    pub avg_rule_usage: f64,
    /// Maximum rule usage
    pub max_rule_usage: usize,
}

impl fmt::Display for GrammarStats {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "Grammar Statistics:")?;
        writeln!(f, "  Rule Count: {}", self.rule_count)?;
        writeln!(f, "  Final Sequence Length: {}", self.final_sequence_length)?;
        writeln!(f, "  Original Sequence Length: {}", self.original_sequence_length)?;
        writeln!(f, "  Compression Ratio: {:.2}x", self.compression_ratio)?;
        writeln!(f, "  Maximum Rule Depth: {}", self.max_rule_depth)?;
        writeln!(f, "  Average Rule Usage: {:.2}", self.avg_rule_usage)?;
        writeln!(f, "  Maximum Rule Usage: {}", self.max_rule_usage)?;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::encode::dna_2bit::encode_dna;

    #[test]
    fn test_sequence_stats_from_string() {
        let seq = "ACGTACGTNNACGT";
        let stats = SequenceStats::from_string(seq, None);
        
        assert_eq!(stats.length, 14);
        assert_eq!(*stats.base_counts.get(&'A').unwrap(), 3);
        assert_eq!(*stats.base_counts.get(&'C').unwrap(), 3);
        assert_eq!(*stats.base_counts.get(&'G').unwrap(), 3);
        assert_eq!(*stats.base_counts.get(&'T').unwrap(), 3);
        assert_eq!(stats.ambiguous_count, 2);
        assert_eq!(stats.gc_percentage, 50.0);
        assert_eq!(stats.max_homopolymer, 2); // NN
    }
    
    #[test]
    fn test_sequence_stats_from_encoded() {
        let seq = b"ACGTACGTACGT";
        let encoded = encode_dna(seq);
        let stats = SequenceStats::from_encoded_bases(&encoded, None);
        
        assert_eq!(stats.length, 12);
        assert_eq!(*stats.base_counts.get(&'A').unwrap(), 3);
        assert_eq!(*stats.base_counts.get(&'C').unwrap(), 3);
        assert_eq!(*stats.base_counts.get(&'G').unwrap(), 3);
        assert_eq!(*stats.base_counts.get(&'T').unwrap(), 3);
        assert_eq!(stats.gc_percentage, 50.0);
        assert_eq!(stats.ambiguous_count, 0);
        assert_eq!(stats.max_homopolymer, 1);
    }
    
    #[test]
    fn test_most_common_kmer() {
        let seq = "ACGTACGTACGT";
        let stats = SequenceStats::from_string(seq, Some(3));
        
        assert!(stats.most_common_kmer.is_some());
        let (kmer, count) = stats.most_common_kmer.unwrap();
        assert_eq!(kmer, "ACG");
        assert_eq!(count, 3);
    }
    
    #[test]
    fn test_display() {
        let seq = "ACGTACGT";
        let stats = SequenceStats::from_string(seq, None);
        
        let display = format!("{}", stats);
        assert!(display.contains("Sequence Length: 8"));
        assert!(display.contains("GC Content: 50.00%"));
    }
} 