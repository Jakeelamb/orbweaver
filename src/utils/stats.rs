use crate::encode::dna_2bit::EncodedBase;
use std::collections::HashMap;
use std::fmt;
use crate::encode::kmer::KMer;

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
        let mut max_homopolymer = if bases.is_empty() { 0 } else { 1 }; // Initialize to 1 if not empty
        let mut current_base = None;
        let mut current_run = 0;

        for base in bases {
            let ch = base.to_char();
            if Some(ch) == current_base {
                current_run += 1;
                max_homopolymer = max_homopolymer.max(current_run); // Use max() instead of if
            } else {
                current_base = Some(ch);
                current_run = 1;
            }
        }

        // Calculate k-mer frequencies if requested
        let most_common_kmer = if let Some(k_size) = k {
            if k_size > 0 && k_size <= bases.len() {
                let mut kmer_counts: HashMap<KMer, usize> = HashMap::new();
                
                for i in 0..=bases.len().saturating_sub(k_size) {
                    if i + k_size <= bases.len() {
                        let kmer_slice = bases[i..i+k_size].to_vec();
                        let kmer = KMer::new(kmer_slice);
                        let canonical_kmer = kmer.canonical();
                        
                        *kmer_counts.entry(canonical_kmer).or_insert(0) += 1;
                    }
                }
                
                kmer_counts.into_iter()
                    .max_by_key(|(_, count)| *count)
                    .map(|(kmer, count)| (kmer.to_string(), count))
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
        // Initialize base counts
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
        let mut max_homopolymer = if sequence.is_empty() { 0 } else { 1 }; // Initialize to 1 if not empty
        let mut current_base = None;
        let mut current_run = 0;

        for ch in sequence.chars().map(|c| c.to_ascii_uppercase()) {
            if Some(ch) == current_base {
                current_run += 1;
                max_homopolymer = max_homopolymer.max(current_run); // Use max() instead of if
            } else {
                current_base = Some(ch);
                current_run = 1;
            }
        }

        // Calculate k-mer frequencies if requested
        let most_common_kmer = if let Some(k_size) = k {
            if k_size > 0 && k_size <= sequence.len() {
                let mut kmer_counts: HashMap<KMer, usize> = HashMap::new();
                
                for i in 0..=sequence.len().saturating_sub(k_size) {
                    if i + k_size <= sequence.len() {
                        let kmer_str = &sequence[i..i+k_size];
                        if let Ok(kmer) = KMer::from_string(kmer_str) {
                            let canonical_kmer = kmer.canonical();
                            *kmer_counts.entry(canonical_kmer).or_insert(0) += 1;
                        } else {
                            // Skip invalid k-mers (containing non-ACGT)
                        }
                    }
                }
                
                kmer_counts.into_iter()
                    .max_by_key(|(_, count)| *count)
                    .map(|(kmer, count)| (kmer.to_string(), count))
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
        // Format GC percentage, handling potential NaN/inf
        if self.gc_percentage.is_finite() {
            writeln!(f, "GC Content: {:.2}%", self.gc_percentage)?;
        } else {
            writeln!(f, "GC Content: N/A")?;
        }
        writeln!(f, "Ambiguous Bases: {}", self.ambiguous_count)?;
        writeln!(f, "Longest Homopolymer: {}", self.max_homopolymer)?;
        
        if let Some((kmer, count)) = &self.most_common_kmer {
            let k_info = if !kmer.is_empty() { 
                format!("{}-mer", kmer.len())
            } else {
                "Kmer".to_string()
            };
            writeln!(f, "Most Common {}: {} (count: {})", k_info, kmer, count)?;
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
    use crate::encode::dna_2bit::EncodedBase;

    #[test]
    fn test_sequence_stats_from_string() {
        let stats = SequenceStats::from_string("ACGTACGT", Some(2));
        assert_eq!(stats.length, 8);
        assert_eq!(stats.base_counts[&'A'], 2);
        assert_eq!(stats.base_counts[&'C'], 2);
        assert_eq!(stats.base_counts[&'G'], 2);
        assert_eq!(stats.base_counts[&'T'], 2);
        assert!((stats.gc_percentage - 50.0).abs() < f64::EPSILON);
        assert_eq!(stats.ambiguous_count, 0);
        assert_eq!(stats.max_homopolymer, 1);
        // Canonical kmers: AC, CG, AC, TA, AC, CG, AC -> AC:4, CG:2, TA:1
        assert_eq!(stats.most_common_kmer, Some(("AC".to_string(), 4))); // Expect AC count 4
    }

    #[test]
    fn test_sequence_stats_from_encoded() {
        let seq = vec![EncodedBase(0), EncodedBase(1), EncodedBase(2), EncodedBase(3), EncodedBase(0), EncodedBase(1)]; // ACGTAC
        let stats = SequenceStats::from_encoded_bases(&seq, Some(2));
        assert_eq!(stats.length, 6);
        assert_eq!(stats.base_counts[&'A'], 2);
        assert_eq!(stats.base_counts[&'C'], 2);
        assert_eq!(stats.base_counts[&'G'], 1);
        assert_eq!(stats.base_counts[&'T'], 1);
        assert!((stats.gc_percentage - 50.0).abs() < f64::EPSILON);
        assert_eq!(stats.ambiguous_count, 0);
        assert_eq!(stats.max_homopolymer, 1);
        // Canonical kmers: AC, CG, AC, TA, AC -> AC:3, CG:1, TA:1
        assert_eq!(stats.most_common_kmer, Some(("AC".to_string(), 3))); // Expect AC count 3
    }

    #[test]
    fn test_most_common_kmer() {
        let seq = "ACGTACGT";
        let stats = SequenceStats::from_string(seq, Some(2));
        // Kmers: AC, CG, GT, TA, AC, CG, GT
        // Canonical: AC, CG, AC, TA, AC, CG, AC
        // Counts: AC=4, CG=2, TA=1
        assert_eq!(stats.most_common_kmer, Some(("AC".to_string(), 4))); // Expect AC, count 4
    }

    #[test]
    fn test_display() {
        let stats = SequenceStats::from_string("ACGT", None);
        let s = format!("{}", stats);
        assert!(s.contains("Sequence Length: 4"));
        assert!(s.contains("A: 1"));
        assert!(s.contains("C: 1"));
        assert!(s.contains("G: 1"));
        assert!(s.contains("T: 1"));
        // Check GC content with a small epsilon for floating point comparison
        let gc_line = s.lines().find(|line| line.starts_with("GC Content:")).unwrap();
        let gc_value_str = gc_line.split(":").nth(1).unwrap().trim().trim_end_matches('%');
        let gc_value: f64 = gc_value_str.parse().unwrap();
        assert!((gc_value - 50.0).abs() < 0.001);
        assert!(s.contains("Ambiguous Bases: 0"));
        assert!(s.contains("Longest Homopolymer: 1"));
    }
} 