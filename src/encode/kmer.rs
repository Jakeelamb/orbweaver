use anyhow::{bail, Result};
use super::dna_2bit::EncodedBase;
use std::collections::HashMap;
use std::hash::{Hash, Hasher};
use std::cmp::Ordering;

/// Represents a fixed-length sequence of DNA bases (kmer)
#[derive(Debug, Clone, Eq)]
pub struct KMer {
    /// Vector of encoded bases
    pub bases: Vec<EncodedBase>,
    /// Length of the kmer
    pub k: usize,
}

impl KMer {
    /// Create a new kmer from a vector of encoded bases
    pub fn new(bases: Vec<EncodedBase>) -> Self {
        let k = bases.len();
        Self { bases, k }
    }

    /// Create a kmer from a DNA string
    pub fn from_string(s: &str) -> Result<Self> {
        if s.is_empty() {
            bail!("Cannot create kmer from empty string");
        }

        let bytes = s.as_bytes();
        let bases: Vec<EncodedBase> = bytes
            .iter()
            .filter_map(|&b| EncodedBase::from_base(b))
            .collect();

        if bases.len() != bytes.len() {
            bail!("Invalid DNA sequence for kmer: {}", s);
        }

        Ok(Self::new(bases))
    }

    /// Get the reverse complement of this kmer
    pub fn reverse_complement(&self) -> Self {
        let rev_comp_bases: Vec<EncodedBase> = self.bases
            .iter()
            .rev()
            .map(|base| base.revcomp())
            .collect();

        Self::new(rev_comp_bases)
    }

    /// Convert to string representation
    pub fn to_string(&self) -> String {
        self.bases.iter().map(|base| base.to_char()).collect()
    }

    /// Check if this kmer is canonical (lexicographically smaller than its reverse complement)
    pub fn is_canonical(&self) -> bool {
        let rev_comp = self.reverse_complement();
        self <= &rev_comp
    }

    /// Get the canonical representation (lexicographically smaller of self or reverse complement)
    pub fn canonical(&self) -> Self {
        let rev_comp = self.reverse_complement();
        if self <= &rev_comp {
            self.clone()
        } else {
            rev_comp
        }
    }
}

impl PartialEq for KMer {
    fn eq(&self, other: &Self) -> bool {
        self.bases == other.bases
    }
}

impl Hash for KMer {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.bases.hash(state);
    }
}

impl PartialOrd for KMer {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for KMer {
    fn cmp(&self, other: &Self) -> Ordering {
        // Compare elements pairwise
        for (a, b) in self.bases.iter().zip(other.bases.iter()) {
            match a.0.cmp(&b.0) {
                Ordering::Equal => continue,
                ordering => return ordering,
            }
        }
        
        // If common prefix is equal, shorter sequence is smaller
        self.bases.len().cmp(&other.bases.len())
    }
}

/// Extracts kmers from a DNA sequence
pub fn extract_kmers(
    sequence: &[EncodedBase], 
    k: usize, 
    canonical: bool
) -> Vec<KMer> {
    if sequence.len() < k {
        return vec![];
    }

    let mut kmers = Vec::with_capacity(sequence.len() - k + 1);
    
    for i in 0..=sequence.len() - k {
        let kmer_bases = sequence[i..i+k].to_vec();
        let mut kmer = KMer::new(kmer_bases);
        
        if canonical {
            kmer = kmer.canonical();
        }
        
        kmers.push(kmer);
    }
    
    kmers
}

/// Counts kmer occurrences in a sequence
pub fn count_kmers(
    sequence: &[EncodedBase], 
    k: usize, 
    canonical: bool
) -> HashMap<KMer, usize> {
    let mut counts = HashMap::new();
    
    for i in 0..=sequence.len().saturating_sub(k) {
        if i + k <= sequence.len() {
            let kmer_bases = sequence[i..i+k].to_vec();
            let mut kmer = KMer::new(kmer_bases);
            
            if canonical {
                kmer = kmer.canonical();
            }
            
            *counts.entry(kmer).or_insert(0) += 1;
        }
    }
    
    counts
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::encode::dna_2bit::EncodedBase;

    #[test]
    fn test_kmer_new() {
        let bases = vec![EncodedBase(0), EncodedBase(1), EncodedBase(2)];
        let kmer = KMer::new(bases.clone());
        assert_eq!(kmer.bases, bases);
        assert_eq!(kmer.k, 3);
    }

    #[test]
    fn test_kmer_from_string() {
        let kmer = KMer::from_string("ACG").unwrap();
        assert_eq!(kmer.bases, vec![EncodedBase(0), EncodedBase(1), EncodedBase(2)]);
        assert_eq!(kmer.k, 3);
    }

    #[test]
    fn test_reverse_complement() {
        let kmer = KMer::from_string("AAGT").unwrap();
        let rev = kmer.reverse_complement();
        assert_eq!(rev.bases, vec![EncodedBase(0), EncodedBase(1), EncodedBase(3), EncodedBase(3)]);
    }

    #[test]
    fn test_canonical() {
        let kmer = KMer::from_string("ACG").unwrap();
        let rev = kmer.reverse_complement();
        let can = kmer.canonical();
        let rev_can = rev.canonical();
        assert_eq!(can, rev_can);
    }

    #[test]
    fn test_cmp() {
        let kmer1 = KMer::from_string("ACG").unwrap();
        let kmer2 = KMer::from_string("CGT").unwrap();
        assert!(kmer1 < kmer2);
    }

    #[test]
    fn test_extract_kmers() {
        let seq = vec![EncodedBase(0), EncodedBase(1), EncodedBase(2), EncodedBase(3)];
        let kmers = extract_kmers(&seq, 2, true);
        // With canonicalization, "AC" and "GT" are canonicalized to the same k-mer
        assert_eq!(kmers.len(), 3);
        assert!(kmers.iter().any(|k| k.bases == vec![EncodedBase(0), EncodedBase(1)]));
    }

    #[test]
    fn test_count_kmers() {
        let seq = vec![EncodedBase(0), EncodedBase(1), EncodedBase(2), EncodedBase(3), EncodedBase(0), EncodedBase(1)]; // ACGTAC
        let counts = count_kmers(&seq, 2, true); // k=2, canonical=true
        
        // Kmers: AC, CG, GT, TA, AC
        // RevComps: GT, CG, AC, TA, GT
        // Canonical: AC, CG, AC, TA, AC
        // Expected counts: AC: 3, CG: 1, TA: 1
        
        let kmer_ac = KMer::from_string("AC").unwrap();
        let kmer_cg = KMer::from_string("CG").unwrap();
        let kmer_ta = KMer::from_string("TA").unwrap();

        assert_eq!(counts.len(), 3, "Expected 3 distinct canonical k-mers");
        assert_eq!(counts.get(&kmer_ac.canonical()), Some(&3), "Canonical AC count should be 3"); // Updated assertion from 2 to 3
        assert_eq!(counts.get(&kmer_cg.canonical()), Some(&1));
        assert_eq!(counts.get(&kmer_ta.canonical()), Some(&1));
    }
} 