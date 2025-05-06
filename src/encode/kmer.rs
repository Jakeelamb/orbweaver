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

    #[test]
    fn test_kmer_new() {
        let bases = vec![
            EncodedBase(0b00), // A
            EncodedBase(0b01), // C
            EncodedBase(0b10), // G
        ];
        let kmer = KMer::new(bases);
        assert_eq!(kmer.k, 3);
        assert_eq!(kmer.to_string(), "ACG");
    }

    #[test]
    fn test_kmer_from_string() {
        let kmer = KMer::from_string("ACGT").unwrap();
        assert_eq!(kmer.k, 4);
        assert_eq!(kmer.to_string(), "ACGT");
        
        // Test invalid sequence
        assert!(KMer::from_string("ACNGT").is_err());
    }

    #[test]
    fn test_reverse_complement() {
        let kmer = KMer::from_string("ACGT").unwrap();
        let revcomp = kmer.reverse_complement();
        assert_eq!(revcomp.to_string(), "ACGT");
        
        let kmer2 = KMer::from_string("AAGT").unwrap();
        let revcomp2 = kmer2.reverse_complement();
        assert_eq!(revcomp2.to_string(), "ACTT");
    }

    #[test]
    fn test_canonical() {
        // ACGT is palindromic when reverse complemented
        let kmer1 = KMer::from_string("ACGT").unwrap();
        assert!(kmer1.is_canonical());
        assert_eq!(kmer1.canonical().to_string(), "ACGT");
        
        // AAGT < ACTT (reverse complement)
        let kmer2 = KMer::from_string("AAGT").unwrap();
        assert!(kmer2.is_canonical());
        assert_eq!(kmer2.canonical().to_string(), "AAGT");
        
        // CCGT > ACGG (reverse complement)
        let kmer3 = KMer::from_string("CCGT").unwrap();
        assert!(!kmer3.is_canonical());
        assert_eq!(kmer3.canonical().to_string(), "ACGG");
    }
    
    #[test]
    fn test_cmp() {
        let kmer1 = KMer::from_string("ACGT").unwrap();
        let kmer2 = KMer::from_string("ACGA").unwrap();
        let kmer3 = KMer::from_string("ACG").unwrap();
        
        // T > A, so ACGT > ACGA
        assert!(kmer1 > kmer2);
        
        // Longer sequence > shorter sequence when common prefix is equal
        assert!(kmer1 > kmer3);
        assert!(kmer2 > kmer3);
    }

    #[test]
    fn test_extract_kmers() {
        let sequence = vec![
            EncodedBase(0b00), // A
            EncodedBase(0b01), // C
            EncodedBase(0b10), // G
            EncodedBase(0b11), // T
            EncodedBase(0b00), // A
        ];
        
        let kmers = extract_kmers(&sequence, 3, false);
        assert_eq!(kmers.len(), 3);
        assert_eq!(kmers[0].to_string(), "ACG");
        assert_eq!(kmers[1].to_string(), "CGT");
        assert_eq!(kmers[2].to_string(), "GTA");
        
        // Test canonical mode
        let canonical_kmers = extract_kmers(&sequence, 3, true);
        assert_eq!(canonical_kmers.len(), 3);
        // These should be canonicalized if their reverse complement is smaller
        assert_eq!(canonical_kmers[0].to_string(), "ACG"); // ACG < CGT
        assert_eq!(canonical_kmers[1].to_string(), "ACG"); // ACG < CGT 
        assert_eq!(canonical_kmers[2].to_string(), "TAC"); // TAC < GTA
    }

    #[test]
    fn test_count_kmers() {
        let sequence = vec![
            EncodedBase(0b00), // A
            EncodedBase(0b01), // C
            EncodedBase(0b10), // G
            EncodedBase(0b11), // T
            EncodedBase(0b00), // A
            EncodedBase(0b01), // C
            EncodedBase(0b10), // G
        ];
        
        let counts = count_kmers(&sequence, 3, false);
        assert_eq!(counts.len(), 5);
        assert_eq!(counts[&KMer::from_string("ACG").unwrap()], 2);
        assert_eq!(counts[&KMer::from_string("CGT").unwrap()], 1);
        assert_eq!(counts[&KMer::from_string("GTA").unwrap()], 1);
        assert_eq!(counts[&KMer::from_string("TAC").unwrap()], 1);
        assert_eq!(counts.get(&KMer::from_string("CGA").unwrap()), None); // Not present
    }
} 