use anyhow::{Context, Result};
use log::{debug, info};
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::Path;
use std::rc::Rc;
use tempfile;

use crate::genome::Chromosome;

/// Represents a motif in the SLP
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct Motif {
    /// The sequence of the motif
    pub sequence: String,
    
    /// Frequency of the motif in the chromosome
    pub frequency: usize,
    
    /// Length of the motif
    pub length: usize,
    
    /// Iteration when this motif was created (0 for base nucleotides)
    pub creation_iteration: usize,
}

/// Represents a rule in the SLP
#[derive(Debug, Clone)]
pub struct Rule {
    /// The left-hand side of the rule (new motif)
    pub lhs: Rc<Motif>,
    
    /// The first motif on the right-hand side
    pub rhs1: Rc<Motif>,
    
    /// The second motif on the right-hand side
    pub rhs2: Rc<Motif>,
}

/// Represents a Straight-Line Program
#[derive(Debug, Clone)]
pub struct SLP {
    /// The rules of the SLP
    pub rules: Vec<Rule>,
    
    /// The dictionary of motifs, keyed by their sequence
    pub motifs: HashMap<String, Rc<Motif>>,
    
    /// The compression ratio achieved by the SLP
    pub compression_ratio: f64,
    
    /// The assembly index (number of unique motifs - 4)
    pub assembly_index: usize,
}

/// Builds an SLP from a sequence
pub struct SLPBuilder {
    /// The current rules of the SLP
    rules: Vec<Rule>,
    
    /// The current dictionary of motifs
    motifs: HashMap<String, Rc<Motif>>,
    
    /// The current iteration count
    iteration: usize,
    
    /// Normalized frequency threshold for stopping
    threshold_factor: f64,
}

impl SLPBuilder {
    /// Create a new SLPBuilder with a default threshold factor
    pub fn new() -> Self {
        Self {
            rules: Vec::new(),
            motifs: HashMap::new(),
            iteration: 0,
            threshold_factor: 1_000_000.0,
        }
    }
    
    /// Set the threshold factor for stopping the SLP construction
    pub fn with_threshold_factor(mut self, factor: f64) -> Self {
        self.threshold_factor = factor;
        self
    }
    
    /// Add a base motif (nucleotide) to the dictionary
    fn add_base_motif(&mut self, base: char, frequency: usize) {
        let sequence = base.to_string();
        let motif = Rc::new(Motif {
            sequence: sequence.clone(),
            frequency,
            length: 1,
            creation_iteration: 0,
        });
        self.motifs.insert(sequence, motif);
    }
    
    /// Initialize the SLP with base nucleotides
    fn initialize(&mut self, sequence: &str) {
        // Count frequencies of base nucleotides
        let mut base_freqs = HashMap::new();
        for base in sequence.chars() {
            *base_freqs.entry(base).or_insert(0) += 1;
        }
        
        // Add base nucleotides to dictionary
        for (base, freq) in base_freqs {
            self.add_base_motif(base, freq);
        }
        
        debug!("Initialized SLP with {} base motifs", self.motifs.len());
    }
    
    /// Find the most frequent adjacent pair of motifs in the sequence representation
    fn find_most_frequent_pair(&self, representation: &[Rc<Motif>]) -> Option<(Rc<Motif>, Rc<Motif>, usize)> {
        let mut pair_frequencies = HashMap::new();
        
        // Count frequencies of adjacent pairs
        for window in representation.windows(2) {
            if let [left, right] = window {
                *pair_frequencies.entry((Rc::clone(left), Rc::clone(right))).or_insert(0) += 1;
            }
        }
        
        // Find the most frequent pair
        pair_frequencies
            .into_iter()
            .max_by(|(pair1, freq1), (pair2, freq2)| {
                // Compare by frequency
                let freq_cmp = freq1.cmp(freq2);
                if freq_cmp != std::cmp::Ordering::Equal {
                    return freq_cmp;
                }
                
                // Tie-breaking: lexicographically by combined sequence
                let seq1 = format!("{}{}", pair1.0.sequence, pair1.1.sequence);
                let seq2 = format!("{}{}", pair2.0.sequence, pair2.1.sequence);
                seq1.cmp(&seq2)
            })
            .map(|((left, right), freq)| (left, right, freq))
    }
    
    /// Build an SLP from a sequence
    pub fn build_from_sequence(&mut self, sequence: &str) -> Result<SLP> {
        self.initialize(sequence);
        
        // Create initial representation as individual characters
        let mut representation: Vec<Rc<Motif>> = sequence
            .chars()
            .map(|c| Rc::clone(self.motifs.get(&c.to_string()).unwrap()))
            .collect();
        
        let original_length = sequence.len();
        let stop_threshold = (original_length as f64 / self.threshold_factor).max(2.0) as usize;
        
        info!("Starting SLP construction with stop threshold: {}", stop_threshold);
        
        // Main loop: find and replace the most frequent pair
        loop {
            self.iteration += 1;
            
            // Find the most frequent pair
            if let Some((left, right, frequency)) = self.find_most_frequent_pair(&representation) {
                // Check if we should stop
                if frequency < stop_threshold {
                    info!(
                        "Stopping SLP construction at iteration {} (frequency {} < threshold {})",
                        self.iteration, frequency, stop_threshold
                    );
                    break;
                }
                
                // Create new motif for the pair
                let new_sequence = format!("{}{}", left.sequence, right.sequence);
                let new_length = left.length + right.length;
                
                let new_motif = Rc::new(Motif {
                    sequence: new_sequence.clone(),
                    frequency,
                    length: new_length,
                    creation_iteration: self.iteration,
                });
                
                // Add rule
                let rule = Rule {
                    lhs: Rc::clone(&new_motif),
                    rhs1: Rc::clone(&left),
                    rhs2: Rc::clone(&right),
                };
                
                self.rules.push(rule);
                self.motifs.insert(new_sequence.clone(), Rc::clone(&new_motif));
                
                // Update representation
                let mut i = 0;
                while i < representation.len() - 1 {
                    if Rc::ptr_eq(&representation[i], &left) && Rc::ptr_eq(&representation[i + 1], &right) {
                        // Replace pair with new motif
                        representation.remove(i + 1);
                        representation[i] = Rc::clone(&new_motif);
                    } else {
                        i += 1;
                    }
                }
                
                debug!(
                    "Iteration {}: Replaced pair '{}/{}' with '{}', frequency: {}, rules: {}",
                    self.iteration, left.sequence, right.sequence, new_sequence, frequency, self.rules.len()
                );
            } else {
                // No more pairs to replace
                break;
            }
        }
        
        // Calculate assembly index and compression ratio
        let assembly_index = self.motifs.len() - 4; // Subtract the 4 base nucleotides
        let final_length = representation.len();
        let compression_ratio = original_length as f64 / final_length as f64;
        
        info!(
            "SLP construction completed: {} iterations, {} motifs, assembly index: {}, compression ratio: {:.2}",
            self.iteration, self.motifs.len(), assembly_index, compression_ratio
        );
        
        Ok(SLP {
            rules: self.rules.clone(),
            motifs: self.motifs.clone(),
            compression_ratio,
            assembly_index,
        })
    }
    
    /// Build an SLP from a VBQ file (using bqtools to decode)
    pub fn build_from_vbq(&mut self, chromosome: &Chromosome) -> Result<SLP> {
        // Use temporary file for decoding
        let temp_file = tempfile::NamedTempFile::new()
            .context("Failed to create temporary file for VBQ decoding")?;
        
        // Run bqtools decode to temporary file
        crate::utils::run_bqtools_decode(&chromosome.vbq_path, temp_file.path())?;
        
        // Read the decoded sequence
        let file = File::open(temp_file.path())
            .context("Failed to open decoded VBQ file")?;
        
        let reader = BufReader::new(file);
        let mut sequence = String::new();
        
        // Skip header line
        let mut lines = reader.lines();
        lines.next(); // Skip FASTA header
        
        // Read sequence lines
        for line in lines {
            let line = line?;
            if !line.starts_with('>') {
                sequence.push_str(line.trim());
            }
        }
        
        // Build SLP
        self.build_from_sequence(&sequence)
    }
}

impl Default for SLPBuilder {
    fn default() -> Self {
        Self::new()
    }
} 