use anyhow::{Context, Result};
use bio::io::fasta;
use log::{debug, info, warn};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use std::fs::File;
use std::path::{Path, PathBuf};
use std::collections::HashMap;

use crate::utils::{
    calculate_n_percentage, clean_chromosome_id, ensure_dir_exists, generate_chromosome_vbq_filename,
    run_bqtools_encode, run_bqtools_decode, validate_sequence_length,
};

/// Represents a chromosome or scaffold from a genome
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Chromosome {
    /// Unique identifier for the chromosome
    pub id: String,
    
    /// Length of the chromosome in base pairs after VBQ processing
    pub processed_length: usize,
    
    /// Original length of the chromosome in base pairs
    pub original_length: usize,
    
    /// Percentage of 'N' bases in the original sequence
    pub original_n_percent: f64,
    
    /// Path to the VBQ file containing the processed chromosome sequence
    pub vbq_path: PathBuf,
    
    /// Assembly index calculated for this chromosome (filled later)
    pub assembly_index: Option<usize>,
    
    /// Graph depth calculated for this chromosome (filled later)
    pub graph_depth: Option<usize>,
    
    /// Average motif length for this chromosome (filled later)
    pub avg_motif_length: Option<f64>,
}

/// Represents a complete genome with multiple chromosomes/scaffolds
#[derive(Debug, Serialize, Deserialize)]
pub struct Genome {
    /// Unique identifier for the genome
    pub id: String,
    
    /// Scientific name of the species
    pub species_name: String,
    
    /// Related group (taxonomic group, clade, etc.)
    pub related_group: String,
    
    /// Average percentage of 'N' bases across all chromosomes
    pub original_n_percent_avg: f64,
    
    /// Total genome size in base pairs after processing
    pub total_processed_size: usize,
    
    /// Average assembly index across all chromosomes (weighted by size)
    pub assembly_index_avg: Option<f64>,
    
    /// Chromosomes contained in this genome
    pub chromosomes: Vec<Chromosome>,
}

/// Processes genomes from FASTA files
pub struct GenomeProcessor {
    /// Path to the input FASTA file
    input_path: PathBuf,
    
    /// Directory for storing output files
    output_dir: PathBuf,
}

impl GenomeProcessor {
    /// Create a new GenomeProcessor
    pub fn new(input_path: &str, output_dir: &str) -> Result<Self> {
        let input_path = PathBuf::from(input_path);
        let output_dir = PathBuf::from(output_dir);
        
        ensure_dir_exists(&output_dir)?;
        
        Ok(Self {
            input_path,
            output_dir,
        })
    }
    
    /// Process the genome from FASTA file
    pub fn process(&self) -> Result<Genome> {
        // Create VBQ directory
        let vbq_dir = self.output_dir.join("vbq");
        ensure_dir_exists(&vbq_dir)?;
        
        // Create temp directory for processing
        let temp_dir = self.output_dir.join("temp");
        ensure_dir_exists(&temp_dir)?;
        
        // Read FASTA file
        let reader = fasta::Reader::from_file(&self.input_path)
            .context(format!("Failed to read FASTA file: {}", self.input_path.display()))?;
        
        // Extract genome ID from file name (will be overridden by the CLI argument)
        let genome_id = self.input_path
            .file_stem()
            .and_then(|s| s.to_str())
            .unwrap_or("unknown")
            .to_string();
        
        // Process each chromosome in parallel
        let chromosomes_result: Result<Vec<_>> = reader
            .records()
            .collect::<Result<Vec<_>, _>>()
            .context("Failed to read FASTA records")?
            .par_iter()
            .map(|record| self.process_chromosome(&vbq_dir, &temp_dir, &genome_id, record))
            .collect();
        
        let chromosomes = chromosomes_result?;
        
        // Calculate genome-level metrics
        let total_processed_size: usize = chromosomes.iter().map(|c| c.processed_length).sum();
        let original_n_percent_avg = if !chromosomes.is_empty() {
            let weighted_sum: f64 = chromosomes
                .iter()
                .map(|c| c.original_n_percent * c.original_length as f64)
                .sum();
            let total_length: usize = chromosomes.iter().map(|c| c.original_length).sum();
            weighted_sum / total_length as f64
        } else {
            0.0
        };
        
        // Create genome object
        let genome = Genome {
            id: genome_id,
            species_name: "Unknown".to_string(), // Will be replaced from CLI args
            related_group: "Unknown".to_string(), // Will be replaced from CLI args
            original_n_percent_avg,
            total_processed_size,
            assembly_index_avg: None, // Will be calculated after SLP construction
            chromosomes,
        };
        
        Ok(genome)
    }
    
    /// Process a single chromosome from a FASTA record
    fn process_chromosome(
        &self,
        vbq_dir: &Path,
        temp_dir: &Path,
        genome_id: &str,
        record: &fasta::Record,
    ) -> Result<Chromosome> {
        let raw_id = record.id();
        let chromosome_id = clean_chromosome_id(raw_id);
        let sequence = std::str::from_utf8(record.seq())?;
        
        // Calculate original N percentage
        let original_n_percent = calculate_n_percentage(sequence);
        let original_length = sequence.len();
        
        info!("Processing chromosome {}: {} bp", chromosome_id, original_length);
        
        // Create temporary FASTA file for this chromosome
        let temp_fasta = temp_dir.join(format!("{}_{}.fa", genome_id, chromosome_id));
        let mut temp_file = File::create(&temp_fasta)
            .context(format!("Failed to create temporary FASTA file: {}", temp_fasta.display()))?;
        
        // Write chromosome to temporary FASTA file
        writeln!(temp_file, ">{}", raw_id)?;
        writeln!(temp_file, "{}", sequence)?;
        
        // Run bqtools encode
        let vbq_path = generate_chromosome_vbq_filename(vbq_dir, genome_id, &chromosome_id);
        run_bqtools_encode(&temp_fasta, &vbq_path, "r")?;
        
        // Validate VBQ file
        let decoded_fasta = temp_dir.join(format!("{}_{}_decoded.fa", genome_id, chromosome_id));
        run_bqtools_decode(&vbq_path, &decoded_fasta)?;
        
        let (_, processed_length) = validate_sequence_length(&temp_fasta, &decoded_fasta)?;
        
        // Create chromosome object
        let chromosome = Chromosome {
            id: chromosome_id,
            processed_length,
            original_length,
            original_n_percent,
            vbq_path,
            assembly_index: None,
            graph_depth: None,
            avg_motif_length: None,
        };
        
        Ok(chromosome)
    }
} 