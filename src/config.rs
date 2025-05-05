use anyhow::Result;
use std::path::PathBuf;

/// Configuration settings derived from CLI arguments.
#[derive(Debug)]
pub struct Config {
    pub input_path: PathBuf,
    pub output_json_path: Option<PathBuf>,
    pub output_text_path: Option<PathBuf>,
    pub output_gfa_path: Option<PathBuf>,
    pub output_dot_path: Option<PathBuf>,
    pub output_fasta_path: Option<PathBuf>,
    pub kmer_size: usize,
    pub min_rule_usage: usize,
    pub max_rule_count: Option<usize>,
    pub reverse_complement_aware: bool,
    pub chunk_size: Option<usize>,
    pub skip_ns: bool,
    pub stats: bool,
    pub visualize: bool,
}

impl Config {
    /// Creates a new Config instance from validated arguments.
    /// Validation should happen before calling this.
    pub fn new(
        input_path: PathBuf,
        output_json_path: Option<PathBuf>,
        output_text_path: Option<PathBuf>,
        output_gfa_path: Option<PathBuf>,
        output_dot_path: Option<PathBuf>,
        output_fasta_path: Option<PathBuf>,
        kmer_size: usize,
        min_rule_usage: usize,
        max_rule_count: Option<usize>,
        reverse_complement_aware: bool,
        chunk_size: Option<usize>,
        skip_ns: bool,
        stats: bool,
        visualize: bool,
    ) -> Self {
        Config {
            input_path,
            output_json_path,
            output_text_path,
            output_gfa_path,
            output_dot_path,
            output_fasta_path,
            kmer_size,
            min_rule_usage,
            max_rule_count,
            reverse_complement_aware,
            chunk_size,
            skip_ns,
            stats,
            visualize,
        }
    }

    /// Validates the configuration settings.
    pub fn validate(&self) -> Result<()> {
        if !self.input_path.exists() {
            anyhow::bail!("Input file not found: {}", self.input_path.display());
        }
        
        if self.kmer_size == 0 { 
            anyhow::bail!("k-mer size must be greater than 0.");
        }
        
        if self.kmer_size > 255 {
            anyhow::bail!("k-mer size must be less than or equal to 255.");
        }
        
        if let Some(chunk_size) = self.chunk_size {
            if chunk_size < 1000 {
                anyhow::bail!("Chunk size must be at least 1000 bases.");
            }
        }
        
        Ok(())
    }
} 