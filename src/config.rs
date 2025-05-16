use anyhow::Result;
use std::path::PathBuf;

/// Configuration settings derived from CLI arguments.
#[derive(Debug)]
pub struct Config {
    /// Path to the input FASTA or FASTQ file.
    pub input_path: PathBuf,
    /// Optional path for JSON output of the grammar.
    pub output_json_path: Option<PathBuf>,
    /// Optional path for plain text output of the grammar.
    pub output_text_path: Option<PathBuf>,
    /// Optional path for GFA (Graphical Fragment Assembly) output.
    pub output_gfa_path: Option<PathBuf>,
    /// Optional path for DOT (Graphviz) output for visualization.
    pub output_dot_path: Option<PathBuf>,
    /// Optional path for FASTA output of reconstructed sequences or blocks.
    pub output_fasta_path: Option<PathBuf>,
    /// The k-mer size to be used for analysis.
    pub kmer_size: usize,
    /// Minimum number of times a rule must be used to be kept.
    pub min_rule_usage: usize,
    /// Optional maximum number of rules to generate.
    pub max_rule_count: Option<usize>,
    /// Whether to consider reverse complements of k-mers.
    pub reverse_complement_aware: bool,
    /// Optional chunk size for processing input files, in bytes.
    pub chunk_size: Option<usize>,
    /// Whether to skip sequences containing 'N' characters.
    pub skip_ns: bool,
    /// Whether to output statistics about the compression and grammar.
    pub stats: bool,
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