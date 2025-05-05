use anyhow::Result;
use std::path::PathBuf;

/// Configuration settings derived from CLI arguments.
#[derive(Debug)]
pub struct Config {
    pub input_path: PathBuf,
    pub output_json_path: Option<PathBuf>,
    pub kmer_size: usize,
    // Add other configuration fields as needed by later tasks
    // pub skip_ns: bool,
    // pub chunk_size: Option<usize>,
    // ...
}

impl Config {
    /// Creates a new Config instance from validated arguments.
    /// Validation should happen before calling this.
    pub fn new(input_path: PathBuf, output_json_path: Option<PathBuf>, kmer_size: usize) -> Self {
        Config {
            input_path,
            output_json_path,
            kmer_size,
        }
    }

    // Potential validation function (or could be done in main.rs)
    // pub fn validate(&self) -> Result<()> {
    //     if !self.input_path.exists() {
    //         anyhow::bail!("Input file not found: {}", self.input_path.display());
    //     }
    //     if self.kmer_size == 0 { 
    //         anyhow::bail!("k-mer size must be greater than 0.");
    //     }
    //     // Add more validation rules
    //     Ok(())
    // }
} 