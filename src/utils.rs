use anyhow::{Context, Result};
use log::{debug, warn};
use std::fs::{self, create_dir_all, File};
use std::io::{self, BufRead, BufReader, Write};
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

/// Ensure that the specified directory exists, creating it if necessary
pub fn ensure_dir_exists(path: &Path) -> Result<()> {
    if !path.exists() {
        create_dir_all(path).context(format!("Failed to create directory: {}", path.display()))?;
    }
    Ok(())
}

/// Count the number of 'N' bases in a sequence
pub fn count_n_bases(sequence: &str) -> usize {
    sequence.chars().filter(|&c| c == 'N' || c == 'n').count()
}

/// Calculate the percentage of 'N' bases in a sequence
pub fn calculate_n_percentage(sequence: &str) -> f64 {
    let n_count = count_n_bases(sequence);
    let total_length = sequence.len();
    if total_length > 0 {
        (n_count as f64 / total_length as f64) * 100.0
    } else {
        0.0
    }
}

/// Run bqtools encode command via CLI
pub fn run_bqtools_encode(
    input_file: &Path,
    output_file: &Path,
    policy: &str,
) -> Result<Output> {
    debug!("Running bqtools encode on {}", input_file.display());
    
    let output = Command::new("bqtools")
        .arg("encode")
        .arg("-i")
        .arg(input_file)
        .arg("-o")
        .arg(output_file)
        .arg("-p")
        .arg(policy)
        .output()
        .context("Failed to execute bqtools encode command")?;
    
    if !output.status.success() {
        let stderr = String::from_utf8_lossy(&output.stderr);
        warn!("bqtools encode error: {}", stderr);
        return Err(anyhow::anyhow!("bqtools encode failed: {}", stderr));
    }
    
    Ok(output)
}

/// Run bqtools decode command via CLI
pub fn run_bqtools_decode(
    input_file: &Path,
    output_file: &Path,
) -> Result<Output> {
    debug!("Running bqtools decode on {}", input_file.display());
    
    let output = Command::new("bqtools")
        .arg("decode")
        .arg("-i")
        .arg(input_file)
        .arg("-o")
        .arg(output_file)
        .output()
        .context("Failed to execute bqtools decode command")?;
    
    if !output.status.success() {
        let stderr = String::from_utf8_lossy(&output.stderr);
        warn!("bqtools decode error: {}", stderr);
        return Err(anyhow::anyhow!("bqtools decode failed: {}", stderr));
    }
    
    Ok(output)
}

/// Compare the length of two files to ensure they match
pub fn validate_sequence_length(original_file: &Path, processed_file: &Path) -> Result<(usize, usize)> {
    let original_len = count_sequence_len(original_file)?;
    let processed_len = count_sequence_len(processed_file)?;
    
    debug!(
        "Length validation: Original={}, Processed={}",
        original_len, processed_len
    );
    
    Ok((original_len, processed_len))
}

/// Count the total sequence length in a FASTA file (ignoring headers)
fn count_sequence_len(file_path: &Path) -> Result<usize> {
    let file = File::open(file_path).context(format!("Failed to open file: {}", file_path.display()))?;
    let reader = BufReader::new(file);
    
    let mut total_len = 0;
    let mut in_header = false;
    
    for line in reader.lines() {
        let line = line?;
        if line.starts_with('>') {
            in_header = true;
        } else if !line.trim().is_empty() {
            if !in_header {
                total_len += line.trim().len();
            }
            in_header = false;
        }
    }
    
    Ok(total_len)
}

/// Generate a unique filename for a chromosome VBQ file
pub fn generate_chromosome_vbq_filename(
    output_dir: &Path,
    genome_id: &str,
    chromosome_id: &str,
) -> PathBuf {
    output_dir.join(format!("{}_{}.vbq", genome_id, chromosome_id))
}

/// Clean chromosome ID by removing problematic characters
pub fn clean_chromosome_id(raw_id: &str) -> String {
    raw_id
        .trim_start_matches('>')
        .split_whitespace()
        .next()
        .unwrap_or("unknown")
        .replace(|c: char| !c.is_alphanumeric() && c != '_', "_")
} 