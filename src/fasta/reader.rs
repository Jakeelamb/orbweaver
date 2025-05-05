use anyhow::{Context, Result};
use bio::io::fasta;
use std::fs::File;
use std::io::BufReader;
use std::path::Path;

/// Reads sequences from a FASTA file, optionally skipping 'N' bases.
/// 
/// Args:
///     fasta_path: Path to the input FASTA file.
///     skip_ns: If true, 'N' and 'n' bases are excluded from the output sequence.
///
/// Returns:
///     A vector of tuples, where each tuple contains the record ID (String) 
///     and the processed sequence data (Vec<u8>).
pub fn read_fasta_sequences(fasta_path: &Path, skip_ns: bool) -> Result<Vec<(String, Vec<u8>)>> {
    println!("Reading FASTA file: {} (skip_ns: {})", fasta_path.display(), skip_ns);

    let file = File::open(fasta_path)
        .with_context(|| format!("Failed to open FASTA file: {}", fasta_path.display()))?;
    let reader = fasta::Reader::new(BufReader::new(file));

    let mut sequences = Vec::new();
    let mut record_count = 0;
    let mut base_count = 0;

    for result in reader.records() {
        let record = result.with_context(|| format!("Failed to read FASTA record from {}", fasta_path.display()))?;
        record_count += 1;
        let record_id = record.id().to_string();
        let original_len = record.seq().len();
        
        let sequence_data: Vec<u8> = if skip_ns {
            record.seq().iter().cloned().filter(|&base| base != b'N' && base != b'n').collect()
        } else {
            record.seq().to_vec()
        };

        base_count += sequence_data.len();
        let skipped_count = original_len - sequence_data.len();

        if skipped_count > 0 {
            println!("  Record '{}': read {} bases (skipped {} 'N' bases).", record_id, sequence_data.len(), skipped_count);
        } else {
             println!("  Record '{}': read {} bases.", record_id, sequence_data.len());
        }

        sequences.push((record_id, sequence_data));
    }

    if record_count == 0 {
        anyhow::bail!("No records found in FASTA file: {}", fasta_path.display());
    }

    println!("Finished reading {} records, total {} bases (after potential N skipping).", record_count, base_count);
    Ok(sequences)
} 