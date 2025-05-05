use anyhow::{Context, Result};
use bio::io::fasta;
use std::fs::File;
use std::io::BufReader;
use std::path::Path;

/// Trait for reading FASTA sequence data in chunks.
/// Implementations allow memory-efficient access to large sequences.
pub trait FastaReader {
    /// Get the total length of a sequence.
    fn get_sequence_length(&self, sequence_id: usize) -> Result<usize>;
    
    /// Read a range of bases from a sequence.
    fn read_sequence_range(&mut self, sequence_id: usize, start: usize, length: usize) -> Result<Vec<u8>>;
}

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

/// An iterator that yields overlapping chunks of a FASTA sequence.
/// Useful for processing large genomes without loading the entire sequence into memory.
pub struct ChunkedSequenceIterator<'a> {
    reader: &'a mut dyn FastaReader,
    chunk_size: usize,
    overlap_size: usize,
    sequence_id: usize,
    pos: usize,
    total_length: Option<usize>,
    exhausted: bool,
}

impl<'a> ChunkedSequenceIterator<'a> {
    pub fn new(reader: &'a mut dyn FastaReader, chunk_size: usize, overlap_size: usize, sequence_id: usize) -> Self {
        // Ensure overlap is smaller than chunk size
        let overlap = overlap_size.min(chunk_size / 2);
        
        ChunkedSequenceIterator {
            reader,
            chunk_size,
            overlap_size: overlap,
            sequence_id,
            pos: 0,
            total_length: None,
            exhausted: false,
        }
    }
}

/// A chunk of sequence with metadata about its position in the original sequence
pub struct SequenceChunk {
    pub data: Vec<u8>,
    pub start_pos: usize,
    pub end_pos: usize,
    pub is_last: bool,
}

impl<'a> Iterator for ChunkedSequenceIterator<'a> {
    type Item = SequenceChunk;
    
    fn next(&mut self) -> Option<Self::Item> {
        if self.exhausted {
            return None;
        }
        
        // Get the total sequence length if not already known
        if self.total_length.is_none() {
            self.total_length = Some(self.reader.get_sequence_length(self.sequence_id).unwrap_or(0));
        }
        
        let total_length = self.total_length.unwrap();
        
        // If we've reached the end, stop
        if self.pos >= total_length {
            return None;
        }
        
        // Calculate how many bases to read in this chunk
        let to_read = self.chunk_size.min(total_length - self.pos);
        
        // Read the chunk
        let data = self.reader.read_sequence_range(self.sequence_id, self.pos, to_read)
            .unwrap_or_else(|_| Vec::new());
        
        if data.is_empty() {
            self.exhausted = true;
            return None;
        }
        
        let start_pos = self.pos;
        let end_pos = self.pos + data.len();
        let is_last = end_pos >= total_length;
        
        // Advance position for next chunk, considering overlap
        if !is_last {
            self.pos += to_read - self.overlap_size;
        } else {
            self.pos = total_length; // Ensure we stop after this chunk
        }
        
        Some(SequenceChunk {
            data,
            start_pos,
            end_pos,
            is_last,
        })
    }
}

/// Implementation of FastaReader for in-memory sequences.
pub struct InMemoryFastaReader {
    sequences: Vec<(String, Vec<u8>)>,
}

impl InMemoryFastaReader {
    pub fn new(sequences: Vec<(String, Vec<u8>)>) -> Self {
        Self { sequences }
    }
    
    pub fn from_fasta(fasta_path: &Path, skip_ns: bool) -> Result<Self> {
        let sequences = read_fasta_sequences(fasta_path, skip_ns)?;
        Ok(Self { sequences })
    }
    
    pub fn get_sequences(&self) -> &Vec<(String, Vec<u8>)> {
        &self.sequences
    }
    
    pub fn get_sequence_count(&self) -> usize {
        self.sequences.len()
    }
    
    pub fn get_sequence_id(&self, index: usize) -> Option<&str> {
        self.sequences.get(index).map(|(id, _)| id.as_str())
    }
}

impl FastaReader for InMemoryFastaReader {
    fn get_sequence_length(&self, sequence_id: usize) -> Result<usize> {
        self.sequences.get(sequence_id)
            .map(|(_, seq)| seq.len())
            .ok_or_else(|| anyhow::anyhow!("Sequence ID {} out of range", sequence_id))
    }
    
    fn read_sequence_range(&mut self, sequence_id: usize, start: usize, length: usize) -> Result<Vec<u8>> {
        if let Some((_, seq)) = self.sequences.get(sequence_id) {
            if start >= seq.len() {
                return Ok(Vec::new());
            }
            
            let end = (start + length).min(seq.len());
            Ok(seq[start..end].to_vec())
        } else {
            Err(anyhow::anyhow!("Sequence ID {} out of range", sequence_id))
        }
    }
} 