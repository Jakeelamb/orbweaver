use anyhow::{Context, Result};
use bio::io::fasta;
use std::fs::File;
use std::io::BufReader;
use std::path::Path;
use memmap2::Mmap;
use rayon::prelude::*;
use std::sync::Arc;

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

// Basic structure to hold FASTA record information parsed from memory-mapped file
#[derive(Clone)]
struct MmapFastaRecord {
    id: String,
    header_offset: usize,
    sequence_offset: usize,
    sequence_length: usize,
}

/// A memory-mapped FASTA reader for efficient, parallel access to large genomes
pub struct MemoryMappedFastaReader {
    // The memory-mapped file
    mmap: Arc<Mmap>,
    // Parsed records information
    records: Vec<MmapFastaRecord>,
    // Whether to skip N bases
    skip_ns: bool,
}

impl MemoryMappedFastaReader {
    /// Create a new memory-mapped FASTA reader
    pub fn new(path: &Path, skip_ns: bool) -> Result<Self> {
        println!("Memory-mapping FASTA file: {}", path.display());
        
        // Open the file and memory map it
        let file = File::open(path)
            .with_context(|| format!("Failed to open FASTA file: {}", path.display()))?;
        
        let mmap = unsafe { Mmap::map(&file) }
            .with_context(|| format!("Failed to memory-map file: {}", path.display()))?;
        
        let mmap = Arc::new(mmap);
        
        // Parse the FASTA header and record information
        let records = Self::parse_fasta_headers(&mmap)
            .with_context(|| format!("Failed to parse FASTA headers from file: {}", path.display()))?;
        
        if records.is_empty() {
            anyhow::bail!("No records found in FASTA file: {}", path.display());
        }
        
        println!("Memory-mapped {} FASTA records", records.len());
        
        Ok(Self {
            mmap,
            records,
            skip_ns,
        })
    }
    
    /// Process FASTA file in memory-mapped mode with parallel chunks
    pub fn process_in_parallel_chunks<F, T>(&self, chunk_size: usize, sequence_id: usize, processor: F) -> Result<Vec<T>>
    where
        F: Fn(&[u8]) -> T + Send + Sync,
        T: Send,
    {
        if sequence_id >= self.records.len() {
            anyhow::bail!("Sequence ID {} out of range (max: {})", sequence_id, self.records.len() - 1);
        }
        
        let record = &self.records[sequence_id];
        println!("Processing record '{}' in parallel chunks", record.id);
        
        // Calculate how many chunks we need
        let num_chunks = (record.sequence_length + chunk_size - 1) / chunk_size;
        println!("Splitting into {} chunks of size {}", num_chunks, chunk_size);
        
        // Process chunks in parallel
        let results: Vec<T> = (0..num_chunks)
            .into_par_iter()
            .map(|chunk_idx| {
                let start = chunk_idx * chunk_size;
                let end = std::cmp::min(start + chunk_size, record.sequence_length);
                
                let data = self.read_sequence_range_internal(sequence_id, start, end - start)
                    .expect("Failed to read sequence range");
                
                processor(&data)
            })
            .collect();
        
        Ok(results)
    }
    
    /// Parse FASTA headers and record information from memory-mapped file
    fn parse_fasta_headers(mmap: &Mmap) -> Result<Vec<MmapFastaRecord>> {
        let mut records = Vec::new();
        let mut pos = 0;
        let data = &mmap[..];
        
        while pos < data.len() {
            // Find the start of a record (>)
            if data[pos] == b'>' {
                let header_offset = pos;
                
                // Find the end of the header line
                let mut header_end = pos;
                while header_end < data.len() && data[header_end] != b'\n' {
                    header_end += 1;
                }
                
                if header_end >= data.len() {
                    break;
                }
                
                // Extract the ID from the header
                let id_end = if let Some(idx) = data[pos + 1..header_end].iter().position(|&b| b == b' ') {
                    pos + 1 + idx
                } else {
                    header_end
                };
                
                let id = String::from_utf8_lossy(&data[pos + 1..id_end]).to_string();
                
                // The sequence starts after the header line
                let sequence_offset = header_end + 1;
                
                // Find the end of the sequence (next '>' or end of file)
                let mut sequence_end = sequence_offset;
                while sequence_end < data.len() && data[sequence_end] != b'>' {
                    sequence_end += 1;
                }
                
                // Calculate the sequence length, excluding whitespace
                let mut seq_len = 0;
                for i in sequence_offset..sequence_end {
                    if !data[i].is_ascii_whitespace() {
                        seq_len += 1;
                    }
                }
                
                records.push(MmapFastaRecord {
                    id,
                    header_offset,
                    sequence_offset,
                    sequence_length: seq_len,
                });
                
                pos = sequence_end;
            } else {
                pos += 1;
            }
        }
        
        Ok(records)
    }
    
    /// Internal method to read a range of sequence data
    fn read_sequence_range_internal(&self, sequence_id: usize, start: usize, length: usize) -> Result<Vec<u8>> {
        if sequence_id >= self.records.len() {
            anyhow::bail!("Sequence ID {} out of range", sequence_id);
        }
        
        let record = &self.records[sequence_id];
        if start >= record.sequence_length {
            return Ok(Vec::new());
        }
        
        let end = std::cmp::min(start + length, record.sequence_length);
        let mut result = Vec::with_capacity(end - start);
        
        // Iterate through the file, skipping whitespace and collecting sequence data
        let mut seq_pos = 0;
        let mut data_pos = record.sequence_offset;
        
        while seq_pos < end && data_pos < self.mmap.len() {
            let byte = self.mmap[data_pos];
            
            if !byte.is_ascii_whitespace() {
                if seq_pos >= start {
                    // Skip 'N' bases if configured to do so
                    if !self.skip_ns || (byte != b'N' && byte != b'n') {
                        result.push(byte);
                    }
                }
                seq_pos += 1;
            }
            
            data_pos += 1;
            
            // Check for end of record
            if data_pos < self.mmap.len() && self.mmap[data_pos] == b'>' {
                break;
            }
        }
        
        Ok(result)
    }
}

impl FastaReader for MemoryMappedFastaReader {
    fn get_sequence_length(&self, sequence_id: usize) -> Result<usize> {
        if sequence_id >= self.records.len() {
            anyhow::bail!("Sequence ID {} out of range", sequence_id);
        }
        
        // Return the cached sequence length
        Ok(self.records[sequence_id].sequence_length)
    }
    
    fn read_sequence_range(&mut self, sequence_id: usize, start: usize, length: usize) -> Result<Vec<u8>> {
        self.read_sequence_range_internal(sequence_id, start, length)
    }
}

/// A stream interface for efficiently processing FASTA sequences
pub struct FastaStream {
    reader: Box<dyn FastaReader>,
    sequence_id: usize,
    chunk_size: usize,
    skip_ns: bool,
    position: usize,
    total_length: Option<usize>,
}

impl FastaStream {
    /// Create a new FASTA stream
    pub fn new(path: &Path) -> Result<Self> {
        // Try memory-mapped reader first
        let result = MemoryMappedFastaReader::new(path, false);
        
        let reader = if let Ok(mmap_reader) = result {
            Box::new(mmap_reader) as Box<dyn FastaReader>
        } else {
            // Fall back to in-memory reader
            println!("Falling back to in-memory FASTA reader");
            let sequences = read_fasta_sequences(path, false)?;
            Box::new(InMemoryFastaReader::new(sequences)) as Box<dyn FastaReader>
        };
        
        Ok(Self {
            reader,
            sequence_id: 0,
            chunk_size: 1024 * 1024, // 1MB default chunk size
            skip_ns: false,
            position: 0,
            total_length: None,
        })
    }
    
    /// Set whether to skip N bases
    pub fn skip_ns(mut self) -> Self {
        self.skip_ns = true;
        self
    }
    
    /// Set the chunk size for reading
    pub fn with_chunk_size(mut self, chunk_size: usize) -> Self {
        self.chunk_size = chunk_size;
        self
    }
    
    /// Set which sequence to read (by index)
    pub fn with_sequence_id(mut self, sequence_id: usize) -> Self {
        self.sequence_id = sequence_id;
        self
    }
}

impl Iterator for FastaStream {
    type Item = Vec<u8>;
    
    fn next(&mut self) -> Option<Self::Item> {
        // Get total length if not already known
        if self.total_length.is_none() {
            self.total_length = self.reader.get_sequence_length(self.sequence_id).ok();
            if self.total_length.is_none() {
                return None; // Error or empty sequence
            }
        }
        
        let total_length = self.total_length.unwrap();
        if self.position >= total_length {
            return None; // End of sequence
        }
        
        // Calculate length to read
        let to_read = std::cmp::min(self.chunk_size, total_length - self.position);
        
        // Read chunk
        match self.reader.read_sequence_range(self.sequence_id, self.position, to_read) {
            Ok(mut data) => {
                // Skip N bases if requested
                if self.skip_ns {
                    data.retain(|&b| b != b'N' && b != b'n');
                }
                
                // Update position
                self.position += to_read;
                
                if !data.is_empty() {
                    Some(data)
                } else {
                    self.next() // Skip empty chunks
                }
            }
            Err(_) => None, // Error reading
        }
    }
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
    pub record_id: String,
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
            record_id: String::new(), // Empty for now, could be populated from reader
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