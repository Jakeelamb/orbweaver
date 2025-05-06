use anyhow::{Result, Context};
use crate::encode::dna_2bit::EncodedBase;
use crate::utils::progress::ProgressTracker;
use std::sync::{Arc, Mutex};
use std::time::Duration;
use rayon::prelude::*;
use log::{debug, info};

/// Configuration for chunking a sequence
#[derive(Debug, Clone)]
pub struct ChunkingConfig {
    /// Size of each chunk
    pub chunk_size: usize,
    /// Amount of overlap between adjacent chunks
    pub overlap_size: usize,
    /// Minimum rule usage required to create a rule
    pub min_rule_usage: usize,
    /// Whether to consider reverse complements
    pub reverse_aware: bool,
    /// Number of threads to use
    pub num_threads: usize,
    /// Whether to show progress
    pub show_progress: bool,
}

impl Default for ChunkingConfig {
    fn default() -> Self {
        Self {
            chunk_size: 100_000,
            overlap_size: 1_000,
            min_rule_usage: 2,
            reverse_aware: true,
            num_threads: num_cpus::get(),
            show_progress: true,
        }
    }
}

/// A chunk of a larger sequence for parallel processing
#[derive(Debug, Clone)]
pub struct Chunk {
    /// Index of this chunk in the original sequence
    pub index: usize,
    /// The actual data for this chunk
    pub data: Vec<EncodedBase>,
    /// Start position in the original sequence
    pub start: usize,
    /// End position in the original sequence
    pub end: usize,
    /// Is this the first chunk?
    pub is_first: bool,
    /// Is this the last chunk?
    pub is_last: bool,
}

/// Split a sequence into chunks for parallel processing
pub fn split_into_chunks(
    sequence: &[EncodedBase], 
    config: &ChunkingConfig
) -> Vec<Chunk> {
    let mut chunks = Vec::new();
    let seq_len = sequence.len();
    
    if seq_len <= config.chunk_size {
        // If the sequence is small enough, just use one chunk
        chunks.push(Chunk {
            index: 0,
            data: sequence.to_vec(),
            start: 0,
            end: seq_len,
            is_first: true,
            is_last: true,
        });
        
        return chunks;
    }
    
    // Calculate how many chunks we'll need
    let effective_chunk_size = config.chunk_size.saturating_sub(config.overlap_size);
    let num_chunks = if effective_chunk_size > 0 {
        (seq_len + effective_chunk_size - 1) / effective_chunk_size
    } else {
        1
    };
    
    debug!("Splitting sequence of length {} into {} chunks", seq_len, num_chunks);
    chunks.reserve(num_chunks);
    
    for i in 0..num_chunks {
        let start = if i == 0 {
            0
        } else {
            i * effective_chunk_size
        };
        
        let end = if i == num_chunks - 1 {
            seq_len
        } else {
            (start + config.chunk_size).min(seq_len)
        };
        
        let data = sequence[start..end].to_vec();
        
        chunks.push(Chunk {
            index: i,
            data,
            start,
            end,
            is_first: i == 0,
            is_last: i == num_chunks - 1,
        });
    }
    
    debug!("Created {} chunks", chunks.len());
    chunks
}

/// Process chunks in parallel using the provided function
pub fn process_chunks_parallel<F, T>(
    chunks: Vec<Chunk>,
    config: &ChunkingConfig,
    process_fn: F,
) -> Result<Vec<T>>
where
    F: Fn(&Chunk) -> Result<T> + Send + Sync,
    T: Send,
{
    // Configure thread pool
    rayon::ThreadPoolBuilder::new()
        .num_threads(config.num_threads)
        .build_global()
        .context("Failed to build thread pool")?;
    
    let num_chunks = chunks.len();
    debug!("Processing {} chunks using {} threads", num_chunks, config.num_threads);
    
    let progress = if config.show_progress {
        let tracker = ProgressTracker::new(num_chunks, "Chunk processing")
            .with_update_interval(Duration::from_millis(500));
        Some(Arc::new(Mutex::new(tracker)))
    } else {
        None
    };
    
    let results: Result<Vec<_>> = chunks
        .into_par_iter()
        .map(|chunk| {
            let result = process_fn(&chunk);
            
            if let Some(tracker) = &progress {
                let mut tracker = tracker.lock().unwrap();
                if tracker.increment() {
                    info!("{}", tracker.status());
                }
            }
            
            result
        })
        .collect();
    
    // Finalize progress tracking
    if let Some(tracker) = progress {
        let mut tracker = tracker.lock().unwrap();
        tracker.finish();
        info!("{}", tracker.status());
    }
    
    results
}

/// Merge a processed chunk with its neighbors, resolving overlaps
pub fn merge_chunks<T: Clone>(chunks: &[T], merge_fn: impl Fn(&T, &T) -> T) -> Vec<T> {
    if chunks.is_empty() {
        return Vec::new();
    }
    
    if chunks.len() == 1 {
        return vec![chunks[0].clone()];
    }
    
    let mut result = Vec::with_capacity(chunks.len());
    
    // First chunk is added as-is
    result.push(chunks[0].clone());
    
    // Merge intermediate chunks
    for i in 1..chunks.len() {
        let merged = merge_fn(&result[i-1], &chunks[i]);
        result.push(merged);
    }
    
    result
}

/// Joins overlapping sequences, preferring the right sequence for overlaps
pub fn join_sequences(left: &[EncodedBase], right: &[EncodedBase], overlap: usize) -> Vec<EncodedBase> {
    if left.is_empty() {
        return right.to_vec();
    }
    
    if right.is_empty() {
        return left.to_vec();
    }
    
    let actual_overlap = overlap.min(left.len()).min(right.len());
    
    if actual_overlap == 0 {
        // No overlap, just concatenate
        let mut result = left.to_vec();
        result.extend_from_slice(right);
        return result;
    }
    
    // Take the left sequence except for the overlapping portion
    let mut result = left[0..left.len() - actual_overlap].to_vec();
    
    // Then add the entire right sequence
    result.extend_from_slice(right);
    
    result
}

/// Find the best overlap between two sequences
pub fn find_best_overlap(
    left: &[EncodedBase], 
    right: &[EncodedBase], 
    max_overlap: usize
) -> usize {
    let max_possible = left.len().min(right.len()).min(max_overlap);
    
    if max_possible == 0 {
        return 0;
    }
    
    // Try different overlap sizes starting from the largest
    for overlap in (1..=max_possible).rev() {
        let left_suffix = &left[left.len() - overlap..];
        let right_prefix = &right[..overlap];
        
        if left_suffix == right_prefix {
            return overlap;
        }
    }
    
    0 // No matching overlap found
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::encode::dna_2bit::EncodedBase;
    
    fn create_test_sequence(len: usize) -> Vec<EncodedBase> {
        (0..len)
            .map(|i| EncodedBase(i as u8 % 4))
            .collect()
    }
    
    #[test]
    fn test_split_into_chunks() {
        let seq = create_test_sequence(1000);
        
        let config = ChunkingConfig {
            chunk_size: 300,
            overlap_size: 50,
            ..Default::default()
        };
        
        let chunks = split_into_chunks(&seq, &config);
        
        // Expected number of chunks with chunk_size=300, overlap_size=50:
        // Effective chunk size = 300 - 50 = 250
        // Number of chunks = ceil(1000 / 250) = 4
        assert_eq!(chunks.len(), 4);
        
        // Check first chunk
        assert_eq!(chunks[0].index, 0);
        assert_eq!(chunks[0].start, 0);
        assert_eq!(chunks[0].end, 300);
        assert_eq!(chunks[0].data.len(), 300);
        assert!(chunks[0].is_first);
        assert!(!chunks[0].is_last);
        
        // Check middle chunk
        assert_eq!(chunks[1].index, 1);
        assert_eq!(chunks[1].start, 250);
        assert_eq!(chunks[1].end, 550);
        assert_eq!(chunks[1].data.len(), 300);
        assert!(!chunks[1].is_first);
        assert!(!chunks[1].is_last);
        
        // Check last chunk
        assert_eq!(chunks[3].index, 3);
        assert_eq!(chunks[3].is_last, true);
    }
    
    #[test]
    fn test_process_chunks_parallel() {
        let seq = create_test_sequence(1000);
        
        let config = ChunkingConfig {
            chunk_size: 300,
            overlap_size: 50,
            show_progress: false,
            ..Default::default()
        };
        
        let chunks = split_into_chunks(&seq, &config);
        
        // Process function that counts distinct values in a chunk
        let process_fn = |chunk: &Chunk| -> Result<usize> {
            let distinct_values = chunk.data
                .iter()
                .map(|base| base.0)
                .collect::<std::collections::HashSet<_>>()
                .len();
            
            Ok(distinct_values)
        };
        
        let results = process_chunks_parallel(chunks, &config, process_fn).unwrap();
        
        // Since we only have 4 possible values (0-3), each chunk should have at most 4 distinct values
        for result in &results {
            assert!(*result <= 4);
        }
        
        assert_eq!(results.len(), 4); // Same number of results as chunks
    }
    
    #[test]
    fn test_join_sequences() {
        // Test with full overlap
        let left = vec![EncodedBase(0), EncodedBase(1), EncodedBase(2)];
        let right = vec![EncodedBase(2), EncodedBase(3), EncodedBase(0)];
        
        let joined = join_sequences(&left, &right, 1);
        assert_eq!(joined.len(), 5);
        assert_eq!(joined, vec![EncodedBase(0), EncodedBase(1), EncodedBase(2), EncodedBase(3), EncodedBase(0)]);
        
        // Test with no overlap
        let joined2 = join_sequences(&left, &right, 0);
        assert_eq!(joined2.len(), 6);
        
        // Test with empty sequences
        let empty: Vec<EncodedBase> = Vec::new();
        assert_eq!(join_sequences(&empty, &right, 0), right);
        assert_eq!(join_sequences(&left, &empty, 0), left);
    }
    
    #[test]
    fn test_find_best_overlap() {
        // Example: ACGT and GTAC
        // Best overlap is GT (2 bases)
        let left = vec![EncodedBase(0), EncodedBase(1), EncodedBase(2), EncodedBase(3)]; // ACGT
        let right = vec![EncodedBase(2), EncodedBase(3), EncodedBase(0), EncodedBase(1)]; // GTAC
        
        let overlap = find_best_overlap(&left, &right, 3);
        assert_eq!(overlap, 2); // GT overlaps
        
        // No overlap example
        let left2 = vec![EncodedBase(0), EncodedBase(0)]; // AA
        let right2 = vec![EncodedBase(1), EncodedBase(1)]; // CC
        
        let overlap2 = find_best_overlap(&left2, &right2, 2);
        assert_eq!(overlap2, 0); // No overlap
    }
} 