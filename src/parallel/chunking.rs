use anyhow::Result;
use crate::encode::dna_2bit::EncodedBase;
use crate::utils::progress::ProgressTracker;
use std::sync::{Arc, Mutex};
use rayon::prelude::*;
use log::debug;
use std::hash::{Hash, Hasher};
use std::collections::hash_map::DefaultHasher;

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
    /// Enable adaptive chunk sizing
    pub adaptive_chunking: bool,
    /// Maximum memory usage per chunk (bytes)
    pub max_memory_per_chunk: Option<usize>,
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
            adaptive_chunking: false,
            max_memory_per_chunk: None,
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

/// Analyzes sequence entropy to determine appropriate chunk sizes
fn analyze_sequence_entropy(sequence: &[EncodedBase], sample_size: usize) -> f64 {
    // Take samples from the sequence to analyze entropy
    let mut samples = Vec::new();
    let step = sequence.len() / sample_size.min(sequence.len()).max(1);
    
    for i in (0..sequence.len()).step_by(step) {
        if samples.len() >= sample_size {
            break;
        }
        samples.push(sequence[i]);
    }
    
    // Calculate Shannon entropy for k-mers in sample
    let k = 3; // Use 3-mers for entropy calculation
    if samples.len() < k {
        return 0.5; // Not enough data, return medium entropy
    }
    
    let mut kmers = std::collections::HashMap::new();
    for i in 0..=(samples.len() - k) {
        let mut hasher = DefaultHasher::new();
        for j in 0..k {
            samples[i + j].hash(&mut hasher);
        }
        let kmer_hash = hasher.finish();
        *kmers.entry(kmer_hash).or_insert(0) += 1;
    }
    
    // Calculate Shannon entropy
    let total_kmers = (samples.len() - k + 1) as f64;
    let mut entropy = 0.0;
    
    for &count in kmers.values() {
        let probability = count as f64 / total_kmers;
        entropy -= probability * probability.log2();
    }
    
    // Normalize to [0,1]
    entropy / 2.0
}

/// Determine optimal chunk size based on sequence properties and available memory
fn calculate_adaptive_chunk_size(
    sequence: &[EncodedBase],
    config: &ChunkingConfig,
    available_memory: usize,
) -> usize {
    // Analyze sequence entropy (higher entropy = more complex sequence = smaller chunks)
    let entropy = analyze_sequence_entropy(sequence, 1000);
    
    // Base chunk size on available memory, capped by config's max_memory_per_chunk if specified
    let memory_based_size = if let Some(max_per_chunk) = config.max_memory_per_chunk {
        (available_memory / config.num_threads).min(max_per_chunk)
    } else {
        available_memory / config.num_threads
    };
    
    // Convert from memory bytes to base count (each EncodedBase is typically 1 byte)
    // Apply entropy factor: lower entropy = larger chunks (more repetitive = more compressible)
    let entropy_factor = 0.5 + (1.0 - entropy) * 0.5; // Range [0.5, 1.0]
    let adaptive_size = (memory_based_size as f64 * entropy_factor) as usize;
    
    // Apply reasonable bounds
    let min_size = 50_000;              // Minimum 50KB
    let max_size = 100_000_000;         // Maximum 100MB
    
    let bounded_size = adaptive_size.clamp(min_size, max_size);
    
    debug!(
        "Adaptive chunk sizing: entropy={:.2}, memory={}, factor={:.2}, size={}",
        entropy, memory_based_size, entropy_factor, bounded_size
    );
    
    bounded_size
}

/// Split a sequence into chunks for parallel processing with adaptive sizing
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
    
    // Determine the chunk size to use
    let chunk_size = if config.adaptive_chunking {
        // Get available memory (fallback to a reasonable default if sysinfo fails)
        let mut sys = sysinfo::System::new_all();
        sys.refresh_memory();
        let available_memory = sys.available_memory() as usize;
        
        calculate_adaptive_chunk_size(sequence, config, available_memory)
    } else {
        config.chunk_size
    };
    
    debug!("Using chunk size: {} for sequence of length {}", chunk_size, seq_len);
    
    // Calculate how many chunks we'll need
    let effective_chunk_size = chunk_size.saturating_sub(config.overlap_size);
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
            (start + chunk_size).min(seq_len)
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
    process_fn: F
) -> Result<Vec<T>>
where
    F: Fn(&Chunk) -> Result<T> + Send + Sync,
    T: Send,
{
    let chunks_count = chunks.len();
    if chunks_count == 0 {
        return Ok(Vec::new());
    }
    
    let progress = if config.show_progress {
        Some(Arc::new(Mutex::new(ProgressTracker::new(
            chunks_count,
            "Processing chunks",
        ))))
    } else {
        None
    };
    
    let results: Vec<Result<T>> = chunks
        .into_par_iter()
        .with_max_len(1) // Each chunk gets its own task
        .map(|chunk| {
            // Process this chunk
            let result = process_fn(&chunk);
            
            // Update progress if tracking
            if let Some(ref progress) = progress {
                let mut progress = progress.lock().unwrap();
                progress.increment();
            }
            
            result
        })
        .collect();
    
    // Finalize progress
    if let Some(progress) = progress {
        let mut progress = progress.lock().unwrap();
        progress.finish();
    }
    
    // Collect and transform results
    let mut processed_results = Vec::with_capacity(chunks_count);
    for (i, result) in results.into_iter().enumerate() {
        match result {
            Ok(value) => processed_results.push(value),
            Err(e) => return Err(anyhow::anyhow!("Error processing chunk {}: {}", i, e)),
        }
    }
    
    Ok(processed_results)
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
    
    fn create_test_sequence(length: usize) -> Vec<EncodedBase> {
        let mut seq = Vec::with_capacity(length);
        for i in 0..length {
            // Create a simple pattern (0,1,2,3,0,1,2,3,...)
            seq.push(EncodedBase((i % 4) as u8));
        }
        seq
    }
    
    #[test]
    fn test_split_into_chunks() {
        let seq = create_test_sequence(1000);
        
        let config = ChunkingConfig {
            chunk_size: 300,
            overlap_size: 50,
            show_progress: false,
            adaptive_chunking: false, 
            max_memory_per_chunk: None,
            ..Default::default()
        };
        
        let chunks = split_into_chunks(&seq, &config);
        
        assert_eq!(chunks.len(), 4);
        assert_eq!(chunks[0].start, 0);
        assert!(chunks[0].is_first);
        assert!(!chunks[0].is_last);
        
        assert_eq!(chunks[3].end, 1000);
        assert!(!chunks[3].is_first);
        assert!(chunks[3].is_last);
        
        // Verify chunk data length
        assert_eq!(chunks[0].data.len(), 300);
        
        // Check overlap
        assert_eq!(chunks[0].end, 300);
        assert_eq!(chunks[1].start, 250); // 300 - 50
    }
    
    #[test]
    fn test_process_chunks_parallel() {
        let seq = create_test_sequence(1000);
        
        let config = ChunkingConfig {
            chunk_size: 300,
            overlap_size: 50,
            show_progress: false,
            adaptive_chunking: false,
            max_memory_per_chunk: None,
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
        let seq1 = vec![EncodedBase(0), EncodedBase(1), EncodedBase(2), EncodedBase(3)];
        let seq2 = vec![EncodedBase(2), EncodedBase(3), EncodedBase(0), EncodedBase(1)];
        
        let overlap = find_best_overlap(&seq1, &seq2, 4);
        assert_eq!(overlap, 2); // "23" matches
        
        let seq3 = vec![EncodedBase(9), EncodedBase(8), EncodedBase(7)];
        let overlap2 = find_best_overlap(&seq1, &seq3, 4);
        assert_eq!(overlap2, 0); // No match
    }
    
    #[test]
    fn test_adaptive_chunk_sizing() {
        let seq = create_test_sequence(10000);
        
        let mut config = ChunkingConfig {
            chunk_size: 1000,
            adaptive_chunking: true,
            max_memory_per_chunk: None,
            ..Default::default()
        };
        
        let available_memory = 10_000_000; // 10MB for testing
        let adaptive_size = calculate_adaptive_chunk_size(&seq, &config, available_memory);
        
        // Should return a reasonable size based on the pattern
        assert!(adaptive_size >= 50_000);
        
        // Test with config max memory per chunk
        config.max_memory_per_chunk = Some(1_000_000); // 1MB max
        let bounded_size = calculate_adaptive_chunk_size(&seq, &config, available_memory);
        
        // Should be bounded by the config
        assert!(bounded_size <= 1_000_000);
    }
} 