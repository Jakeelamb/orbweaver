//! SIMD-optimized DNA encoding and processing operations.
//!
//! This module provides high-performance vectorized operations for:
//! - DNA base encoding (ASCII -> 2-bit)
//! - DNA base decoding (2-bit -> ASCII)
//! - Bulk digram counting
//! - Fingerprint computation
//!
//! Uses lookup tables and batch processing for efficient auto-vectorization.

use super::dna_2bit::EncodedBase;

/// Lookup table for ASCII -> 2-bit encoding.
/// Index by ASCII value (0-255), returns encoded value or 0xFF for invalid.
/// This approach allows for efficient batch processing and auto-vectorization.
#[rustfmt::skip]
static ASCII_TO_2BIT: [u8; 256] = {
    let mut table = [0xFFu8; 256];
    table[b'A' as usize] = 0b00;
    table[b'a' as usize] = 0b00;
    table[b'C' as usize] = 0b01;
    table[b'c' as usize] = 0b01;
    table[b'G' as usize] = 0b10;
    table[b'g' as usize] = 0b10;
    table[b'T' as usize] = 0b11;
    table[b't' as usize] = 0b11;
    table
};

/// Lookup table for 2-bit -> ASCII decoding.
static TWOBIT_TO_ASCII: [u8; 4] = [b'A', b'C', b'G', b'T'];

/// SIMD-friendly batch encoding of DNA bases.
/// Uses lookup table for efficient vectorization.
/// Returns (encoded_bases, count of valid bases).
///
/// This is significantly faster than the iterator-based approach for large sequences
/// because it avoids branching and enables auto-vectorization.
#[inline]
pub fn encode_dna_batch(input: &[u8], output: &mut Vec<EncodedBase>) {
    output.clear();
    output.reserve(input.len());

    // Process in chunks for better cache utilization and vectorization
    const CHUNK_SIZE: usize = 64;
    let mut temp = [0u8; CHUNK_SIZE];

    for chunk in input.chunks(CHUNK_SIZE) {
        let chunk_len = chunk.len();

        // Lookup all values in the chunk
        for (i, &byte) in chunk.iter().enumerate() {
            temp[i] = ASCII_TO_2BIT[byte as usize];
        }

        // Filter valid bases (not 0xFF) and push to output
        for i in 0..chunk_len {
            if temp[i] != 0xFF {
                output.push(EncodedBase(temp[i]));
            }
        }
    }
}

/// Encode DNA directly to a pre-allocated buffer.
/// Returns the number of valid bases encoded.
/// This version is useful when you know the output size won't exceed input size.
#[inline]
pub fn encode_dna_into(input: &[u8], output: &mut [EncodedBase]) -> usize {
    let mut out_idx = 0;

    for &byte in input {
        let encoded = ASCII_TO_2BIT[byte as usize];
        if encoded != 0xFF && out_idx < output.len() {
            output[out_idx] = EncodedBase(encoded);
            out_idx += 1;
        }
    }

    out_idx
}

/// Count valid DNA bases in a sequence without allocating.
/// Useful for pre-allocating output buffers.
#[inline]
pub fn count_valid_bases(input: &[u8]) -> usize {
    input
        .iter()
        .filter(|&&b| ASCII_TO_2BIT[b as usize] != 0xFF)
        .count()
}

/// Batch decode 2-bit encoded bases to ASCII.
#[inline]
pub fn decode_dna_batch(input: &[EncodedBase], output: &mut Vec<u8>) {
    output.clear();
    output.reserve(input.len());

    for &base in input {
        if (base.0 as usize) < 4 {
            output.push(TWOBIT_TO_ASCII[base.0 as usize]);
        }
    }
}

/// Compute digram IDs for an entire sequence at once.
/// Each digram ID is computed as: (base1 << 2) | base2, giving values 0-15.
/// Returns a vector of (position, digram_id) pairs.
///
/// This is optimized for counting digrams in bulk.
#[inline]
pub fn compute_digram_ids(sequence: &[EncodedBase]) -> Vec<(usize, u8)> {
    if sequence.len() < 2 {
        return Vec::new();
    }

    let mut result = Vec::with_capacity(sequence.len() - 1);

    for i in 0..sequence.len() - 1 {
        let digram_id = (sequence[i].0 << 2) | sequence[i + 1].0;
        result.push((i, digram_id));
    }

    result
}

/// Count occurrences of each possible digram (0-15 for 2-bit encoded DNA).
/// Returns an array of 16 counts.
#[inline]
pub fn count_digrams(sequence: &[EncodedBase]) -> [usize; 16] {
    let mut counts = [0usize; 16];

    if sequence.len() < 2 {
        return counts;
    }

    // Process sequence to count digrams
    for i in 0..sequence.len() - 1 {
        let digram_id = ((sequence[i].0 as usize) << 2) | (sequence[i + 1].0 as usize);
        counts[digram_id] += 1;
    }

    counts
}

/// Find positions of a specific digram in the sequence.
/// digram_id is (base1 << 2) | base2.
#[inline]
pub fn find_digram_positions(sequence: &[EncodedBase], digram_id: u8) -> Vec<usize> {
    if sequence.len() < 2 {
        return Vec::new();
    }

    let mut positions = Vec::new();

    for i in 0..sequence.len() - 1 {
        let current_digram = (sequence[i].0 << 2) | sequence[i + 1].0;
        if current_digram == digram_id {
            positions.push(i);
        }
    }

    positions
}

/// Batch compute digram counts for multiple chunks in parallel.
/// Uses Rayon for parallel processing.
#[cfg(feature = "parallel")]
pub fn count_digrams_parallel(sequences: &[&[EncodedBase]]) -> [usize; 16] {
    use rayon::prelude::*;

    sequences
        .par_iter()
        .map(|seq| count_digrams(seq))
        .reduce(
            || [0usize; 16],
            |mut acc, counts| {
                for i in 0..16 {
                    acc[i] += counts[i];
                }
                acc
            },
        )
}

/// Parallel digram counting for a single large sequence by chunking.
pub fn count_digrams_parallel_chunked(sequence: &[EncodedBase], chunk_size: usize) -> [usize; 16] {
    use rayon::prelude::*;

    if sequence.len() < 2 {
        return [0usize; 16];
    }

    // Split into overlapping chunks to handle boundary digrams
    let chunks: Vec<&[EncodedBase]> = sequence
        .windows(chunk_size.max(2))
        .step_by(chunk_size.saturating_sub(1).max(1))
        .collect();

    if chunks.is_empty() {
        return count_digrams(sequence);
    }

    // Count in parallel
    let partial_counts: Vec<[usize; 16]> = chunks
        .par_iter()
        .map(|chunk| count_digrams(chunk))
        .collect();

    // Merge counts
    let mut total = [0usize; 16];
    for counts in partial_counts {
        for i in 0..16 {
            total[i] += counts[i];
        }
    }

    total
}

/// Fast fingerprint computation using FxHash.
/// This is used by the LCG algorithm for content-based rule identification.
#[inline]
pub fn fingerprint_sequence(sequence: &[EncodedBase]) -> u64 {
    use std::hash::{Hash, Hasher};

    let mut hasher = rustc_hash::FxHasher::default();
    for base in sequence {
        base.0.hash(&mut hasher);
    }
    hasher.finish()
}

/// Compute fingerprints for all windows of a given size.
/// Uses a rolling approach for efficiency.
pub fn compute_window_fingerprints(sequence: &[EncodedBase], window_size: usize) -> Vec<u64> {
    if sequence.len() < window_size {
        return Vec::new();
    }

    let mut fingerprints = Vec::with_capacity(sequence.len() - window_size + 1);

    for window in sequence.windows(window_size) {
        fingerprints.push(fingerprint_sequence(window));
    }

    fingerprints
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_encode_dna_batch() {
        let input = b"ACGTACGT";
        let mut output = Vec::new();
        encode_dna_batch(input, &mut output);

        assert_eq!(output.len(), 8);
        assert_eq!(output[0], EncodedBase(0b00)); // A
        assert_eq!(output[1], EncodedBase(0b01)); // C
        assert_eq!(output[2], EncodedBase(0b10)); // G
        assert_eq!(output[3], EncodedBase(0b11)); // T
    }

    #[test]
    fn test_encode_skips_invalid() {
        let input = b"ACNGTXACGT";
        let mut output = Vec::new();
        encode_dna_batch(input, &mut output);

        // N and X should be skipped
        assert_eq!(output.len(), 8);
    }

    #[test]
    fn test_count_digrams() {
        let sequence = vec![
            EncodedBase(0), // A
            EncodedBase(1), // C
            EncodedBase(0), // A
            EncodedBase(1), // C
            EncodedBase(2), // G
        ];

        let counts = count_digrams(&sequence);

        // AC appears twice (positions 0 and 2)
        assert_eq!(counts[(0 << 2) | 1], 2); // AC

        // CA appears once (position 1)
        assert_eq!(counts[(1 << 2) | 0], 1); // CA

        // CG appears once (position 3)
        assert_eq!(counts[(1 << 2) | 2], 1); // CG
    }

    #[test]
    fn test_find_digram_positions() {
        let sequence = vec![
            EncodedBase(0), // A
            EncodedBase(1), // C
            EncodedBase(0), // A
            EncodedBase(1), // C
        ];

        let ac_id = (0 << 2) | 1; // AC
        let positions = find_digram_positions(&sequence, ac_id);

        assert_eq!(positions, vec![0, 2]);
    }

    #[test]
    fn test_fingerprint_deterministic() {
        let seq1 = vec![EncodedBase(0), EncodedBase(1), EncodedBase(2)];
        let seq2 = vec![EncodedBase(0), EncodedBase(1), EncodedBase(2)];
        let seq3 = vec![EncodedBase(0), EncodedBase(1), EncodedBase(3)];

        assert_eq!(fingerprint_sequence(&seq1), fingerprint_sequence(&seq2));
        assert_ne!(fingerprint_sequence(&seq1), fingerprint_sequence(&seq3));
    }

    #[test]
    fn test_decode_batch() {
        let input = vec![
            EncodedBase(0),
            EncodedBase(1),
            EncodedBase(2),
            EncodedBase(3),
        ];
        let mut output = Vec::new();
        decode_dna_batch(&input, &mut output);

        assert_eq!(output, b"ACGT");
    }
}
