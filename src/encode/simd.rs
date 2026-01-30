//! SIMD-optimized DNA encoding and processing operations.
//!
//! This module provides high-performance vectorized operations for:
//! - DNA base encoding (ASCII -> 2-bit)
//! - DNA base decoding (2-bit -> ASCII)
//! - Bulk digram counting
//! - Fingerprint computation
//!
//! Uses lookup tables and batch processing for efficient auto-vectorization.
//! Also provides explicit SIMD intrinsics for x86_64 (AVX2) and aarch64 (NEON)
//! when available.

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

// =============================================================================
// Explicit SIMD Intrinsics
// =============================================================================

/// Check if AVX2 is available at runtime
#[cfg(target_arch = "x86_64")]
pub fn has_avx2() -> bool {
    is_x86_feature_detected!("avx2")
}

#[cfg(not(target_arch = "x86_64"))]
pub fn has_avx2() -> bool {
    false
}

/// Check if NEON is available (always true on aarch64)
#[cfg(target_arch = "aarch64")]
pub fn has_neon() -> bool {
    true
}

#[cfg(not(target_arch = "aarch64"))]
pub fn has_neon() -> bool {
    false
}

/// AVX2-optimized DNA encoding (x86_64 only)
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx2")]
unsafe fn encode_dna_avx2(input: &[u8], output: &mut Vec<EncodedBase>) {
    use std::arch::x86_64::*;

    output.clear();
    output.reserve(input.len());

    // Create lookup vectors for ASCII to 2-bit conversion
    // We'll process in 32-byte chunks
    let chunks = input.chunks_exact(32);
    let remainder = chunks.remainder();

    // Character constants
    let a_upper = _mm256_set1_epi8(b'A' as i8);
    let c_upper = _mm256_set1_epi8(b'C' as i8);
    let g_upper = _mm256_set1_epi8(b'G' as i8);
    let t_upper = _mm256_set1_epi8(b'T' as i8);
    let a_lower = _mm256_set1_epi8(b'a' as i8);
    let c_lower = _mm256_set1_epi8(b'c' as i8);
    let g_lower = _mm256_set1_epi8(b'g' as i8);
    let t_lower = _mm256_set1_epi8(b't' as i8);

    // Encoding values
    let val_0 = _mm256_set1_epi8(0);
    let val_1 = _mm256_set1_epi8(1);
    let val_2 = _mm256_set1_epi8(2);
    let val_3 = _mm256_set1_epi8(3);
    let val_invalid = _mm256_set1_epi8(-1); // 0xFF

    for chunk in chunks {
        let data = _mm256_loadu_si256(chunk.as_ptr() as *const __m256i);

        // Compare with each valid character (upper and lower case)
        let is_a = _mm256_or_si256(
            _mm256_cmpeq_epi8(data, a_upper),
            _mm256_cmpeq_epi8(data, a_lower),
        );
        let is_c = _mm256_or_si256(
            _mm256_cmpeq_epi8(data, c_upper),
            _mm256_cmpeq_epi8(data, c_lower),
        );
        let is_g = _mm256_or_si256(
            _mm256_cmpeq_epi8(data, g_upper),
            _mm256_cmpeq_epi8(data, g_lower),
        );
        let is_t = _mm256_or_si256(
            _mm256_cmpeq_epi8(data, t_upper),
            _mm256_cmpeq_epi8(data, t_lower),
        );

        // Create result vector with encoded values
        let result = _mm256_blendv_epi8(
            _mm256_blendv_epi8(
                _mm256_blendv_epi8(
                    _mm256_blendv_epi8(val_invalid, val_0, is_a),
                    val_1, is_c,
                ),
                val_2, is_g,
            ),
            val_3, is_t,
        );

        // Extract and filter valid bases
        let mut result_arr = [0i8; 32];
        _mm256_storeu_si256(result_arr.as_mut_ptr() as *mut __m256i, result);

        for &val in &result_arr {
            if val != -1 {
                output.push(EncodedBase(val as u8));
            }
        }
    }

    // Handle remainder with scalar code
    for &byte in remainder {
        let encoded = ASCII_TO_2BIT[byte as usize];
        if encoded != 0xFF {
            output.push(EncodedBase(encoded));
        }
    }
}

/// NEON-optimized DNA encoding (aarch64 only)
#[cfg(target_arch = "aarch64")]
unsafe fn encode_dna_neon(input: &[u8], output: &mut Vec<EncodedBase>) {
    use std::arch::aarch64::*;

    output.clear();
    output.reserve(input.len());

    // Process in 16-byte chunks
    let chunks = input.chunks_exact(16);
    let remainder = chunks.remainder();

    // Character constants
    let a_upper = vdupq_n_u8(b'A');
    let c_upper = vdupq_n_u8(b'C');
    let g_upper = vdupq_n_u8(b'G');
    let t_upper = vdupq_n_u8(b'T');
    let a_lower = vdupq_n_u8(b'a');
    let c_lower = vdupq_n_u8(b'c');
    let g_lower = vdupq_n_u8(b'g');
    let t_lower = vdupq_n_u8(b't');

    // Encoding values
    let val_0 = vdupq_n_u8(0);
    let val_1 = vdupq_n_u8(1);
    let val_2 = vdupq_n_u8(2);
    let val_3 = vdupq_n_u8(3);
    let val_invalid = vdupq_n_u8(0xFF);

    for chunk in chunks {
        let data = vld1q_u8(chunk.as_ptr());

        // Compare with each valid character
        let is_a = vorrq_u8(vceqq_u8(data, a_upper), vceqq_u8(data, a_lower));
        let is_c = vorrq_u8(vceqq_u8(data, c_upper), vceqq_u8(data, c_lower));
        let is_g = vorrq_u8(vceqq_u8(data, g_upper), vceqq_u8(data, g_lower));
        let is_t = vorrq_u8(vceqq_u8(data, t_upper), vceqq_u8(data, t_lower));

        // Build result using bitwise select
        let mut result = val_invalid;
        result = vbslq_u8(is_a, val_0, result);
        result = vbslq_u8(is_c, val_1, result);
        result = vbslq_u8(is_g, val_2, result);
        result = vbslq_u8(is_t, val_3, result);

        // Extract and filter valid bases
        let mut result_arr = [0u8; 16];
        vst1q_u8(result_arr.as_mut_ptr(), result);

        for &val in &result_arr {
            if val != 0xFF {
                output.push(EncodedBase(val));
            }
        }
    }

    // Handle remainder with scalar code
    for &byte in remainder {
        let encoded = ASCII_TO_2BIT[byte as usize];
        if encoded != 0xFF {
            output.push(EncodedBase(encoded));
        }
    }
}

/// Best-effort SIMD DNA encoding that uses the optimal implementation for the current platform.
///
/// Falls back to scalar code if SIMD is not available.
pub fn encode_dna_simd(input: &[u8], output: &mut Vec<EncodedBase>) {
    #[cfg(target_arch = "x86_64")]
    {
        if has_avx2() {
            // SAFETY: We checked that AVX2 is available
            unsafe {
                encode_dna_avx2(input, output);
                return;
            }
        }
    }

    #[cfg(target_arch = "aarch64")]
    {
        // NEON is always available on aarch64
        unsafe {
            encode_dna_neon(input, output);
            return;
        }
    }

    // Fallback to batch encoding
    encode_dna_batch(input, output);
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

    #[test]
    fn test_encode_dna_simd() {
        let input = b"ACGTACGT";
        let mut output = Vec::new();
        super::encode_dna_simd(input, &mut output);

        assert_eq!(output.len(), 8);
        assert_eq!(output[0], EncodedBase(0b00)); // A
        assert_eq!(output[1], EncodedBase(0b01)); // C
        assert_eq!(output[2], EncodedBase(0b10)); // G
        assert_eq!(output[3], EncodedBase(0b11)); // T
    }

    #[test]
    fn test_encode_dna_simd_with_invalid() {
        let input = b"ACNGTXACGT";
        let mut output = Vec::new();
        super::encode_dna_simd(input, &mut output);

        // N and X should be skipped
        assert_eq!(output.len(), 8);
    }

    #[test]
    fn test_encode_dna_simd_long_sequence() {
        // Test with a sequence longer than SIMD register size
        let input: Vec<u8> = (0..100).map(|i| match i % 4 {
            0 => b'A',
            1 => b'C',
            2 => b'G',
            _ => b'T',
        }).collect();

        let mut output = Vec::new();
        super::encode_dna_simd(&input, &mut output);

        assert_eq!(output.len(), 100);
        for (i, base) in output.iter().enumerate() {
            assert_eq!(base.0, (i % 4) as u8);
        }
    }

    #[test]
    fn test_simd_detection() {
        // Just verify the detection functions compile and run
        let _ = super::has_avx2();
        let _ = super::has_neon();
    }
}
