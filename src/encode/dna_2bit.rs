use anyhow::{bail, Result};
use serde::{Serialize, Deserialize};

/// Represents a 2-bit encoded DNA base.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord, Serialize, Deserialize)]
pub struct EncodedBase(pub u8);

impl EncodedBase {
    /// Creates an EncodedBase from a standard DNA character (A, C, G, T, case-insensitive).
    /// Returns None for any other character.
    pub fn from_base(base: u8) -> Option<Self> {
        match base {
            b'A' | b'a' => Some(Self(0b00)),
            b'C' | b'c' => Some(Self(0b01)),
            b'G' | b'g' => Some(Self(0b10)),
            b'T' | b't' => Some(Self(0b11)),
            _ => None, // Non ACGT bases are ignored
        }
    }

    /// Converts the encoded base back to its character representation ('A', 'C', 'G', 'T').
    /// Returns 'N' for any invalid encoded value (should not happen with valid EncodedBase).
    pub fn to_char(&self) -> char {
        match self.0 {
            0b00 => 'A',
            0b01 => 'C',
            0b10 => 'G',
            0b11 => 'T',
            _ => 'N', // Fallback for unexpected values
        }
    }

    /// Returns the reverse complement of the encoded base.
    /// A (00) <-> T (11)
    /// C (01) <-> G (10)
    pub fn revcomp(&self) -> Self {
        // XOR with 11 flips the bits correctly for ACGT complements
        Self(self.0 ^ 0b11)
    }
}

/// Encodes a slice of DNA bytes into a vector of EncodedBase.
/// Skips any bytes that are not valid ACGT bases.
pub fn encode_dna(seq: &[u8]) -> Vec<EncodedBase> {
    seq.iter().filter_map(|&b| EncodedBase::from_base(b)).collect()
}

/// Maps a DNA base (A, C, G, T, case-insensitive) to its 2-bit representation.
/// Returns None for non-ACGT characters.
pub fn base_to_2bit(base: u8) -> Option<u8> {
    EncodedBase::from_base(base).map(|eb| eb.0)
}

/// Maps a 2-bit representation back to a DNA base character ('A', 'C', 'G', 'T').
/// Panics if the input is not a valid 2-bit value (0-3).
pub fn twobit_to_base(twobit: u8) -> u8 {
    match twobit {
        0b00 => b'A',
        0b01 => b'C',
        0b10 => b'G',
        0b11 => b'T',
        _ => panic!("Invalid 2-bit value: {}", twobit),
    }
}

/// Encodes a DNA sequence (standard u8 slice) into a 2-bit packed representation.
/// Each byte in the output vector holds 4 bases.
/// Non-ACGT characters result in an error.
/// The length of the original sequence is needed for decoding later.
pub fn encode_dna_packed(sequence: &[u8]) -> Result<Vec<u8>> {
    let mut encoded = Vec::with_capacity((sequence.len() + 3) / 4);
    let mut current_byte = 0u8;
    let mut bit_offset = 0;

    for base in sequence {
        match base_to_2bit(*base) {
            Some(twobit) => {
                current_byte |= twobit << (6 - bit_offset);
                bit_offset += 2;
                if bit_offset == 8 {
                    encoded.push(current_byte);
                    current_byte = 0;
                    bit_offset = 0;
                }
            }
            None => {
                bail!("Invalid DNA base '{}' found in sequence for 2-bit packed encoding.", *base as char);
            }
        }
    }

    // Push the last byte if it contains any data
    if bit_offset > 0 {
        encoded.push(current_byte);
    }

    Ok(encoded)
}

/// Decodes a 2-bit packed DNA sequence back into a standard u8 vector.
/// Requires the original sequence length to handle padding correctly.
pub fn decode_dna_packed(encoded: &[u8], original_length: usize) -> Vec<u8> {
    let mut decoded = Vec::with_capacity(original_length);
    let mut bases_decoded = 0;

    for &byte in encoded {
        for i in 0..4 {
            if bases_decoded >= original_length {
                break; // Stop if we've decoded the original length
            }
            let shift = 6 - (i * 2);
            let twobit = (byte >> shift) & 0b11; // Extract 2 bits
            decoded.push(twobit_to_base(twobit));
            bases_decoded += 1;
        }
    }

    decoded
}

/// Decodes a specific segment from a 2-bit packed DNA sequence.
/// This is useful for chunked processing of large genomes.
/// 
/// Args:
///     encoded: The 2-bit encoded sequence.
///     start_pos: The starting position in the original sequence.
///     length: The number of bases to decode.
///     total_length: The total length of the original sequence.
///
/// Returns:
///     The decoded segment as a Vec<u8>.
pub fn decode_dna_packed_segment(encoded: &[u8], start_pos: usize, length: usize, total_length: usize) -> Vec<u8> {
    if start_pos >= total_length {
        return Vec::new();
    }

    let actual_length = length.min(total_length - start_pos);
    let mut decoded = Vec::with_capacity(actual_length);
    
    // Calculate byte index and bit offset for start_pos
    let start_byte = start_pos / 4;
    let start_offset = (start_pos % 4) * 2;
    
    let mut bases_decoded = 0;
    let mut current_byte_idx = start_byte;
    let mut current_offset = start_offset;

    while bases_decoded < actual_length && current_byte_idx < encoded.len() {
        let byte = encoded[current_byte_idx];
        
        // Extract the 2-bit value at the current offset
        let shift = 6 - current_offset;
        let twobit = (byte >> shift) & 0b11;
        decoded.push(twobit_to_base(twobit));
        
        bases_decoded += 1;
        current_offset += 2;
        
        if current_offset >= 8 {
            current_byte_idx += 1;
            current_offset = 0;
        }
    }

    decoded
}

/// Computes the reverse complement of a DNA sequence.
/// Handles case-insensitivity (output is uppercase).
/// Non-ACGT characters are skipped/ignored.
pub fn reverse_complement(sequence: &[u8]) -> Vec<u8> {
    sequence
        .iter()
        .rev() // Iterate in reverse
        .filter_map(|&base| match base.to_ascii_uppercase() {
            b'A' => Some(b'T'),
            b'T' => Some(b'A'),
            b'C' => Some(b'G'),
            b'G' => Some(b'C'),
            _ => None, // Skip non-ACGT bases
        })
        .collect()
}


#[cfg(test)]
mod tests {
    use super::*; // Imports items from the parent module (encoder.rs)

    #[test]
    fn test_from_base() {
        assert_eq!(EncodedBase::from_base(b'A'), Some(EncodedBase(0b00)));
        assert_eq!(EncodedBase::from_base(b'c'), Some(EncodedBase(0b01)));
        assert_eq!(EncodedBase::from_base(b'G'), Some(EncodedBase(0b10)));
        assert_eq!(EncodedBase::from_base(b't'), Some(EncodedBase(0b11)));
        assert_eq!(EncodedBase::from_base(b'N'), None);
        assert_eq!(EncodedBase::from_base(b'x'), None);
    }

    #[test]
    fn test_to_char() {
        assert_eq!(EncodedBase(0b00).to_char(), 'A');
        assert_eq!(EncodedBase(0b01).to_char(), 'C');
        assert_eq!(EncodedBase(0b10).to_char(), 'G');
        assert_eq!(EncodedBase(0b11).to_char(), 'T');
        // Test invalid internal value
        assert_eq!(EncodedBase(0b100).to_char(), 'N');
    }

    #[test]
    fn test_revcomp() {
        assert_eq!(EncodedBase(0b00).revcomp(), EncodedBase(0b11)); // A -> T
        assert_eq!(EncodedBase(0b11).revcomp(), EncodedBase(0b00)); // T -> A
        assert_eq!(EncodedBase(0b01).revcomp(), EncodedBase(0b10)); // C -> G
        assert_eq!(EncodedBase(0b10).revcomp(), EncodedBase(0b01)); // G -> C
    }

    #[test]
    fn test_encode_dna_function() {
        let seq = b"ACGTNacgtX";
        let encoded = encode_dna(seq);
        let expected = vec![
            EncodedBase(0b00), // A
            EncodedBase(0b01), // C
            EncodedBase(0b10), // G
            EncodedBase(0b11), // T
            EncodedBase(0b00), // a
            EncodedBase(0b01), // c
            EncodedBase(0b10), // g
            EncodedBase(0b11), // t
        ];
        assert_eq!(encoded, expected);
    }

    #[test]
    fn test_encode_dna() {
        let seq1 = b"ACGT";
        let encoded1 = encode_dna(seq1);
        assert_eq!(encoded1, vec![EncodedBase(0), EncodedBase(1), EncodedBase(2), EncodedBase(3)]);

        let seq2 = b"acgTNA"; // Skips N
        let encoded2 = encode_dna(seq2);
        assert_eq!(encoded2, vec![EncodedBase(0), EncodedBase(1), EncodedBase(2), EncodedBase(3), EncodedBase(0)]);
        
        let seq3 = b"";
        let encoded3 = encode_dna(seq3);
        assert!(encoded3.is_empty());
    }

    #[test]
    fn test_twobit_to_base_conversion() {
        assert_eq!(twobit_to_base(0b00), b'A');
        assert_eq!(twobit_to_base(0b01), b'C');
        assert_eq!(twobit_to_base(0b10), b'G');
        assert_eq!(twobit_to_base(0b11), b'T');
    }

    #[test]
    #[should_panic(expected = "Invalid 2-bit value")]
    fn test_twobit_to_base_invalid() {
        twobit_to_base(0b100); // Should panic
    }

    #[test]
    fn test_encode_decode_dna_packed() {
        let original = b"ACGTACGTACGT";
        
        // Encode
        let packed = encode_dna_packed(original).unwrap();
        
        // Each byte holds 4 bases, so 12 bases need 3 bytes
        assert_eq!(packed.len(), 3);
        
        // First byte should be: 00 01 10 11 = A C G T = 00011011 = 27
        assert_eq!(packed[0], 0b00011011);
        
        // Second byte should be: 00 01 10 11 = A C G T = 00011011 = 27
        assert_eq!(packed[1], 0b00011011);
        
        // Third byte should be: 00 01 10 11 = A C G T = 00011011 = 27
        assert_eq!(packed[2], 0b00011011);
        
        // Decode
        let decoded = decode_dna_packed(&packed, original.len());
        
        // Compare
        assert_eq!(decoded, original);
    }

    #[test]
    fn test_encode_dna_packed_invalid_base() {
        let result = encode_dna_packed(b"ACGNT");
        assert!(result.is_err());
    }

    #[test]
    fn test_reverse_complement() {
        let seq = b"ACGT";
        let revcomp = reverse_complement(seq);
        assert_eq!(revcomp, b"ACGT");
        
        let seq2 = b"AACGTT";
        let revcomp2 = reverse_complement(seq2);
        assert_eq!(revcomp2, b"AACGTT");
        
        let seq3 = b"AACGNT";
        let revcomp3 = reverse_complement(seq3);
        assert_eq!(revcomp3, b"ACGTT"); // N is skipped
    }
} 