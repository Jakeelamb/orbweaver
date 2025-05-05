use anyhow::{bail, Result};

/// Maps a DNA base (A, C, G, T, case-insensitive) to its 2-bit representation.
/// Returns None for non-ACGT characters.
pub fn base_to_2bit(base: u8) -> Option<u8> {
    match base.to_ascii_uppercase() {
        b'A' => Some(0b00),
        b'C' => Some(0b01),
        b'G' => Some(0b10),
        b'T' => Some(0b11),
        _ => None,
    }
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
pub fn encode_dna_2bit(sequence: &[u8]) -> Result<Vec<u8>> {
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
                bail!("Invalid DNA base '{}' found in sequence for 2-bit encoding.", *base as char);
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
pub fn decode_dna_2bit(encoded: &[u8], original_length: usize) -> Vec<u8> {
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
    fn test_base_to_2bit_conversion() {
        assert_eq!(base_to_2bit(b'A'), Some(0b00));
        assert_eq!(base_to_2bit(b'c'), Some(0b01));
        assert_eq!(base_to_2bit(b'G'), Some(0b10));
        assert_eq!(base_to_2bit(b't'), Some(0b11));
        assert_eq!(base_to_2bit(b'N'), None);
        assert_eq!(base_to_2bit(b'x'), None);
    }

    #[test]
    fn test_twobit_to_base_conversion() {
        assert_eq!(twobit_to_base(0b00), b'A');
        assert_eq!(twobit_to_base(0b01), b'C');
        assert_eq!(twobit_to_base(0b10), b'G');
        assert_eq!(twobit_to_base(0b11), b'T');
    }

    #[test]
    #[should_panic]
    fn test_twobit_to_base_invalid() {
        twobit_to_base(4);
    }

    #[test]
    fn test_encode_decode_dna_2bit() {
        let seq1 = b"ACGT"; // 1 byte
        let encoded1 = encode_dna_2bit(seq1).unwrap();
        assert_eq!(encoded1.len(), 1);
        // A=00, C=01, G=10, T=11 -> 00011011 = 0x1B
        assert_eq!(encoded1[0], 0b00011011);
        assert_eq!(decode_dna_2bit(&encoded1, seq1.len()), seq1);

        let seq2 = b"A"; // 1 byte, padded
        let encoded2 = encode_dna_2bit(seq2).unwrap();
        assert_eq!(encoded2.len(), 1);
        // A=00 -> 00000000 (bits 6-7 are 00, rest are 0)
        assert_eq!(encoded2[0], 0b00000000);
        assert_eq!(decode_dna_2bit(&encoded2, seq2.len()), seq2);

        let seq3 = b"ACG"; // 1 byte, padded
        let encoded3 = encode_dna_2bit(seq3).unwrap();
        assert_eq!(encoded3.len(), 1);
        // A=00, C=01, G=10 -> 00011000 = 0x18
        assert_eq!(encoded3[0], 0b00011000);
        assert_eq!(decode_dna_2bit(&encoded3, seq3.len()), seq3);

        let seq4 = b"ACGTACGT"; // 2 bytes
        let encoded4 = encode_dna_2bit(seq4).unwrap();
        assert_eq!(encoded4.len(), 2);
        assert_eq!(encoded4[0], 0b00011011);
        assert_eq!(encoded4[1], 0b00011011);
        assert_eq!(decode_dna_2bit(&encoded4, seq4.len()), seq4);
        
        let seq5 = b"ACGTACGTA"; // 3 bytes, padded
        let encoded5 = encode_dna_2bit(seq5).unwrap();
        assert_eq!(encoded5.len(), 3);
        assert_eq!(encoded5[0], 0b00011011);
        assert_eq!(encoded5[1], 0b00011011);
        assert_eq!(encoded5[2], 0b00000000); // A -> 00xxxxxx
        assert_eq!(decode_dna_2bit(&encoded5, seq5.len()), seq5);
    }

    #[test]
    fn test_encode_dna_2bit_invalid_base() {
        let seq = b"ACGTN";
        assert!(encode_dna_2bit(seq).is_err());
        let seq2 = b"ACGT ";
         assert!(encode_dna_2bit(seq2).is_err());
    }

    #[test]
    fn test_reverse_complement() {
        assert_eq!(reverse_complement(b"ACGT"), b"ACGT");
        assert_eq!(reverse_complement(b"AAAA"), b"TTTT");
        assert_eq!(reverse_complement(b"GATTACA"), b"TGTAATC");
        assert_eq!(reverse_complement(b"acgt"), b"ACGT");
        assert_eq!(reverse_complement(b"AGCT"), b"AGCT");
        assert_eq!(reverse_complement(b""), b"");
        assert_eq!(reverse_complement(b"ACGNT"), b"ACGT"); // Skips N
    }
} 