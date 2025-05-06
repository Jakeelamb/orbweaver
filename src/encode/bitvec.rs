use anyhow::{Result, bail};
use super::dna_2bit::EncodedBase;
use std::fmt;
use std::iter::FromIterator;

/// BitVector represents a DNA sequence as a compact bit array.
/// Each base is encoded using 2 bits:
/// A: 00, C: 01, G: 10, T: 11
#[derive(Clone, PartialEq, Eq)]
pub struct BitVector {
    /// The actual bit data storage
    data: Vec<u8>,
    /// Number of bases in the sequence
    length: usize,
    /// Current position for streaming operations
    position: usize,
}

impl BitVector {
    /// Create a new empty BitVector
    pub fn new() -> Self {
        Self {
            data: Vec::new(),
            length: 0,
            position: 0,
        }
    }

    /// Create a BitVector with a preallocated capacity
    pub fn with_capacity(capacity: usize) -> Self {
        // Each byte stores 4 bases
        let byte_capacity = (capacity + 3) / 4;
        Self {
            data: Vec::with_capacity(byte_capacity),
            length: 0,
            position: 0,
        }
    }

    /// Create a BitVector from a slice of encoded bases
    pub fn from_encoded_bases(bases: &[EncodedBase]) -> Self {
        if bases.is_empty() {
            return Self::new();
        }
        
        let mut bitvec = Self::with_capacity(bases.len());
        
        // Process in chunks of 4 bases for better efficiency
        let mut i = 0;
        while i + 4 <= bases.len() {
            // Pack 4 bases into a single byte
            let mut byte = 0u8;
            byte |= (bases[i].0 & 0b11) << 6;
            byte |= (bases[i+1].0 & 0b11) << 4;
            byte |= (bases[i+2].0 & 0b11) << 2;
            byte |= bases[i+3].0 & 0b11;
            
            bitvec.data.push(byte);
            bitvec.length += 4;
            i += 4;
        }
        
        // Handle the remaining bases (0-3)
        if i < bases.len() {
            let mut byte = 0u8;
            let mut shift = 6;
            
            while i < bases.len() {
                byte |= (bases[i].0 & 0b11) << shift;
                shift -= 2;
                i += 1;
                bitvec.length += 1;
            }
            
            bitvec.data.push(byte);
        }
        
        bitvec
    }

    /// Create a BitVector from a DNA string
    pub fn from_string(s: &str) -> Result<Self> {
        let mut bitvec = Self::with_capacity(s.len());
        
        // Process in chunks of 4 characters for better efficiency
        let bytes = s.as_bytes();
        let mut i = 0;
        
        while i + 4 <= bytes.len() {
            let mut byte = 0u8;
            
            // Convert and pack 4 bases
            for j in 0..4 {
                match EncodedBase::from_base(bytes[i+j]) {
                    Some(encoded) => {
                        byte |= (encoded.0 & 0b11) << (6 - j*2);
                    },
                    None => bail!("Invalid DNA base '{}' found in sequence", bytes[i+j] as char),
                }
            }
            
            bitvec.data.push(byte);
            bitvec.length += 4;
            i += 4;
        }
        
        // Handle remaining characters (0-3)
        if i < bytes.len() {
            let mut byte = 0u8;
            let mut shift = 6;
            
            while i < bytes.len() {
                match EncodedBase::from_base(bytes[i]) {
                    Some(encoded) => {
                        byte |= (encoded.0 & 0b11) << shift;
                        shift -= 2;
                        bitvec.length += 1;
                    },
                    None => bail!("Invalid DNA base '{}' found in sequence", bytes[i] as char),
                }
                i += 1;
            }
            
            bitvec.data.push(byte);
        }
        
        Ok(bitvec)
    }

    /// Add a base to the end of the vector
    pub fn push(&mut self, base: EncodedBase) {
        let byte_index = self.length / 4;
        let bit_offset = (self.length % 4) * 2;
        
        // Ensure we have enough space
        if byte_index >= self.data.len() {
            self.data.push(0);
        }
        
        // Clear the bits and then set them
        self.data[byte_index] &= !(0b11 << (6 - bit_offset));
        self.data[byte_index] |= (base.0 & 0b11) << (6 - bit_offset);
        
        self.length += 1;
    }

    /// Get the base at the specified index
    pub fn get(&self, index: usize) -> Option<EncodedBase> {
        if index >= self.length {
            return None;
        }
        
        let byte_index = index / 4;
        let bit_offset = (index % 4) * 2;
        
        let byte = self.data[byte_index];
        let base_bits = (byte >> (6 - bit_offset)) & 0b11;
        
        Some(EncodedBase(base_bits))
    }

    /// Get length (number of bases)
    pub fn len(&self) -> usize {
        self.length
    }
    
    /// Check if the BitVector is empty
    pub fn is_empty(&self) -> bool {
        self.length == 0
    }
    
    /// Extract a slice of the BitVector
    pub fn slice(&self, start: usize, end: usize) -> Result<Self> {
        if start > end || end > self.length {
            bail!("Invalid slice range: {}..{} (length: {})", start, end, self.length);
        }
        
        // For small slices or slices that are not byte-aligned, use the simple approach
        if end - start < 16 || start % 4 != 0 || end % 4 != 0 {
            let mut result = Self::with_capacity(end - start);
            for i in start..end {
                if let Some(base) = self.get(i) {
                    result.push(base);
                }
            }
            return Ok(result);
        }
        
        // For large byte-aligned slices, copy the bytes directly
        let start_byte = start / 4;
        let end_byte = (end + 3) / 4;
        
        let result = Self {
            data: self.data[start_byte..end_byte].to_vec(),
            length: end - start,
            position: 0,
        };
        
        Ok(result)
    }
    
    /// Convert the BitVector to a String
    pub fn to_string(&self) -> String {
        let mut result = String::with_capacity(self.length);
        
        // Process in chunks for better efficiency
        for byte_idx in 0..self.data.len() {
            let byte = self.data[byte_idx];
            let bases_in_byte = std::cmp::min(4, self.length - byte_idx * 4);
            
            for i in 0..bases_in_byte {
                let shift = 6 - i * 2;
                let base_bits = (byte >> shift) & 0b11;
                result.push(EncodedBase(base_bits).to_char());
            }
        }
        
        result
    }
    
    /// Get the reverse complement of this sequence
    pub fn reverse_complement(&self) -> Self {
        let mut result = Self::with_capacity(self.length);
        
        // For very short sequences, use the simple approach
        if self.length < 8 {
            for i in (0..self.length).rev() {
                if let Some(base) = self.get(i) {
                    result.push(base.revcomp());
                }
            }
            return result;
        }
        
        // For longer sequences, process in chunks for better efficiency
        let mut bases = Vec::with_capacity(self.length);
        
        // First extract all bases
        for i in 0..self.length {
            if let Some(base) = self.get(i) {
                bases.push(base);
            }
        }
        
        // Then push them in reverse order with complement
        for base in bases.iter().rev() {
            result.push(base.revcomp());
        }
        
        result
    }
    
    /// Iterate through the bases in the BitVector
    pub fn iter(&self) -> BitVectorIterator {
        BitVectorIterator {
            bitvec: self,
            position: 0,
        }
    }
    
    /// Read the next chunk of bases for streaming operations
    pub fn read_chunk(&mut self, max_size: usize) -> Vec<EncodedBase> {
        let mut chunk = Vec::with_capacity(max_size);
        let end_pos = std::cmp::min(self.position + max_size, self.length);
        
        for i in self.position..end_pos {
            if let Some(base) = self.get(i) {
                chunk.push(base);
            }
        }
        
        self.position = end_pos;
        chunk
    }
    
    /// Reset the position for streaming operations
    pub fn reset_position(&mut self) {
        self.position = 0;
    }
    
    /// Check if we've reached the end of the stream
    pub fn eof(&self) -> bool {
        self.position >= self.length
    }
    
    /// Calculate memory usage in bytes
    pub fn memory_usage(&self) -> usize {
        // Size of the data vector + overhead
        self.data.capacity() + 
        // Size of the Vec struct itself (3 pointers)
        3 * std::mem::size_of::<usize>() +
        // Size of the BitVector struct fields
        2 * std::mem::size_of::<usize>()
    }
}

/// Iterator for BitVector
pub struct BitVectorIterator<'a> {
    bitvec: &'a BitVector,
    position: usize,
}

impl<'a> Iterator for BitVectorIterator<'a> {
    type Item = EncodedBase;
    
    fn next(&mut self) -> Option<Self::Item> {
        if self.position < self.bitvec.len() {
            let base = self.bitvec.get(self.position);
            self.position += 1;
            base
        } else {
            None
        }
    }
}

impl FromIterator<EncodedBase> for BitVector {
    fn from_iter<I: IntoIterator<Item = EncodedBase>>(iter: I) -> Self {
        let mut bitvec = Self::new();
        for base in iter {
            bitvec.push(base);
        }
        bitvec
    }
}

impl fmt::Debug for BitVector {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "BitVector [len={}] '{}'", self.length, self.to_string())
    }
}

impl From<&[EncodedBase]> for BitVector {
    fn from(bases: &[EncodedBase]) -> Self {
        Self::from_encoded_bases(bases)
    }
}

// Add a convenient method to estimate memory savings
pub fn estimate_memory_savings(original_size: usize, two_bit_size: usize) -> (f64, String) {
    let original_bytes = original_size;
    let two_bit_bytes = (two_bit_size + 3) / 4; // Round up to nearest byte
    
    let percentage = if original_bytes > 0 {
        100.0 * (original_bytes - two_bit_bytes) as f64 / original_bytes as f64
    } else {
        0.0
    };
    
    let description = if percentage >= 70.0 {
        "Excellent (>70%)"
    } else if percentage >= 50.0 {
        "Good (50-70%)"
    } else if percentage >= 30.0 {
        "Moderate (30-50%)"
    } else {
        "Minimal (<30%)"
    };
    
    (percentage, description.to_string())
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_new_bitvector() {
        let bv = BitVector::new();
        assert_eq!(bv.len(), 0);
        assert!(bv.is_empty());
    }
    
    #[test]
    fn test_push_and_get() {
        let mut bv = BitVector::new();
        bv.push(EncodedBase(0b00)); // A
        bv.push(EncodedBase(0b01)); // C
        bv.push(EncodedBase(0b10)); // G
        bv.push(EncodedBase(0b11)); // T
        
        assert_eq!(bv.len(), 4);
        assert_eq!(bv.get(0), Some(EncodedBase(0b00)));
        assert_eq!(bv.get(1), Some(EncodedBase(0b01)));
        assert_eq!(bv.get(2), Some(EncodedBase(0b10)));
        assert_eq!(bv.get(3), Some(EncodedBase(0b11)));
        assert_eq!(bv.get(4), None);
    }
    
    #[test]
    fn test_from_string() {
        let bv = BitVector::from_string("ACGT").unwrap();
        assert_eq!(bv.len(), 4);
        assert_eq!(bv.to_string(), "ACGT");
        
        // Test with invalid characters
        let result = BitVector::from_string("ACNGT");
        assert!(result.is_err());
    }
    
    #[test]
    fn test_to_string() {
        let bv = BitVector::from_string("ACGTACGT").unwrap();
        assert_eq!(bv.to_string(), "ACGTACGT");
    }
    
    #[test]
    fn test_slice() {
        let bv = BitVector::from_string("ACGTACGT").unwrap();
        let slice = bv.slice(2, 6).unwrap();
        assert_eq!(slice.to_string(), "GTAC");
        
        // Test invalid slice
        let invalid_slice = bv.slice(6, 2);
        assert!(invalid_slice.is_err());
        
        let out_of_bounds = bv.slice(0, 9);
        assert!(out_of_bounds.is_err());
    }
    
    #[test]
    fn test_reverse_complement() {
        let bv = BitVector::from_string("ACGT").unwrap();
        let revcomp = bv.reverse_complement();
        assert_eq!(revcomp.to_string(), "ACGT"); // ACGT is its own reverse complement
        
        let bv2 = BitVector::from_string("AAGT").unwrap();
        let revcomp2 = bv2.reverse_complement();
        assert_eq!(revcomp2.to_string(), "ACTT");
    }
    
    #[test]
    fn test_iter() {
        let bv = BitVector::from_string("ACGT").unwrap();
        let bases: Vec<EncodedBase> = bv.iter().collect();
        
        assert_eq!(bases.len(), 4);
        assert_eq!(bases[0], EncodedBase(0b00));
        assert_eq!(bases[1], EncodedBase(0b01));
        assert_eq!(bases[2], EncodedBase(0b10));
        assert_eq!(bases[3], EncodedBase(0b11));
    }
    
    #[test]
    fn test_storage_efficiency() {
        // Verify we're using 2 bits per base
        let dna = "ACGTACGTACGTACGTACGTACGT"; // 24 bases
        let bv = BitVector::from_string(dna).unwrap();
        
        // 24 bases * 2 bits per base = 48 bits = 6 bytes
        assert_eq!(bv.data.len(), 6);
    }
} 