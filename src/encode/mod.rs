pub mod dna_2bit;
pub mod kmer;
pub mod bitvec;
pub mod simd;

pub use dna_2bit::EncodedBase;
pub use simd::{encode_dna_batch, count_digrams, fingerprint_sequence}; 