use anyhow::{Context, Result};
use bio::io::fasta;
// Use types from the vbinseq 0.1.7 crate
use vbinseq::header::VBinseqHeader;
use vbinseq::writer::VBinseqWriterBuilder;
use std::fs::File;
use std::io::BufWriter;
use std::path::PathBuf;

/// Encodes a FASTA file into VBQ format.
///
/// Args:
///     input: Path to the input FASTA file.
///     output: Path to the output VBQ file.
///     level: Compression level (0-9).
pub fn run_fasta_to_vbq(input: PathBuf, output: PathBuf, level: u8) -> Result<()> {
    // Validate compression level
    if level > 9 {
        // vbinseq 0.1.7 might only support level 0 (uncompressed) or >0 (compressed)
        // We'll keep the check but note that only 0 or 1 might be meaningful internally.
        println!("Warning: Compression level > 0 specified. vbinseq 0.1.7 might only distinguish between 0 and >0.");
    }

    println!("Reading FASTA from: {}", input.display());
    println!("Writing VBQ to: {}", output.display());
    println!("Using compression level: {} (interpreted by vbinseq 0.1.7)", level);

    // --- Input FASTA Reader ---
    let reader = fasta::Reader::from_file(&input)
        .with_context(|| format!("Failed to open input FASTA file: {}", input.display()))?;

    // --- Output VBQ Writer ---
    let output_file = File::create(&output)
        .with_context(|| format!("Failed to create output VBQ file: {}", output.display()))?;
    let writer = BufWriter::new(output_file);

    // --- VBQ Header (using vbinseq 0.1.7 header) ---
    let qual = false; // FASTA files typically don't have quality scores
    let compressed = level > 0;
    let paired = false; // Assuming single-end sequences from FASTA
    // Reverting to positional arguments based on compiler error
    let header = VBinseqHeader::new(qual, compressed, paired);
    // header.qual = qual;
    // header.compressed = compressed;
    // header.paired = paired;
    // block_size might be set automatically by the builder/writer

    // --- VBQ Writer (using vbinseq 0.1.7 builder) ---
    let mut writer = VBinseqWriterBuilder::default()
        .header(header)
        // If compression level needs explicit setting, check vbinseq docs/API
        // .compression_level(if compressed { Some(level) } else { None })
        .build(writer)
        .context("Failed to build VBinseqWriter (vbinseq 0.1.7)")?;

    // --- Process Records ---
    let mut record_count = 0;
    let mut total_bases = 0;
    for result in reader.records() {
        let record = result.with_context(|| format!("Failed to read record from {}", input.display()))?;
        record_count += 1;
        total_bases += record.seq().len();

        let seq = record.seq();

        // Modify to use write_nucleotides_quality, passing empty quality slice
        let quality: &[u8] = b""; // Empty slice for quality
        writer.write_nucleotides_quality((record_count - 1) as u64, seq, quality)
            .with_context(|| format!("Failed to write record '{}' (using write_nucleotides_quality) to VBQ file (vbinseq 0.1.7)", record.id()))?;
    }

    // Finalize the VBQ file
    writer.finish().context("Failed to finalize VBQ file (vbinseq 0.1.7)")?;

    println!(
        "Successfully encoded {} records ({} bases) to {} using vbinseq 0.1.7",
        record_count,
        total_bases,
        output.display()
    );

    Ok(())
}

// --- New Custom Binary Format Encoder --- 

use byteorder::{LittleEndian, WriteBytesExt};
use bitnuc::encode_alloc; // Use bitnuc directly
use std::io::Write; // Need write_all and flush

struct EncodedSequenceData {
    original_len: u64,
    packed_len: u64, // Now length in u64 elements
    packed_data: Vec<u64>, // Corrected type
}

/// Encodes a FASTA file into a custom Orbweaver binary format (.orb).
/// Format:
///   - num_sequences: u64
///   - for each sequence:
///     - original_length: u64
///     - packed_length (in u64 elements): u64 
///     - packed_data: [u64] (Little Endian)
pub fn run_fasta_to_orb_bin(input: PathBuf, output: PathBuf) -> Result<()> {
    println!("Reading FASTA from: {}", input.display());
    println!("Writing custom ORB format to: {}", output.display());

    // --- Input FASTA Reader ---
    let reader = fasta::Reader::from_file(&input)
        .with_context(|| format!("Failed to open input FASTA file: {}", input.display()))?;

    // --- Collect encoded data first ---
    let mut encoded_sequences: Vec<EncodedSequenceData> = Vec::new();
    let mut total_bases = 0;

    for result in reader.records() {
        let record = result.with_context(|| format!("Failed to read record from {}", input.display()))?;
        let seq_bytes = record.seq();
        let original_len = seq_bytes.len() as u64;
        total_bases += original_len;

        // Pre-process sequence: Replace non-ACGT with 'A'
        let mut processed_seq = seq_bytes.to_vec();
        let mut replacements = 0;
        for base in processed_seq.iter_mut() {
            match base.to_ascii_uppercase() { // Handle upper/lower case
                b'A' | b'C' | b'G' | b'T' => { /* Keep valid */ }
                _ => { // Replace anything else
                    *base = b'A';
                    replacements += 1;
                }
            }
        }
        if replacements > 0 {
            println!("    [Info] Replaced {} non-ACGT characters with 'A' in record {}", 
                     replacements, record.id());
        }

        // Encode using bitnuc -> Vec<u64>
        let packed_data: Vec<u64> = encode_alloc(&processed_seq) // Use processed sequence
            .context(format!("Failed to encode sequence for record {} using bitnuc", 
                             record.id()))?;
        let packed_len = packed_data.len() as u64; // Length is number of u64s

        encoded_sequences.push(EncodedSequenceData {
            original_len,
            packed_len,
            packed_data,
        });
    }

    let num_sequences = encoded_sequences.len() as u64;
    println!(
        "Read {} records ({} bases), preparing to write to {}",
        num_sequences,
        total_bases,
        output.display()
    );

    // --- Output ORB Writer ---
    let output_file = File::create(&output)
        .with_context(|| format!("Failed to create output ORB file: {}", output.display()))?;
    let mut writer = BufWriter::new(output_file);

    // Write header (number of sequences)
    writer.write_u64::<LittleEndian>(num_sequences)
        .context("Failed to write number of sequences to ORB file")?;

    // Write each sequence's metadata and data
    let mut written_count = 0;
    for encoded_data in encoded_sequences {
        writer.write_u64::<LittleEndian>(encoded_data.original_len)
            .with_context(|| format!("Failed to write original length for sequence {}", written_count))?;
        // Write packed length (number of u64 elements)
        writer.write_u64::<LittleEndian>(encoded_data.packed_len)
            .with_context(|| format!("Failed to write packed length for sequence {}", written_count))?;
        // Write the actual packed u64 data elements
        for packed_u64 in encoded_data.packed_data {
             writer.write_u64::<LittleEndian>(packed_u64)
                 .with_context(|| format!("Failed to write packed data element for sequence {}", written_count))?;
        }
        written_count += 1;
    }

    // Ensure buffer is flushed
    writer.flush().context("Failed to flush ORB writer buffer")?;

    println!(
        "Successfully encoded {} records ({} bases) to {}",
        written_count,
        total_bases,
        output.display()
    );

    Ok(())
} 