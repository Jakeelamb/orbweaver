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

// --- New Custom Binary Format Encoder (Decoded Data) --- 

use byteorder::{LittleEndian, WriteBytesExt};
// Removed bitnuc::encode_alloc as we write decoded data
use std::io::{Write, Seek, SeekFrom}; // Need Seek for header update

// Removed EncodedSequenceData struct

/// Encodes a FASTA file into a custom Orbweaver binary format (.orb) storing DECODED data.
/// Format:
///   - num_sequences: u64 (Little Endian)
///   - total_decoded_length: u64 (Little Endian)
///   - concatenated_decoded_data: [u8] (Raw sequence bytes)
pub fn run_fasta_to_orb_bin(input: PathBuf, output: PathBuf) -> Result<()> {
    println!("Reading FASTA from: {}", input.display());
    println!("Writing custom ORB format (decoded data) to: {}", output.display());

    // --- Input FASTA Reader ---
    let reader = fasta::Reader::from_file(&input)
        .with_context(|| format!("Failed to open input FASTA file: {}", input.display()))?;

    // --- Output ORB Writer ---
    let mut output_file = File::create(&output)
        .with_context(|| format!("Failed to create output ORB file: {}", output.display()))?;
    // Write placeholder for header (2 * u64 = 16 bytes)
    let header_placeholder: [u8; 16] = [0; 16];
    output_file.write_all(&header_placeholder)
        .context("Failed to write header placeholder")?;
    
    // Use BufWriter for the sequence data part
    let mut writer = BufWriter::new(&mut output_file); // Pass mutable reference

    let mut num_sequences: u64 = 0;
    let mut total_decoded_length: u64 = 0;

    for result in reader.records() {
        let record = result.with_context(|| format!("Failed to read record from {}", input.display()))?;
        let seq_bytes = record.seq();
        num_sequences += 1;

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

        // Write processed sequence directly
        writer.write_all(&processed_seq)
             .with_context(|| format!("Failed to write sequence data for record {}", record.id()))?;
        total_decoded_length += processed_seq.len() as u64;

    }

    // Ensure BufWriter buffer is flushed before writing header
    writer.flush().context("Failed to flush ORB writer buffer")?;
    
    // Drop the BufWriter to regain ownership of the File for seeking
    drop(writer);

    println!(
        "Finished writing sequence data. Read {} records ({} bases). Writing header...",
        num_sequences,
        total_decoded_length
    );

    // Seek back to the beginning to write the actual header
    // output_file is already available after BufWriter is dropped
    output_file.seek(SeekFrom::Start(0))
        .context("Failed to seek to beginning of ORB file for header write")?;

    output_file.write_u64::<LittleEndian>(num_sequences)
        .context("Failed to write final number of sequences to ORB header")?;
    output_file.write_u64::<LittleEndian>(total_decoded_length)
        .context("Failed to write final total decoded length to ORB header")?;

    println!(
        "Successfully encoded {} records ({} bases) to {}",
        num_sequences,
        total_decoded_length,
        output.display()
    );

    Ok(())
} 