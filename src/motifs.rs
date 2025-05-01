use anyhow::{Context, Result};
use bio::data_structures::suffix_array::{lcp, suffix_array};
use bitnuc::decode;
use byteorder::{LittleEndian, ReadBytesExt};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufReader, Read};
use std::path::Path;
use std::str;
use sysinfo::System;

/// Analyzes a single chunk of sequence data for repeats.
///
/// Args:
///    chunk_data: The byte slice representing the sequence chunk.
///    min_len: The minimum length of repeated substrings to report.
///
/// Returns:
///    A HashMap where keys are the repeated substrings (of length >= min_len)
///    found within this chunk, and values are their occurrence counts within the chunk.
fn analyze_chunk(chunk_data: &[u8], min_len: usize) -> Result<HashMap<String, usize>> {
    if chunk_data.is_empty() || min_len == 0 {
        return Ok(HashMap::new());
    }

    // Append the sentinel byte for suffix array construction
    let mut seq_bytes_with_sentinel = chunk_data.to_vec();
    seq_bytes_with_sentinel.push(b'\x00'); // Use a sentinel not in the alphabet

    let seq_bytes = &seq_bytes_with_sentinel;
    let n = chunk_data.len(); // Use original chunk length for boundary checks

    let sa = suffix_array(seq_bytes);
    let lcp = lcp(seq_bytes, &sa);

    let mut repeat_counts: HashMap<String, usize> = HashMap::new();
    let mut stack: Vec<(usize, usize)> = Vec::new(); // (lcp_value, index_in_sa)

    for i in 1..=sa.len() {
        let current_lcp = if i < sa.len() {
            // The closure takes the value directly (isize), not a reference.
            lcp.get(i).map_or(0, |val| if val < 0 { 0 } else { val as usize })
        } else {
            0 // LCP is 0 for the fictitious suffix after the last one
        };

        let mut last_popped_index = i; // Tracks the original SA index for the current LCP interval start
        while !stack.is_empty() && stack.last().unwrap().0 > current_lcp {
            let (prev_lcp, prev_index) = stack.pop().unwrap();
            last_popped_index = prev_index; // Update start index as we pop

            // Interval starts after the index of the element below the popped one on the stack
            let interval_start_sa_idx = stack.last().map_or(0, |&(_, idx)| idx);
            // Interval count is the number of suffixes covered by this LCP value
            let interval_count = i - interval_start_sa_idx;

            if prev_lcp >= min_len && interval_count >= 2 {
                // Use the SA index from before the pop started (last_popped_index) for representative suffix
                let suffix_start = sa[last_popped_index] as usize; // Use last_popped_index here

                // Ensure the repeat slice doesn't exceed the original chunk bounds (n)
                if suffix_start + prev_lcp <= n {
                    let repeat_slice = &chunk_data[suffix_start..suffix_start + prev_lcp];
                    match str::from_utf8(repeat_slice) {
                        Ok(s) => {
                            let count_entry = repeat_counts.entry(s.to_string()).or_insert(0);
                             // The interval_count represents the number of suffixes sharing this prefix (frequency)
                            if interval_count > *count_entry {
                                 *count_entry = interval_count;
                            }
                        }
                        Err(_) => {
                             // This should be rare if input is DNA, but handle non-UTF8 possibility
                            eprintln!(
                                "Warning: Found non-UTF8 repeat candidate (length {}), skipping.",
                                prev_lcp
                            );
                        }
                    }
                } else {
                     // This can happen if the repeat extends into the sentinel byte, ignore it.
                     // eprintln!("Debug: Repeat candidate length {} starting at {} extends beyond chunk boundary {}", prev_lcp, suffix_start, n);
                }
            }
        }

        // If the stack is empty or the new LCP is greater than the top, push it.
        // We use last_popped_index to associate the LCP value with the correct starting suffix index
        // for the interval it represents.
        if stack.is_empty() || stack.last().unwrap().0 < current_lcp {
            stack.push((current_lcp, last_popped_index)); // Push current LCP and its starting SA index
        }
         // If current_lcp == stack top, we extend the interval, no stack change needed.
         // If current_lcp < stack top, the while loop handled it.
    }

    // Filter repeats that don't meet the minimum count threshold (already handled by interval_count >= 2 check)
    // repeat_counts.retain(|_, count| *count >= 2); // Not strictly needed anymore

    Ok(repeat_counts)
}


/// Finds repeats in sequences from a custom Orbweaver binary format (.orb) using chunking.
///
/// This version reads the sequence data in chunks to handle potentially large files
/// that might exceed available memory if processed whole. It uses overlapping chunks
/// to catch repeats that might span chunk boundaries. The reported counts are heuristic,
/// representing the maximum count observed for a repeat in any single chunk where it appeared.
///
/// Args:
///    input_path: Path to the input ORB file.
///    min_len: The minimum length of repeated substrings to report.
///
/// Returns:
///    A HashMap where keys are the repeated substrings (of length >= min_len)
///    and values are their *maximum observed* occurrence counts within any single chunk.
pub fn find_repeats(input_path: &Path, min_len: usize) -> Result<HashMap<String, usize>> {
    if min_len == 0 {
        return Ok(HashMap::new());
    }

     // --- Determine Chunk Size and Overlap ---
    let mut sys = System::new_all();
    sys.refresh_memory();
    let total_memory = sys.total_memory(); // Total physical memory in bytes

    // Heuristic: Aim for chunks requiring roughly 1/10th of total memory for SA/LCP arrays.
    // Suffix array memory is approx 5N to 10N bytes. Let's be conservative (15N).
    // So, N = total_memory / 10 / 15
    let chunk_size = (total_memory / 150).max(1024 * 1024) as usize; // Min 1MB chunk size
    // Overlap should be larger than min_len, maybe related to chunk size.
    let overlap = (chunk_size / 20).max(min_len * 2 + 1) as usize; // Ensure overlap > min_len

    println!(
        "System total memory: {} bytes. Using chunk size: {} bytes, overlap: {} bytes",
        total_memory, chunk_size, overlap
    );
    println!(
        "Reading sequence(s) from ORB file: {} (min_len: {})",
        input_path.display(),
        min_len
    );

    let input_file = File::open(input_path)
        .with_context(|| format!("Failed to open input ORB file: {}", input_path.display()))?;
    let mut reader = BufReader::new(input_file);

    let num_sequences = reader
        .read_u64::<LittleEndian>()
        .context("Failed to read number of sequences")?;
    println!("  Found {} sequences in file.", num_sequences);

    let mut buffer: Vec<u8> = Vec::with_capacity(chunk_size + overlap); // Buffer for incoming sequence data
    let mut global_repeat_counts: HashMap<String, usize> = HashMap::new();
    let mut total_bases_processed: u64 = 0;
    let mut chunk_index = 0;

    for i in 0..num_sequences {
        let original_len = reader
            .read_u64::<LittleEndian>()
            .with_context(|| format!("Failed to read original length for sequence {}", i))?;
        let packed_len = reader
            .read_u64::<LittleEndian>()
            .with_context(|| format!("Failed to read packed length for sequence {}", i))?;

        // println!("    Reading sequence {}: original_len={}, packed_len={} (u64s)",
        //          i, original_len, packed_len);

        if original_len == 0 {
            continue; // Skip empty sequences
        }

        let mut packed_data_u64: Vec<u64> = vec![0; packed_len as usize];
        reader
            .read_u64_into::<LittleEndian>(&mut packed_data_u64)
            .with_context(|| format!("Failed to read packed u64 data for sequence {}", i))?;

        // Decode directly into a temporary buffer or extend the main buffer
        // Decoding requires exact size, so decode then extend
        let mut decoded_sequence: Vec<u8> = vec![0; original_len as usize];
        decode(&packed_data_u64, original_len as usize, &mut decoded_sequence)
            .context(format!("Failed to decode sequence {} using bitnuc", i))?;

        buffer.extend_from_slice(&decoded_sequence);
        total_bases_processed += original_len;

        // Process chunks as buffer fills
        while buffer.len() >= chunk_size {
            println!(
                "Processing chunk {} ({} bytes)...",
                chunk_index, chunk_size
            );
            let chunk_to_process = &buffer[0..chunk_size];

            match analyze_chunk(chunk_to_process, min_len) {
                Ok(local_counts) => {
                    for (repeat, local_count) in local_counts {
                        global_repeat_counts
                            .entry(repeat)
                            .and_modify(|global_count| {
                                *global_count = (*global_count).max(local_count);
                            })
                            .or_insert(local_count);
                    }
                }
                Err(e) => {
                    eprintln!(
                        "Error analyzing chunk {}: {}. Skipping chunk.",
                        chunk_index, e
                    );
                    // Optionally: return Err(e) or continue processing other chunks
                }
            }

            // Prepare buffer for the next chunk: keep only the overlap
            let next_start = chunk_size.saturating_sub(overlap);
            buffer.drain(0..next_start); // More efficient than creating a new vec
            chunk_index += 1;
        }
    }

    // Process the final remaining data in the buffer
    if !buffer.is_empty() {
        println!(
            "Processing final chunk {} ({} bytes)...",
            chunk_index,
            buffer.len()
        );
        match analyze_chunk(&buffer, min_len) {
            Ok(local_counts) => {
                for (repeat, local_count) in local_counts {
                    global_repeat_counts
                        .entry(repeat)
                        .and_modify(|global_count| {
                            *global_count = (*global_count).max(local_count);
                        })
                        .or_insert(local_count);
                }
            }
            Err(e) => {
                eprintln!(
                    "Error analyzing final chunk {}: {}. Results might be incomplete.",
                    chunk_index, e
                );
                // Decide if this is a fatal error
            }
        }
    }

    println!(
        "Chunked repeat analysis complete. Total bases processed: {}. Found {} unique potential repeats (using max count heuristic).",
        total_bases_processed,
        global_repeat_counts.len()
    );

     // Final filter based on the aggregated max counts
    global_repeat_counts.retain(|_, count| *count >= 2);

    println!("Filtered down to {} repeats occurring at least twice in at least one chunk.", global_repeat_counts.len());


    Ok(global_repeat_counts)
}

// --- Tests (Optional: Need to be adapted for chunking or use smaller test files) ---
// #[cfg(test)]
// mod tests {
//     use super::*;
//     use std::io::Write;
//     use tempfile::NamedTempFile;
//     use byteorder::{LittleEndian, WriteBytesExt};
//     use bitnuc::encode; // Assuming encode function exists and works complementarily to decode

//     // Helper to create a dummy ORB file for testing
//     fn create_test_orb(sequences: &[&[u8]], path: &Path) -> Result<()> {
//         let mut file = File::create(path)?;
//         file.write_u64::<LittleEndian>(sequences.len() as u64)?;

//         for seq in sequences {
//             let original_len = seq.len() as u64;
//             // Simulate packing with bitnuc::encode (replace with actual if available)
//              let (packed_data, packed_len_u64) = encode(seq)?; // Assuming encode returns (Vec<u64>, u64)


//             file.write_u64::<LittleEndian>(original_len)?;
//             file.write_u64::<LittleEndian>(packed_len_u64)?;
//              for &word in &packed_data {
//                 file.write_u64::<LittleEndian>(word)?;
//             }
//         }
//         Ok(())
//     }

//     // Test case needs adjustment for chunking logic and bitnuc::encode availability
//     #[test]
//     fn test_find_repeats_chunked() {
//         // This test needs bitnuc::encode to create the file,
//         // and the expected results might change based on chunking boundaries.
//         // It's complex to set up without the exact encode function.
//         // Skipping detailed implementation for now.
//         assert!(true); // Placeholder
//     }
// }