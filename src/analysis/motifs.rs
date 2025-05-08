use anyhow::{Context, Result};
use bio::data_structures::suffix_array::{lcp, suffix_array};
use byteorder::{LittleEndian, ReadBytesExt};
use std::collections::{HashMap, BinaryHeap, HashSet};
use std::fs::{File};
use std::path::Path;
use sysinfo::{System, SystemExt};
// use rayon::prelude::*; // Removed unused rayon
use crossbeam_channel::{bounded, Receiver, Sender};
use memmap2::Mmap;
use std::sync::Arc;
use std::hash::{Hasher};
use twox_hash::XxHash64;
use std::cmp::Reverse;

/// Analyzes a single chunk of sequence data for repeats.
///
/// Args:
///    chunk_data: The byte slice representing the sequence chunk.
///    min_len: The minimum length of repeated substrings to report.
///
/// Returns:
///    A HashMap where keys are the repeated substrings (as Vec<u8>) found within this chunk,
///    and values are their occurrence counts within the chunk.
fn analyze_chunk(chunk_data: &[u8], min_len: usize) -> Result<HashMap<Vec<u8>, usize>> {
    if chunk_data.is_empty() || min_len == 0 {
        return Ok(HashMap::new());
    }

    let mut seq_bytes_with_sentinel = chunk_data.to_vec();
    seq_bytes_with_sentinel.push(b'\x00');

    let seq_bytes = &seq_bytes_with_sentinel;
    let n = chunk_data.len();

    let sa = suffix_array(seq_bytes);
    let lcp = lcp(seq_bytes, &sa);

    let mut repeat_counts: HashMap<Vec<u8>, usize> = HashMap::new(); // Key is Vec<u8>
    let mut stack: Vec<(usize, usize)> = Vec::new();

    for i in 1..=sa.len() {
        let current_lcp = if i < sa.len() {
            lcp.get(i).map_or(0, |val| if val < 0 { 0 } else { val as usize })
        } else {
            0
        };

        let mut last_popped_index = i;
        while !stack.is_empty() && stack.last().unwrap().0 > current_lcp {
            let (prev_lcp, prev_index) = stack.pop().unwrap();
            last_popped_index = prev_index;

            let interval_start_sa_idx = stack.last().map_or(0, |&(_, idx)| idx);
            let interval_count = i - interval_start_sa_idx;

            if prev_lcp >= min_len && interval_count >= 2 {
                let suffix_start = sa[last_popped_index] as usize;

                if suffix_start + prev_lcp <= n {
                    let repeat_slice = &chunk_data[suffix_start..suffix_start + prev_lcp];
                    // Directly use the byte slice, converting to Vec<u8> for the key
                    let repeat_key = repeat_slice.to_vec();
                    let count_entry = repeat_counts.entry(repeat_key).or_insert(0);
                    if interval_count > *count_entry {
                        *count_entry = interval_count;
                    }
                    // Removed the str::from_utf8 check here, handle potential non-utf8 during final merge
                } else {
                    // Repeat candidate extends beyond chunk boundary (into sentinel), ignore.
                }
            }
        }

        if stack.is_empty() || stack.last().unwrap().0 < current_lcp {
            stack.push((current_lcp, last_popped_index));
        }
    }

    Ok(repeat_counts)
}

/// Finds repeats in sequences from a custom Orbweaver binary format (.orb) using
/// parallel chunking on a memory-mapped file and suffix arrays.
///
/// Memory-maps the ORB file (which must contain decoded data), analyzes chunks
/// in parallel using Rayon, and merges results.
/// Uses Vec<u8> internally for repeat keys, converting to String for final output.
/// Reported counts are heuristic (max count observed in any chunk).
///
/// Args:
///    input_path: Path to the input ORB file (decoded data format).
///    min_len: The minimum length of repeated substrings to report.
///
/// Returns:
///    A HashMap where keys are the repeated substrings (String, length >= min_len)
///    and values are their *maximum observed* occurrence counts within any single chunk.
pub fn find_repeats(input_path: &Path, min_len: usize) -> Result<HashMap<String, usize>> {
    if min_len == 0 {
        return Ok(HashMap::new());
    }

    let mut sys = System::new_all();
    sys.refresh_memory();
    let total_memory = sys.total_memory();
    let chunk_size = (total_memory / 150).max(1024 * 1024) as usize;
    let overlap = (chunk_size / 20).max(min_len * 2 + 1) as usize;

    println!(
        "System total memory: {} bytes. Using chunk size: {} bytes, overlap: {} bytes",
        total_memory, chunk_size, overlap
    );
    println!(
        "Processing memory-mapped ORB file: {} (min_len: {}) with Suffix Array in parallel...",
        input_path.display(), min_len
    );

    // --- Memory Map the File --- 
    let input_file = File::open(input_path)
        .with_context(|| format!("Failed to open input ORB file: {}", input_path.display()))?;

    // Read header first (cannot read from Mmap directly easily)
    // Need a way to read first 16 bytes without consuming the reader for Mmap
    // Re-open or use a cursor/reader on the file before mapping
    let mut header_reader = std::io::BufReader::new(&input_file); 
    let num_sequences = header_reader
        .read_u64::<LittleEndian>()
        .context("Failed to read number of sequences from ORB header")?;
    let total_decoded_length = header_reader
        .read_u64::<LittleEndian>()
        .context("Failed to read total decoded length from ORB header")?;
    let header_size: usize = 16;
    drop(header_reader); // Drop reader before mapping

    println!("  Found {} sequences, total decoded length: {}. Mapping file...", num_sequences, total_decoded_length);

    // Map the entire file and wrap in Arc for sharing
    let mmap = Arc::new(unsafe { Mmap::map(&input_file)? });

    // Get a slice to the sequence data part of the map
    // Need to check length against Arc<Mmap>
    if (header_size + total_decoded_length as usize) > mmap.len() {
        anyhow::bail!(
            "ORB file header indicates total length {} + {} which exceeds file size {}",
            header_size, total_decoded_length, mmap.len()
        );
    }
    // We still create the slice here, but tasks will borrow from the Arc<Mmap>
    let sequence_data: &[u8] = &mmap[header_size..header_size + total_decoded_length as usize];
    println!("  Successfully memory-mapped sequence data ({} bytes).", sequence_data.len());

    // --- Parallel Chunk Processing --- 
    let (sender, receiver): (Sender<Result<HashMap<Vec<u8>, usize>>>, Receiver<Result<HashMap<Vec<u8>, usize>>>) = bounded(rayon::current_num_threads() * 2);
    let mut tasks_spawned = 0;

    // Calculate number of chunks based on total length, chunk size, and overlap
    let step_size = chunk_size.saturating_sub(overlap);
    if step_size == 0 {
        anyhow::bail!("Overlap size {} must be smaller than chunk size {}", overlap, chunk_size);
    }
    let num_chunks = if sequence_data.is_empty() { 
        0 
    } else { 
        // Calculate how many full steps we can take
        let full_steps = sequence_data.len().saturating_sub(chunk_size) / step_size;
        // Total chunks = initial chunk + full steps + potentially one final partial chunk
        1 + full_steps + if sequence_data.len() > chunk_size + full_steps * step_size { 1 } else { 0 }
    };

    println!("Spawning {} analysis tasks for {} bytes...", num_chunks, sequence_data.len());

    // Remove mmap_ref, tasks will capture Arc
    //let mmap_ref = &mmap; 

    for i in 0..num_chunks {
        let start = i * step_size;
        let end = (start + chunk_size).min(sequence_data.len());
        if start >= end { continue; }
        
        // Get a slice for the current chunk directly from the memory map
        // Slice is temporary for the spawn call setup
        //let chunk_slice = &sequence_data[start..end]; 

        let sender_clone = sender.clone();
        let current_chunk_index = i;
        let mmap_clone = Arc::clone(&mmap); // Clone the Arc for the closure
        
        rayon::spawn(move || {
            // Re-slice inside the closure using the cloned Arc
            // This ensures the slice is valid for the lifetime of the closure
            let sequence_data_in_task: &[u8] = &mmap_clone[header_size..header_size + total_decoded_length as usize];
            let chunk_slice_in_task = &sequence_data_in_task[start..end];

            let result = analyze_chunk(chunk_slice_in_task, min_len);
            if let Err(e) = &result {
                eprintln!("Error analyzing chunk {}: {}", current_chunk_index, e);
            }
            if sender_clone.send(result).is_err() {
                 eprintln!("Warning: Failed to send result for chunk {} (receiver dropped).", current_chunk_index);
             }
        });
        tasks_spawned += 1;
    }

    // Drop the original sender 
    drop(sender);

    // --- Collect and Merge Results (remains the same) --- 
    println!("Waiting for {} analysis tasks to complete and merging results...", tasks_spawned);
    let mut global_repeat_counts: HashMap<String, usize> = HashMap::new(); // Final result map

    for result in receiver {
        match result {
            Ok(local_counts) => {
                 for (repeat_vec, local_count) in local_counts {
                     match String::from_utf8(repeat_vec) { 
                         Ok(repeat_str) => {
                             global_repeat_counts
                                 .entry(repeat_str)
                                 .and_modify(|global_count| {
                                     *global_count = (*global_count).max(local_count);
                                 })
                                 .or_insert(local_count);
                         }
                         Err(_) => {
                             eprintln!("Warning: Skipping non-UTF8 repeat candidate of length {}", local_count);
                         }
                     }
                 }
            }
            Err(_) => {}
        }
    }

    println!(
        "Parallel analysis complete. Found {} unique potential repeats (using max count heuristic).",
        global_repeat_counts.len()
    );

    global_repeat_counts.retain(|_, count| *count >= 2);
    println!("Filtered down to {} repeats occurring at least twice in at least one chunk.", global_repeat_counts.len());

    // No longer need explicit drop(mmap_ref)
    //drop(mmap_ref); 
    Ok(global_repeat_counts)
}

// --- Rolling Hash Implementation ---

const BASE: u64 = 101; // Prime base for hashing
const MOD: u64 = 1_000_000_007; // Large prime modulus

/// Calculates rolling hashes for a given window size `k` over a sequence chunk.
///
/// Args:
///    chunk_data: The byte slice representing the sequence chunk.
///    k: The window size (minimum length of substring to hash).
///
/// Returns:
///    A HashMap where keys are the hash values of substrings of length `k`,
///    and values are their occurrence counts within the chunk.
fn rolling_hash_chunk(chunk_data: &[u8], k: usize) -> Result<HashMap<u64, usize>> {
    if k == 0 || k > chunk_data.len() {
        return Ok(HashMap::new()); // Cannot hash if k is 0 or larger than data
    }

    let mut hash_counts = HashMap::new();
    let mut current_hash = 0u64;
    let mut base_k = 1u64; // BASE^(k-1) % MOD

    // Calculate BASE^k % MOD for removing the leading character
    for _ in 1..k {
        base_k = (base_k * BASE) % MOD;
    }

    // Calculate the hash of the first window
    for i in 0..k {
        current_hash = (current_hash * BASE + chunk_data[i] as u64) % MOD;
    }
    *hash_counts.entry(current_hash).or_insert(0) += 1;

    // Slide the window across the rest of the chunk
    for i in k..chunk_data.len() {
        // Remove leading character's contribution
        let leading_char_val = chunk_data[i - k] as u64;
        current_hash = (current_hash + MOD - (leading_char_val * base_k % MOD)) % MOD;
        
        // Add trailing character's contribution
        let trailing_char_val = chunk_data[i] as u64;
        current_hash = (current_hash * BASE + trailing_char_val) % MOD;

        // Increment count for the new hash
        *hash_counts.entry(current_hash).or_insert(0) += 1;
    }

    // Filter out hashes that only occurred once (optional, depends on definition)
    // hash_counts.retain(|_, count| *count >= 2);

    Ok(hash_counts)
}

/// Finds frequent substring hashes in sequences from a memory-mapped Orbweaver binary file (.orb)
/// using parallel chunking and rolling hashes (Rabin-Karp style).
///
/// Memory-maps the ORB file, computes rolling hashes for a fixed window size (`min_len`)
/// in parallel using Rayon on slices of the mapped data, and merges the hash counts.
///
/// Args:
///    input_path: Path to the input ORB file (decoded data format).
///    min_len: The window size (k) for rolling hash calculation.
///
/// Returns:
///    A HashMap where keys are the u64 hash values of substrings of length `min_len`
///    and values are their *maximum observed* occurrence counts within any single chunk.
pub fn find_repeats_hashed(input_path: &Path, min_len: usize) -> Result<HashMap<u64, usize>> {
    if min_len == 0 {
        println!("Warning: min_len must be > 0 for rolling hash. Returning empty results.");
        return Ok(HashMap::new());
    }

    let mut sys = System::new_all();
    sys.refresh_memory();
    let total_memory = sys.total_memory();
    let chunk_size = (total_memory / 100).max(1024 * 1024) as usize;
    let overlap = (chunk_size / 20).max(min_len) as usize;

    println!(
        "System total memory: {} bytes. Using chunk size: {} bytes, overlap: {} bytes",
        total_memory, chunk_size, overlap
    );
    println!(
        "Processing memory-mapped ORB file: {} (window_size k={}) with Rolling Hash in parallel...",
        input_path.display(), min_len
    );

    // --- Memory Map the File --- 
    let input_file = File::open(input_path)
        .with_context(|| format!("Failed to open input ORB file: {}", input_path.display()))?;
    let mut header_reader = std::io::BufReader::new(&input_file); 
    let num_sequences = header_reader
        .read_u64::<LittleEndian>()
        .context("Failed to read number of sequences from ORB header")?;
    let total_decoded_length = header_reader
        .read_u64::<LittleEndian>()
        .context("Failed to read total decoded length from ORB header")?;
    let header_size: usize = 16;
    drop(header_reader); 

    println!("  Found {} sequences, total decoded length: {}. Mapping file...", num_sequences, total_decoded_length);
    // Map the entire file and wrap in Arc
    let mmap = Arc::new(unsafe { Mmap::map(&input_file)? });

    // Check length against Arc<Mmap>
    if (header_size + total_decoded_length as usize) > mmap.len() {
        anyhow::bail!(
            "ORB file header indicates total length {} + {} which exceeds file size {}",
            header_size, total_decoded_length, mmap.len()
        );
    }
    // Create slice temporarily for setup
    let sequence_data: &[u8] = &mmap[header_size..header_size + total_decoded_length as usize];
    println!("  Successfully memory-mapped sequence data ({} bytes).", sequence_data.len());

    // --- Parallel Chunk Processing --- 
    let (sender, receiver): (Sender<Result<HashMap<u64, usize>>>, Receiver<Result<HashMap<u64, usize>>>) = bounded(rayon::current_num_threads() * 2);
    let mut tasks_spawned = 0;

    let step_size = chunk_size.saturating_sub(overlap);
    if step_size == 0 {
        anyhow::bail!("Overlap size {} must be smaller than chunk size {}", overlap, chunk_size);
    }
    let num_chunks = if sequence_data.is_empty() { 
        0 
    } else {
        let full_steps = sequence_data.len().saturating_sub(chunk_size) / step_size;
        1 + full_steps + if sequence_data.len() > chunk_size + full_steps * step_size { 1 } else { 0 }
    };

    println!("Spawning {} hash analysis tasks for {} bytes...", num_chunks, sequence_data.len());

    // Remove mmap_ref
    //let mmap_ref = &mmap; 
    let k = min_len; // Rolling hash window size

    for i in 0..num_chunks {
        let start = i * step_size;
        let end = (start + chunk_size).min(sequence_data.len());
         if start >= end || end.saturating_sub(start) < k { continue; }

        //let chunk_slice = &sequence_data[start..end];

        let sender_clone = sender.clone();
        let current_chunk_index = i;
        let mmap_clone = Arc::clone(&mmap); // Clone Arc for closure

        rayon::spawn(move || {
            // Re-slice inside closure
            let sequence_data_in_task: &[u8] = &mmap_clone[header_size..header_size + total_decoded_length as usize];
            let chunk_slice_in_task = &sequence_data_in_task[start..end];

            let result = rolling_hash_chunk(chunk_slice_in_task, k);
            if let Err(e) = &result {
                eprintln!("Error hashing chunk {}: {}", current_chunk_index, e);
            }
            if sender_clone.send(result).is_err() {
                eprintln!("Warning: Failed to send hash result for chunk {} (receiver dropped).", current_chunk_index);
            }
        });
        tasks_spawned += 1;
    }

    drop(sender);

    // --- Collect and Merge Hash Results (remains the same) --- 
    println!("Waiting for {} hash analysis tasks to complete and merging results...", tasks_spawned);
    let mut global_hash_counts: HashMap<u64, usize> = HashMap::new(); // Final result map (hash -> max_count)

    for result in receiver {
        match result {
            Ok(local_hash_counts) => {
                for (hash_val, local_count) in local_hash_counts {
                    global_hash_counts
                        .entry(hash_val)
                        .and_modify(|global_count| {
                            *global_count = (*global_count).max(local_count);
                        })
                        .or_insert(local_count);
                }
            }
            Err(_) => {}
        }
    }

    println!(
        "Parallel rolling hash analysis complete. Found {} unique potential hashes (using max count heuristic).",
        global_hash_counts.len()
    );

    global_hash_counts.retain(|_, count| *count >= 2);
    println!("Filtered down to {} hashes occurring at least twice in at least one chunk.", global_hash_counts.len());
    println!("Note: These are hash values. Collisions are possible.");

    // No longer need explicit drop
    //drop(mmap_ref);
    Ok(global_hash_counts)
}

// --- Minimizer Sampling ---

/// Computes minimizers for a given sequence chunk.
/// Returns a list of (minimizer_hash, position) tuples.
/// Avoids adding consecutive identical minimizers.
fn get_minimizers_for_chunk(seq: &[u8], k: usize, w: usize) -> Result<Vec<(u64, usize)>> {
    if k == 0 || w == 0 || k > seq.len() || w > seq.len() || k + w - 1 > seq.len() {
        return Ok(Vec::new()); // Invalid parameters or sequence too short
    }

    let num_windows = seq.len() - (k + w - 1) + 1; // Total number of windows

    let mut minimizers = Vec::new();
    let mut last_min: Option<(u64, usize)> = None; // Store the last added minimizer (hash, pos)

    // Use a deque to efficiently find the minimum in the sliding window of k-mer hashes
    use std::collections::VecDeque;
    let mut window_hashes: VecDeque<(u64, usize)> = VecDeque::with_capacity(w); // Stores (hash, index_in_seq)

    // Process the first window to initialize
    for j in 0..w {
        let kmer_start = j;
        let kmer_end = kmer_start + k;
        if kmer_end > seq.len() { break; } // Ensure k-mer is within bounds
        let kmer = &seq[kmer_start..kmer_end];

        let mut hasher = XxHash64::default();
        hasher.write(kmer);
        let hash = hasher.finish();
        let kmer_pos = kmer_start; // Position of the k-mer start

        // Maintain deque property: monotonically increasing hashes (indices don't matter here)
        while !window_hashes.is_empty() && window_hashes.back().unwrap().0 >= hash {
            window_hashes.pop_back();
        }
        window_hashes.push_back((hash, kmer_pos));
    }

    // Add the minimizer for the first window if valid
    if let Some(&(min_hash, min_pos)) = window_hashes.front() {
         let current_min = (min_hash, min_pos);
         minimizers.push(current_min);
         last_min = Some(current_min);
    }

    // Slide the window across the rest of the sequence
    for i in 1..num_windows {
        // Index of the k-mer entering the window (relative to seq start)
        let entering_kmer_idx = i + w - 1;
        // Index of the k-mer leaving the window (relative to seq start)
        let leaving_kmer_idx = i - 1;

        // Remove the hash of the k-mer that slid out of the window from the front of the deque
        if !window_hashes.is_empty() && window_hashes.front().unwrap().1 == leaving_kmer_idx {
            window_hashes.pop_front();
        }

        // Calculate hash of the new k-mer entering the window
        let kmer_start = entering_kmer_idx;
        let kmer_end = kmer_start + k;
         if kmer_end > seq.len() { break; } // Ensure k-mer is within bounds
        let kmer = &seq[kmer_start..kmer_end];
        let mut hasher = XxHash64::default();
        hasher.write(kmer);
        let hash = hasher.finish();
        let kmer_pos = kmer_start;

        // Maintain deque property
        while !window_hashes.is_empty() && window_hashes.back().unwrap().0 >= hash {
            window_hashes.pop_back();
        }
        window_hashes.push_back((hash, kmer_pos));

        // The minimizer for the current window is at the front of the deque
         if let Some(&(min_hash, min_pos)) = window_hashes.front() {
             let current_min = (min_hash, min_pos);
             // Add if it's different from the last one added
             if Some(current_min) != last_min {
                 minimizers.push(current_min);
                 last_min = Some(current_min);
             }
         }
    }

    Ok(minimizers)
}


/// Finds frequent minimizer hashes using parallel chunking.
/// Returns a map of minimizer_hash -> total_count.
pub fn find_minimizers(input_path: &Path, k: usize, w: usize) -> Result<HashMap<u64, usize>> {
    if k == 0 || w == 0 {
        println!("Warning: k and w must be > 0 for minimizer sampling.");
        return Ok(HashMap::new());
    }

    // Use the generic processor
    let chunk_results: Vec<Vec<(u64, usize)>> = process_orb_file_in_chunks(
        input_path,
        k + w -1, // Effective "length" for overlap calculation
        "Minimizer Sampling",
        move |chunk| get_minimizers_for_chunk(chunk, k, w), // Pass k and w
    )?;

    // --- Merge Results ---
    println!("Merging minimizer results from {} chunks...", chunk_results.len());
    let mut global_minimizer_counts: HashMap<u64, usize> = HashMap::new();

    for local_minimizers in chunk_results {
        for (hash, _pos) in local_minimizers {
            *global_minimizer_counts.entry(hash).or_insert(0) += 1;
        }
    }

    println!(
        "Parallel minimizer sampling complete. Found {} unique minimizer hashes.",
        global_minimizer_counts.len()
    );

    // Optionally filter by count if needed
    // global_minimizer_counts.retain(|_, count| *count >= min_occurrence_threshold);

    Ok(global_minimizer_counts)
}


// --- MinHash Sketching ---

/// Generates a MinHash sketch (N smallest unique hashes) for a sequence chunk.
fn generate_chunk_sketch(chunk_data: &[u8], k: usize, sketch_size_n: usize) -> Result<Vec<u64>> {
    if k == 0 || sketch_size_n == 0 || k > chunk_data.len() {
        return Ok(Vec::new());
    }

    let mut heap: BinaryHeap<Reverse<u64>> = BinaryHeap::with_capacity(sketch_size_n + 1);
    let mut seen_hashes = HashSet::new(); // Track hashes seen *within this chunk*

    for i in 0..=chunk_data.len().saturating_sub(k) {
        let kmer = &chunk_data[i..i + k];
        let mut hasher = XxHash64::default();
        hasher.write(kmer);
        let hash = hasher.finish();

        // Process only if hash is new within this chunk
        if seen_hashes.insert(hash) {
            // Add the hash wrapped in Reverse to treat BinaryHeap as a min-heap for pushes,
            // but it stores the actual hash values.
            heap.push(Reverse(hash));
            // If heap exceeds N, remove the largest hash (due to Reverse wrapper)
            if heap.len() > sketch_size_n {
                heap.pop();
            }
        }
    }

    // Extract the hashes from the heap
    let mut sketch: Vec<u64> = heap.into_iter().map(|reversed_hash| reversed_hash.0).collect();
    // Sort the sketch for consistent output and to ensure the assertion check passes
    sketch.sort_unstable();
    Ok(sketch) // Return the chunk's sketch (N smallest unique hashes found in it, sorted)
}


/// Generates a MinHash sketch for the entire sequence using parallel chunking.
pub fn generate_sketch(input_path: &Path, k: usize, sketch_size_n: usize) -> Result<Vec<u64>> {
    if k == 0 || sketch_size_n == 0 {
         println!("Warning: k and sketch_size_n must be > 0 for MinHash sketching.");
        return Ok(Vec::new());
    }

    // Use the generic processor to get sketches from each chunk
    let chunk_sketches: Vec<Vec<u64>> = process_orb_file_in_chunks(
        input_path,
        k, // k-mer size for overlap
        "MinHash Sketching",
        move |chunk| generate_chunk_sketch(chunk, k, sketch_size_n), // Pass k and N
    )?;

    // --- Merge Sketches ---
    println!("Merging sketches from {} chunks...", chunk_sketches.len());
    // Collect all unique hashes from all chunk sketches
    let mut global_hashes = HashSet::new();
    for sketch in chunk_sketches {
        for hash in sketch {
            global_hashes.insert(hash);
        }
    }

    // Now, select the overall N smallest unique hashes from the global set
    let mut final_heap: BinaryHeap<Reverse<u64>> = BinaryHeap::with_capacity(sketch_size_n + 1);
    for hash in global_hashes {
        final_heap.push(Reverse(hash));
        if final_heap.len() > sketch_size_n {
            final_heap.pop();
        }
    }

    let mut final_sketch: Vec<u64> = final_heap.into_iter().map(|r| r.0).collect();
    final_sketch.sort_unstable(); // Sort for consistent output

    println!(
        "Parallel MinHash sketching complete. Generated final sketch of size {}.",
        final_sketch.len()
    );

    Ok(final_sketch)
}


// --- Tests ---
#[cfg(test)]
mod tests {
    use super::*;
    use std::path::{Path, PathBuf};
    use std::fs::{self, File};
    use byteorder::{LittleEndian, WriteBytesExt};
    use std::io::Write;
    use std::collections::HashSet;
    use anyhow::Result;

    // Helper to create a dummy ORB file for testing
    fn create_test_orb(path: &Path, data: &[u8]) -> Result<()> {
        let mut file = File::create(path)?;
        file.write_u64::<LittleEndian>(1)?; // num_sequences
        file.write_u64::<LittleEndian>(data.len() as u64)?; // total_decoded_length
        file.write_all(data)?;
        Ok(())
    }

    #[test]
    fn test_get_minimizers_simple() {
        // k=3, w=4. Sequence: ACGTACGTAC
        // Windows (k=3):
        // ACG TAC GTA CGT -> Window 1 (ACG, CGT, GTA, TAC): Hashes -> Min?
        // CGT GTA TAC GTA -> Window 2 (CGT, GTA, TAC, GTA): Hashes -> Min?
        // GTA TAC GTA ACA -> Window 3 (GTA, TAC, GTA, ACA): Hashes -> Min?
        // ...etc
        let seq = b"ACGTACGTAC";
        let k = 3;
        let w = 4; // Window of 4 k-mers
        // Expected windows of k-mers:
        // [ACG, CGT, GTA, TAC] -> k-mer start indices 0, 1, 2, 3
        // [CGT, GTA, TAC, GTA] -> k-mer start indices 1, 2, 3, 4
        // [GTA, TAC, GTA, ACA] -> k-mer start indices 2, 3, 4, 5
        // [TAC, GTA, ACA, CGT] -> k-mer start indices 3, 4, 5, 6 -> seq end is 9, last kmer starts at 7 (TAC)
        let minimizers = get_minimizers_for_chunk(seq, k, w).unwrap();
        // The exact hashes depend on the hasher, but we expect roughly len(seq) - (k+w-1) results,
        // potentially fewer due to duplicate removal.
        // n = 10, k = 3, w = 4. n - (k+w-1) + 1 = 10 - (3+4-1) + 1 = 10 - 6 + 1 = 5 windows
        // The number of minimizers depends on hash values and deduplication.
        // Let's just check it runs without error for now.
        assert!(!minimizers.is_empty());
        // Example manual check (hashes vary, logic check):
        // Win 1: k-mers at 0,1,2,3.
        // Win 2: k-mers at 1,2,3,4.
        // Win 3: k-mers at 2,3,4,5.
        // Win 4: k-mers at 3,4,5,6.
        // Win 5: k-mers at 4,5,6,7. <- Last window
    }

    #[test]
    fn test_get_minimizers_deduplication() {
        // Sequence designed so adjacent windows might have the same minimizer
        let seq = b"AGTAGTAGT"; // k=2, w=3
        // Windows (k=2):
        // [AG, GT, TA] -> k-mers at 0, 1, 2
        // [GT, TA, AG] -> k-mers at 1, 2, 3
        // [TA, AG, GT] -> k-mers at 2, 3, 4
        // [AG, GT, TA] -> k-mers at 3, 4, 5
        // [GT, TA, AG] -> k-mers at 4, 5, 6
        // [TA, AG, GT] -> k-mers at 5, 6, 7 <-- Last window (kmer GT starts at 7)
        // n=9, k=2, w=3. num_windows = 9 - (2+3-1) + 1 = 9 - 4 + 1 = 6
        let minimizers = get_minimizers_for_chunk(seq, 2, 3).unwrap();
        // We expect fewer than 6 minimizers if deduplication works.
        // Exact count depends on hash function.
        let mut unique_hashes = HashSet::new();
        let mut last_hash = None;
        for (h, _) in minimizers.iter() {
            if last_hash != Some(*h) {
                 unique_hashes.insert(*h);
                 last_hash = Some(*h);
            }
        }
         println!("Minimizers found (k=2, w=3): {:?}", minimizers);
         assert!(minimizers.len() <= 6); // Should be less than or equal to num windows
         // Check if consecutive entries are different
         for i in 0..minimizers.len().saturating_sub(1) {
             assert_ne!(minimizers[i], minimizers[i+1], "Consecutive minimizers should be different");
         }
    }

    #[test]
    fn test_generate_chunk_sketch_simple() {
        let seq = b"ACGTACGTACGT";
        let k = 3;
        let n = 5; // Get 5 smallest unique hashes
        let sketch = generate_chunk_sketch(seq, k, n).unwrap();
        assert!(sketch.len() <= n); // Sketch size should be at most N
        // Check if sorted
        assert!(sketch.windows(2).all(|w| w[0] <= w[1]));
        // Check uniqueness
        let unique_count = sketch.iter().collect::<HashSet<_>>().len();
        assert_eq!(sketch.len(), unique_count);
    }

     #[test]
    fn test_generate_chunk_sketch_duplicates() {
        let seq = b"AAAAAAAAAA"; // All k-mers are AAA
        let k = 3;
        let n = 5;
        let sketch = generate_chunk_sketch(seq, k, n).unwrap();
        assert_eq!(sketch.len(), 1); // Only one unique hash (AAA)
    }

    #[test]
    fn test_generate_chunk_sketch_n_larger_than_unique() {
        let seq = b"ACGTACGT"; // Fewer unique k-mers than N
        let k = 3;
        let n = 10;
        let sketch = generate_chunk_sketch(seq, k, n).unwrap();
        let unique_kmers = (0..=seq.len()-k).map(|i| &seq[i..i+k]).collect::<HashSet<_>>();
        assert_eq!(sketch.len(), unique_kmers.len()); // Should contain all unique hashes
    }

    // Integration test (requires creating a file)
    #[test]
    fn test_find_minimizers_integration() -> Result<()> {
        let test_file = PathBuf::from("test_minimizers.orb");
        let seq_data = b"GTGTGTACACACGTGTGTACACAC"; // Some repeating patterns
        create_test_orb(&test_file, seq_data)?;

        let k = 4;
        let w = 5;
        let minimizer_counts = find_minimizers(&test_file, k, w)?;

        assert!(!minimizer_counts.is_empty());
        // Check if counts seem reasonable (e.g., some hashes should have count > 1)
        assert!(minimizer_counts.values().any(|&count| count > 1));

        fs::remove_file(test_file)?; // Clean up
        Ok(())
    }

    #[test]
    fn test_generate_sketch_integration() -> Result<()> {
        let test_file = PathBuf::from("test_sketch.orb");
        let seq_data = b"ACGTACGTACGTACGTACGTACGTACGTACGT";
        create_test_orb(&test_file, seq_data)?;

        let k = 5;
        let n = 10; // Sketch size
        let sketch = generate_sketch(&test_file, k, n)?;

        assert!(sketch.len() <= n);
        // Check uniqueness and sorted order
         let unique_count = sketch.iter().collect::<HashSet<_>>().len();
         assert_eq!(sketch.len(), unique_count);
         assert!(sketch.windows(2).all(|w| w[0] <= w[1]));


        fs::remove_file(test_file)?; // Clean up
        Ok(())
    }
}


// --- Generic Parallel Chunk Processing ---

/// Generic helper to process a memory-mapped ORB file in parallel chunks.
///
/// Args:
///    input_path: Path to the input ORB file.
///    min_len_or_k: Minimum length for repeats or k-mer size for hashing/sketching.
///    context_str: A string describing the operation (for logging).
///    chunk_processor: A closure that takes a chunk slice and returns a Result<T>.
///
/// Returns:
///    A vector containing the results from each successfully processed chunk.
fn process_orb_file_in_chunks<T: Send + 'static>(
    input_path: &Path,
    min_len_or_k: usize, // Used for overlap calculation primarily
    context_str: &str,
    chunk_processor: impl Fn(&[u8]) -> Result<T> + Send + Sync + Copy + 'static,
) -> Result<Vec<T>> {

    let mut sys = System::new_all();
    sys.refresh_memory();
    let total_memory = sys.total_memory();
    // Adjust memory divisor based on typical needs of the algorithms
    let memory_divisor = match context_str {
        "Suffix Array" => 150, // Suffix array can be memory intensive
        _ => 100, // Default for hashing/sketching
    };
    let chunk_size = (total_memory / memory_divisor).max(1024 * 1024) as usize;
    // Ensure overlap is sufficient for k-mers or min_len repeats
    let overlap = (chunk_size / 20).max(min_len_or_k * 2 + 1) as usize;

    println!(
        "System total memory: {} bytes. Using chunk size: {} bytes, overlap: {} bytes for {}",
        total_memory, chunk_size, overlap, context_str
    );
    println!(
        "Processing memory-mapped ORB file: {} ({}) in parallel...",
        input_path.display(), context_str
    );

    // --- Memory Map the File ---
    let input_file = File::open(input_path)
        .with_context(|| format!("Failed to open input ORB file: {}", input_path.display()))?;
    let mut header_reader = std::io::BufReader::new(&input_file);
    let num_sequences = header_reader
        .read_u64::<LittleEndian>()
        .context("Failed to read number of sequences from ORB header")?;
    let total_decoded_length = header_reader
        .read_u64::<LittleEndian>()
        .context("Failed to read total decoded length from ORB header")?;
    let header_size: usize = 16;
    drop(header_reader);

    println!("  Found {} sequences, total decoded length: {}. Mapping file...", num_sequences, total_decoded_length);

    let mmap = Arc::new(unsafe { Mmap::map(&input_file)? });

    if (header_size + total_decoded_length as usize) > mmap.len() {
        anyhow::bail!(
            "ORB file header indicates total length {} + {} which exceeds file size {}",
            header_size, total_decoded_length, mmap.len()
        );
    }
    // Check if there's any data to process
    if total_decoded_length == 0 {
        println!("  Warning: ORB file contains no sequence data (total_decoded_length is 0). Returning empty results.");
        return Ok(Vec::new());
    }

    let sequence_data_len = total_decoded_length as usize;
    println!("  Successfully memory-mapped sequence data ({} bytes).", sequence_data_len);

    // --- Parallel Chunk Processing Setup ---
    let (sender, receiver): (Sender<Result<T>>, Receiver<Result<T>>) = bounded(rayon::current_num_threads() * 2);
    let mut tasks_spawned = 0;

    let step_size = chunk_size.saturating_sub(overlap);
    if step_size == 0 {
        anyhow::bail!("Overlap size {} must be smaller than chunk size {}", overlap, chunk_size);
    }

    // Correct calculation for number of chunks
     let num_chunks = if sequence_data_len == 0 {
         0
     } else if sequence_data_len <= chunk_size {
         1 // Only one chunk if data fits entirely
     } else {
         // Calculate how many full steps we can take *after* the first chunk
         let remaining_data = sequence_data_len - chunk_size;
         let full_steps = remaining_data / step_size;
         // Total chunks = initial chunk + full steps + potentially one final partial chunk
         1 + full_steps + if remaining_data % step_size > 0 { 1 } else { 0 }
     };


    println!("Spawning {} analysis tasks for {} bytes using {}...", num_chunks, sequence_data_len, context_str);

    // --- Spawn Tasks ---
    for i in 0..num_chunks {
        let start = i * step_size;
        // Ensure end does not exceed the actual mapped data length
        let end = (start + chunk_size).min(sequence_data_len);

        // Skip invalid chunks (e.g., if start >= end due to calculation edge cases)
         if start >= end {
             // This might happen if sequence_data_len is very small or zero
             println!("  Skipping invalid chunk index {}: start={}, end={}", i, start, end);
             continue;
         }

        let sender_clone = sender.clone();
        let mmap_clone = Arc::clone(&mmap); // Clone Arc for the closure

        rayon::spawn(move || {
            // Re-slice inside the closure using the cloned Arc and calculated indices
            let sequence_data_in_task: &[u8] = &mmap_clone[header_size..header_size + sequence_data_len];
            // Boundary check: Ensure end does not exceed sequence_data_in_task length
             let actual_end = end.min(sequence_data_in_task.len());
             if start >= actual_end {
                 // Should not happen if previous check passed, but good for robustness
                 eprintln!("Error: Invalid slice range in task: start={}, end={}", start, actual_end);
                 let _ = sender_clone.send(Err(anyhow::anyhow!("Invalid slice range in task")));
                 return;
             }
            let chunk_slice_in_task = &sequence_data_in_task[start..actual_end];

            let result = chunk_processor(chunk_slice_in_task);
            if let Err(e) = &result {
                eprintln!("Error processing chunk {}: {}", i, e);
            }
            if sender_clone.send(result).is_err() {
                 eprintln!("Warning: Failed to send result for chunk {} (receiver dropped).", i);
             }
        });
        tasks_spawned += 1;
    }

    // Drop the original sender so the receiver knows when all tasks are done
    drop(sender);

    // --- Collect Results ---
    println!("Waiting for {} analysis tasks to complete and collecting results...", tasks_spawned);
    let mut collected_results = Vec::with_capacity(tasks_spawned);
    for result in receiver {
        match result {
            Ok(chunk_result) => collected_results.push(chunk_result),
            Err(e) => eprintln!("Received error from worker thread: {}", e), // Log errors, but don't push them
        }
    }

    println!("Collected results from {} successful tasks.", collected_results.len());
    Ok(collected_results)
}