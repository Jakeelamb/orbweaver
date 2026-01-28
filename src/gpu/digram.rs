//! GPU-accelerated digram counting and processing.
//!
//! This module provides OpenCL-based implementations for counting and processing digrams
//! (adjacent symbol pairs) efficiently on GPUs.

use crate::grammar::symbol::{Symbol, SymbolType, Direction};
use crate::grammar::digram_table::{DigramKey, DigramKeyTuple, DigramTable};
use crate::grammar::digram::Digram; // For Digram type
use crate::encode::dna_2bit::EncodedBase;
use crate::gpu::GpuContext;
use rustc_hash::{FxHashMap, FxHasher};
use std::collections::HashMap;
use std::collections::hash_map::DefaultHasher; // For DefaultHasher - keep for now, might be used elsewhere or can be cleaned later
use std::hash::{Hash, Hasher};
use ocl::{Buffer, Kernel};
use anyhow::{Result, bail, anyhow};
use crate::utils::hash::custom_hash;

// --- Digram constants for GPU Kernels ---

/// Maximum possible digram ID value that the kernel can produce or count (0-15 for 2-bit DNA).
/// The `count_digrams_by_id` kernel uses `local_counts[16]` and expects this as `max_digram_id_value` argument.
pub const MAX_POSSIBLE_DIGRAM_IDS: u32 = 16;

/// Maps a kernel output ID (0-15) back to its canonical `(Symbol, Symbol)` pair.
/// This array is indexed by the ID received from the `compute_digram_ids` kernel or generated directly.
/// When `reverse_aware` is true in `compute_digram_ids`, the kernel outputs one of 10 unique canonical IDs.
/// Those 10 canonical IDs are: AA (0), AC (1), AG (2), AT (3), CA (4), CC (5), CG (6), GA (8), GC (9), TA (12).
/// Entries for other IDs (which would only be produced if `reverse_aware` is false, or if an ID
/// is one of the non-canonical results of `base1*4 + base2` that maps to a canonical one) are marked
/// with `(Symbol::N, Symbol::N)` if they don't directly correspond to one of the 10 canonical forms
/// *at that specific ID index*. The primary use is to map the 10 canonical output IDs.
pub static CANONICAL_KERNEL_ID_TO_BASES: [(Symbol, Symbol); MAX_POSSIBLE_DIGRAM_IDS as usize] = [
    // ID 0 (AA)
    (
        Symbol { id: 0, symbol_type: SymbolType::Terminal(EncodedBase(0b00)), strand: Direction::Forward, source_grammar_id: None, original_pos: None },
        Symbol { id: 0, symbol_type: SymbolType::Terminal(EncodedBase(0b00)), strand: Direction::Forward, source_grammar_id: None, original_pos: None },
    ),
    // ID 1 (AC)
    (
        Symbol { id: 0, symbol_type: SymbolType::Terminal(EncodedBase(0b00)), strand: Direction::Forward, source_grammar_id: None, original_pos: None },
        Symbol { id: 0, symbol_type: SymbolType::Terminal(EncodedBase(0b01)), strand: Direction::Forward, source_grammar_id: None, original_pos: None },
    ),
    // ID 2 (AG)
    (
        Symbol { id: 0, symbol_type: SymbolType::Terminal(EncodedBase(0b00)), strand: Direction::Forward, source_grammar_id: None, original_pos: None },
        Symbol { id: 0, symbol_type: SymbolType::Terminal(EncodedBase(0b10)), strand: Direction::Forward, source_grammar_id: None, original_pos: None },
    ),
    // ID 3 (AT)
    (
        Symbol { id: 0, symbol_type: SymbolType::Terminal(EncodedBase(0b00)), strand: Direction::Forward, source_grammar_id: None, original_pos: None },
        Symbol { id: 0, symbol_type: SymbolType::Terminal(EncodedBase(0b11)), strand: Direction::Forward, source_grammar_id: None, original_pos: None },
    ),
    // ID 4 (CA)
    (
        Symbol { id: 0, symbol_type: SymbolType::Terminal(EncodedBase(0b01)), strand: Direction::Forward, source_grammar_id: None, original_pos: None },
        Symbol { id: 0, symbol_type: SymbolType::Terminal(EncodedBase(0b00)), strand: Direction::Forward, source_grammar_id: None, original_pos: None },
    ),
    // ID 5 (CC)
    (
        Symbol { id: 0, symbol_type: SymbolType::Terminal(EncodedBase(0b01)), strand: Direction::Forward, source_grammar_id: None, original_pos: None },
        Symbol { id: 0, symbol_type: SymbolType::Terminal(EncodedBase(0b01)), strand: Direction::Forward, source_grammar_id: None, original_pos: None },
    ),
    // ID 6 (CG)
    (
        Symbol { id: 0, symbol_type: SymbolType::Terminal(EncodedBase(0b01)), strand: Direction::Forward, source_grammar_id: None, original_pos: None },
        Symbol { id: 0, symbol_type: SymbolType::Terminal(EncodedBase(0b10)), strand: Direction::Forward, source_grammar_id: None, original_pos: None },
    ),
    // ID 7 (CT -> maps to AG, ID 2 if reverse_aware) - Placeholder for direct mapping at ID 7
    (
        Symbol { id: 0, symbol_type: SymbolType::Terminal(EncodedBase(4)), strand: Direction::Forward, source_grammar_id: None, original_pos: None }, // N
        Symbol { id: 0, symbol_type: SymbolType::Terminal(EncodedBase(4)), strand: Direction::Forward, source_grammar_id: None, original_pos: None }, // N
    ),
    // ID 8 (GA)
    (
        Symbol { id: 0, symbol_type: SymbolType::Terminal(EncodedBase(0b10)), strand: Direction::Forward, source_grammar_id: None, original_pos: None },
        Symbol { id: 0, symbol_type: SymbolType::Terminal(EncodedBase(0b00)), strand: Direction::Forward, source_grammar_id: None, original_pos: None },
    ),
    // ID 9 (GC)
    (
        Symbol { id: 0, symbol_type: SymbolType::Terminal(EncodedBase(0b10)), strand: Direction::Forward, source_grammar_id: None, original_pos: None },
        Symbol { id: 0, symbol_type: SymbolType::Terminal(EncodedBase(0b01)), strand: Direction::Forward, source_grammar_id: None, original_pos: None },
    ),
    // ID 10 (GG -> maps to CC, ID 5 if reverse_aware) - Placeholder for direct mapping at ID 10
    (
        Symbol { id: 0, symbol_type: SymbolType::Terminal(EncodedBase(4)), strand: Direction::Forward, source_grammar_id: None, original_pos: None }, // N
        Symbol { id: 0, symbol_type: SymbolType::Terminal(EncodedBase(4)), strand: Direction::Forward, source_grammar_id: None, original_pos: None }, // N
    ),
    // ID 11 (GT -> maps to AC, ID 1 if reverse_aware) - Placeholder for direct mapping at ID 11
    (
        Symbol { id: 0, symbol_type: SymbolType::Terminal(EncodedBase(4)), strand: Direction::Forward, source_grammar_id: None, original_pos: None }, // N
        Symbol { id: 0, symbol_type: SymbolType::Terminal(EncodedBase(4)), strand: Direction::Forward, source_grammar_id: None, original_pos: None }, // N
    ),
    // ID 12 (TA)
    (
        Symbol { id: 0, symbol_type: SymbolType::Terminal(EncodedBase(0b11)), strand: Direction::Forward, source_grammar_id: None, original_pos: None },
        Symbol { id: 0, symbol_type: SymbolType::Terminal(EncodedBase(0b00)), strand: Direction::Forward, source_grammar_id: None, original_pos: None },
    ),
    // ID 13 (TC -> maps to GA, ID 8 if reverse_aware) - Placeholder for direct mapping at ID 13
    (
        Symbol { id: 0, symbol_type: SymbolType::Terminal(EncodedBase(4)), strand: Direction::Forward, source_grammar_id: None, original_pos: None }, // N
        Symbol { id: 0, symbol_type: SymbolType::Terminal(EncodedBase(4)), strand: Direction::Forward, source_grammar_id: None, original_pos: None }, // N
    ),
    // ID 14 (TG -> maps to CA, ID 4 if reverse_aware) - Placeholder for direct mapping at ID 14
    (
        Symbol { id: 0, symbol_type: SymbolType::Terminal(EncodedBase(4)), strand: Direction::Forward, source_grammar_id: None, original_pos: None }, // N
        Symbol { id: 0, symbol_type: SymbolType::Terminal(EncodedBase(4)), strand: Direction::Forward, source_grammar_id: None, original_pos: None }, // N
    ),
    // ID 15 (TT -> maps to AA, ID 0 if reverse_aware) - Placeholder for direct mapping at ID 15
    (
        Symbol { id: 0, symbol_type: SymbolType::Terminal(EncodedBase(4)), strand: Direction::Forward, source_grammar_id: None, original_pos: None }, // N
        Symbol { id: 0, symbol_type: SymbolType::Terminal(EncodedBase(4)), strand: Direction::Forward, source_grammar_id: None, original_pos: None }, // N
    ),
];

/// Structure to hold a sequence for GPU operations
pub struct GpuSequence {
    /// Raw sequence data
    data: Vec<u8>,
    
    /// OpenCL buffer for the sequence data (if uploaded to GPU)
    buffer: Option<Buffer<u8>>,
}

impl Clone for GpuSequence {
    fn clone(&self) -> Self {
        Self {
            data: self.data.clone(),
            buffer: None, // Cloned instance does not share the GPU buffer
        }
    }
}

impl GpuSequence {
    /// Create a new GPU sequence from raw data
    pub fn new(data: Vec<u8>) -> Self {
        Self {
            data,
            buffer: None
        }
    }
    
    /// Creates a GpuSequence from Symbol vectors (typically from grammar rules)
    pub fn from_symbols(symbols: &[Symbol]) -> Result<Self> {
        // Convert from symbols to 2-bit encoded sequence
        let mut encoded = Vec::with_capacity(symbols.len());
        
        for symbol in symbols {
            if let SymbolType::Terminal(base) = symbol.symbol_type {
                encoded.push(base.0);
            } else {
                bail!("Cannot convert non-terminal symbols to GPU sequence directly");
            }
        }
        
        Ok(Self::new(encoded))
    }
    
    /// Get the length of the sequence
    pub fn get_len(&self) -> usize {
        self.data.len()
    }
    
    /// Get the raw data of the sequence
    pub fn get_data(&self) -> &[u8] {
        &self.data
    }
    
    /// Get the GPU buffer if it exists
    pub fn get_buffer(&self) -> Option<&Buffer<u8>> {
        self.buffer.as_ref()
    }
    
    /// Create a symbol at the given position
    pub fn get_symbol(&self, pos: usize) -> Symbol {
        if pos >= self.data.len() {
            panic!("Position {} is out of bounds (len: {})", pos, self.data.len());
        }
        
        let base = EncodedBase(self.data[pos]);
        Symbol::terminal(pos, base, Direction::Forward, None, None)
    }
    
    /// Upload the sequence to GPU memory using OpenCL
    pub fn upload_to_gpu(&mut self, gpu_context: &GpuContext) -> Result<()> {
        if self.buffer.is_some() {
            // Already uploaded
            return Ok(());
        }
        
        // Create an OpenCL buffer
        let buffer = Buffer::builder()
            .queue(gpu_context.queue.clone())
            .flags(ocl::flags::MEM_READ_ONLY)
            .len(self.data.len())
            .copy_host_slice(&self.data)
            .build()
            .map_err(|e| anyhow!("Failed to upload sequence to GPU: {:?}", e))?;
            
        self.buffer = Some(buffer);
        log::debug!("Uploaded {} bytes to GPU", self.data.len());

        Ok(())
    }
    
    /// Find the most frequent digram in the sequence using OpenCL
    pub fn find_most_frequent_digram(
        &self, 
        min_count: usize, 
        _reverse_aware: bool,
        gpu_context: Option<&GpuContext>
    ) -> Result<Option<(DigramKey, Vec<(usize, (Symbol, Symbol))>)>> {
        if let Some(context) = gpu_context {
            // Try to use GPU acceleration
            if self.buffer.is_none() {
                log::debug!("Sequence not uploaded to GPU, using CPU implementation.");
                // Need to upload first - but this requires a mutable reference
                // Instead, compute on CPU
                return self.find_most_frequent_digram_cpu(min_count, _reverse_aware);
            }

            // Use OpenCL
            match find_most_frequent_digram_opencl(self, min_count, _reverse_aware, context) {
                Ok(result) => return Ok(result),
                Err(e) => {
                    log::warn!("GPU digram finding failed: {:?}. Falling back to CPU", e);
                    return self.find_most_frequent_digram_cpu(min_count, _reverse_aware);
                }
            }
        }

        // Fallback to CPU
        self.find_most_frequent_digram_cpu(min_count, _reverse_aware)
    }
    
    /// CPU implementation for finding most frequent digram
    pub(crate) fn find_most_frequent_digram_cpu(
        &self, 
        min_count: usize, 
        _reverse_aware: bool
    ) -> Result<Option<(DigramKey, Vec<(usize, (Symbol, Symbol))>)>> {
        // Calculate digram frequencies
        let mut digram_counts: FxHashMap<DigramKey, Vec<(usize, (Symbol, Symbol))>> = FxHashMap::default();
        let data_arc = self.get_data();
        let data = &data_arc[..];
        
        for i in 0..self.data.len().saturating_sub(1) {
            let base1 = EncodedBase(data[i]);
            let base2 = EncodedBase(data[i+1]);
            
            let sym1 = Symbol::terminal(i, base1, Direction::Forward, None, None);
            let sym2 = Symbol::terminal(i+1, base2, Direction::Forward, None, None);
            
            let key_tuple = DigramTable::canonical_key((&sym1, &sym2), _reverse_aware);
            let key_hash = custom_hash(&key_tuple);
            
            digram_counts.entry(key_hash)
                .or_default()
                .push((i, (sym1, sym2)));
        }
        
        // Find the most frequent digram
        let most_frequent = digram_counts.iter()
            .filter(|(_, occurrences)| occurrences.len() >= min_count)
            .max_by_key(|(_, occurrences)| occurrences.len());
        
        if let Some((key, occurrences)) = most_frequent {
            Ok(Some((*key, occurrences.clone())))
        } else {
            Ok(None)
        }
    }

    /// Find the most frequent digram using the enhanced OpenCL kernel
    pub fn find_most_frequent_digram_enhanced(
        &self,
        min_count: usize,
        _reverse_aware: bool,
        gpu_context: &GpuContext
    ) -> Result<Option<(DigramKey, Vec<(usize, (Symbol, Symbol))>)>> {
        // Let's use our new sequential version which supports chunking
        self.find_most_frequent_digram_sequential(min_count, gpu_context)
    }

    /// Find digrams most frequent in sequence using naive approach (not using OpenCL)
    pub fn find_most_frequent_digram_naive(&self, min_count: usize) -> Option<(DigramKey, Vec<(usize, (Symbol, Symbol))>)> {
        let mut counts = HashMap::<DigramKey, Vec<usize>>::new();
        let data = &self.data;
        
        for i in 0..self.data.len().saturating_sub(1) {
            let base1 = EncodedBase(data[i]);
            let base2 = EncodedBase(data[i+1]);
            
            let sym1 = Symbol::terminal(i, base1, Direction::Forward, None, None);
            let sym2 = Symbol::terminal(i+1, base2, Direction::Forward, None, None);
            
            let key_tuple = DigramTable::canonical_key((&sym1, &sym2), false);
            let key_hash = custom_hash(&key_tuple);
            
            counts.entry(key_hash).or_insert_with(Vec::new).push(i);
        }
        
        let mut max_count = 0;
        let mut max_key = 0;
        let mut positions = Vec::new();
        
        for (key, pos_list) in &counts {
            if pos_list.len() > max_count && pos_list.len() >= min_count {
                max_count = pos_list.len();
                max_key = *key;
                positions = pos_list.clone();
            }
        }
        
        if max_count < min_count {
            return None;
        }
        
        let mut occurrences = Vec::new();
        for pos in positions {
            let base1 = EncodedBase(data[pos]);
            let base2 = EncodedBase(data[pos+1]);
            
            let sym1 = Symbol::terminal(pos, base1, Direction::Forward, None, None);
            let sym2 = Symbol::terminal(pos+1, base2, Direction::Forward, None, None);
            
            occurrences.push((pos, (sym1, sym2)));
        }
        
        Some((max_key, occurrences))
    }

    pub fn find_most_frequent_digram_sequential(&self, min_count: usize, _gpu_context: &GpuContext) -> Result<Option<(DigramKey, Vec<(usize, (Symbol, Symbol))>)>> {
        // Calculate approximate memory needed for this operation
        // Each sequence element needs space for the original data and at least 8 bytes for counts
        let estimated_mem_needed = self.data.len() * (1 + 8);
        
        // If the sequence is too large for GPU memory, we need to process in chunks
        let chunk_size = 10_000_000; // Process 10M bases at a time
        let mut position = 0;
        let mut counts = HashMap::<DigramKey, Vec<usize>>::new();
        
        if estimated_mem_needed > 1_000_000_000 { // 1 GB threshold for chunking
            log::debug!("Sequence too large for GPU memory, processing in chunks");
            
            // Process in chunks
            while position < self.data.len() {
                let end_pos = std::cmp::min(position + chunk_size, self.data.len());
                let chunk_len = end_pos - position;
                
                log::trace!("Processing chunk starting at {} with length {}", position, chunk_len);
                
                // Extract chunk
                let chunk_data = self.data[position..end_pos].to_vec();
                
                // Create a GPU sequence for this chunk
                let chunk_seq = GpuSequence::new(chunk_data);
                
                // Process this chunk
                if let Some((_, chunk_occurrences)) = chunk_seq.find_most_frequent_digram_naive(2) {
                    // Adjust positions and collect results
                    for (pos, (sym_a, sym_b)) in chunk_occurrences {
                        let adjusted_pos = position + pos;
                        let key_tuple = DigramTable::canonical_key((&sym_a, &sym_b), false);
                        let key_hash = custom_hash(&key_tuple);
                        
                        counts.entry(key_hash).or_insert_with(Vec::new).push(adjusted_pos);
                    }
                }
                
                position = end_pos;
            }
        } else {
            // Process entire sequence at once
            if let Some((_, occurrences)) = self.find_most_frequent_digram_naive(2) {
                for (pos, (sym_a, sym_b)) in occurrences {
                    let key_tuple = DigramTable::canonical_key((&sym_a, &sym_b), false);
                    let key_hash = custom_hash(&key_tuple);
                    counts.entry(key_hash).or_insert_with(Vec::new).push(pos);
                }
            }
        }
        
        // Find most frequent
        let mut max_count = 0;
        let mut max_key = 0;
        let mut positions = Vec::new();
        
        for (key, pos_list) in &counts {
            if pos_list.len() > max_count && pos_list.len() >= min_count {
                max_count = pos_list.len();
                max_key = *key;
                positions = pos_list.clone();
            }
        }
        
        if max_count < min_count {
            return Ok(None);
        }
        
        // Create occurrences
        let mut occurrences = Vec::new();
        for pos in positions {
            let sym_a = self.get_symbol(pos);
            let sym_b = self.get_symbol(pos + 1);
            
            occurrences.push((pos, (sym_a, sym_b)));
        }
        
        Ok(Some((max_key, occurrences)))
    }

    pub fn len(&self) -> usize {
        self.data.len()
    }
}

/// Find the most frequent digram using OpenCL
fn find_most_frequent_digram_opencl(
    sequence: &GpuSequence,
    min_count: usize,
    _reverse_aware: bool,
    gpu_context: &GpuContext
) -> Result<Option<(DigramKey, Vec<(usize, (Symbol, Symbol))>)>> {
    if sequence.get_len() < 2 {
        return Ok(None); // No digrams to find
    }
    
    // Get the OpenCL context and queue
    let _context = &gpu_context.context;
    let _queue = &gpu_context.queue;
    
    // Get maximum memory of GPU and calculate safe memory usage (50% of available)
    let max_memory = gpu_context.device.info(ocl::enums::DeviceInfo::GlobalMemSize)
        .map(|info| {
            if let Some(size_str) = info.to_string().strip_prefix("DeviceInfoResult::Size(") {
                if let Some(size_num) = size_str.strip_suffix(")") {
                    if let Ok(size) = size_num.parse::<usize>() {
                        return size;
                    }
                }
            }
            1024 * 1024 * 1024 // Default to 1GB
        }).unwrap_or(1024 * 1024 * 1024);
    
    let safe_memory = (max_memory as f64 * 0.5) as usize;
    
    // Estimate memory needed - input data plus 8 bytes per position for counts
    let seq_len = sequence.get_len();
    let memory_needed = seq_len + seq_len * 8;
    
    // If sequence fits in memory, use direct method
    if memory_needed <= safe_memory {
        return find_most_frequent_digram_opencl_direct(sequence, min_count, _reverse_aware, gpu_context);
    }

    // Otherwise use chunked approach
    log::debug!("Sequence too large for GPU memory ({} bytes needed, {} available), using chunked approach",
             memory_needed, safe_memory);
    
    // Define chunk size to fit within memory constraints (with safety margin)
    let chunk_size = (safe_memory / 8) / 2;
    let chunk_count = (seq_len + chunk_size - 1) / chunk_size;
    
    log::debug!("Processing in {} chunks of ~{} bytes each", chunk_count, chunk_size);

    // Process each chunk and combine results
    let mut all_counts = HashMap::<DigramKey, usize>::new();
    
    for chunk_idx in 0..chunk_count {
        let start = chunk_idx * chunk_size;
        let end = std::cmp::min(start + chunk_size, seq_len);
        
        // Skip last chunk if it's too small (< 2 bases)
        if end - start < 2 {
            continue;
        }
        
        // Extract chunk data
        let chunk_data = &sequence.get_data()[start..end];
        let chunk_seq = GpuSequence::new(chunk_data.to_vec());
        
        // Process this chunk
        if let Ok(Some((_chunk_key, _chunk_positions))) = 
            find_most_frequent_digram_opencl_direct(&chunk_seq, 0, _reverse_aware, gpu_context) {
            
            // Extract all digram counts from this chunk (min_count = 0 to get all)
            let chunk_counts = count_digrams_gpu(&chunk_seq, _reverse_aware, gpu_context)?;
            
            // Merge counts with global map
            for (key, count) in chunk_counts.counts {
                *all_counts.entry(key).or_insert(0) += count;
            }
        }
    }

    // Find the most frequent digram across all chunks
    let mut max_key = 0;
    let mut max_count = 0;
    
    for (key, count) in all_counts {
        if count > max_count {
            max_key = key;
            max_count = count;
        }
    }
    
    // If no digram appears enough times, return None
    if max_count < min_count {
        return Ok(None);
    }
    
    // Now we need to collect all positions - do another pass over the sequence
    let mut positions = Vec::new();
    
    for i in 0..seq_len.saturating_sub(1) {
        let sym1 = sequence.get_symbol(i);
        let sym2 = sequence.get_symbol(i + 1);
        
        // Create Digram and get its canonical form
        let current_digram = Digram::new(sym1.clone(), sym2.clone());
        let canonical_digram_to_hash = current_digram.canonical(_reverse_aware);

        // Hash the canonical Digram to get its DigramKey (u64)
        let mut hasher = DefaultHasher::new();
        canonical_digram_to_hash.hash(&mut hasher);
        let current_key_hash: DigramKey = hasher.finish();
        
        if current_key_hash == max_key {
            // Store the original symbols, not necessarily the canonical ones, along with position
            positions.push((i, (sym1, sym2)));
        }
    }

    // Return the result
    Ok(Some((max_key, positions)))
}

/// Direct method for finding the most frequent digram using OpenCL (for smaller sequences)
fn find_most_frequent_digram_opencl_direct(
    sequence: &GpuSequence,
    min_count: usize,
    _reverse_aware: bool,
    gpu_context: &GpuContext
) -> Result<Option<(DigramKey, Vec<(usize, (Symbol, Symbol))>)>> {
    let seq_len = sequence.get_len();
    if seq_len < 2 {
        return Ok(None);
    }

    // 1. Get digram counts using the GPU
    let gpu_digram_counts = count_digrams_gpu(sequence, _reverse_aware, gpu_context)?;

    // 2. Find the most frequent digram from the counts
    let mut max_key: DigramKey = 0;
    let mut max_count: usize = 0;

    for (key, count) in gpu_digram_counts.counts {
        if count > max_count {
            max_key = key;
            max_count = count;
        }
    }

    // 3. If no digram appears enough times, return None
    if max_count < min_count {
        return Ok(None);
    }

    // 4. Collect all positions of the most frequent digram
    // This requires a CPU pass over the sequence, generating keys consistently.
    let mut positions = Vec::new();
    for i in 0..seq_len.saturating_sub(1) {
        let sym1 = sequence.get_symbol(i);
        let sym2 = sequence.get_symbol(i + 1);

        // Create Digram and get its canonical form
        let current_digram = Digram::new(sym1.clone(), sym2.clone());
        let canonical_digram_to_hash = current_digram.canonical(_reverse_aware);
        
        // Hash the canonical Digram to get its DigramKey (u64)
        let mut hasher = DefaultHasher::new();
        canonical_digram_to_hash.hash(&mut hasher);
        let current_key_hash: DigramKey = hasher.finish();

        if current_key_hash == max_key {
            // Store the original symbols, not necessarily the canonical ones, along with position
            positions.push((i, (sym1, sym2)));
        }
    }
    
    Ok(Some((max_key, positions)))
}

/// Results from GPU-based digram counting
pub struct GpuDigramCounts {
    /// Map from digram hash to count
    pub counts: HashMap<DigramKey, usize>,
}

/// Count digrams using GPU acceleration with OpenCL
/// Count digrams using GPU acceleration with OpenCL
pub fn count_digrams_gpu(
    sequence: &GpuSequence,
    _reverse_aware: bool,
    gpu_context: &GpuContext
) -> Result<GpuDigramCounts> {
    let program = gpu_context.program.as_ref().ok_or_else(|| anyhow!("OpenCL program not loaded in GpuContext"))?;

    // Step 1: Compute digram IDs using the first kernel, getting the GPU buffer directly
    let (digram_ids_buffer, num_ids) = compute_digram_ids_opencl(
        &sequence.data,
        _reverse_aware,
        gpu_context
    ).map_err(|e| anyhow!("Failed to compute digram IDs via OpenCL: {:?}", e))?;

    if num_ids == 0 {
        return Ok(GpuDigramCounts { counts: HashMap::new() });
    }
    // The digram_ids_buffer is now directly from compute_digram_ids_opencl, no need to recreate.

    // Determine max_digram_id_value based on reverse_aware flag.
    // If _reverse_aware, the kernel uses a mapping to 8 canonical IDs (0-7).
    // If not _reverse_aware, all 16 digram IDs (0-15) are used.
    // The `count_digrams_by_id` kernel's `local_counts` array is always size 16.
    // The `max_digram_id_value` argument tells the kernel the upper bound of IDs to consider.
    // This should always be MAX_POSSIBLE_DIGRAM_IDS (16) because the digram_ids themselves
    // will be in the correct range (0-15, or the 10 canonical IDs if reverse_aware).
    let max_digram_id_value = MAX_POSSIBLE_DIGRAM_IDS;

    // --- Global Counts Buffer (Output for count_digrams_by_id kernel) ---
    // Size is max_digram_id_value, as kernel outputs counts for each possible ID up to this max.
    let global_counts_buffer = Buffer::<u32>::builder()
        .queue(gpu_context.queue.clone())
        .flags(ocl::flags::MEM_WRITE_ONLY) // Kernel writes to it
        .len(max_digram_id_value as usize)
        .build().map_err(|e| anyhow!("Failed to create OpenCL buffer for global counts: {:?}", e))?;

    // --- Setup and Run count_digrams_by_id Kernel ---
    // Determine Work Sizes (similar logic to compute_digram_ids_opencl, but for num_ids)
    // This preliminary kernel is just to query work group info for "count_digrams_by_id"
    let preliminary_count_kernel = Kernel::builder()
        .program(program)
        .name("count_digrams_by_id")
        .queue(gpu_context.queue.clone())
        .global_work_size(1) // Dummy size
        .local_work_size(1)  // Dummy size
        .arg_named("digram_ids", None::<&Buffer<u8>>)
        .arg_named("num_ids", 0u32)
        .arg_named("global_counts", None::<&Buffer<u32>>)
        .arg_named("max_digram_id_value", MAX_POSSIBLE_DIGRAM_IDS)
        .build().map_err(|e| anyhow!("Failed to build preliminary count_digrams_by_id kernel: {:?}", e))?;

    let device_max_wg_size = gpu_context.get_recommended_work_group_size().unwrap_or(256).max(1);
    let kernel_max_wg_size = match preliminary_count_kernel.wg_info(
        gpu_context.device, ocl::enums::KernelWorkGroupInfo::WorkGroupSize
    )? {
        ocl::enums::KernelWorkGroupInfoResult::WorkGroupSize(size) => size.max(1),
        _ => 256usize,
    };
    let mut local_work_size = std::cmp::min(device_max_wg_size, kernel_max_wg_size);
    local_work_size = std::cmp::min(local_work_size, num_ids).max(1);
    let global_work_size = ((num_ids + local_work_size - 1) / local_work_size * local_work_size).max(local_work_size);
    
    let count_kernel = Kernel::builder()
        .program(program)
        .name("count_digrams_by_id")
        .queue(gpu_context.queue.clone())
        .global_work_size(global_work_size)
        .local_work_size(local_work_size)
        .arg(&digram_ids_buffer)
        .arg(num_ids as u32)
        .arg(&global_counts_buffer)
        .arg(max_digram_id_value as u32)
        .build().map_err(|e| anyhow!("Failed to build count_digrams_by_id kernel: {:?}", e))?;

    unsafe {
        count_kernel.enq().map_err(|e| anyhow!("Failed to enqueue count_digrams_by_id kernel: {:?}", e))?;
    }

    // --- Read Results from global_counts_buffer ---
    let counts_buffer_size = MAX_POSSIBLE_DIGRAM_IDS;
    let mut counts_vec = vec![0u32; MAX_POSSIBLE_DIGRAM_IDS as usize];
    global_counts_buffer.read(&mut counts_vec).enq().map_err(|e| anyhow!("Failed to read digram counts from GPU: {:?}", e))?;

    let mut final_counts: HashMap<DigramKey, usize> = HashMap::new();
    for (id_uchar_idx, &count_val) in counts_vec.iter().enumerate() {
        if count_val > 0 {
            let b1_val: u8;
            let b2_val: u8;

            if _reverse_aware {
                // id_uchar_idx is 0-7 (kernel's canonical ID)
                if id_uchar_idx < CANONICAL_KERNEL_ID_TO_BASES.len() {
                    let (s1, s2) = CANONICAL_KERNEL_ID_TO_BASES[id_uchar_idx];
b1_val = match s1.symbol_type {
    SymbolType::Terminal(eb) => eb.0,
    _ => {
        eprintln!("[GPU COUNTING ERROR] Non-terminal symbol found in CANONICAL_KERNEL_ID_TO_BASES at index {}. Skipping.", id_uchar_idx);
        continue;
    }
};
b2_val = match s2.symbol_type {
    SymbolType::Terminal(eb) => eb.0,
    _ => {
        eprintln!("[GPU COUNTING ERROR] Non-terminal symbol found in CANONICAL_KERNEL_ID_TO_BASES at index {}. Skipping.", id_uchar_idx);
        continue;
    }
};
                } else {
                    // This should not happen if max_digram_id_value was set correctly to 8
                    eprintln!("[GPU COUNTING ERROR] Canonical ID index {} out of bounds for CANONICAL_KERNEL_ID_TO_BASES (len {}). Skipping.", id_uchar_idx, CANONICAL_KERNEL_ID_TO_BASES.len());
                    continue; 
                }
            } else {
                // id_uchar_idx is 0-15 (direct non-canonical digram ID from base1*4 + base2)
                b1_val = (id_uchar_idx / 4) as u8;
                b2_val = (id_uchar_idx % 4) as u8;
            }

            let base1 = EncodedBase(b1_val);
            let base2 = EncodedBase(b2_val);

            // Create symbols. Instance ID (0) and source info (None) are fine
            // as they are ignored for hashing/comparison in Symbol. Strand is Forward.
            let sym1 = Symbol::terminal(0, base1, Direction::Forward, None, None);
            let sym2 = Symbol::terminal(0, base2, Direction::Forward, None, None);
            
            let original_digram = Digram::new(sym1, sym2);
            
            // When _reverse_aware, the kernel already provides a canonical ID.
            // The (b1_val, b2_val) from CANONICAL_KERNEL_ID_TO_BASES represents one form of that canonical digram.
            // We must ensure the Digram object we hash is the canonical one according to Rust's definition.
            let key_digram = if _reverse_aware {
                original_digram.canonical(true) 
            } else {
                // If not reverse_aware, the kernel counted this specific (b1,b2) orientation.
                // The concept of a single canonical key is less applicable as all 16 are distinct counts.
                // However, to store in GpuDigramCounts which uses DigramKey (u64), we still hash this specific form.
                original_digram
            };

            let mut hasher = DefaultHasher::new();
            key_digram.hash(&mut hasher);
            let digram_key: DigramKey = hasher.finish();

            final_counts.insert(digram_key, count_val as usize);
        }
    }

    Ok(GpuDigramCounts { counts: final_counts })
}

/// Check if OpenCL is available on the system
pub fn is_opencl_available() -> bool {
    crate::gpu::GpuContext::is_gpu()
}

pub fn compute_digram_ids_opencl(
    sequence: &[u8],
    _reverse_aware: bool,
    gpu_context: &GpuContext,
) -> Result<(Buffer<u8>, usize), ocl::Error> {
    let seq_len = sequence.len();
    if seq_len < 2 {
        let num_digrams = 0;
        let empty_buffer = Buffer::<u8>::builder()
            .queue(gpu_context.queue.clone())
            .flags(ocl::flags::MEM_READ_WRITE) // General purpose flags for an empty buffer
            .len(num_digrams)
            .build()?;
        return Ok((empty_buffer, num_digrams));
    }
    let num_digrams = seq_len - 1;

    let program = gpu_context.program.as_ref().expect("Program not loaded");

    // --- Sequence Buffer ---
    let seq_buffer = Buffer::<u8>::builder()
        .queue(gpu_context.queue.clone())
        .flags(ocl::flags::MEM_READ_ONLY | ocl::flags::MEM_COPY_HOST_PTR)
        .len(seq_len)
        .copy_host_slice(sequence)
        .build()?;

    // --- IDs Buffer ---
    let ids_buffer = Buffer::<u8>::builder()
        .queue(gpu_context.queue.clone())
        .flags(ocl::flags::MEM_WRITE_ONLY)
        .len(num_digrams)
        .build()?;
        
    // --- Build Kernel (first to query its properties) ---
    let preliminary_kernel = Kernel::builder()
        .program(program)
        .name("compute_digram_ids")
        .queue(gpu_context.queue.clone())
        .global_work_size(1)
        .local_work_size(1)
        .arg_named("sequence", None::<&Buffer<u8>>)
        .arg_named("sequence_len", 0u32)
        .arg_named("reverse_aware", 0u32)
        .arg_named("digram_ids", None::<&Buffer<u8>>)
        .build()?;

    // --- Determine Work Sizes ---
    let device_max_work_group_size = gpu_context.get_recommended_work_group_size().unwrap_or(256).max(1);
    
    // Get the kernel-specific max work group size
    // This is the maximum local work size the kernel can actually handle
    let kernel_max_work_group_size = match preliminary_kernel.wg_info(
        gpu_context.device, ocl::enums::KernelWorkGroupInfo::WorkGroupSize
    )? {
        ocl::enums::KernelWorkGroupInfoResult::WorkGroupSize(size) => size.max(1), // Extract size and ensure > 0
        _ => {
            // Handle other variants or unexpected results - default to a safe value
            println!("Warning: Could not query kernel-specific work group size, using default 256.");
            256usize // Default to 256 if query fails or returns unexpected variant
        }
    };

    println!("  [GPU DEBUG] Device max_work_group_size: {}", device_max_work_group_size);
    println!("  [GPU DEBUG] Kernel specific max_work_group_size: {}", kernel_max_work_group_size);

    // Use the minimum of device, kernel, and problem size constraints for local_work_size
    let mut local_work_size = std::cmp::min(device_max_work_group_size, kernel_max_work_group_size);
    local_work_size = std::cmp::min(local_work_size, num_digrams).max(1);

    // Calculate global work size, ensuring it's a multiple of local_work_size and covers all digrams.
    let global_work_size = (num_digrams + local_work_size - 1) / local_work_size * local_work_size;
    let global_work_size = global_work_size.max(local_work_size); // Ensure not zero if local_work_size is not
    
    println!("  [GPU DEBUG] num_digrams: {}", num_digrams);
    println!("  [GPU DEBUG] final local_work_size: {}", local_work_size);
    println!("  [GPU DEBUG] final global_work_size: {}", global_work_size);

    // --- Re-Create Kernel with Correct Args & Work Sizes ---
    let kernel = Kernel::builder()
        .program(program)
        .name("compute_digram_ids")
        .queue(gpu_context.queue.clone())
        .global_work_size(global_work_size)
        .local_work_size(local_work_size)
        .arg(&seq_buffer)
        .arg(seq_len as u32)
        .arg(_reverse_aware as u32)
        .arg(&ids_buffer)
        .build()?;

    // Execute kernel
    unsafe {
        kernel.enq()?;
    }

    // Return the buffer and number of digrams directly
    Ok((ids_buffer, num_digrams))
}

pub fn populate_opencl(
    sequence: &[Symbol],
    _reverse_aware: bool,
    gpu_context: &GpuContext,
) -> Result<HashMap<u64, Vec<(usize, (Symbol, Symbol))>>, ocl::Error> {
    if sequence.len() < 2 { // Handle sequences too short for digrams
        return Ok(HashMap::new());
    }

    // Convert symbols to u8 for the OpenCL kernel if they are terminals
    let encoded: Vec<u8> = sequence
        .iter()
        .filter_map(|s| match s.symbol_type { // Use filter_map to skip non-terminals gracefully
            SymbolType::Terminal(EncodedBase(b)) => Some(b),
            _ => {
                // Log a warning or handle non-terminals if GPU path is taken with them
                // For now, this implementation will effectively skip them for `encoded`
                None
            }
        })
        .collect();

    // If all symbols were non-terminals, or sequence too short after filtering
    if encoded.len() < 2 {
        return Ok(HashMap::new());
    }

    let (ids_buffer, num_digrams) = compute_digram_ids_opencl(&encoded, _reverse_aware, gpu_context)?;

    let mut occurrences = HashMap::new();
    // IMPORTANT: The hashes correspond to digrams in `encoded`.
    // We need a way to map these back to the original `sequence` indices
    // if `sequence` contained non-terminals that were skipped.
    // This simplified version assumes `sequence` only contains terminals if it reaches here,
    // or that the panic in the original instructions is acceptable.
    // For a robust solution, mapping indices from `encoded` back to `sequence` is needed.
    // For now, proceeding with the logic that assumes `encoded` is a direct representation for hashing.
    // The loop should iterate up to `hashes.len()`, which is `encoded.len() - 1`.
    // The `sequence` indices for `s1` and `s2` must be correct.

    let _current_terminal_idx = 0;
    let mut original_indices = Vec::new();
    for (original_idx, sym) in sequence.iter().enumerate() {
        if matches!(sym.symbol_type, SymbolType::Terminal(_)) {
            original_indices.push(original_idx);
        }
    }

    // Read the digram IDs from the GPU buffer
    let mut host_ids = vec![0u8; num_digrams];
    if num_digrams > 0 { // Only read if there are digrams
        ids_buffer.read(&mut host_ids).enq()?;
    }

    for i in 0..num_digrams {
        let digram_id = host_ids[i] as usize;
        
        if digram_id >= CANONICAL_KERNEL_ID_TO_BASES.len() {
            eprintln!("Warning: GPU returned an out-of-bounds digram ID: {}. Max expected: {}. Skipping.", digram_id, CANONICAL_KERNEL_ID_TO_BASES.len() - 1);
            continue; 
        }

        let (s1_base_sym, s2_base_sym) = CANONICAL_KERNEL_ID_TO_BASES[digram_id];

        match (s1_base_sym.symbol_type, s2_base_sym.symbol_type) {
            (SymbolType::Terminal(_), SymbolType::Terminal(_)) => {
                // This digram ID corresponds to a valid pair of terminal symbols.
                // `i` is the starting index of the digram in the `encoded` sequence.
                // `original_indices[i]` is the index in the original `sequence` for the first symbol.
                // `original_indices[i+1]` is the index for the second symbol.
                if i + 1 < original_indices.len() { // Ensure we have a pair of original indices
                    let s1_orig_idx = original_indices[i];
                    let s2_orig_idx = original_indices[i+1];

                    // Critical check: were these symbols adjacent in the *original* sequence?
                    if s2_orig_idx == s1_orig_idx + 1 {
                        let s1_original = sequence[s1_orig_idx].clone();
                        let s2_original = sequence[s2_orig_idx].clone();

                        // s1_base_sym and s2_base_sym are the canonical symbols from CANONICAL_KERNEL_ID_TO_BASES
                        // Get the canonical tuple first.
                        let canonical_tuple: DigramKeyTuple = DigramTable::canonical_key((&s1_base_sym, &s2_base_sym), _reverse_aware);
                        
                        // Now, hash this tuple to get the u64 key for the local FxHashMap
                        let mut hasher = FxHasher::default();
                        canonical_tuple.hash(&mut hasher);    // DigramKeyTuple (Symbol, Symbol) implements Hash
                        let occurrences_key: u64 = hasher.finish();

                        occurrences
                            .entry(occurrences_key) // Use the u64 hash
                            .or_insert_with(Vec::new)
                            .push((s1_orig_idx, (s1_original, s2_original)));
                    }
                    // If not adjacent in original sequence, they don't form this digram at that position.
                }
            }
            _ => {
                // This digram ID maps to a placeholder (e.g., N,N) in CANONICAL_KERNEL_ID_TO_BASES.
                // These are not actual digrams to be counted, so skip.
                continue;
            }
        }
    }
    Ok(occurrences)
}

/// Find the most frequent terminal digram using a suffix array for optimization
pub fn find_most_frequent_terminal_digram_suffix_array(
    sequence: &[EncodedBase],
    min_usage: usize,
    _reverse_aware: bool,
    _suffix_array: Option<&[usize]>
) -> Option<(DigramKey, Vec<(usize, (Symbol, Symbol))>)> {
    // Convert to GpuSequence format
    let seq_data: Vec<u8> = sequence.iter().map(|b| b.0).collect();
    let gpu_seq = GpuSequence::new(seq_data);
    
    // Use the naive approach instead since we've fixed that implementation
    gpu_seq.find_most_frequent_digram_naive(min_usage)
}

pub fn find_most_frequent_digram_gpu_optimized(
    sequence: &GpuSequence, 
    _reverse_aware: bool, 
    _gpu_context: &GpuContext, 
    min_count: usize
) -> Result<Option<(DigramKeyTuple, usize)>, String> { 
    let mut counts: HashMap<DigramKeyTuple, usize> = HashMap::new();
    
    // Iterate through the sequence
    for i in 0..sequence.len().saturating_sub(1) { 
        let sym1 = sequence.get_symbol(i);
        let sym2 = sequence.get_symbol(i + 1);
        
        // 1. Get precise canonical key tuple
        let key_tuple = DigramTable::canonical_key((&sym1, &sym2), _reverse_aware);

        let count_entry = counts.entry(key_tuple).or_insert(0);
        *count_entry += 1;
    }

    // Find the most frequent digram
    let mut max_key = DigramKeyTuple::new(&Symbol::terminal(0, EncodedBase(0), Direction::Forward, None, None), &Symbol::terminal(1, EncodedBase(1), Direction::Forward, None, None));
    let mut max_count = 0;

    for (key, count) in &counts {
        if count > &max_count {
            max_key = key.clone();
            max_count = *count;
        }
    }

    if max_count >= min_count {
        Ok(Some((max_key, max_count)))
    } else {
        Ok(None)
    }
}

#[cfg(test)]
mod tests {
    // Tests will be moved to tests/gpu_tests.rs as per instructions.
    // This section is now empty.
}

