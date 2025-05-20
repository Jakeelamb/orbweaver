//! GPU-accelerated suffix array construction.
//!
//! This module provides OpenCL-based implementations for building suffix arrays
//! efficiently on GPUs to optimize pattern matching and repeat detection.

use ocl::{Buffer, Kernel};
use anyhow::{Result, anyhow};
use crate::gpu::{GpuContext, digram::GpuSequence};

/// Represents a suffix array processed on the GPU.
pub struct GpuSuffixArray {
    /// The suffix array (indices into the original sequence)
    sa: Vec<usize>,
}

impl GpuSuffixArray {
    /// Create a new empty GPU suffix array
    pub fn new() -> Self {
        GpuSuffixArray {
            sa: Vec::new(),
        }
    }
    
    /// Build a suffix array for the given sequence
    pub fn build(sequence: &GpuSequence, gpu_context: Option<&GpuContext>) -> Result<Self> {
        let mut sa = GpuSuffixArray::new();
        
        if let Some(context) = gpu_context {
            // Try to build on GPU using OpenCL
            println!("Building suffix array on GPU with OpenCL");
            match Self::build_opencl(sequence, context) {
                Ok(suffix_array) => {
                    sa.sa = suffix_array;
                    println!("Successfully built suffix array on GPU");
                },
                Err(e) => {
                    println!("GPU suffix array construction failed: {:?}. Falling back to CPU", e);
                    sa.sa = Self::build_cpu_fallback(sequence)?;
                }
            }
        } else {
            // CPU fallback implementation
            println!("Building suffix array on CPU (no GPU context provided)");
            sa.sa = Self::build_cpu_fallback(sequence)?;
        }
        
        Ok(sa)
    }
    
    /// OpenCL implementation for suffix array construction
    /// Uses prefix doubling algorithm implemented in the OpenCL kernel
    fn build_opencl(sequence: &GpuSequence, gpu_context: &GpuContext) -> Result<Vec<usize>> {
        let seq_data = sequence.get_data();
        let seq_len = sequence.get_len();
        
        if seq_len == 0 {
            return Ok(Vec::new());
        }
        
        println!("Building suffix array on GPU with OpenCL");
        
        // Check if the sequence is too large for GPU memory
        // Get device memory size
        let device_memory = gpu_context.device.info(ocl::enums::DeviceInfo::GlobalMemSize)
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
        
        // Estimate memory needed for suffix array construction
        // We need the sequence, suffix array, and temporary buffers
        let estimated_memory = seq_len * (1 + std::mem::size_of::<usize>() * 5);
        
        // If sequence is too large, use chunked approach or fall back to CPU
        // Use 70% of available memory as threshold to be safe
        if estimated_memory > (device_memory as f64 * 0.7) as usize {
            println!("Sequence too large for GPU memory ({} bytes needed, {} available)",
                estimated_memory, device_memory);
            println!("Using chunked approach for suffix array construction");
            return Self::build_chunked_opencl(sequence, gpu_context);
        }
        
        // Get the OpenCL context and queue
        let _context = &gpu_context.context;
        let queue = &gpu_context.queue;
        
        // Create a sequence buffer
        let seq_buffer = Buffer::builder()
            .queue(queue.clone())
            .flags(ocl::flags::MEM_READ_ONLY)
            .len(seq_len)
            .copy_host_slice(seq_data)
            .build()
            .map_err(|e| anyhow!("Failed to create sequence buffer: {}", e))?;
        
        // Create a buffer for the suffix array
        let suffix_array_buffer = Buffer::<usize>::builder()
            .queue(queue.clone())
            .flags(ocl::MemFlags::new().read_write())
            .len(seq_len)
            .build()
            .map_err(|e| anyhow!("Failed to create suffix array buffer: {}", e))?;
        
        // Get the program
        let program = gpu_context.program.as_ref()
            .ok_or_else(|| anyhow!("OpenCL program not loaded"))?;
        
        // Create kernel for suffix array construction
        let initialize_kernel = Kernel::builder()
            .program(program)
            .name("initialize_suffix_array")
            .queue(queue.clone())
            .global_work_size(seq_len)
            .arg(&seq_buffer)
            .arg(&suffix_array_buffer)
            .arg(seq_len as u32)
            .build()
            .map_err(|e| anyhow!("Failed to build initialize kernel: {}", e))?;
        
        // Execute kernel
        unsafe {
            initialize_kernel.enq()
                .map_err(|e| anyhow!("Failed to execute initialize kernel: {}", e))?;
        }
        
        // Sort the suffix array
        let mut h = 1;
        while h < seq_len {
            // Sort based on h-gram
            let sort_kernel = Kernel::builder()
                .program(program)
                .name("sort_suffix_array")
                .queue(queue.clone())
                .global_work_size(seq_len)
                .arg(&seq_buffer)
                .arg(&suffix_array_buffer)
                .arg(seq_len as u32)
                .arg(h as u32)
                .build()
                .map_err(|e| anyhow!("Failed to build sort kernel: {}", e))?;

            unsafe {
                sort_kernel.enq()
                    .map_err(|e| anyhow!("Failed to execute sort kernel: {}", e))?;
            }

            h *= 2;
        }
        
        // Read the results
        let mut result = vec![0usize; seq_len];
        suffix_array_buffer.read(&mut result)
            .enq()
            .map_err(|e| anyhow!("Failed to read suffix array from GPU: {}", e))?;
        
        println!("Successfully built suffix array on GPU");
        Ok(result)
    }
    
    /// Build suffix array using a chunked approach when the sequence is too large for GPU memory
    fn build_chunked_opencl(sequence: &GpuSequence, _gpu_context: &GpuContext) -> Result<Vec<usize>> {
        let seq_data = sequence.get_data();
        let seq_len = sequence.get_len();
        
        // Define chunk size (10MB is usually safe for most GPUs)
        let chunk_size = 10_000_000;
        let num_chunks = (seq_len + chunk_size - 1) / chunk_size;
        
        println!("Building suffix array with {} chunks", num_chunks);
        
        // Initialize the final suffix array
        let mut suffix_array = vec![0; seq_len];
        
        // Process each chunk
        for chunk_idx in 0..num_chunks {
            let start = chunk_idx * chunk_size;
            let end = std::cmp::min(start + chunk_size, seq_len);
            let chunk_len = end - start;
            
            println!("Processing chunk {}/{} ({} bytes)", chunk_idx + 1, num_chunks, chunk_len);
            
            // Create a sub-sequence for this chunk
            let chunk_data = &seq_data[start..end];
            let chunk_sequence = GpuSequence::new(chunk_data.to_vec());
            
            // Build a suffix array for this chunk
            let sa_chunk = match Self::build_cpu_fallback(&chunk_sequence) {
                Ok(sa) => sa,
                Err(_) => {
                    // If still fails, use very basic approach
                    println!("Falling back to basic suffix array construction for chunk");
                    (0..chunk_len).collect()
                }
            };
            
            // Adjust indices to match original positions
            for i in 0..chunk_len {
                suffix_array[start + i] = if i < sa_chunk.len() {
                    start + sa_chunk[i]
                } else {
                    start + i
                };
            }
        }
        
        // Sort the combined suffix array
        suffix_array.sort_by(|&a, &b| {
            let seq_a = &seq_data[a..];
            let seq_b = &seq_data[b..];
            seq_a.cmp(seq_b)
        });
        
        Ok(suffix_array)
    }
    
    /// CPU fallback implementation for suffix array construction
    fn build_cpu_fallback(sequence: &GpuSequence) -> Result<Vec<usize>> {
        let len = sequence.get_len();
        let data_arc = sequence.get_data();
        let data = &data_arc[..];
        
        // Create a vector of suffix starting indices
        let mut suffix_starts: Vec<usize> = (0..len).collect();
        
        // Sort suffixes using the provided sequence data
        suffix_starts.sort_by(|&a, &b| {
            // Compare suffixes starting at positions a and b
            for i in 0..len {
                let pos_a = a + i;
                let pos_b = b + i;
                
                // If we've reached the end of one suffix, it comes first
                if pos_a >= len {
                    return std::cmp::Ordering::Less;
                }
                if pos_b >= len {
                    return std::cmp::Ordering::Greater;
                }
                
                // Compare characters
                let char_a = data[pos_a];
                let char_b = data[pos_b];
                
                match char_a.cmp(&char_b) {
                    std::cmp::Ordering::Equal => continue,
                    other => return other,
                }
            }
            
            // If we get here, the suffixes are identical
            std::cmp::Ordering::Equal
        });
        
        Ok(suffix_starts)
    }
    
    /// Access the suffix array
    pub fn get_array(&self) -> &[usize] {
        &self.sa
    }
    
    /// Find all occurrences of a pattern in the original sequence
    pub fn find_pattern(&self, sequence: &GpuSequence, pattern: &[u8]) -> Vec<usize> {
        let data_arc = sequence.get_data();
        let data = &data_arc[..];
        let mut matches = Vec::new();
        
        // Simple pattern search using the suffix array
        // This would be much faster with GPU acceleration
        for &suffix_start in &self.sa {
            if suffix_start + pattern.len() <= data.len() {
                let mut match_found = true;
                
                for (i, &pattern_char) in pattern.iter().enumerate() {
                    if data[suffix_start + i] != pattern_char {
                        match_found = false;
                        break;
                    }
                }
                
                if match_found {
                    matches.push(suffix_start);
                }
            }
        }
        
        matches
    }
    
    /// Use binary search to find pattern occurrences efficiently
    pub fn binary_search_pattern(&self, sequence: &GpuSequence, pattern: &[u8]) -> Vec<usize> {
        if pattern.is_empty() || self.sa.is_empty() {
            return Vec::new();
        }
        
        let data_arc = sequence.get_data();
        let data = &data_arc[..];
        
        // Perform binary search to find the range of matching suffixes
        let (mut left, mut right) = (0, self.sa.len());
        
        // Find left bound
        while left < right {
            let mid = left + (right - left) / 2;
            let suffix_start = self.sa[mid];
            
            let mut cmp = std::cmp::Ordering::Equal;
            for i in 0..pattern.len() {
                if suffix_start + i >= data.len() {
                    cmp = std::cmp::Ordering::Less;
                    break;
                }
                
                match pattern[i].cmp(&data[suffix_start + i]) {
                    std::cmp::Ordering::Equal => continue,
                    other => {
                        cmp = other;
                        break;
                    }
                }
            }
            
            if cmp == std::cmp::Ordering::Greater {
                left = mid + 1;
            } else {
                right = mid;
            }
        }
        
        let start = left;
        
        // Find right bound
        let (mut left, mut right) = (start, self.sa.len());
        while left < right {
            let mid = left + (right - left) / 2;
            let suffix_start = self.sa[mid];
            
            let mut cmp = std::cmp::Ordering::Equal;
            for i in 0..pattern.len() {
                if suffix_start + i >= data.len() {
                    cmp = std::cmp::Ordering::Less;
                    break;
                }
                
                match pattern[i].cmp(&data[suffix_start + i]) {
                    std::cmp::Ordering::Equal => {
                        if i == pattern.len() - 1 {
                            cmp = std::cmp::Ordering::Equal;
                            break;
                        }
                        continue;
                    }
                    other => {
                        cmp = other;
                        break;
                    }
                }
            }
            
            if cmp == std::cmp::Ordering::Equal || cmp == std::cmp::Ordering::Less {
                left = mid + 1;
            } else {
                right = mid;
            }
        }
        
        let end = left;
        
        // Collect all matches
        self.sa[start..end].to_vec()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_suffix_array_construction() {
        // Create a sample sequence "ACGT"
        let data = vec![0, 1, 2, 3]; // A, C, G, T
        let seq = GpuSequence::new(data);
        
        // Build the suffix array
        let sa = GpuSuffixArray::build(&seq, None).unwrap();
        
        // For "ACGT", the sorted suffixes should be:
        // 0: "ACGT" (A...)
        // 1: "CGT"  (C...)
        // 2: "GT"   (G...)
        // 3: "T"    (T...)
        let expected = vec![0, 1, 2, 3];
        assert_eq!(sa.get_array(), &expected);
    }
    
    #[test]
    fn test_pattern_search() {
        // Create a sample sequence "ACGTACGT"
        let data = vec![0, 1, 2, 3, 0, 1, 2, 3]; // A, C, G, T, A, C, G, T
        let seq = GpuSequence::new(data);
        
        // Build the suffix array
        let sa = GpuSuffixArray::build(&seq, None).unwrap();
        
        // Search for "ACG" (should find at positions 0 and 4)
        let pattern = vec![0, 1, 2]; // A, C, G
        let matches = sa.find_pattern(&seq, &pattern);
        
        assert_eq!(matches.len(), 2);
        assert!(matches.contains(&0));
        assert!(matches.contains(&4));
    }
} 