//! GPU-accelerated LCG (Locally Consistent Grammar) operations.
//!
//! This module provides GPU-accelerated implementations of the most
//! computationally intensive parts of LCG parsing:
//! - Position fingerprint computation
//! - Local minimum detection for cut points
//! - Phrase fingerprint computation

use super::GpuContext;
use crate::encode::dna_2bit::EncodedBase;
use crate::grammar::lcg::Fingerprint;
use anyhow::{Context, Result};
use ocl::{Buffer, Kernel, MemFlags};
use std::sync::Arc;

/// Minimum sequence length to use GPU acceleration
/// Below this threshold, CPU is more efficient due to GPU overhead
pub const GPU_THRESHOLD: usize = 10_000;

/// GPU-accelerated LCG parser
pub struct GpuLcgParser {
    gpu_context: Arc<GpuContext>,
    work_group_size: usize,
}

impl GpuLcgParser {
    /// Create a new GPU LCG parser
    pub fn new(gpu_context: Arc<GpuContext>) -> Result<Self> {
        let work_group_size = gpu_context.get_recommended_work_group_size()?;
        Ok(Self {
            gpu_context,
            work_group_size,
        })
    }

    /// Check if GPU acceleration is available and beneficial for the given sequence length
    pub fn should_use_gpu(&self, sequence_len: usize) -> bool {
        sequence_len >= GPU_THRESHOLD && self.gpu_context.program.is_some()
    }

    /// Compute position fingerprints on the GPU.
    ///
    /// Each position gets a fingerprint based on the digram at that position.
    /// This is the parallel version of LCG's fingerprint computation.
    pub fn compute_fingerprints(&self, sequence: &[EncodedBase]) -> Result<Vec<Fingerprint>> {
        let program = self
            .gpu_context
            .program
            .as_ref()
            .context("GPU kernels not loaded")?;

        let sequence_len = sequence.len();
        if sequence_len == 0 {
            return Ok(Vec::new());
        }

        // Convert to raw bytes for GPU
        let sequence_bytes: Vec<u8> = sequence.iter().map(|b| b.0).collect();

        // Create GPU buffers
        let sequence_buffer = Buffer::<u8>::builder()
            .queue(self.gpu_context.queue.clone())
            .flags(MemFlags::new().read_only())
            .len(sequence_len)
            .copy_host_slice(&sequence_bytes)
            .build()
            .context("Failed to create sequence buffer")?;

        let fingerprints_buffer = Buffer::<u64>::builder()
            .queue(self.gpu_context.queue.clone())
            .flags(MemFlags::new().write_only())
            .len(sequence_len)
            .build()
            .context("Failed to create fingerprints buffer")?;

        // Create and execute kernel
        let kernel = Kernel::builder()
            .program(program)
            .name("compute_position_fingerprints")
            .queue(self.gpu_context.queue.clone())
            .global_work_size(round_up_to_work_group(sequence_len, self.work_group_size))
            .arg(&sequence_buffer)
            .arg(&(sequence_len as u32))
            .arg(&fingerprints_buffer)
            .build()
            .context("Failed to build compute_position_fingerprints kernel")?;

        unsafe {
            kernel
                .enq()
                .context("Failed to enqueue compute_position_fingerprints kernel")?;
        }

        // Read results
        let mut fingerprints = vec![0u64; sequence_len];
        fingerprints_buffer
            .read(&mut fingerprints)
            .enq()
            .context("Failed to read fingerprints from GPU")?;

        Ok(fingerprints)
    }

    /// Find cut points using GPU-accelerated local minimum detection.
    ///
    /// Returns a vector of positions where cuts should be made.
    pub fn find_cut_points(
        &self,
        fingerprints: &[Fingerprint],
        half_window: usize,
        min_phrase_len: usize,
        max_phrase_len: usize,
    ) -> Result<Vec<usize>> {
        let program = self
            .gpu_context
            .program
            .as_ref()
            .context("GPU kernels not loaded")?;

        let num_fingerprints = fingerprints.len();
        if num_fingerprints == 0 {
            return Ok(Vec::new());
        }

        // Create GPU buffers
        let fingerprints_buffer = Buffer::<u64>::builder()
            .queue(self.gpu_context.queue.clone())
            .flags(MemFlags::new().read_only())
            .len(num_fingerprints)
            .copy_host_slice(fingerprints)
            .build()
            .context("Failed to create fingerprints buffer")?;

        let is_cut_point_buffer = Buffer::<u8>::builder()
            .queue(self.gpu_context.queue.clone())
            .flags(MemFlags::new().read_write())
            .len(num_fingerprints)
            .build()
            .context("Failed to create is_cut_point buffer")?;

        let cut_point_count_buffer = Buffer::<u32>::builder()
            .queue(self.gpu_context.queue.clone())
            .flags(MemFlags::new().read_write())
            .len(1)
            .copy_host_slice(&[0u32])
            .build()
            .context("Failed to create cut_point_count buffer")?;

        // Create and execute find_local_minima kernel
        let kernel = Kernel::builder()
            .program(program)
            .name("find_local_minima")
            .queue(self.gpu_context.queue.clone())
            .global_work_size(round_up_to_work_group(num_fingerprints, self.work_group_size))
            .arg(&fingerprints_buffer)
            .arg(&(num_fingerprints as u32))
            .arg(&(half_window as u32))
            .arg(&(min_phrase_len as u32))
            .arg(&(max_phrase_len as u32))
            .arg(&is_cut_point_buffer)
            .arg(&cut_point_count_buffer)
            .build()
            .context("Failed to build find_local_minima kernel")?;

        unsafe {
            kernel
                .enq()
                .context("Failed to enqueue find_local_minima kernel")?;
        }

        // Read is_cut_point array
        let mut is_cut_point = vec![0u8; num_fingerprints];
        is_cut_point_buffer
            .read(&mut is_cut_point)
            .enq()
            .context("Failed to read is_cut_point from GPU")?;

        // Collect cut points (filter on CPU since it's a simple operation)
        // We also apply min_phrase_len and max_phrase_len constraints here
        let mut cut_points = Vec::new();
        let mut last_cut = 0;

        for (pos, &is_cut) in is_cut_point.iter().enumerate() {
            if is_cut == 1 {
                let phrase_len = pos - last_cut;
                // Include if meets min_phrase_len or is first position
                if pos == 0 || phrase_len >= min_phrase_len {
                    cut_points.push(pos);
                    last_cut = pos;
                }
            }
            // Force cut at max_phrase_len
            else if pos > 0 && pos - last_cut >= max_phrase_len {
                cut_points.push(pos);
                last_cut = pos;
            }
        }

        Ok(cut_points)
    }

    /// Compute phrase fingerprints given cut points.
    ///
    /// Returns a fingerprint for each phrase defined by the cut points.
    pub fn compute_phrase_fingerprints(
        &self,
        sequence: &[EncodedBase],
        cut_points: &[usize],
    ) -> Result<Vec<Fingerprint>> {
        let program = self
            .gpu_context
            .program
            .as_ref()
            .context("GPU kernels not loaded")?;

        let sequence_len = sequence.len();
        let num_phrases = cut_points.len();

        if num_phrases == 0 || sequence_len == 0 {
            return Ok(Vec::new());
        }

        // Convert sequence to bytes
        let sequence_bytes: Vec<u8> = sequence.iter().map(|b| b.0).collect();

        // Convert cut points to u32
        let cut_points_u32: Vec<u32> = cut_points.iter().map(|&p| p as u32).collect();

        // Create GPU buffers
        let sequence_buffer = Buffer::<u8>::builder()
            .queue(self.gpu_context.queue.clone())
            .flags(MemFlags::new().read_only())
            .len(sequence_len)
            .copy_host_slice(&sequence_bytes)
            .build()
            .context("Failed to create sequence buffer")?;

        let cut_points_buffer = Buffer::<u32>::builder()
            .queue(self.gpu_context.queue.clone())
            .flags(MemFlags::new().read_only())
            .len(num_phrases)
            .copy_host_slice(&cut_points_u32)
            .build()
            .context("Failed to create cut_points buffer")?;

        let phrase_fingerprints_buffer = Buffer::<u64>::builder()
            .queue(self.gpu_context.queue.clone())
            .flags(MemFlags::new().write_only())
            .len(num_phrases)
            .build()
            .context("Failed to create phrase_fingerprints buffer")?;

        // Create and execute kernel
        let kernel = Kernel::builder()
            .program(program)
            .name("compute_phrase_fingerprints")
            .queue(self.gpu_context.queue.clone())
            .global_work_size(round_up_to_work_group(num_phrases, self.work_group_size))
            .arg(&sequence_buffer)
            .arg(&(sequence_len as u32))
            .arg(&cut_points_buffer)
            .arg(&(num_phrases as u32))
            .arg(&phrase_fingerprints_buffer)
            .build()
            .context("Failed to build compute_phrase_fingerprints kernel")?;

        unsafe {
            kernel
                .enq()
                .context("Failed to enqueue compute_phrase_fingerprints kernel")?;
        }

        // Read results
        let mut phrase_fingerprints = vec![0u64; num_phrases];
        phrase_fingerprints_buffer
            .read(&mut phrase_fingerprints)
            .enq()
            .context("Failed to read phrase_fingerprints from GPU")?;

        Ok(phrase_fingerprints)
    }
}

/// Round up a value to the nearest multiple of work_group_size
fn round_up_to_work_group(value: usize, work_group_size: usize) -> usize {
    ((value + work_group_size - 1) / work_group_size) * work_group_size
}

#[cfg(test)]
mod tests {
    use super::*;

    fn encode_seq(s: &[u8]) -> Vec<EncodedBase> {
        s.iter()
            .filter_map(|&b| EncodedBase::from_base(b))
            .collect()
    }

    #[test]
    fn test_gpu_threshold() {
        // Just verify the threshold constant is reasonable
        assert!(GPU_THRESHOLD >= 1000);
        assert!(GPU_THRESHOLD <= 100_000);
    }

    // Note: GPU tests require actual GPU hardware, so they're marked as ignored by default
    #[test]
    #[ignore]
    fn test_gpu_fingerprints() {
        let gpu_context = match GpuContext::new() {
            Ok(ctx) => Arc::new(ctx),
            Err(_) => {
                println!("GPU not available, skipping test");
                return;
            }
        };

        let parser = GpuLcgParser::new(gpu_context).expect("Failed to create GPU parser");

        let sequence = encode_seq(b"ACGTACGTACGTACGT");
        let fingerprints = parser
            .compute_fingerprints(&sequence)
            .expect("Failed to compute fingerprints");

        assert_eq!(fingerprints.len(), sequence.len());
    }

    #[test]
    #[ignore]
    fn test_gpu_cut_points() {
        let gpu_context = match GpuContext::new() {
            Ok(ctx) => Arc::new(ctx),
            Err(_) => {
                println!("GPU not available, skipping test");
                return;
            }
        };

        let parser = GpuLcgParser::new(gpu_context).expect("Failed to create GPU parser");

        // Create sample fingerprints
        let fingerprints: Vec<u64> = (0..100).map(|i| (i * 37 + 17) as u64).collect();

        let cut_points = parser
            .find_cut_points(&fingerprints, 4, 2, 256)
            .expect("Failed to find cut points");

        // Should always have at least one cut point (position 0)
        assert!(!cut_points.is_empty());
        assert_eq!(cut_points[0], 0);
    }
}
