use anyhow::{Result, Context, anyhow};
use ocl::{Buffer, Kernel};
use std::collections::HashMap;
use crate::encode::dna_2bit::EncodedBase;
use crate::fasta::reader::FastaStream;
use crate::gpu::GpuContext;

/// Counts the occurrences of each base (A, C, G, T) in a FASTA stream using GPU acceleration.
///
/// Args:
///     fasta_stream: A stream providing chunks of DNA sequence data.
///     gpu_context: The GpuContext for interacting with the OpenCL device.
///
/// Returns:
///     A Result containing a HashMap mapping each EncodedBase to its count,
///     or an error if GPU processing fails.
pub fn count_bases_gpu_streaming(
    mut fasta_stream: FastaStream, // Take ownership and make mutable
    gpu_context: &GpuContext,
) -> Result<HashMap<EncodedBase, usize>> {
    println!("Starting GPU-accelerated base counting via streaming...");

    let queue = &gpu_context.queue;
    let program = gpu_context.program.as_ref()
        .ok_or_else(|| anyhow!("OpenCL program not loaded in GpuContext"))?;

    // Create a GPU buffer for the 4 base counts (A, C, G, T), initialized to zero.
    let counts_buffer = Buffer::<u32>::builder()
        .queue(queue.clone())
        .flags(ocl::flags::MEM_READ_WRITE)
        .len(4)
        .fill_val(0u32) // Initialize counts to 0
        .build()
        .with_context(|| "Failed to create GPU buffer for base counts")?;

    let mut total_bases_processed: usize = 0;

    // Process the FASTA stream chunk by chunk
    while let Some(chunk_data_u8) = fasta_stream.next() { // FastaStream yields Vec<u8>
        if chunk_data_u8.is_empty() {
            continue;
        }
        total_bases_processed += chunk_data_u8.len();

        // The FastaStream already gives Vec<u8> which might be ASCII or other encodings.
        // We need to ensure it's 2-bit encoded for the kernel.
        // For this, we first convert to EncodedBase, then to raw u8 for the kernel.
        let encoded_bases: Vec<EncodedBase> = chunk_data_u8.iter().filter_map(|&b| EncodedBase::from_base(b)).collect();
        let kernel_input_data: Vec<u8> = encoded_bases.iter().map(|eb| eb.0).collect();

        if kernel_input_data.is_empty() {
            println!("Skipping empty or non-DNA chunk after encoding.");
            continue;
        }
        
        let chunk_len = kernel_input_data.len();

        // Create a GPU buffer for the current sequence chunk.
        let sequence_chunk_buffer = Buffer::<u8>::builder()
            .queue(queue.clone())
            .flags(ocl::flags::MEM_READ_ONLY | ocl::flags::MEM_COPY_HOST_PTR)
            .len(chunk_len)
            .copy_host_slice(&kernel_input_data)
            .build()
            .with_context(|| format!("Failed to create GPU buffer for sequence chunk of size {}", chunk_len))?;

        // Build the kernel for counting bases in the current chunk.
        let count_kernel = Kernel::builder()
            .program(program)
            .name("count_bases_kernel")
            .queue(queue.clone())
            .global_work_size(chunk_len) // One work-item per base in the chunk
            .arg(&sequence_chunk_buffer)
            .arg(&counts_buffer)
            .arg(chunk_len as u32)
            .build()
            .with_context(|| "Failed to build count_bases_kernel")?;

        // Execute the kernel.
        unsafe {
            count_kernel.enq().with_context(|| "Failed to enqueue count_bases_kernel")?;
        }
        queue.finish().with_context(|| "Failed to finish queue after count_bases_kernel")?; // Ensure kernel execution completes
    }

    // Read the final counts back from the GPU buffer.
    let mut cpu_counts = vec![0u32; 4];
    counts_buffer.read(&mut cpu_counts).enq().with_context(|| "Failed to read base counts from GPU")?;
    queue.finish().with_context(|| "Failed to finish queue after reading counts")?;


    // Convert the raw counts (u32) into the HashMap<EncodedBase, usize>.
    let mut base_counts_map = HashMap::new();
    base_counts_map.insert(EncodedBase(0), cpu_counts[0] as usize);
    base_counts_map.insert(EncodedBase(1), cpu_counts[1] as usize);
    base_counts_map.insert(EncodedBase(2), cpu_counts[2] as usize);
    base_counts_map.insert(EncodedBase(3), cpu_counts[3] as usize);

    println!("Finished GPU-accelerated base counting. Processed {} bases.", total_bases_processed);
    println!("Raw counts from GPU: A:{}, C:{}, G:{}, T:{}", cpu_counts[0], cpu_counts[1], cpu_counts[2], cpu_counts[3]);

    Ok(base_counts_map)
}


#[cfg(test)]
mod tests {
    use super::*;
    use crate::fasta::reader::FastaStream;
    use crate::gpu::GpuContext;
    use std::path::PathBuf;
    use std::io::Write;
    use tempfile::NamedTempFile;

    fn create_test_fasta_file(content: &str) -> NamedTempFile {
        let mut file = NamedTempFile::new().unwrap();
        file.write_all(content.as_bytes()).unwrap();
        file
    }

    #[test]
    #[ignore] // Ignores test by default, requires a GPU and OpenCL environment.
              // Remove ignore to run test manually.
    fn test_count_bases_gpu_streaming_basic() -> Result<()> {
        let gpu_context = GpuContext::new()
            .with_context(|| "Failed to create GpuContext for test")?;

        let fasta_content = ">test_seq\nACGTACGTNNNACGT\n";
        let test_file = create_test_fasta_file(fasta_content);
        let fasta_path = PathBuf::from(test_file.path());

        let stream = FastaStream::new(&fasta_path)?
            .skip_ns() // Ensure Ns are skipped before encoding
            .with_chunk_size(10); // Small chunk size for testing

        let counts = count_bases_gpu_streaming(stream, &gpu_context)?;

        assert_eq!(counts.get(&EncodedBase(0)), Some(&3));
        assert_eq!(counts.get(&EncodedBase(1)), Some(&3));
        assert_eq!(counts.get(&EncodedBase(2)), Some(&3));
        assert_eq!(counts.get(&EncodedBase(3)), Some(&3));
        
        Ok(())
    }

    #[test]
    #[ignore]
    fn test_count_bases_gpu_streaming_large_chunks() -> Result<()> {
        let gpu_context = GpuContext::new()?;
        let fasta_content = ">test_seq\n".to_string() + &"ACGT".repeat(1000) + "\n"; // 4000 bases
        let test_file = create_test_fasta_file(&fasta_content);
        let fasta_path = PathBuf::from(test_file.path());

        let stream = FastaStream::new(&fasta_path)?.with_chunk_size(2048);
        let counts = count_bases_gpu_streaming(stream, &gpu_context)?;

        assert_eq!(counts.get(&EncodedBase(0)), Some(&1000));
        assert_eq!(counts.get(&EncodedBase(1)), Some(&1000));
        assert_eq!(counts.get(&EncodedBase(2)), Some(&1000));
        assert_eq!(counts.get(&EncodedBase(3)), Some(&1000));
        Ok(())
    }

    #[test]
    #[ignore]
    fn test_count_bases_gpu_streaming_empty_file() -> Result<()> {
        let gpu_context = GpuContext::new()?;
        let fasta_content = ">empty_seq\n";
        let test_file = create_test_fasta_file(fasta_content);
        let fasta_path = PathBuf::from(test_file.path());

        let stream = FastaStream::new(&fasta_path)?;
        let counts = count_bases_gpu_streaming(stream, &gpu_context)?;

        assert_eq!(counts.get(&EncodedBase(0)), Some(&0));
        assert_eq!(counts.get(&EncodedBase(1)), Some(&0));
        assert_eq!(counts.get(&EncodedBase(2)), Some(&0));
        assert_eq!(counts.get(&EncodedBase(3)), Some(&0));
        Ok(())
    }
    
    #[test]
    #[ignore]
    fn test_count_bases_gpu_streaming_all_one_base() -> Result<()> {
        let gpu_context = GpuContext::new()?;
        let fasta_content = ">all_a\nAAAA\nAAAA\nAA"; // 10 'A's
        let test_file = create_test_fasta_file(fasta_content);
        let fasta_path = PathBuf::from(test_file.path());

        let stream = FastaStream::new(&fasta_path)?.with_chunk_size(5);
        let counts = count_bases_gpu_streaming(stream, &gpu_context)?;
        
        assert_eq!(counts.get(&EncodedBase(0)), Some(&10));
        assert_eq!(counts.get(&EncodedBase(1)), Some(&0));
        assert_eq!(counts.get(&EncodedBase(2)), Some(&0));
        assert_eq!(counts.get(&EncodedBase(3)), Some(&0));
        Ok(())
    }
} 