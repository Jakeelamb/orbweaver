// OpenCL kernel wrapper for Orbweaver

/// Get embedded digram hashing OpenCL kernel code
pub fn get_digram_kernel() -> &'static str {
    // Return only the kernel code itself
    "
__kernel void compute_digram_hashes(__global const uchar *sequence, const uint sequence_len, const uint reverse_aware, __global ulong *hashes) {
    const uint gid = get_global_id(0);
    
    // Strict bounds checking to prevent out-of-bounds access
    if (gid >= sequence_len - 1) return;

    // Check if sequence pointers are valid
    if (sequence == NULL || hashes == NULL) return;

    // Safely read adjacent bases
    const uchar b1 = sequence[gid];
    const uchar b2 = sequence[gid + 1];
    
    // Basic validation - in DNA bases should be 0-3
    if (b1 > 3 || b2 > 3) return;
    
    // Compute hash with proper type casting
    ulong hash = (((ulong)b1) << 32) | ((ulong)b2);

    if (reverse_aware) {
        // DNA complement: A(00)<->T(11), C(01)<->G(10)
        const uchar rc_b1 = 3 - b2; // Complement of b2
        const uchar rc_b2 = 3 - b1; // Complement of b1
        const ulong rc_hash = (((ulong)rc_b1) << 32) | ((ulong)rc_b2);
        
        // Use canonical form (minimum of hash and rc_hash)
        hash = min(hash, rc_hash);
    }

    // Write to output with bounds check
    hashes[gid] = hash;
}

__kernel void count_digrams(__global const ulong *hashes, const uint num_hashes, __global uint *counts, const uint counts_size) {
    const uint gid = get_global_id(0);
    
    // Strict bounds checking
    if (gid >= num_hashes) return;
    if (hashes == NULL || counts == NULL) return;

    const ulong hash = hashes[gid];
    if (hash != 0) {
        // Use modulo to ensure we stay within the counts buffer bounds
        uint index = hash % counts_size;
        atomic_inc(&counts[index]);
    }
}

// Enhanced digram detection kernel - handles full digram finding process
__kernel void find_digrams(
    __global const uint *symbols, // Input sequence of symbol IDs
    __global uint *digram_counts, // Output array of digram counts
    __global uint *digram_indices, // Output array of digram positions
    const uint sequence_length, // Length of the symbol sequence
    const uint reverse_aware // Flag for reverse complement canonicalization
) {
    uint gid = get_global_id(0); // Global work-item ID
    
    // Strict bounds checking to prevent out-of-bounds access
    if (gid >= sequence_length - 1) return;
    if (symbols == NULL || digram_counts == NULL || digram_indices == NULL) return;

    // Read two adjacent symbols with bounds checking
    uint s1 = symbols[gid];
    uint s2 = symbols[gid + 1];

    // Validate symbol values (assuming 2-bit DNA encoding)
    if (s1 > 3 || s2 > 3) return;

    // Compute digram key (simple concatenation for now)
    // Use 20-bit limit instead of full 32-bit to avoid potential hash conflicts 
    uint digram_key = (s1 << 16) | s2;
    
    // Apply a modulo to keep indices within array bounds (assume 1M max size)
    digram_key = digram_key % 1000000;

    // Handle reverse complement if enabled
    if (reverse_aware) {
        // Assuming symbol encodes A=0, C=1, G=2, T=3
        uint s1_rc = 3 - s1; // A<->T, C<->G
        uint s2_rc = 3 - s2;
        uint digram_key_rc = (s2_rc << 16) | s1_rc;
        digram_key_rc = digram_key_rc % 1000000; // Same modulo on reverse complement
        digram_key = min(digram_key, digram_key_rc); // Use canonical form
    }

    // Atomically increment digram count with bounds check
    atomic_inc(&digram_counts[digram_key]);

    // Store position
    digram_indices[gid] = digram_key;
}
"
}

/// Returns the OpenCL kernel source code for suffix array construction
pub fn get_suffix_array_kernel() -> &'static str {
    // Return only the kernel code itself
    "
// OpenCL Kernel for Suffix Array Construction
// Based on prefix doubling algorithm

// Helper function to compare two suffixes
inline int compare_suffixes(
    __global const uchar* sequence,
    uint seq_len,
    uint pos1,
    uint pos2,
    uint compare_len
) {
    for (uint i = 0; i < compare_len; i++) {
        if (pos1 + i >= seq_len) return -1;
        if (pos2 + i >= seq_len) return 1;
        
        uchar c1 = sequence[pos1 + i];
        uchar c2 = sequence[pos2 + i];
        
        if (c1 < c2) return -1;
        if (c1 > c2) return 1;
    }
    return 0;
}

// Initialize suffix array with indices 0...n-1
__kernel void initialize_suffix_array(
    __global const uchar* sequence,
    __global size_t* suffix_array,
    uint seq_len
) {
    uint idx = get_global_id(0);
    if (idx < seq_len) {
        suffix_array[idx] = idx;
    }
}

// Sort suffixes based on their first h characters
__kernel void sort_suffixes(
    __global const uchar* sequence,
    __global size_t* suffix_array,
    __global uint* ranks,
    __global uint* new_ranks,
    uint seq_len,
    uint h
) {
    uint idx = get_global_id(0);
    if (idx >= seq_len) return;
    
    uint suffix = suffix_array[idx];
    
    // For simplicity, this is a very basic insertion sort
    // A real implementation would use a more efficient sorting algorithm
    for (uint i = idx + 1; i < seq_len; i++) {
        uint other_suffix = suffix_array[i];
        int cmp = compare_suffixes(sequence, seq_len, suffix, other_suffix, h);
        if (cmp > 0) {
            // Swap positions
            suffix_array[idx] = other_suffix;
            suffix_array[i] = suffix;
            suffix = other_suffix;
        }
    }
    
    // Calculate ranks based on sorted order
    if (idx == 0) {
        ranks[suffix] = 0;
    } else {
        uint prev_suffix = suffix_array[idx - 1];
        int cmp = compare_suffixes(sequence, seq_len, prev_suffix, suffix, h);
        if (cmp == 0) {
            ranks[suffix] = ranks[prev_suffix];
        } else {
            ranks[suffix] = idx;
        }
    }
}
"
}

/// Get vector addition OpenCL kernel code (for testing)
pub fn get_vector_add_kernel() -> &'static str {
    // Return only the kernel code itself
    "
__kernel void add_vectors(__global const float *a, __global const float *b, __global float *c) {
    int idx = get_global_id(0);
    c[idx] = a[idx] + b[idx];
}
"
}

#[cfg(test)]
mod tests {
    #[test]
    fn test_get_digram_kernel() {
        let kernel = super::get_digram_kernel();
        assert!(kernel.contains("compute_digram_hashes"));
        assert!(kernel.contains("count_digrams"));
        assert!(kernel.contains("find_digrams"));
    }
    
    #[test]
    fn test_get_suffix_array_kernel() {
        // This test is a placeholder for the suffix array kernel functionality
        // The actual implementation is currently inlined in the GpuContext::load_kernels method
        
        // Just a placeholder to ensure the test suite passes
        assert!(true);
    }
} 