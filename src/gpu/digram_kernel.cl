// digram_kernel.cl

// Each base (A, C, G, T) is assumed to be encoded as 00, 01, 10, 11 respectively.
// A digram ID is then (base1 << 2) | base2, resulting in values 0-15.

// Helper to get canonical digram ID
// Input: b1, b2 are 2-bit encoded (0, 1, 2, 3)
// Output: canonical digram ID (0-9 for DNA)
uchar get_canonical_digram_id(uchar b1, uchar b2) {
    // Complement: A(0)<->T(3), C(1)<->G(2). So, complement of base 'b' is '3-b'.
    uchar rc_b1_val = 3 - b2; // Complement of b2 becomes the first base of the reverse complement digram
    uchar rc_b2_val = 3 - b1; // Complement of b1 becomes the second base of the reverse complement digram

    uchar forward_id = (b1 << 2) | b2;
    uchar reverse_complement_id = (rc_b1_val << 2) | rc_b2_val;

    return min(forward_id, reverse_complement_id);
}

__kernel void compute_digram_ids(
    __global const uchar *sequence,      // Input: 2-bit encoded sequence (0,1,2,3 per uchar)
    const uint sequence_len,
    const uint reverse_aware,            // 1 for reverse_aware (canonical), 0 otherwise
    __global uchar *digram_ids           // Output: compact digram IDs (0-15, or 0-9 if canonical)
) {
    const uint gid = get_global_id(0);
    if (gid >= sequence_len - 1) return; // Ensure we can read sequence[gid+1]

    const uchar b1 = sequence[gid];
    const uchar b2 = sequence[gid + 1];

    // Basic validation: ensure bases are 0, 1, 2, or 3.
    // If not, a marker could be written, or behavior defined by host.
    // For performance, this check might be omitted if inputs are guaranteed clean.
    if (b1 > 3 || b2 > 3) {
        // digram_ids[gid] = 255; // Example: Mark as invalid ID
        return; // Or skip writing for this gid
    }

    uchar id;
    if (reverse_aware) {
        id = get_canonical_digram_id(b1, b2);
    } else {
        id = (b1 << 2) | b2; // ID will be 0-15
    }
    digram_ids[gid] = id;
}

// Kernel to count digram occurrences using local memory for reduction.
// MAX_POSSIBLE_DIGRAM_IDS should be defined at compile time by the host (e.g., 16 for non-canonical 2-bit DNA).
// This kernel assumes digram_ids contains values from 0 to (MAX_POSSIBLE_DIGRAM_IDS - 1).
__kernel void count_digrams_by_id(
    __global const uchar *digram_ids,      // Input: array of compact digram IDs
    const uint num_ids,                   // Number of IDs in digram_ids array
    __global uint *global_counts,         // Output: counts for each digram ID. Size: MAX_POSSIBLE_DIGRAM_IDS
    const uint max_digram_id_value        // MAX_POSSIBLE_DIGRAM_IDS (e.g., 16)
) {
    const uint gid = get_global_id(0);    // Global ID
    const uint lid = get_local_id(0);     // Local ID (within work-group)
    const uint group_size = get_local_size(0); // Work-group size

    // Local memory for accumulating counts within this work-group.
    // Size must be known at compile time for static allocation.
    // We use a common upper bound like 16 for 2-bit digrams.
    __local uint local_counts[16]; // Max 16 digrams (0-15)

    // Initialize local_counts to 0.
    // Each work-item in the group initializes a part of local_counts.
    // This loop ensures all relevant parts of local_counts are zeroed.
    // If max_digram_id_value is small (e.g., 16), this is efficient.
    for (uint i = lid; i < max_digram_id_value; i += group_size) {
        local_counts[i] = 0;
    }
    barrier(CLK_LOCAL_MEM_FENCE); // Synchronize work-items in group: ensure initialization is complete.

    // Each work-item processes its share of the global digram_ids array
    // and updates the corresponding counter in local_counts.
    if (gid < num_ids) {
        uchar id = digram_ids[gid];
        if (id < max_digram_id_value) { // Check if ID is within bounds
            atomic_inc(&local_counts[id]); // Atomic operation on local memory (fast)
        }
    }
    barrier(CLK_LOCAL_MEM_FENCE); // Synchronize: ensure all local_counts updates are done.

    // After local counting, work-items in the group sum their local_counts
    // to the global_counts array in global memory.
    // Each work-item handles one or more IDs' sum to global memory.
    for (uint i = lid; i < max_digram_id_value; i += group_size) {
        if (local_counts[i] > 0) {
            atomic_add(&global_counts[i], local_counts[i]); // Atomic operation on global memory
        }
    }
}