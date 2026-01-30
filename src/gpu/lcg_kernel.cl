// lcg_kernel.cl
// GPU kernels for Locally Consistent Grammar (LCG) acceleration
//
// These kernels accelerate the two most computationally intensive parts of LCG:
// 1. Computing position fingerprints (digram-based hashes)
// 2. Finding local minima for cut point detection

// FxHash-compatible constants (matches rustc_hash)
#define FX_HASH_K 0x517cc1b727220a95UL

// Compute FxHash-compatible fingerprint for a digram (two 2-bit encoded bases)
inline ulong compute_digram_fingerprint(uchar b1, uchar b2) {
    ulong hash = 0UL;
    // First base
    hash = (hash ^ (ulong)b1) * FX_HASH_K;
    // Second base
    hash = (hash ^ (ulong)b2) * FX_HASH_K;
    return hash;
}

// Kernel 1: Compute position fingerprints for all positions in the sequence
// Each position gets a fingerprint based on the digram at that position (current + next base)
// This is the parallel version of LCG's fingerprint computation
__kernel void compute_position_fingerprints(
    __global const uchar *sequence,      // Input: 2-bit encoded sequence (0,1,2,3 per uchar)
    const uint sequence_len,             // Length of the sequence
    __global ulong *fingerprints         // Output: fingerprint for each position
) {
    const uint gid = get_global_id(0);

    if (gid >= sequence_len) return;

    if (gid + 1 < sequence_len) {
        // Compute digram fingerprint
        const uchar b1 = sequence[gid];
        const uchar b2 = sequence[gid + 1];
        fingerprints[gid] = compute_digram_fingerprint(b1, b2);
    } else {
        // Last position - just hash the single base
        fingerprints[gid] = (ulong)sequence[gid];
    }
}

// Kernel 2: Find local minima in the fingerprint array
// A position is a local minimum if its fingerprint is <= all fingerprints in the window
// This identifies potential cut points for LCG parsing
__kernel void find_local_minima(
    __global const ulong *fingerprints,  // Input: fingerprint array
    const uint num_fingerprints,         // Number of fingerprints
    const uint half_window,              // Half of the window size for local minimum detection
    const uint min_phrase_len,           // Minimum phrase length
    const uint max_phrase_len,           // Maximum phrase length
    __global uchar *is_cut_point,        // Output: 1 if this position is a cut point, 0 otherwise
    __global uint *cut_point_count       // Output: atomic counter for number of cut points
) {
    const uint gid = get_global_id(0);

    // Position 0 is always a cut point
    if (gid == 0) {
        is_cut_point[gid] = 1;
        atomic_inc(cut_point_count);
        return;
    }

    if (gid >= num_fingerprints) {
        is_cut_point[gid] = 0;
        return;
    }

    const ulong fp = fingerprints[gid];

    // Calculate window bounds
    const uint window_start = (gid > half_window) ? (gid - half_window) : 0;
    const uint window_end = min(gid + half_window + 1, num_fingerprints);

    // Check if this is a local minimum
    uchar is_local_min = 1;
    for (uint i = window_start; i < window_end && is_local_min; i++) {
        if (fingerprints[i] < fp) {
            is_local_min = 0;
        }
    }

    // Determine if this should be a cut point
    // We mark it as a cut point if:
    // 1. It's a local minimum and we've accumulated enough for min_phrase_len
    // 2. We've hit the max_phrase_len limit
    // Note: The actual phrase length validation is done on the host side
    // This kernel just marks potential cut points

    is_cut_point[gid] = is_local_min;

    if (is_local_min) {
        atomic_inc(cut_point_count);
    }
}

// Kernel 3: Compact cut points into a dense array
// After finding all cut points, this kernel compacts them for efficient host-side processing
__kernel void compact_cut_points(
    __global const uchar *is_cut_point,  // Input: boolean array of cut points
    const uint num_positions,            // Total number of positions
    __global uint *cut_points,           // Output: compact array of cut point positions
    __global uint *cut_point_index       // Atomic index for compact array
) {
    const uint gid = get_global_id(0);

    if (gid >= num_positions) return;

    if (is_cut_point[gid]) {
        const uint idx = atomic_inc(cut_point_index);
        cut_points[idx] = gid;
    }
}

// Kernel 4: Compute phrase fingerprints
// Given cut points, compute the fingerprint for each phrase
// This is used for rule identification and deduplication
__kernel void compute_phrase_fingerprints(
    __global const uchar *sequence,      // Input: 2-bit encoded sequence
    const uint sequence_len,             // Length of the sequence
    __global const uint *cut_points,     // Input: array of cut point positions
    const uint num_phrases,              // Number of phrases (= num_cut_points for most cases)
    __global ulong *phrase_fingerprints  // Output: fingerprint for each phrase
) {
    const uint gid = get_global_id(0);

    if (gid >= num_phrases) return;

    const uint start = cut_points[gid];
    const uint end = (gid + 1 < num_phrases) ? cut_points[gid + 1] : sequence_len;

    if (start >= end) {
        phrase_fingerprints[gid] = 0UL;
        return;
    }

    // Compute FxHash-compatible fingerprint for the phrase
    ulong hash = 0UL;
    for (uint i = start; i < end; i++) {
        hash = (hash ^ (ulong)sequence[i]) * FX_HASH_K;
    }
    phrase_fingerprints[gid] = hash;
}
