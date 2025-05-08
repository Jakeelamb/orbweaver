__kernel void compute_digram_hashes(__global const uchar *sequence, const uint sequence_len, const uint reverse_aware, __global ulong *hashes) {
    const uint gid = get_global_id(0);
    if (gid >= sequence_len - 1) return;

    const uchar b1 = sequence[gid];
    const uchar b2 = sequence[gid + 1];
    ulong hash = ((ulong)b1 << 32) | b2;

    if (reverse_aware) {
        // DNA complement: A(00)<->T(11), C(01)<->G(10)
        const uchar rc_b1 = 3 - b2; // Complement of b2
        const uchar rc_b2 = 3 - b1; // Complement of b1
        const ulong rc_hash = ((ulong)rc_b1 << 32) | rc_b2;
        
        // Use canonical form (minimum of hash and rc_hash)
        hash = min(hash, rc_hash);
    }

    hashes[gid] = hash;
}

__kernel void count_digrams(__global const ulong *hashes, const uint num_hashes, __global uint *counts, const uint counts_size) {
    const uint gid = get_global_id(0);
    if (gid >= num_hashes) return;

    const ulong hash = hashes[gid];
    if (hash != 0) {
        atomic_inc(&counts[hash % counts_size]);
    }
} 