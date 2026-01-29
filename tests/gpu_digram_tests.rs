//! tests/gpu_digram_tests.rs

use anyhow::Result;
use std::collections::HashMap;
use std::hash::{Hash, Hasher};
use rustc_hash::FxHasher;

use orbweaver::{
    grammar::digram_table::{DigramKey, DigramTable},
    grammar::symbol::{Symbol, SymbolType, Direction},
    encode::dna_2bit::EncodedBase,
    gpu::{
        digram::{count_digrams_gpu, GpuSequence, CANONICAL_KERNEL_ID_TO_BASES, MAX_POSSIBLE_DIGRAM_IDS},
        GpuContext,
    },
};


// Helper to create a DigramKey (hash) from two chars and reverse_aware flag
fn get_key_hash(s1_char: char, s2_char: char, reverse_aware: bool) -> DigramKey {
    let s1_base = EncodedBase::from_base(s1_char as u8).expect("Invalid char for s1_base");
    let s2_base = EncodedBase::from_base(s2_char as u8).expect("Invalid char for s2_base");
    // Create dummy symbols for canonical_key. ID and direction don't affect the key for terminals normally.
    let sym1 = Symbol::terminal(0, s1_base, Direction::Forward, None, None);
    let sym2 = Symbol::terminal(0, s2_base, Direction::Forward, None, None);
    
    let tuple_key = DigramTable::canonical_key((&sym1, &sym2), reverse_aware);
    let mut hasher = FxHasher::default();
    tuple_key.hash(&mut hasher);
    hasher.finish()
}

// Helper function to initialize GpuContext for tests
// This might need adjustment based on how GpuContext is typically created and if it requires specific setup.
fn setup_gpu_context() -> Result<GpuContext> {
    GpuContext::new() // Assuming default context is okay for these tests
}

#[test]
fn test_canonical_kernel_id_to_bases_mapping() {
    fn term(base_char: char) -> Symbol {
    let base = EncodedBase::from_base(base_char as u8).expect("Invalid char for EncodedBase in test helper");
        Symbol::terminal(0, base, Direction::Forward, None, None)
    }

    // Expected mappings for the 10 canonical digrams based on the kernel logic
    // (ID, (Base1, Base2))
    let expected_canonical_mappings: HashMap<u32, (Symbol, Symbol)> = [
        (0,  (term('A'), term('A'))), // AA
        (1,  (term('A'), term('C'))), // AC
        (2,  (term('A'), term('G'))), // AG
        (3,  (term('A'), term('T'))), // AT
        (4,  (term('C'), term('A'))), // CA
        (5,  (term('C'), term('C'))), // CC
        (6,  (term('C'), term('G'))), // CG
        // ID 7 (CT) maps to AG (ID 2), so CANONICAL_KERNEL_ID_TO_BASES[7] should be N,N
        (8,  (term('G'), term('A'))), // GA
        (9,  (term('G'), term('C'))), // GC
        // ID 10 (GG) maps to CC (ID 5), so CANONICAL_KERNEL_ID_TO_BASES[10] should be N,N
        // ID 11 (GT) maps to AC (ID 1), so CANONICAL_KERNEL_ID_TO_BASES[11] should be N,N
        (12, (term('T'), term('A'))), // TA
        // ID 13 (TC) maps to GA (ID 8), so CANONICAL_KERNEL_ID_TO_BASES[13] should be N,N
        // ID 14 (TG) maps to CA (ID 4), so CANONICAL_KERNEL_ID_TO_BASES[14] should be N,N
        // ID 15 (TT) maps to AA (ID 0), so CANONICAL_KERNEL_ID_TO_BASES[15] should be N,N
    ].iter().map(|(id, (s1, s2))| (*id, (s1.clone(), s2.clone()))).collect();

    assert_eq!(CANONICAL_KERNEL_ID_TO_BASES.len(), MAX_POSSIBLE_DIGRAM_IDS as usize, "Mapping array size should match MAX_POSSIBLE_DIGRAM_IDS");

    for id in 0..(MAX_POSSIBLE_DIGRAM_IDS as usize) {
        let (s1_map, s2_map) = CANONICAL_KERNEL_ID_TO_BASES[id];
        if let Some((s1_exp, s2_exp)) = expected_canonical_mappings.get(&(id as u32)) {
            let s1_map_char = match s1_map.symbol_type { SymbolType::Terminal(eb) => Some(eb.to_char()), _ => None }.expect("s1_map not terminal or not char");
                let s1_exp_char = match s1_exp.symbol_type { SymbolType::Terminal(eb) => Some(eb.to_char()), _ => None }.expect("s1_exp not terminal or not char");
                assert_eq!(s1_map_char, s1_exp_char, "Mismatch for ID {} base1", id);
            let s2_map_char = match s2_map.symbol_type { SymbolType::Terminal(eb) => Some(eb.to_char()), _ => None }.expect("s2_map not terminal or not char");
                let s2_exp_char = match s2_exp.symbol_type { SymbolType::Terminal(eb) => Some(eb.to_char()), _ => None }.expect("s2_exp not terminal or not char");
                assert_eq!(s2_map_char, s2_exp_char, "Mismatch for ID {} base2", id);
        } else {
            // For IDs not in our expected_canonical_mappings, they should be (N, N)
            assert_eq!(s1_map.symbol_type, SymbolType::Terminal(EncodedBase(0b11)), "Expected placeholder (T) for base1 of ID {}", id);
            assert_eq!(s2_map.symbol_type, SymbolType::Terminal(EncodedBase(0b11)), "Expected placeholder (T) for base2 of ID {}", id);
        }
    }
}

// Helper to create a GpuSequence from a string of bases like "ACGT"
fn sequence_from_str(s: &str) -> Result<GpuSequence> {
    let symbols: Vec<Symbol> = s.chars().enumerate().map(|(i, c)| {
        let base = match c {
            'A' => EncodedBase(0b00),
            'C' => EncodedBase(0b01),
            'G' => EncodedBase(0b10),
            'T' => EncodedBase(0b11),
            _ => panic!("Invalid base in test sequence string"),
        };
        Symbol::terminal(i, base, Direction::Forward, None, None)
    }).collect();
    GpuSequence::from_symbols(&symbols)
}

// Helper to compare GpuDigramCounts.counts with an expected HashMap<DigramKey, usize>
fn assert_counts_match(actual_counts: &HashMap<orbweaver::grammar::digram_table::DigramKey, usize>, expected_counts: &HashMap<orbweaver::grammar::digram_table::DigramKey, usize>) {
    assert_eq!(actual_counts, expected_counts, "Digram counts do not match expected counts.");
}

#[test]
fn test_count_digrams_gpu_empty_sequence() -> Result<()> {
    let gpu_context = setup_gpu_context()?;
    let sequence = sequence_from_str("")?;
    let result = count_digrams_gpu(&sequence, false, &gpu_context)?;
    assert!(result.counts.is_empty(), "Counts should be empty for an empty sequence");
    Ok(())
}

#[test]
fn test_count_digrams_gpu_simple_no_reverse_aware() -> Result<()> {
    let gpu_context = setup_gpu_context()?;
    let sequence = sequence_from_str("ACGTAC")?;
    // Expected: AC:2, CG:1, GT:1, TA:1
    let result = count_digrams_gpu(&sequence, false, &gpu_context)?;
    
    let mut expected: HashMap<DigramKey, usize> = HashMap::new();
    expected.insert(get_key_hash('A', 'C', false), 2);
    expected.insert(get_key_hash('C', 'G', false), 1);
    expected.insert(get_key_hash('G', 'T', false), 1);
    expected.insert(get_key_hash('T', 'A', false), 1);
    assert_counts_match(&result.counts, &expected);
    Ok(())
}

#[test]
fn test_count_digrams_gpu_simple_reverse_aware() -> Result<()> {
    let gpu_context = setup_gpu_context()?;
    let sequence = sequence_from_str("ACGTAC")?;
    // Sequence: AC, CG, GT, TA, AC
    // AC (original) -> AC key
    // CG (original) -> CG key
    // GT (rev_comp of AC) -> AC key (because reverse_aware=true)
    // TA (original) -> TA key
    // AC (original) -> AC key
    // Expected: AC:3, CG:1, TA:1
    let result = count_digrams_gpu(&sequence, true, &gpu_context)?;
    
    let mut expected: HashMap<DigramKey, usize> = HashMap::new();
    expected.insert(get_key_hash('A', 'C', true), 3);
    expected.insert(get_key_hash('C', 'G', true), 1);
    expected.insert(get_key_hash('T', 'A', true), 1);
    assert_counts_match(&result.counts, &expected);
    Ok(())
}

#[test]
fn test_count_digrams_gpu_single_digram_type() -> Result<()> {
    let gpu_context = setup_gpu_context()?;
    let sequence = sequence_from_str("AAAAA")?;
    // Expected: AA:4
    let result_no_reverse = count_digrams_gpu(&sequence, false, &gpu_context)?;
    let mut expected_aa: HashMap<DigramKey, usize> = HashMap::new();
    expected_aa.insert(get_key_hash('A', 'A', false), 4);
    assert_counts_match(&result_no_reverse.counts, &expected_aa);

    // With reverse_aware=true, AA is canonical to TT. Kernel output is ID 0 (AA).
    // So, count should still be for AA.
    let result_reverse = count_digrams_gpu(&sequence, true, &gpu_context)?;
    // With reverse_aware=true, AA is canonical to TT. Kernel output is ID 0 (AA).
    // So, count should still be for AA's hash when reverse_aware is true for AA itself.
    let mut expected_aa_rev: HashMap<DigramKey, usize> = HashMap::new();
    expected_aa_rev.insert(get_key_hash('A', 'A', true), 4);
    assert_counts_match(&result_reverse.counts, &expected_aa_rev);
    Ok(())
}

#[test]
fn test_count_digrams_gpu_reverse_complements() -> Result<()> {
    let gpu_context = setup_gpu_context()?;
    let sequence = sequence_from_str("ACGTGTCA")?;
    // Sequence: AC, CG, GT, TG, GT, TC, CA
    // No reverse aware:
    // AC:1, CG:1, GT:2, TG:1, TC:1, CA:1
    let result_no_reverse = count_digrams_gpu(&sequence, false, &gpu_context)?;
    let mut expected_no_reverse: HashMap<DigramKey, usize> = HashMap::new();
    expected_no_reverse.insert(get_key_hash('A', 'C', false), 1);
    expected_no_reverse.insert(get_key_hash('C', 'G', false), 1);
    expected_no_reverse.insert(get_key_hash('G', 'T', false), 2);
    expected_no_reverse.insert(get_key_hash('T', 'G', false), 1);
    expected_no_reverse.insert(get_key_hash('T', 'C', false), 1);
    expected_no_reverse.insert(get_key_hash('C', 'A', false), 1);
    assert_counts_match(&result_no_reverse.counts, &expected_no_reverse);

    // Reverse aware:
    // AC (AC, GT, GT) -> AC:3 (AC is ID 1, GT is ID 11 which maps to AC)
    // CG (CG) -> CG:1 (CG is ID 6)
    // TG (TG) -> CA:1 (TG is ID 14 which maps to CA, CA is ID 4)
    // TC (TC) -> GA:1 (TC is ID 13 which maps to GA, GA is ID 8)
    // CA (CA) -> CA:1 (CA is ID 4)
    // Combined: AC:3, CG:1, CA:2, GA:1
    let result_reverse = count_digrams_gpu(&sequence, true, &gpu_context)?;
    let mut expected_reverse: HashMap<DigramKey, usize> = HashMap::new();
    // AC (AC, GT, GT) -> AC:3 (AC is ID 1, GT is ID 11 which maps to AC)
    expected_reverse.insert(get_key_hash('A', 'C', true), 3);
    // CG (CG) -> CG:1 (CG is ID 6)
    expected_reverse.insert(get_key_hash('C', 'G', true), 1);
    // CA (CA, TG) -> CA:2 (CA is ID 4, TG is ID 14 which maps to CA)
    expected_reverse.insert(get_key_hash('C', 'A', true), 2);
    // GA (TC) -> GA:1 (TC is ID 13 which maps to GA, GA is ID 8)
    expected_reverse.insert(get_key_hash('G', 'A', true), 1);
    assert_counts_match(&result_reverse.counts, &expected_reverse);
    Ok(())
}

#[test]
fn test_count_digrams_gpu_longer_sequence_with_ns_in_mapping() -> Result<()> {
    // Test with digrams that map to N,N in CANONICAL_KERNEL_ID_TO_BASES
    // e.g., CT (ID 7), GG (ID 10), TT (ID 15)
    // When reverse_aware=true, these are canonicalized by the kernel:
    // CT -> AG (ID 2)
    // GG -> CC (ID 5)
    // TT -> AA (ID 0)
    let gpu_context = setup_gpu_context()?;
    let sequence = sequence_from_str("CTGGTT")?;
    // Sequence: CT, TG, GG, GT, TT
    // No reverse aware:
    // CT:1, TG:1, GG:1, GT:1, TT:1
    let result_no_reverse = count_digrams_gpu(&sequence, false, &gpu_context)?;
    let mut expected_no_reverse: HashMap<DigramKey, usize> = HashMap::new();
    expected_no_reverse.insert(get_key_hash('C', 'T', false), 1);
    expected_no_reverse.insert(get_key_hash('T', 'G', false), 1);
    expected_no_reverse.insert(get_key_hash('G', 'G', false), 1);
    expected_no_reverse.insert(get_key_hash('G', 'T', false), 1);
    expected_no_reverse.insert(get_key_hash('T', 'T', false), 1);
    assert_counts_match(&result_no_reverse.counts, &expected_no_reverse);

    // Reverse aware:
    // CT -> AG (ID 2) -> AG:1
    // TG -> CA (ID 4) -> CA:1
    // GG -> CC (ID 5) -> CC:1
    // GT -> AC (ID 1) -> AC:1
    // TT -> AA (ID 0) -> AA:1
    let result_reverse = count_digrams_gpu(&sequence, true, &gpu_context)?;
    let mut expected_reverse: HashMap<DigramKey, usize> = HashMap::new();
    // CT -> AG (ID 2) -> AG:1
    expected_reverse.insert(get_key_hash('A', 'G', true), 1);
    // TG -> CA (ID 4) -> CA:1
    expected_reverse.insert(get_key_hash('C', 'A', true), 1);
    // GG -> CC (ID 5) -> CC:1
    expected_reverse.insert(get_key_hash('C', 'C', true), 1);
    // GT -> AC (ID 1) -> AC:1
    expected_reverse.insert(get_key_hash('A', 'C', true), 1);
    // TT -> AA (ID 0) -> AA:1
    expected_reverse.insert(get_key_hash('A', 'A', true), 1);
    assert_counts_match(&result_reverse.counts, &expected_reverse);
    Ok(())
}
