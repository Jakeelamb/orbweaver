// Placeholder for the parallel grammar engine

use crate::encode::dna_2bit::EncodedBase;
use crate::fasta::reader::SequenceChunk;
use crate::grammar::engine::{Grammar, Sequitur};
use crate::grammar::rule::Rule;
use crate::grammar::symbol::{Symbol, SymbolType, Direction};
use crate::utils::hash::canonical_hash_symbols;
use crate::parallel::chunking::{ChunkingConfig, split_into_chunks, Chunk};
use crate::grammar::builder::GrammarBuilder;
use anyhow::{Result, Context};
use rayon::prelude::*;
use std::collections::{HashMap, BTreeMap};
use std::time::Instant;
use union_find::UnionFind;
use union_find::QuickUnionUf;
use union_find::UnionBySize;

/// Metrics for the parallel Sequitur algorithm
#[derive(Debug, Default)]
pub struct ParallelMetrics {
    /// Time spent processing individual chunks
    pub chunk_processing_time: std::time::Duration,
    /// Time spent merging grammars
    pub merge_time: std::time::Duration,
    /// Time spent deduplicating rules
    pub deduplication_time: std::time::Duration,
    /// Total time spent in the algorithm
    pub total_time: std::time::Duration,
    /// Number of chunks processed
    pub chunk_count: usize,
    /// Number of rules before merging
    pub rules_before_merge: usize,
    /// Number of rules after merging
    pub rules_after_merge: usize,
}

/// Runs the Sequitur algorithm in parallel on chunked input
pub fn parallel_sequitur(
    sequence: &[EncodedBase],
    config: ChunkingConfig,
) -> Result<(Grammar, ParallelMetrics)> {
    let total_start = Instant::now();
    let mut metrics: ParallelMetrics = ParallelMetrics::default();

    // Configure Rayon thread pool
    rayon::ThreadPoolBuilder::new()
        .num_threads(config.num_threads)
        .build_global()?;

    println!("Using {} threads for parallel processing.", config.num_threads);

    // --- Step 1: Split sequence into chunks ---
    // Assuming split_into_chunks is defined elsewhere and returns Vec<Chunk>
    let chunks = split_into_chunks(sequence, &config);
    let num_chunks = chunks.len();
    metrics.chunk_count = num_chunks; // Track chunk count
    println!("Split sequence into {} chunks.", num_chunks);

    // --- Step 2: Process chunks in parallel ---
    let chunk_processing_start = Instant::now();
    let chunk_results: Vec<Result<Grammar>> = chunks
        .into_par_iter()
        .map(|chunk| {
             let mut builder = GrammarBuilder::new(config.min_rule_usage, config.reverse_aware);
             builder.build_grammar(&chunk.data)?; // Pass chunk data to builder
             let (seq, rules) = builder.get_grammar();
             let depth = builder.get_max_rule_depth();
             Ok(Grammar { sequence: seq.clone(), rules: rules.clone(), max_depth: depth })
        })
        .collect();
    metrics.chunk_processing_time = chunk_processing_start.elapsed();

    // Handle potential errors from chunk processing
    let mut chunk_grammars = Vec::with_capacity(num_chunks);
    for result in chunk_results {
        chunk_grammars.push(result.with_context(|| "Error processing sequence chunk")?);
    }

    // --- Step 3: Merge grammars ---
    // Correctly call merge_grammars and destructure the result
    let (merged_grammar, merge_metrics) = merge_grammars(chunk_grammars, &config, sequence.len())?;

    // Combine metrics
    metrics.merge_time = merge_metrics.merge_time;
    metrics.deduplication_time = merge_metrics.deduplication_time;
    metrics.rules_before_merge = merge_metrics.rules_before_merge;
    metrics.rules_after_merge = merge_metrics.rules_after_merge; // Use the value from merge_metrics
    metrics.total_time = total_start.elapsed();

    Ok((merged_grammar, metrics))
}

/// Merges multiple grammars into a single grammar
fn merge_grammars(grammars: Vec<Grammar>, config: &ChunkingConfig, total_sequence_len: usize) -> Result<(Grammar, ParallelMetrics)> {
    let start_merge = Instant::now();
    let mut metrics = ParallelMetrics::default();

    if grammars.is_empty() {
        // Manually create default Grammar and return with metrics
        let default_grammar = Grammar {
            sequence: Vec::new(),
            rules: HashMap::new(),
            max_depth: 0,
        };
        return Ok((default_grammar, metrics));
    }

    // Deduplicate rules using the union-find approach
    let deduplication_start = Instant::now();
    let (rule_map, merged_rules) = deduplicate_rules(&grammars);
    metrics.deduplication_time = deduplication_start.elapsed();
    // Calculate rules_before_merge here using input `grammars`
    metrics.rules_before_merge = grammars.iter().map(|g| g.rules.len()).sum();
    metrics.rules_after_merge = merged_rules.len(); // Corrected access

    // --- Step 5: Reconstruct the final sequence by stitching chunks ---
    let mut final_sequence = Vec::with_capacity(total_sequence_len); // Estimate capacity
    let mut current_original_pos = 0; // Track position in the original, uncompressed sequence

    for (chunk_index, grammar) in grammars.iter().enumerate() {
        let grammar_ptr = grammar as *const Grammar;
        let chunk_start_in_original = grammar.sequence.first().map_or(current_original_pos, |s| s.id); // Approximate start pos
        let chunk_seq = &grammar.sequence;

        // Determine the segment of the chunk sequence to append
        let sequence_segment_to_add = if chunk_index == 0 {
            // For the first chunk, take the whole sequence
            chunk_seq
        } else {
            // --- Simplified Overlap Handling ---
            let target_start_pos_original = current_original_pos;
            let start_index_in_chunk = chunk_seq.iter().position(|symbol| {
                 symbol.id >= target_start_pos_original
            }).unwrap_or(0);

            if start_index_in_chunk < chunk_seq.len() {
                &chunk_seq[start_index_in_chunk..]
            } else {
                &[]
            }
        };

        // Append the relevant segment, remapping rule IDs
        for symbol in sequence_segment_to_add {
            let remapped_symbol = match symbol.symbol_type {
                SymbolType::Terminal(_) => symbol.clone(),
                SymbolType::NonTerminal(old_rule_id) => {
                    if let Some(&new_rule_id) = rule_map.get(&(grammar_ptr, old_rule_id)) {
                        Symbol::non_terminal(symbol.id, new_rule_id, symbol.strand)
                    } else {
                        eprintln!("Warning: Rule ID {} from grammar {:?} not found in map during sequence reconstruction.", old_rule_id, grammar_ptr);
                        symbol.clone()
                    }
                }
            };
            final_sequence.push(remapped_symbol);
        }

        // Update the current position in the original sequence (simplified update)
        if let Some(last_symbol) = sequence_segment_to_add.last() {
             current_original_pos = last_symbol.id + 1;
        }
    }

    println!("Warning: Final sequence reconstruction logic is complex and uses heuristics for overlap handling. Length might deviate slightly. Final length: {}", final_sequence.len());

    // Calculate merge time excluding deduplication
    metrics.merge_time = start_merge.elapsed().saturating_sub(metrics.deduplication_time);

    // Calculate final max depth
    let max_depth = calculate_max_depth(&merged_rules);

    Ok((
        Grammar {
            sequence: final_sequence,
            rules: merged_rules,
            max_depth,
        },
        metrics,
    ))
}

/// Deduplicate rules across multiple grammars using Union-Find
fn deduplicate_rules(grammars: &[Grammar]) -> (HashMap<(*const Grammar, usize), usize>, HashMap<usize, Rule>) {
    // --- Step 1: Collect all rules and map old IDs to a temporary global ID ---\
    let mut old_to_temp_global = HashMap::new(); // Maps (grammar_ptr, old_rule_id) -> temp_global_id
    // Stores (rule_hash, grammar_ptr, old_rule_id, symbols_vec) by temp_global_id
    let mut temp_global_to_info: Vec<(u64, *const Grammar, usize, Vec<Symbol>)> = Vec::new();
    let mut temp_global_counter = 0;

    for grammar in grammars {
        for (old_rule_id, rule) in &grammar.rules {
            let grammar_ptr = grammar as *const Grammar;
            let rule_hash = canonical_hash_symbols(&rule.symbols);

            old_to_temp_global.insert((grammar_ptr, *old_rule_id), temp_global_counter);
            // Store symbols instead of the whole rule clone initially
            temp_global_to_info.push((rule_hash, grammar_ptr, *old_rule_id, rule.symbols.clone()));
            temp_global_counter += 1;
        }
    }

    let num_rules = temp_global_counter;
    if num_rules == 0 {
        return (HashMap::new(), HashMap::new()); // No rules to merge
    }

    // --- Step 2: Initialize Union-Find and map rules by hash ---
    let mut uf: QuickUnionUf<UnionBySize> = QuickUnionUf::new(num_rules);
    let mut hash_to_representative_temp_id = HashMap::new();

    for temp_global_id in 0..num_rules {
        let (rule_hash, _, _, _) = temp_global_to_info[temp_global_id];

        if let Some(representative_temp_id) = hash_to_representative_temp_id.get(&rule_hash) {
            uf.union(temp_global_id, *representative_temp_id);
        } else {
            hash_to_representative_temp_id.insert(rule_hash, temp_global_id);
        }
    }

    // --- Step 3: Determine final representative rules and assign new IDs ---
    let mut final_rule_map = HashMap::new(); // Maps (grammar_ptr, old_rule_id) -> final_new_rule_id
    let mut final_merged_rules_info = BTreeMap::new(); // Stores (final_new_id) -> (representative_grammar_ptr, representative_old_id, representative_symbols)
    let mut representative_to_new_id = HashMap::new(); // Maps temp_representative_id -> final_new_rule_id
    let mut next_final_id_counter = 0;

    for temp_global_id in 0..num_rules {
        let representative_temp_id = uf.find(temp_global_id);
        let (_, grammar_ptr, old_rule_id, _) = temp_global_to_info[temp_global_id].clone(); // Don't need symbols here

        // Assign a final ID if this representative hasn't been seen yet
        let final_new_rule_id = *representative_to_new_id.entry(representative_temp_id).or_insert_with(|| {
            let new_id = next_final_id_counter;
            next_final_id_counter += 1;

            // Get the info corresponding to the representative temp ID
            let (_, rep_grammar_ptr, rep_old_id, rep_symbols) = temp_global_to_info[representative_temp_id].clone();

            // Store the representative's original context and symbols with the new final ID
            final_merged_rules_info.insert(new_id, (rep_grammar_ptr, rep_old_id, rep_symbols));

            new_id
        });

        // Map the original rule to its final new representative ID
        final_rule_map.insert((grammar_ptr, old_rule_id), final_new_rule_id);
    }

    // --- Step 4: Create final rules and update non-terminal references ---
    let mut final_merged_rules_map: HashMap<usize, Rule> = HashMap::new();

    for (final_id, (rep_grammar_ptr, _rep_old_id, rep_symbols)) in final_merged_rules_info {
        let mut final_symbols = Vec::with_capacity(rep_symbols.len());
        let mut needs_update = false; // Track if any symbol was updated

        for original_symbol in rep_symbols {
            if let SymbolType::NonTerminal(old_child_id) = original_symbol.symbol_type {
                // Use the representative's context (rep_grammar_ptr) and the child's old ID
                // to find the final ID in the map.
                if let Some(final_new_child_id) = final_rule_map.get(&(rep_grammar_ptr, old_child_id)) {
                    if *final_new_child_id != old_child_id { // Avoid redundant clones if ID is same
                         // Create a new symbol with the updated ID
                         final_symbols.push(Symbol {
                            id: original_symbol.id, // Keep original instance ID? Or update? Let's keep for now.
                            symbol_type: SymbolType::NonTerminal(*final_new_child_id),
                            strand: original_symbol.strand,
                         });
                         needs_update = true;
                    } else {
                        final_symbols.push(original_symbol); // No change needed
                    }
                } else {
                    // Child ID not found in map - should not happen in a consistent grammar merge.
                    eprintln!("Warning: Could not remap non-terminal ID {} (from grammar {:?}) in rule {}", old_child_id, rep_grammar_ptr, final_id);
                    final_symbols.push(original_symbol); // Keep original symbol
                }
            } else {
                // Terminal symbol, add as is
                final_symbols.push(original_symbol);
            }
        }

        // Create the final rule object
        // TODO: Recalculate usage_count and positions based on merging if needed
        final_merged_rules_map.insert(final_id, Rule {
            id: final_id,
            symbols: final_symbols,
            usage_count: 0, // Placeholder - usage needs recalculation
            positions: Vec::new(), // Placeholder - positions need recalculation/merging
            depth: None, // Placeholder - depth needs recalculation
        });
    }

    (final_rule_map, final_merged_rules_map)
}

/// Calculate the maximum rule depth in a grammar
fn calculate_max_depth(rules: &HashMap<usize, Rule>) -> usize {
    let mut max_depth = 0;
    let mut depths = HashMap::new();
    
    for rule_id in rules.keys() {
        max_depth = max_depth.max(calculate_rule_depth(*rule_id, rules, &mut depths));
    }
    
    max_depth
}

/// Calculate the depth of a specific rule
fn calculate_rule_depth(
    rule_id: usize, 
    rules: &HashMap<usize, Rule>,
    depths: &mut HashMap<usize, usize>
) -> usize {
    // Check if we've already calculated this rule's depth
    if let Some(&depth) = depths.get(&rule_id) {
        return depth;
    }
    
    // Base case: rule not found
    let rule = match rules.get(&rule_id) {
        Some(r) => r,
        None => return 0,
    };
    
    // Find the maximum depth of any non-terminal in this rule
    let mut max_child_depth = 0;
    for symbol in &rule.symbols {
        if let SymbolType::NonTerminal(child_id) = symbol.symbol_type {
            // Skip self-references to avoid infinite recursion
            if child_id != rule_id {
                max_child_depth = max_child_depth.max(calculate_rule_depth(child_id, rules, depths));
            }
        }
    }
    
    // This rule's depth is 1 + the max depth of its children
    let depth = 1 + max_child_depth;
    depths.insert(rule_id, depth);
    
    depth
}

/// Convert a SequenceChunk to our new Chunk type
fn convert_chunk(chunk: &SequenceChunk) -> Chunk {
    // Convert the Vec<u8> to Vec<EncodedBase>
    let encoded_data: Vec<EncodedBase> = chunk.data.iter()
        .filter_map(|&b| EncodedBase::from_base(b))
        .collect();

    Chunk {
        index: 0, // We don't have this information in SequenceChunk
        data: encoded_data,
        start: chunk.start_pos,
        end: chunk.end_pos,
        is_first: chunk.start_pos == 0,
        is_last: chunk.is_last,
    }
}

/// Creates overlapping chunks from a sequence of bases
fn create_chunks(bases: &[EncodedBase], chunk_size: usize, overlap: usize) -> Vec<SequenceChunk> {
    let mut chunks = Vec::new();
    let mut start = 0;
    
    while start < bases.len() {
        let end = (start + chunk_size).min(bases.len());
        let is_last = end == bases.len();
        
        // Convert EncodedBase chunk back to Vec<u8> for SequenceChunk
        let data_u8: Vec<u8> = bases[start..end].iter().map(|b| b.to_char() as u8).collect();

        chunks.push(SequenceChunk {
            data: data_u8,
            start_pos: start,
            end_pos: end,
            record_id: String::new(),
            is_last,
        });
        
        if is_last {
            break;
        }
        
        // Move start position for next chunk, ensuring overlap
        // Ensure we don't subtract past zero with saturating_sub
        start = end.saturating_sub(overlap.min(chunk_size.saturating_sub(1)));
        
        // Prevent infinite loops if overlap is too large or start doesn't advance
        if start >= end.saturating_sub(1) && !is_last { 
             start = end; // Force advance if stuck
        }
    }
    
    chunks
}

#[cfg(test)]
mod tests {
    use super::*;
    
    fn make_bases(seq: &str) -> Vec<EncodedBase> {
        seq.chars()
            .filter_map(|c| match c {
                'A' | 'a' => Some(EncodedBase(0)),
                'C' | 'c' => Some(EncodedBase(1)),
                'G' | 'g' => Some(EncodedBase(2)),
                'T' | 't' => Some(EncodedBase(3)),
                _ => None,
            })
            .collect()
    }
    
    #[test]
    fn test_create_chunks() {
        let bases = make_bases("ACGTACGTACGT");
        
        // Chunk size = 5, overlap = 2
        let chunks = create_chunks(&bases, 5, 2);
        
        assert_eq!(chunks.len(), 4);
        
        // First chunk: 0-5 (ACGTA)
        assert_eq!(chunks[0].data.len(), 5);
        assert_eq!(chunks[0].start_pos, 0);
        assert_eq!(chunks[0].end_pos, 5);
        assert!(!chunks[0].is_last);
        
        // Second chunk: 3-8 (TACGT)
        assert_eq!(chunks[1].data.len(), 5);
        assert_eq!(chunks[1].start_pos, 3);
        assert_eq!(chunks[1].end_pos, 8);
        assert!(!chunks[1].is_last);
        
        // Third chunk: 6-11 (ACGTA)
        assert_eq!(chunks[2].data.len(), 5);
        assert_eq!(chunks[2].start_pos, 6);
        assert_eq!(chunks[2].end_pos, 11);
        assert!(!chunks[2].is_last);
        
        // Fourth chunk: 9-12 (CGT)
        assert_eq!(chunks[3].data.len(), 3);
        assert_eq!(chunks[3].start_pos, 9);
        assert_eq!(chunks[3].end_pos, 12);
        assert!(chunks[3].is_last);
    }
    
    #[test]
    fn test_parallel_sequitur() -> Result<()> {
        let bases = make_bases("ACGTACGTACGT");
        
        let config = ChunkingConfig {
            chunk_size: 5,
            overlap_size: 2,
            min_rule_usage: 2,
            reverse_aware: false,
            num_threads: 2,
            show_progress: false,
        };
        
        let (grammar, metrics) = parallel_sequitur(&bases, config)?;
        
        // Simple test to ensure we get a valid grammar
        assert!(!grammar.sequence.is_empty());
        assert!(metrics.chunk_count > 0);
        
        Ok(())
    }
    
    #[test]
    fn test_merge_grammars() -> Result<()> {
        // TODO: Implement test for merging grammars
        Ok(())
    }
} 