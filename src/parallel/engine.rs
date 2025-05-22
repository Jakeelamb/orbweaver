// Placeholder for the parallel grammar engine

use crate::encode::dna_2bit::EncodedBase;
use crate::grammar::engine::Grammar;
use crate::grammar::rule::Rule;
use crate::grammar::symbol::{Symbol, SymbolType};
use crate::utils::hash::canonical_hash_symbols;
use crate::parallel::chunking::{ChunkingConfig, split_into_chunks};
use crate::grammar::builder::GrammarBuilder;
use anyhow::{Result, Context};
use rayon::prelude::*;
use std::collections::HashMap;
use std::time::Instant;
use union_find::UnionFind;
use union_find::QuickUnionUf;
use union_find::UnionBySize;
use crate::gpu::GpuContext; // Import GpuContext
use crate::analysis::assembly_index::calculate_rule_assembly_indices; // Added for Assembly Index

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
    gpu_context_param: Option<&GpuContext> 
) -> Result<(Grammar, ParallelMetrics)> {
    let total_start = Instant::now();
    let mut metrics: ParallelMetrics = ParallelMetrics::default();

    // Configure Rayon thread pool, ignore error if already initialized
    let _ = rayon::ThreadPoolBuilder::new()
        .num_threads(config.num_threads)
        .build_global();

    println!("Using {} threads for parallel processing.", config.num_threads);
    
    if config.use_gpu {
        println!("GPU acceleration enabled for chunk processing");
    }

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
             
             // Configure the builder
             if let Some(max_rules) = config.max_memory_per_chunk {
                 builder = builder.with_max_rules(max_rules / 100); 
             }
             
             // Build grammar using GPU or CPU based on config and context availability
             if config.use_gpu && gpu_context_param.is_some() {
                 builder = builder.with_gpu(gpu_context_param);
                 builder.build_grammar_with_gpu(&chunk.data, chunk.index)?;
             } else {
                 if config.use_gpu && gpu_context_param.is_none() {
                     // Print warning only once or use logging if available
                     // For now, skipping the print to avoid excessive output
                 }
                 builder.build_grammar(&chunk.data, chunk.index)?;
             }
             
             let (seq, rules) = builder.get_grammar();
             let depth = builder.get_max_rule_depth();
             Ok(Grammar { 
                 sequence: seq.clone(), 
                 rules: rules.clone(), 
                 max_depth: depth,
                 origins: HashMap::new(), // Initialize origins here
             })
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
/// This function is now public.
pub fn merge_grammars(grammars: Vec<Grammar>, _config: &ChunkingConfig, total_sequence_len: usize) -> Result<(Grammar, ParallelMetrics)> {
    let start_merge = Instant::now();
    let mut metrics = ParallelMetrics::default();

    if grammars.is_empty() {
        // Manually create default Grammar and return with metrics
        let default_grammar = Grammar {
            sequence: Vec::new(),
            rules: HashMap::new(),
            max_depth: 0,
            origins: HashMap::new(), // Initialize origins as empty
        };
        return Ok((default_grammar, metrics));
    }

    // Assign unique IDs to input grammars for tracking origins
    let grammars_with_ids: Vec<(usize, Grammar)> = grammars.into_iter().enumerate().collect();

    // Deduplicate rules using the union-find approach, passing grammars with IDs
    let deduplication_start = Instant::now();
    let (rule_map, mut merged_rules, rule_origins) = deduplicate_rules(&grammars_with_ids); // Pass grammars with IDs
    metrics.deduplication_time = deduplication_start.elapsed();
    // Calculate rules_before_merge here using input `grammars`
    metrics.rules_before_merge = grammars_with_ids.iter().map(|(_, g)| g.rules.len()).sum(); // Use grammars_with_ids
    metrics.rules_after_merge = merged_rules.len(); // Corrected access

    // --- Step 5: Reconstruct the final sequence by stitching chunks ---
    let mut final_sequence = Vec::with_capacity(total_sequence_len); // Estimate capacity
    let mut current_original_pos = 0; // Track position in the original, uncompressed sequence

    for (chunk_index, (grammar_id, grammar)) in grammars_with_ids.iter().enumerate() { // Iterate with grammar_id
        let _chunk_start_in_original = grammar.sequence.first().map_or(current_original_pos, |s| s.id); // Approximate start pos
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
                    // Use grammar_id and old_rule_id for lookup
                    if let Some(&new_rule_id) = rule_map.get(&(*grammar_id, old_rule_id)) {
                        Symbol::non_terminal(symbol.id, new_rule_id, symbol.strand)
                    } else {
                        eprintln!("Warning: Rule ID {} from grammar {} not found in map during sequence reconstruction.", old_rule_id, grammar_id);
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

    // Calculate Assembly Indices for the merged rules
    // It's important to operate on a mutable reference if the function modifies the rules in place.
    // Assuming merged_rules needs to be mutable for this call, or the function handles it.
    // If calculate_rule_assembly_indices doesn't modify in place but returns new map, adjust accordingly.
    if let Err(e) = calculate_rule_assembly_indices(&mut merged_rules) {
        eprintln!("Warning: Failed to calculate assembly indices during grammar merge: {}. Proceeding without them.", e);
        // Depending on policy, this could be a hard error: return Err(e.into());
    }

    // Calculate final max depth
    let max_depth = calculate_max_depth(&merged_rules);

    Ok((
        Grammar {
            sequence: final_sequence,
            rules: merged_rules,
            max_depth,
            origins: rule_origins,
        },
        metrics,
    ))
}

/// Deduplicate rules across multiple grammars using Union-Find
/// Takes Vec<(usize, Grammar)> where usize is the assigned grammar ID
/// Returns: (map_from_old_to_new_id, map_of_new_rules, map_from_new_id_to_origins)
fn deduplicate_rules(grammars_with_ids: &[(usize, Grammar)]) 
    -> (HashMap<(usize, usize), usize>, HashMap<usize, Rule>, HashMap<usize, Vec<(usize, usize)>>) 
{
    // --- Step 1: Collect all rules and map old IDs to a temporary global ID ---
    let mut old_to_temp_global = HashMap::new(); // Maps (grammar_id, old_rule_id) -> temp_global_id
    let mut temp_global_to_info: Vec<(u64, usize, usize, Vec<Symbol>)> = Vec::new(); // (hash, grammar_id, old_rule_id, symbols)
    let mut temp_global_counter = 0;

    for (grammar_id, grammar) in grammars_with_ids {
        for (old_rule_id, rule) in &grammar.rules {
            let rule_hash = canonical_hash_symbols(&rule.symbols);
            old_to_temp_global.insert((*grammar_id, *old_rule_id), temp_global_counter);
            temp_global_to_info.push((rule_hash, *grammar_id, *old_rule_id, rule.symbols.clone()));
            temp_global_counter += 1;
        }
    }

    let num_rules = temp_global_counter;
    if num_rules == 0 {
        return (HashMap::new(), HashMap::new(), HashMap::new()); 
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

    // --- Step 3 (Modified): Determine final representatives, origins, and total usage ---
    let mut representative_to_final_id: HashMap<usize, usize> = HashMap::new();
    let mut final_id_to_origins: HashMap<usize, Vec<(usize, usize)>> = HashMap::new(); // (grammar_id, old_rule_id)
    let mut final_id_to_representative_info: HashMap<usize, (Vec<Symbol>, usize)> = HashMap::new(); // (rep_symbols, total_usage)
    let mut next_final_id_counter = 0;

    // Helper map to quickly find grammar by ID
    let grammar_map: HashMap<usize, &Grammar> = grammars_with_ids.iter().map(|(id, g)| (*id, g)).collect();

    for temp_global_id in 0..num_rules {
        let representative_temp_id = uf.find(temp_global_id);
        let (_, grammar_id, old_rule_id, _) = temp_global_to_info[temp_global_id]; 

        // Safely get usage count from original grammar using grammar_id
        let usage_count = grammar_map.get(&grammar_id)
            .and_then(|g| g.rules.get(&old_rule_id))
            .map_or(0, |r| r.usage_count);

        let final_new_id = *representative_to_final_id.entry(representative_temp_id).or_insert_with(|| {
            let new_id = next_final_id_counter;
            next_final_id_counter += 1;
            // Get representative info 
            let (_, _, _, rep_symbols) = temp_global_to_info[representative_temp_id].clone(); 
            final_id_to_representative_info.insert(new_id, (rep_symbols, 0)); 
            new_id
        });

        // Add origin
        final_id_to_origins.entry(final_new_id).or_default().push((grammar_id, old_rule_id));
        // Add usage count
        if let Some((_, total_usage)) = final_id_to_representative_info.get_mut(&final_new_id) {
            *total_usage += usage_count;
        }
    }
    
    // Create final_rule_map (old -> new)
    let mut final_rule_map = HashMap::new(); // Maps (grammar_id, old_rule_id) -> final_new_id
    for (final_id, origins) in &final_id_to_origins {
        for (grammar_id, old_rule_id) in origins {
             final_rule_map.insert((*grammar_id, *old_rule_id), *final_id);
        }
    }

    // --- Step 4 (Modified): Create final rules with remapped symbols and total usage ---
    let mut final_merged_rules_map: HashMap<usize, Rule> = HashMap::new();
    for (final_id, (rep_symbols, total_usage_count)) in final_id_to_representative_info {
        let mut final_symbols = Vec::with_capacity(rep_symbols.len());

        // Remap non-terminals within the representative symbols
        for original_symbol in rep_symbols {
            if let SymbolType::NonTerminal(old_child_id) = original_symbol.symbol_type {
                // Find the grammar_id associated with the *representative* rule
                // This requires finding an origin for this final_id
                let (origin_grammar_id, _) = final_id_to_origins.get(&final_id)
                                                .and_then(|origins| origins.first()) // Take the first origin as the context for remapping children
                                                .cloned()
                                                .ok_or_else(|| anyhow::anyhow!("No origin found for final rule ID {}", final_id)).unwrap(); 
                                                
                // Lookup using the representative's original grammar_id and the old_child_id
                if let Some(final_new_child_id) = final_rule_map.get(&(origin_grammar_id, old_child_id)) {
                    final_symbols.push(Symbol {
                        id: original_symbol.id, 
                        symbol_type: SymbolType::NonTerminal(*final_new_child_id),
                        strand: original_symbol.strand,
                        source_grammar_id: None, // Not directly from an original sequence after merge
                        original_pos: None,      // Not directly from an original sequence after merge
                    });
                } else {
                    eprintln!("Warning: Could not remap non-terminal ID {} (from grammar {}) in rule {}", old_child_id, origin_grammar_id, final_id);
                    final_symbols.push(original_symbol); // Keep original
                }
            } else {
                final_symbols.push(original_symbol); // Terminal
            }
        }

        // Create the final rule
        final_merged_rules_map.insert(final_id, Rule {
            id: final_id,
            symbols: final_symbols,
            usage_count: total_usage_count, // Use summed usage count
            positions: Vec::new(),     // Positions are lost/invalid after merge
            depth: None,               // Depth needs recalculation
            assembly_index: None,      // Initialize assembly_index
        });
    }

    (final_rule_map, final_merged_rules_map, final_id_to_origins)
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
    fn test_parallel_sequitur() -> Result<()> {
        let bases = make_bases("ACGTACGTACGT");
        
        let config = ChunkingConfig {
            chunk_size: 5,
            overlap_size: 2,
            min_rule_usage: 2,
            reverse_aware: true,
            num_threads: 2,
            show_progress: false,
            adaptive_chunking: false,
            max_memory_per_chunk: None,
            use_gpu: false,
        };
        
        let (grammar, metrics) = parallel_sequitur(&bases, config, None)?;
        
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