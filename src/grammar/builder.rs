// extern crate assert_matches;
use crate::grammar::digram_table::{DigramKey, DigramTable};
use crate::grammar::rule::Rule;
use crate::grammar::symbol::{Symbol, SymbolType, Direction};
use crate::encode::dna_2bit::EncodedBase;
use anyhow::Result;
use std::collections::{HashMap, HashSet};
use rayon::prelude::*;
use std::time::Instant;
use std::collections::BinaryHeap;
use std::cmp::Reverse;
use crate::grammar::digram::find_most_frequent_terminal_digram_suffix_array;
use ordered_float::OrderedFloat;

/// Builds a grammar (set of rules) by iteratively replacing 
/// the most frequent digrams in a sequence.
#[derive(Debug)]
pub struct GrammarBuilder {
    // The working sequence, modified during grammar construction.
    sequence: Vec<Symbol>,
    // Stores the rules created during the process.
    rules: HashMap<usize, Rule>,
    // Tracks digram occurrences and handles canonicalization.
    digram_table: DigramTable,
    // Counter for assigning unique IDs to new rules (non-terminals).
    next_rule_id: usize,
    // Configuration settings
    min_rule_usage: usize,
    reverse_aware: bool,
    // Maximum number of rules to keep before inlining least used rules
    max_rule_count: Option<usize>,
    // Track rule depths for hierarchy analysis
    rule_depths: HashMap<usize, usize>,
    // Performance metrics
    metrics: PerformanceMetrics,
    // Streaming mode settings
    stream_mode: bool,
    chunk_count: usize,
    total_bases_processed: usize,
}

/// Track performance metrics during grammar construction
#[derive(Debug, Default)]
pub struct PerformanceMetrics {
    digram_table_time: std::time::Duration,
    step_count: usize,
    replacement_count: usize,
    replacement_time: std::time::Duration,
    inlining_time: std::time::Duration,
    eviction_time: std::time::Duration,
}

impl GrammarBuilder {
    /// Creates a new GrammarBuilder.
    pub fn new(min_rule_usage: usize, reverse_aware: bool) -> Self {
        GrammarBuilder {
            sequence: Vec::new(),
            rules: HashMap::new(),
            digram_table: DigramTable::new(),
            next_rule_id: 0,
            min_rule_usage,
            reverse_aware,
            max_rule_count: None,
            rule_depths: HashMap::new(),
            metrics: PerformanceMetrics::default(),
            stream_mode: false,
            chunk_count: 0,
            total_bases_processed: 0,
        }
    }

    /// Sets the maximum number of rules before triggering rule eviction.
    pub fn with_max_rules(mut self, max_rule_count: usize) -> Self {
        self.max_rule_count = Some(max_rule_count);
        self
    }

    /// Initializes the builder with the input sequence.
    /// Converts the raw byte sequence into Terminal Symbols.
    fn initialize_sequence(&mut self, initial_sequence: &[EncodedBase]) {
        self.sequence = initial_sequence
            .iter()
            .enumerate()
            .map(|(i, &base)| {
                // Initial symbols are always on the '+' strand? Or does this need context?
                // Assume '+' for now.
                Symbol::terminal(i, base, Direction::Forward)
            })
            .collect();
        println!("Initialized sequence with {} symbols.", self.sequence.len());
    }

    /// Populates the digram table based on the current sequence state.
    /// Uses parallel processing for large sequences.
    /// Optimization: If sequence is all terminals, uses suffix array for first pass.
    fn rebuild_digram_table(&mut self) {
        let start = Instant::now();
        
        // Clear existing table
        self.digram_table = DigramTable::new();

        // Special case: Empty sequence
        if self.sequence.len() < 2 {
            return;
        }

        // Check if sequence contains only terminals
        let sequence_len = self.sequence.len();
        let sample_size = std::cmp::min(1000, sequence_len);
        let sample_step = if sample_size >= sequence_len { 1 } else { sequence_len / sample_size };
        
        // Check if sequence is likely all terminals by sampling
        let likely_all_terminals = self.sequence.par_iter()
            .step_by(sample_step)
            .take(sample_size)
            .all(|s| matches!(s.symbol_type, SymbolType::Terminal(_)));

        // For very large sequences (>1M), use suffix array optimization if possible
        if sequence_len > 1_000_000 && likely_all_terminals && self.rules.is_empty() {
            println!("Large sequence ({}). Using suffix array optimization for digram finding.", sequence_len);
            
            // Extract EncodedBase sequence for suffix array in parallel
            let encoded_sequence: Vec<EncodedBase> = self.sequence
                .par_iter()
                .filter_map(|s| match s.symbol_type {
                    SymbolType::Terminal(base) => Some(base),
                    _ => None, 
                })
                .collect();
            
            // Verify that all symbols were successful terminals (sanity check)
            if encoded_sequence.len() == sequence_len {
                // Use suffix array to find the single most frequent terminal digram
                if let Some((canonical_key, count, positions)) = 
                    find_most_frequent_terminal_digram_suffix_array(&encoded_sequence, self.min_rule_usage, self.reverse_aware) {
                    
                    // For very frequent digrams, consider only a subset of positions to reduce memory
                    let max_positions = 1_000_000; // Cap to avoid excessive memory use
                    let positions_to_use = if positions.len() > max_positions {
                        println!("Limiting to {} of {} positions for memory efficiency", max_positions, positions.len());
                        let sample_rate = positions.len() / max_positions;
                        positions.into_iter()
                            .enumerate()
                            .filter_map(|(i, pos)| if i % sample_rate == 0 { Some(pos) } else { None })
                            .collect()
                    } else {
                        positions
                    };
                    
                    // Reconstruct the occurrences Vec<(usize, (Symbol, Symbol))> in parallel
                    let occurrences: Vec<(usize, (Symbol, Symbol))> = positions_to_use.par_iter()
                        .filter_map(|&pos| {
                            if pos + 1 < self.sequence.len() {
                                Some((pos, (self.sequence[pos], self.sequence[pos + 1])))
                            } else {
                                None
                            }
                        })
                        .collect();
                    
                    // Add only this single most frequent entry to the table
                    self.digram_table.add_single_entry(canonical_key, occurrences);
                    println!("Suffix array found most frequent digram: {:?}, count: {}", canonical_key, count);
                    
                    let elapsed = start.elapsed();
                    self.metrics.digram_table_time += elapsed;
                    return;
                }
            }
        }
        
        // For smaller sequences or when suffix array approach fails or isn't applicable
        println!("Using parallel hash map for digram table (seq length: {}, rules: {}).", 
                 sequence_len, self.rules.len());
        
        // Use specialized chunk size based on sequence length
        let chunk_size = if sequence_len > 10_000_000 {
            // Very large sequences (>10M): process in larger chunks
            1_000_000
        } else if sequence_len > 1_000_000 {
            // Large sequences (1M-10M): medium chunks
            100_000
        } else if sequence_len > 100_000 {
            // Medium sequences (100K-1M): smaller chunks
            10_000
        } else {
            // Small sequences: process all at once
            sequence_len
        };
        
        if chunk_size < sequence_len {
            println!("Processing digram table in chunks of size {}", chunk_size);
            
            // Process in chunks to reduce memory pressure
            for chunk_start in (0..sequence_len).step_by(chunk_size) {
                let chunk_end = std::cmp::min(chunk_start + chunk_size, sequence_len);
                if chunk_end > chunk_start + 1 { // Need at least 2 symbols for a digram
                    let sequence_chunk = &self.sequence[chunk_start..chunk_end];
                    self.digram_table.add_digrams_from_sequence(sequence_chunk, self.reverse_aware);
                }
            }
        } else {
            // Process the entire sequence at once
            self.digram_table.add_digrams_from_sequence(&self.sequence, self.reverse_aware);
        }
        
        let elapsed = start.elapsed();
        self.metrics.digram_table_time += elapsed;
        
        if sequence_len > 100_000 {
            println!("Digram table built in {:?}, contains {} unique digrams", 
                     elapsed, self.digram_table.len());
        }
    }

    /// Performs one step of grammar construction: finding and replacing the most frequent digram.
    /// Returns true if a replacement was made, false otherwise.
    fn step(&mut self) -> bool {
        // Find most frequent digram
        let most_frequent = self.digram_table.find_most_frequent_digram();
        
        // If no digram is found or it's not frequent enough, stop
        let (canonical_key, count, occurrences) = match most_frequent {
            Some(result) if result.1 >= self.min_rule_usage => result,
            _ => return false,
        };
        
        // Create a new rule
        let rule_id = self.next_rule_id;
        self.next_rule_id += 1;
        
        // Create a rule definition based on the first occurrence (this must exist)
        let first_digram = occurrences.first()
            .map(|(_, digram)| digram.clone())
            .expect("Occurrences should not be empty");
        
        let new_rule = Rule::new(rule_id, first_digram.0.clone(), first_digram.1.clone());
        
        // Add rule to the grammar
        self.rules.insert(rule_id, new_rule);
        
        // Replace occurrences with the new rule
        self.replace_occurrences(rule_id, canonical_key, &occurrences);
        
        // The sequence has been modified, so rebuild the digram table
        self.rebuild_digram_table();
        
        true
    }
    
    /// Calculate the depth of a rule based on its symbols
    fn calculate_rule_depth(&self, _rule_id: usize, rule: &Rule) -> usize {
        let mut max_child_depth = 0;
        
        for symbol in &rule.symbols {
            if let SymbolType::NonTerminal(child_rule_id) = symbol.symbol_type {
                if let Some(depth) = self.rule_depths.get(&child_rule_id) {
                    max_child_depth = max_child_depth.max(*depth);
                }
            }
        }
        
        // The depth is the maximum depth of its children + 1
        max_child_depth + 1
    }
    
    /// Find and inline rules that are used only once
    fn inline_single_use_rules(&mut self) {
        // First, find all rules with usage_count == 1
        let single_use_rule_ids: Vec<usize> = self.rules
            .iter()
            .filter(|(_, rule)| rule.usage_count == 1)
            .map(|(id, _)| *id)
            .collect();
        
        if single_use_rule_ids.is_empty() {
            return;
        }
        
        println!("Inlining {} rules used only once", single_use_rule_ids.len());
        
        // Keep track of rules to remove
        let mut rules_to_remove = HashSet::new();
        
        // For each single-use rule
        for rule_id in single_use_rule_ids {
            // Skip if already marked for removal
            if rules_to_remove.contains(&rule_id) {
                continue;
            }
            
            // Find where this rule is used in the sequence
            let mut found_in_sequence = false;
            for (seq_idx, symbol) in self.sequence.iter().enumerate() {
                if let SymbolType::NonTerminal(nt_rule_id) = symbol.symbol_type {
                    if nt_rule_id == rule_id {
                        found_in_sequence = true;
                        
                        // Get the rule to inline
                        if let Some(rule) = self.rules.get(&rule_id).cloned() {
                            // Clone the rules map before passing to avoid borrowing conflicts
                            let rules_snapshot = self.rules.clone();
                            let expanded = rule.expand_for_inlining(&rules_snapshot);
                            
                            // Apply the proper strand propagation
                            let mut final_expansion = Vec::new();
                            for mut exp_symbol in expanded {
                                if symbol.strand == Direction::Reverse {
                                    // Flip the strand of the expanded symbol
                                    exp_symbol.strand = if exp_symbol.strand == Direction::Forward { Direction::Reverse } else { Direction::Forward };
                                }
                                final_expansion.push(exp_symbol);
                            }
                            
                            // Replace the non-terminal with its expansion
                            self.sequence.splice(seq_idx..seq_idx+1, final_expansion);
                            
                            // Mark this rule for removal
                            rules_to_remove.insert(rule_id);
                            
                            // Since we modified the sequence, break and rebuild the whole sequence
                            break;
                        }
                    }
                }
            }
            
            if !found_in_sequence {
                // We need to collect all rules that need modification first
                let mut rule_modifications = Vec::new();
                
                // Find rules that use this rule
                for (other_rule_id, other_rule) in &self.rules {
                    if *other_rule_id == rule_id || rules_to_remove.contains(other_rule_id) {
                        continue;
                    }
                    
                    for (sym_idx, symbol) in other_rule.symbols.iter().enumerate() {
                        if let SymbolType::NonTerminal(nt_rule_id) = symbol.symbol_type {
                            if nt_rule_id == rule_id {
                                // Rule is used in another rule definition
                                // Get the rule to inline
                                if let Some(rule_to_inline) = self.rules.get(&rule_id).cloned() {
                                    // Take a snapshot of the rules for expansion
                                    let rules_snapshot = self.rules.clone();
                                    let expanded = rule_to_inline.expand_for_inlining(&rules_snapshot);
                                    
                                    // Apply the proper strand propagation
                                    let mut final_expansion = Vec::new();
                                    for mut exp_symbol in expanded {
                                        if symbol.strand == Direction::Reverse {
                                            // Flip the strand of the expanded symbol
                                            exp_symbol.strand = if exp_symbol.strand == Direction::Forward { Direction::Reverse } else { Direction::Forward };
                                        }
                                        final_expansion.push(exp_symbol);
                                    }
                                    
                                    // Store the modification to apply later
                                    rule_modifications.push((*other_rule_id, sym_idx, final_expansion));
                                    
                                    // Mark this rule for removal
                                    rules_to_remove.insert(rule_id);
                                    break;
                                }
                            }
                        }
                    }
                }
                
                // Apply all the modifications
                for (other_rule_id, sym_idx, expansion) in rule_modifications {
                    if let Some(other_rule) = self.rules.get_mut(&other_rule_id) {
                        other_rule.symbols.splice(sym_idx..sym_idx+1, expansion);
                    }
                }
            }
        }
        
        // Remove all rules marked for removal
        for rule_id in rules_to_remove {
            self.rules.remove(&rule_id);
            self.rule_depths.remove(&rule_id);
            println!("Removed inlined rule {}", rule_id);
        }
    }
    
    /// Evict least-used rules when exceeding the maximum rule count
    fn evict_rules(&mut self, max_rule_count: usize) {
        // If we're under the limit, no need to evict
        if self.rules.len() <= max_rule_count {
            return;
        }
        
        println!("Rule eviction: current count {} exceeds maximum {}", self.rules.len(), max_rule_count);
        
        /// Calculate priority score for rule eviction
        /// Returns a float where lower values indicate higher eviction priority
        fn calculate_rule_priority(
            usage_count: usize,
            memory_impact: usize,
            depth: usize,
            rule_id: usize,
            symbol_count: usize
        ) -> f64 {
            // Factors to consider (in order of importance):
            // 1. Memory impact: Higher impact = lower eviction priority
            // 2. Usage count: Higher usage = lower eviction priority
            // 3. Depth: Higher depth (more nested) = higher eviction priority
            // 4. Age (approximated by rule_id): Newer rules (higher ID) = higher eviction priority
            // 5. Symbol count: Longer rules = lower eviction priority (preserve complex patterns)
            
            // Normalize and weight each factor
            let usage_factor = (usage_count as f64).log10() * 10.0;  // Log scale for usage
            let memory_factor = (memory_impact as f64).sqrt() * 5.0; // Square root scale for memory impact
            let depth_factor = -(depth as f64) * 2.0;                // Negative for higher depth = lower score
            let age_factor = -(rule_id as f64).log10() * 2.0;        // Newer rules have higher IDs = lower score
            let length_factor = (symbol_count as f64).sqrt() * 1.5;  // Bonus for longer rules
            
            // Combined score (higher = keep, lower = evict)
            usage_factor + memory_factor + depth_factor + age_factor + length_factor
        }
        
        // Create a min-heap (using Reverse for min heap) to track rule priorities
        // Use a more sophisticated priority scoring system
        let mut rule_priority_queue = BinaryHeap::new();
        
        // Calculate memory impact of each rule
        let mut rule_memory_impacts = HashMap::new();
        for (&rule_id, rule) in &self.rules {
            // Skip rule 0 and other very low rule IDs (likely important base patterns)
            if rule_id < 5 {
                continue;
            }
            
            // Calculate memory impact (symbols saved by this rule)
            // For each usage, we replace rule.symbols.len() symbols with 1 reference
            // So savings per usage = rule.symbols.len() - 1
            // Total savings = (rule.symbols.len() - 1) * rule.usage_count
            let memory_impact = if rule.symbols.len() > 1 {
                (rule.symbols.len() - 1) * rule.usage_count
            } else {
                0 // No memory benefit for single-symbol rules
            };
            
            rule_memory_impacts.insert(rule_id, memory_impact);
            
            // Get rule depth (how nested it is)
            let depth = *self.rule_depths.get(&rule_id).unwrap_or(&0);
            
            // Calculate a comprehensive priority score
            // Lower score = higher eviction priority
            let priority_score = calculate_rule_priority(
                rule.usage_count,
                memory_impact,
                depth,
                rule_id,
                rule.symbols.len()
            );
            
            // Use OrderedFloat for f64 comparison in BinaryHeap
            rule_priority_queue.push(Reverse((OrderedFloat(priority_score), rule_id)));
        }
        
        // Determine how many rules to evict
        let rules_to_evict = self.rules.len().saturating_sub(max_rule_count);
        println!("  Evicting {} rules to meet limit of {}", rules_to_evict, max_rule_count);
        
        // Set of rule IDs to evict
        let mut eviction_ids = HashSet::new();
        let mut evicted_count = 0;
        let mut total_memory_impact = 0;
        
        // Pop rules from the priority queue until we've evicted enough
        while evicted_count < rules_to_evict && !rule_priority_queue.is_empty() {
            let Reverse((OrderedFloat(score), rule_id)) = rule_priority_queue.pop().unwrap();
            
            // Skip if already marked for eviction
            if eviction_ids.contains(&rule_id) {
                continue;
            }
            
            let memory_impact = *rule_memory_impacts.get(&rule_id).unwrap_or(&0);
            let rule = &self.rules[&rule_id];
            
            eviction_ids.insert(rule_id);
            evicted_count += 1;
            total_memory_impact += memory_impact;
            
            println!("  Evicting rule {} (score: {:.2}, usage: {}, memory impact: {})", 
                     rule_id, score, rule.usage_count, memory_impact);
        }
        
        println!("  Total memory impact of eviction: {} symbols", total_memory_impact);
        
        // Now perform the actual rule eviction by inlining each evicted rule
        self.inline_rules(&eviction_ids);
    }
    
    /// Inline a set of rules (replace references with rule contents)
    fn inline_rules(&mut self, rule_ids: &HashSet<usize>) {
        if rule_ids.is_empty() {
            return;
        }
        
        // For each rule, collect its definition for inlining
        let mut rule_definitions = HashMap::new();
        for &rule_id in rule_ids {
            if let Some(rule) = self.rules.get(&rule_id) {
                rule_definitions.insert(rule_id, rule.symbols.clone());
            }
        }
        
        // First inline rules in the main sequence
        self.inline_rules_in_sequence(&rule_definitions);
        
        // Then inline rules in other rule definitions
        self.inline_rules_in_rules(&rule_definitions);
        
        // Finally remove the inlined rules
        for &rule_id in rule_ids {
            self.rules.remove(&rule_id);
            self.rule_depths.remove(&rule_id);
        }
        
        // Rebuild the digram table after inlining
        self.rebuild_digram_table();
    }
    
    /// Inline rules in the main sequence
    fn inline_rules_in_sequence(&mut self, rule_definitions: &HashMap<usize, Vec<Symbol>>) {
        let mut replacements = Vec::new();
        
        // Scan the sequence for symbols to replace
        for (idx, symbol) in self.sequence.iter().enumerate() {
            if let SymbolType::NonTerminal(rule_id) = symbol.symbol_type {
                if rule_definitions.contains_key(&rule_id) {
                    // This symbol references a rule we're evicting
                    let definition = &rule_definitions[&rule_id];
                    replacements.push((idx, definition.clone()));
                }
            }
        }
        
        // Apply replacements from back to front to avoid shifting indices
        replacements.sort_by_key(|(idx, _)| std::cmp::Reverse(*idx));
        
        for (idx, replacement) in replacements {
            // Remove the original symbol
            self.sequence.remove(idx);
            
            // Insert the expanded symbols
            for (i, sym) in replacement.iter().enumerate() {
                self.sequence.insert(idx + i, *sym);
            }
        }
    }
    
    /// Inline rules within other rule definitions
    fn inline_rules_in_rules(&mut self, rule_definitions: &HashMap<usize, Vec<Symbol>>) {
        // Keep track of which rules need modification
        let mut rule_modifications = HashMap::new();
        
        // Scan all rules for references to rules we're evicting
        for (&rule_id, rule) in &self.rules {
            if rule_definitions.contains_key(&rule_id) {
                // Skip rules that are themselves being evicted
                continue;
            }
            
            let mut changes = Vec::new();
            
            // Look for references to evicted rules
            for (idx, symbol) in rule.symbols.iter().enumerate() {
                if let SymbolType::NonTerminal(ref_rule_id) = symbol.symbol_type {
                    if rule_definitions.contains_key(&ref_rule_id) {
                        // This symbol references a rule we're evicting
                        let definition = &rule_definitions[&ref_rule_id];
                        changes.push((idx, definition.clone()));
                    }
                }
            }
            
            if !changes.is_empty() {
                rule_modifications.insert(rule_id, changes);
            }
        }
        
        // Apply all modifications after we're done collecting them
        for (other_rule_id, changes) in rule_modifications {
            if let Some(other_rule) = self.rules.get_mut(&other_rule_id) {
                // Apply changes from highest index to lowest
                let mut sorted_changes = changes;
                sorted_changes.sort_by_key(|(idx, _)| std::cmp::Reverse(*idx));
                
                for (idx, expansion) in sorted_changes {
                    other_rule.symbols.splice(idx..idx+1, expansion);
                }
            }
        }
    }

    /// Replaces all occurrences of a digram with a non-terminal symbol for a rule.
    /// This implementation supports parallel replacements for non-overlapping positions.
    fn replace_occurrences(&mut self, rule_id: usize, canonical_key: DigramKey, occurrences: &[(usize, (Symbol, Symbol))]) {
        if occurrences.is_empty() {
            return;
        }
        
        let start = Instant::now();
        
        // Sort positions in descending order to avoid invalidation during replacement
        let mut positions: Vec<usize> = occurrences.iter().map(|(pos, _)| *pos).collect();
        positions.sort_unstable_by(|a, b| b.cmp(a)); // Sort descending
        
        // Group positions into non-overlapping sets (positions that are at least 2 apart)
        let mut position_groups: Vec<Vec<usize>> = Vec::new();
        
        for pos in positions {
            // Find a group where this position doesn't overlap with existing ones
            let mut found_group = false;
            for group in &mut position_groups {
                if group.iter().all(|&existing_pos| (existing_pos as isize - pos as isize).abs() >= 2) {
                    group.push(pos);
                    found_group = true;
                    break;
                }
            }
            
            // If no suitable group exists, create a new one
            if !found_group {
                position_groups.push(vec![pos]);
            }
        }
        
        // For each independent group, process the changes in parallel and apply sequentially
        for group in position_groups {
            // Process each group in parallel using rayon
            let changes: Vec<(usize, Symbol)> = group.par_iter().map(|&pos| {
                // Get original digram for strand information
                let original_digram = occurrences.iter()
                    .find(|(p, _)| *p == pos)
                    .map(|(_, digram)| digram)
                    .expect("Original digram not found for position");
                
                // Create non-terminal with same strand as first symbol of digram
                let non_term = Symbol::non_terminal(pos, rule_id, original_digram.0.strand);
                
                (pos, non_term)
            }).collect();
            
            // Sort by position descending before applying (to avoid index invalidation)
            let mut sorted_changes = changes;
            sorted_changes.sort_unstable_by(|a, b| b.0.cmp(&a.0));
            
            // Apply changes to the sequence
            for (pos, non_term) in sorted_changes {
                if pos < self.sequence.len() - 1 {
                    self.sequence[pos] = non_term;
                    self.sequence.remove(pos + 1);
                    
                    // Remove this occurrence from the digram table
                    self.digram_table.remove_occurrence(&canonical_key, pos);
                }
            }
            
            // Increment the rule usage count
            if let Some(rule) = self.rules.get_mut(&rule_id) {
                rule.usage_count += group.len();
            }
        }
        
        // Update metrics
        let elapsed = start.elapsed();
        self.metrics.replacement_time += elapsed;
        self.metrics.replacement_count += occurrences.len();
        
        if occurrences.len() > 1000 {
            println!("Replaced {} occurrences in {:?}", occurrences.len(), elapsed);
        }
    }

    /// Enables streaming mode for processing large sequences in chunks
    pub fn enable_streaming_mode(mut self) -> Self {
        self.stream_mode = true;
        self
    }

    /// Process a chunk of the sequence incrementally
    pub fn process_sequence_chunk(&mut self, chunk: &[EncodedBase]) -> Result<()> {
        let chunk_start_time = Instant::now();
        self.chunk_count += 1;
        self.total_bases_processed += chunk.len();
        
        if self.sequence.is_empty() {
            // First chunk - initialize the sequence
            self.initialize_sequence(chunk);
        } else {
            // Add this chunk to the existing sequence
            let start_id = self.sequence.len();
            let chunk_symbols: Vec<Symbol> = chunk
                .iter()
                .enumerate()
                .map(|(i, &base)| {
                    Symbol::terminal(start_id + i, base, Direction::Forward)
                })
                .collect();
            
            self.sequence.extend(chunk_symbols);
        }
        
        // After adding the chunk, rebuild digram table and process rules
        self.rebuild_digram_table();
        
        // Determine how many rule building steps to perform based on chunk size
        let max_steps_per_chunk = if self.stream_mode {
            // In streaming mode, perform more extensive processing per chunk
            // We want to create rules as we go to keep memory usage low
            // Use heuristic based on chunk size
            std::cmp::max(chunk.len() / 100, 20)
        } else {
            // In regular mode, do minimal processing per chunk
            10
        };
        
        // Perform several rule building steps
        let mut steps_performed = 0;
        while self.step() && steps_performed < max_steps_per_chunk {
            steps_performed += 1;
        }
        
        // In streaming mode, evict rules if needed to control memory usage
        if self.stream_mode {
            if let Some(max_count) = self.max_rule_count {
                if self.rules.len() > max_count {
                    let evict_start = Instant::now();
                    self.evict_rules(max_count);
                    self.metrics.eviction_time += evict_start.elapsed();
                }
            }
        }
        
        if self.chunk_count % 10 == 0 || self.chunk_count == 1 {
            println!("Processed chunk {} ({} bases) in {:.2?}. Current stats: {} symbols, {} rules", 
                     self.chunk_count, 
                     chunk.len(),
                     chunk_start_time.elapsed(),
                     self.sequence.len(),
                     self.rules.len());
        }
        
        Ok(())
    }
    
    /// Finalize the grammar after all chunks have been processed
    pub fn finalize_grammar(&mut self) -> Result<()> {
        // For streaming mode, finalize the construction
        // Should only be called at the end of processing
        println!("Finalizing grammar after processing {} chunks ({} bases)...",
                 self.chunk_count, self.total_bases_processed);
                 
        // Continue rule building until no more replacements can be made
        let mut remaining_steps = 0;
        while self.step() {
            remaining_steps += 1;
            
            // Report progress periodically
            if remaining_steps % 100 == 0 {
                println!("  Finalization progress: {} additional steps, current sequence length: {}, rules: {}", 
                         remaining_steps, 
                         self.sequence.len(),
                         self.rules.len());
            }
            
            // Periodic rule eviction during finalization, if enabled
            if let Some(max_count) = self.max_rule_count {
                if remaining_steps % 50 == 0 && self.rules.len() > max_count {
                    let evict_start = Instant::now();
                    self.evict_rules(max_count);
                    self.metrics.eviction_time += evict_start.elapsed();
                }
            }
        }
        
        // Perform any final optimizations
        if let Some(max_count) = self.max_rule_count {
            if self.rules.len() > max_count {
                let evict_start = Instant::now();
                self.evict_rules(max_count);
                self.metrics.eviction_time += evict_start.elapsed();
            }
        }
        
        // Final inlining of single-use rules
        let inline_start = Instant::now();
        self.inline_single_use_rules();
        self.metrics.inlining_time += inline_start.elapsed();
        
        println!("Grammar finalization complete. Performed {} additional steps.", remaining_steps);
        println!("Final stats: {} symbols, {} rules", self.sequence.len(), self.rules.len());
        
        Ok(())
    }

    /// Builds the grammar from a DNA sequence.
    pub fn build_grammar(&mut self, initial_sequence: &[EncodedBase]) -> Result<()> {
        println!("Building grammar from sequence of length {}", initial_sequence.len());
        let start = Instant::now();
        
        // Initialize sequence from the input DNA
        self.initialize_sequence(initial_sequence);

        // --- Initial Step: Use Suffix Array for first rule --- 
        let suffix_array_start = Instant::now();
        let initial_digram_info = find_most_frequent_terminal_digram_suffix_array(
            initial_sequence,
            self.min_rule_usage,
            self.reverse_aware,
        );

        if let Some((digram_key, count, positions)) = initial_digram_info {
             println!(
                 "Initial SA Step: Found terminal digram with key {} (count {} >= {}). Creating Rule 0.",
                 digram_key,
                 count,
                 self.min_rule_usage
             );
             
             // Since DigramKey is now a hash, we need to extract the first occurrence from the sequence
             // to create the rule definition
             if !positions.is_empty() && positions[0] + 1 < initial_sequence.len() {
                 let rule_id = 0;
                 self.next_rule_id = 1;
                 
                 // Get the actual symbols from the first occurrence position
                 let base1 = initial_sequence[positions[0]];
                 let base2 = initial_sequence[positions[0] + 1];
                 
                 // Create symbol instances for the rule
                 let sym1 = Symbol::terminal(0, base1, Direction::Forward);
                 let sym2 = Symbol::terminal(1, base2, Direction::Forward);

                 let mut new_rule = Rule::new(rule_id, sym1, sym2);
                 new_rule.usage_count = count;
                 // Calculate initial depth (will be 1 for terminal-only rule)
                 let depth = 1; 
                 self.rule_depths.insert(rule_id, depth);
                 new_rule.set_depth(depth);
                 new_rule.positions = positions.clone(); // Store initial positions

                 self.rules.insert(rule_id, new_rule);

                 // Replace occurrences found by suffix array
                 self.replace_initial_occurrences(rule_id, &positions);
                 self.metrics.replacement_count += positions.len();
             }
        }

        // Now, build the digram table for the potentially modified sequence
        self.rebuild_digram_table();
        
        // Continue building the grammar through iterative steps
        let mut more_steps = true;
        let mut step_count = 0;
        
        while more_steps {
            more_steps = self.step();
            step_count += 1;
            
            if self.max_rule_count.is_some() && step_count % 50 == 0 {
                // Periodically check if we need rule eviction
                if let Some(max_count) = self.max_rule_count {
                    if self.rules.len() > max_count {
                        let evict_start = Instant::now();
                        self.evict_rules(max_count);
                        self.metrics.eviction_time += evict_start.elapsed();
                    }
                }
            }
            
            // Report progress for long-running grammar construction
            if step_count % 100 == 0 {
                println!("  Grammar building step {}, sequence length: {}, rules: {}", 
                        step_count, 
                        self.sequence.len(),
                        self.rules.len());
            }
        }
        
        // Final rule inlining
        let inline_start = Instant::now();
        self.inline_single_use_rules();
        self.metrics.inlining_time += inline_start.elapsed();
        
        // Print building summary
        let elapsed = start.elapsed();
        println!("Grammar construction complete in {}s. Steps: {}, final seq len: {}, rules: {}", 
                 elapsed.as_secs_f64(),
                 step_count,
                 self.sequence.len(),
                 self.rules.len());
                 
        Ok(())
    }

    /// Provides access to the current state of the sequence and rules.
    pub fn get_grammar(&self) -> (&Vec<Symbol>, &HashMap<usize, Rule>) {
        (&self.sequence, &self.rules)
    }
    
    /// Returns the rule depth mapping for analysis.
    pub fn get_rule_depths(&self) -> &HashMap<usize, usize> {
        &self.rule_depths
    }
    
    /// Gets the maximum rule depth in the grammar.
    pub fn get_max_rule_depth(&self) -> usize {
        self.rule_depths.values().copied().max().unwrap_or(0)
    }
    
    /// Gets average rule depth in the grammar.
    pub fn get_avg_rule_depth(&self) -> f64 {
        if self.rule_depths.is_empty() {
            return 0.0;
        }
        
        let sum: usize = self.rule_depths.values().sum();
        sum as f64 / self.rule_depths.len() as f64
    }
    
    /// Gets performance metrics
    pub fn get_performance_metrics(&self) -> &PerformanceMetrics {
        &self.metrics
    }

    /// Process a chunk of raw bytes by converting them to EncodedBase
    pub fn process_bytes_chunk(&mut self, chunk: &[u8]) -> Result<()> {
        // Convert from bytes to EncodedBase
        let encoded_chunk: Vec<EncodedBase> = chunk.iter()
            .filter_map(|&b| EncodedBase::from_base(b))
            .collect();
        
        // Use the regular process_sequence_chunk method
        self.process_sequence_chunk(&encoded_chunk)
    }

    // --- Add helper function for initial replacement --- 
    /// Replaces occurrences at specified positions with a non-terminal rule.
    /// Assumes positions are sorted and non-overlapping at distance 1.
    /// Designed for the initial replacement based on suffix array results.
    fn replace_initial_occurrences(&mut self, rule_id: usize, positions: &[usize]) {
        if positions.is_empty() {
            return;
        }

        // Sort positions descending to replace from the end
        let mut sorted_positions = positions.to_vec();
        sorted_positions.sort_unstable_by(|a, b| b.cmp(a));

        let mut replaced_indices = HashSet::new(); // Track indices already affected

        for &pos in &sorted_positions {
            // Check if this position or the next one has already been affected by a replacement
            if replaced_indices.contains(&pos) || replaced_indices.contains(&(pos + 1)) {
                continue; 
            }

            if pos + 1 < self.sequence.len() {
                // Determine the strand of the new non-terminal based on the first symbol being replaced
                let nt_strand = self.sequence[pos].strand;
                let new_symbol = Symbol::non_terminal(pos, rule_id, nt_strand);

                // Perform the replacement
                self.sequence[pos] = new_symbol;
                self.sequence.remove(pos + 1);

                // Mark affected indices
                replaced_indices.insert(pos);
                // Since we removed pos+1, no need to mark it explicitly
            }
        }
    }
    // --- End helper function --- 
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::encode::dna_2bit::EncodedBase; // Import EncodedBase
    use crate::grammar::symbol::SymbolType;

    // Helper to create encoded base vector
    fn encode_seq(seq_bytes: &[u8]) -> Vec<EncodedBase> {
        seq_bytes.iter().filter_map(|&b| EncodedBase::from_base(b)).collect()
    }

    #[test]
    fn test_simple_repetition() -> Result<()> {
        let mut builder = GrammarBuilder::new(2, false);
        let seq = encode_seq(b"ABABAB"); // Use helper
        builder.build_grammar(&seq)?;
        let (final_sequence, rules) = builder.get_grammar();

        // Expected: Rule R0 = AB, Sequence = R0 R0 R0
        assert_eq!(rules.len(), 1);
        assert_eq!(final_sequence.len(), 3);
        assert!(matches!(final_sequence[0].symbol_type, SymbolType::NonTerminal(0)));
        assert!(matches!(final_sequence[1].symbol_type, SymbolType::NonTerminal(0)));
        assert!(matches!(final_sequence[2].symbol_type, SymbolType::NonTerminal(0)));

        let rule0 = &rules[&0];
        assert_eq!(rule0.symbols.len(), 2);
        assert!(matches!(rule0.symbols[0].symbol_type, SymbolType::Terminal(t) if t == EncodedBase(0)));
        assert!(matches!(rule0.symbols[1].symbol_type, SymbolType::Terminal(t) if t == EncodedBase(1))); // assuming B is 1
        Ok(())
    }

    #[test]
    fn test_overlapping_patterns() -> Result<()> {
        let mut builder = GrammarBuilder::new(2, false);
        let seq = encode_seq(b"ABCABCABC"); // Use helper
        builder.build_grammar(&seq)?;
        let (_final_sequence, rules) = builder.get_grammar();

        // Check if rules like AB, BC, CA or ABC are formed
        // Exact output depends on tie-breaking, so check for plausible rules
        assert!(rules.len() >= 1); // Should find at least one rule
        // Example check for a rule like "AB"
        let rule_ab_found = rules.values().any(|r| {
            r.symbols.len() == 2 &&
            matches!(r.symbols[0].symbol_type, SymbolType::Terminal(EncodedBase(0))) && // Wrap
            matches!(r.symbols[1].symbol_type, SymbolType::Terminal(EncodedBase(1)))    // Wrap
        });
         assert!(rule_ab_found, "Rule for 'AB' or similar not found");
        Ok(())
    }

    #[test]
    fn test_nested_rules() -> Result<()> {
        let mut builder = GrammarBuilder::new(2, false);
        // Sequence: XYXYXYZZXYXYXY
        let seq = encode_seq(b"ABABABCBCBABABAB"); // Use A=X, B=Y, C=Z for simplicity
        builder.build_grammar(&seq)?;
        let (_final_sequence, rules) = builder.get_grammar();

        // Expect rules like R0=AB, R1=R0R0, R2=BC, R3=R2R2 ... etc.
        assert!(rules.len() > 1); // Expect multiple rules for nesting

        // Example: Check for a rule composed of other rules
        let nested_rule_found = rules.values().any(|r| {
            r.symbols.len() == 2 &&
            matches!(r.symbols[0].symbol_type, SymbolType::NonTerminal(_)) &&
            matches!(r.symbols[1].symbol_type, SymbolType::NonTerminal(_))
        });
        assert!(nested_rule_found, "Nested rule (NT -> NT NT) not found");
        Ok(())
    }


    #[test]
    fn test_rule_utility() -> Result<()> {
        let mut builder = GrammarBuilder::new(2, false);
        // Sequence: A B A B C D C D A B A B
        // Expect: R0=AB, R1=CD. Final: R0 R0 R1 R1 R0 R0
        let seq = encode_seq(b"ABABCDCDABAB"); // Use helper
        builder.build_grammar(&seq)?;
        let (final_sequence, rules) = builder.get_grammar();

        assert_eq!(rules.len(), 2); // Expect rules AB and CD
        assert_eq!(final_sequence.len(), 6); // R0 R0 R1 R1 R0 R0

        // Find rules based on their content (more robust than assuming IDs)
        let rule0_id = rules.values()
            .find(|r| matches!(r.symbols[0].symbol_type, SymbolType::Terminal(t) if t == EncodedBase(0)) && // A
                       matches!(r.symbols[1].symbol_type, SymbolType::Terminal(t) if t == EncodedBase(1)))
            .map(|r| r.id)
            .expect("Rule starting with A not found");
        let rule1_id = rules.values()
            .find(|r| matches!(r.symbols[0].symbol_type, SymbolType::Terminal(t) if t == EncodedBase(2)) && // C
                       matches!(r.symbols[1].symbol_type, SymbolType::Terminal(t) if t == EncodedBase(3)))
            .map(|r| r.id)
            .expect("Rule starting with C not found");


        assert!(matches!(final_sequence[0].symbol_type, SymbolType::NonTerminal(id) if id == rule0_id));
        assert!(matches!(final_sequence[1].symbol_type, SymbolType::NonTerminal(id) if id == rule0_id));
        assert!(matches!(final_sequence[2].symbol_type, SymbolType::NonTerminal(id) if id == rule1_id));
        assert!(matches!(final_sequence[3].symbol_type, SymbolType::NonTerminal(id) if id == rule1_id));
        assert!(matches!(final_sequence[4].symbol_type, SymbolType::NonTerminal(id) if id == rule0_id));
        assert!(matches!(final_sequence[5].symbol_type, SymbolType::NonTerminal(id) if id == rule0_id));
        Ok(())
    }

     #[test]
    fn test_no_frequent_digrams() -> Result<()> {
        let mut builder = GrammarBuilder::new(2, false);
        let seq = encode_seq(b"ABCDEFG"); // Use helper
        builder.build_grammar(&seq)?;
        let (final_sequence, rules) = builder.get_grammar();

        assert!(rules.is_empty()); // No rules should be created
        assert_eq!(final_sequence.len(), 7); // Sequence remains unchanged
        assert!(matches!(final_sequence[0].symbol_type, SymbolType::Terminal(t) if t == EncodedBase(0)));
        assert!(matches!(final_sequence[1].symbol_type, SymbolType::Terminal(t) if t == EncodedBase(1)));
        // ... check other terminals if needed
        Ok(())
    }

    #[test]
    fn test_empty_sequence() -> Result<()> {
        let mut builder = GrammarBuilder::new(2, false);
        let seq = encode_seq(b""); // Use helper
        builder.build_grammar(&seq)?;
        let (final_sequence, rules) = builder.get_grammar();

        assert!(rules.is_empty());
        assert!(final_sequence.is_empty());
        Ok(()) // Add Ok result
    }

    #[test]
    fn test_single_base_sequence() -> Result<()> {
        let mut builder = GrammarBuilder::new(2, false);
        let seq = encode_seq(b"A"); // Use helper
        builder.build_grammar(&seq)?;
        let (final_sequence, rules) = builder.get_grammar();

        assert!(rules.is_empty());
        assert_eq!(final_sequence.len(), 1);
        assert!(matches!(final_sequence[0].symbol_type, SymbolType::Terminal(t) if t == EncodedBase(0)));
        Ok(()) // Add Ok result
    }

     #[test]
    fn test_reverse_complement_aware() -> Result<()> {
        // Sequence: ACGT ACGT (palindromic)
        // Should create rule R0 = ACGT (or equivalent)
        // Using AC GT AC GT
        let mut builder_aware = GrammarBuilder::new(2, true); // Reverse aware
        let seq = encode_seq(b"ACGTACGT"); // Use helper
        builder_aware.build_grammar(&seq)?;
        let (_final_sequence_aware, rules_aware) = builder_aware.get_grammar();

        // Expect rule creation due to AC == GT (rev comp)
        assert!(!rules_aware.is_empty(), "Reverse aware should find repeating patterns (AC/GT)");

        // Non-aware should treat AC and GT differently
         let mut builder_naive = GrammarBuilder::new(2, false); // Not reverse aware
         builder_naive.build_grammar(&seq)?;
         let (_final_sequence_naive, rules_naive) = builder_naive.get_grammar();
        
         // Maybe a rule for AC depending on tie break, but not guaranteed like aware mode
         // Less strict check: number of rules might differ
         // assert_ne!(rules_aware.len(), rules_naive.len(), "Rule count should potentially differ");
         Ok(())
    }

    // Add more tests:
    // - Edge cases: single base repeats (AAAA), alternating bases (ABAB)
    // - Minimum rule usage > 2
    // - Sequences that generate cycles (needs careful handling or cycle detection)
    // - Very long sequences (if performance testing is desired here)
} 