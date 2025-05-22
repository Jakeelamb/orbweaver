// extern crate assert_matches;
use crate::grammar::digram_table::{DigramTable, DigramKeyTuple, DigramSource};
use crate::grammar::kmer_table::KmerTable;
use crate::grammar::rule::Rule;
use crate::grammar::symbol::{Symbol, SymbolType, Direction};
use crate::encode::dna_2bit::EncodedBase;
use anyhow::Result;
use std::collections::{HashMap, HashSet};
use std::time::Instant;
use crate::gpu::GpuContext;
use crate::gpu::digram::GpuSequence;
use crate::grammar::engine::Grammar;
use crate::analysis::assembly_index::calculate_rule_assembly_indices; // Added for Assembly Index
use log;

/// Builds a grammar (set of rules) by iteratively replacing 
/// the most frequent patterns (digrams or k-mers) in a sequence.
#[derive(Debug)]
pub struct GrammarBuilder {
    // The working sequence, modified during grammar construction.
    sequence: Vec<Symbol>,
    // Stores the rules created during the process.
    rules: HashMap<usize, Rule>,
    // Tracks digram occurrences and handles canonicalization.
    digram_table: DigramTable,
    // Tracks k-mer occurrences (used when k > 2)
    kmer_table: Option<KmerTable>,
    // Counter for assigning unique IDs to new rules (non-terminals).
    next_rule_id: usize,
    // Configuration settings
    min_rule_usage: usize,
    reverse_aware: bool,
    // Maximum number of rules to keep before inlining least used rules
    max_rule_count: Option<usize>,
    // Size of k-mers to use (k=2 for digrams)
    kmer_size: usize,
    // Track rule depths for hierarchy analysis
    rule_depths: HashMap<usize, usize>,
    // Performance metrics
    metrics: PerformanceMetrics,
    // Streaming mode settings
    stream_mode: bool,
    chunk_count: usize,
    total_bases_processed: usize,
    // Store the GPU context if provided
    gpu_context: Option<GpuContext>,
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
            kmer_table: None,
            next_rule_id: 0,
            min_rule_usage,
            reverse_aware,
            max_rule_count: None,
            kmer_size: 2, // Default to digrams
            rule_depths: HashMap::new(),
            metrics: PerformanceMetrics::default(),
            stream_mode: false,
            chunk_count: 0,
            total_bases_processed: 0,
            // Initialize gpu_context to None
            gpu_context: None,
        }
    }

    /// Sets the maximum number of rules before triggering rule eviction.
    pub fn with_max_rules(mut self, max_rule_count: usize) -> Self {
        self.max_rule_count = Some(max_rule_count);
        self
    }
    
    /// Sets the k-mer size for grammar construction (k ≥ 2).
    pub fn with_kmer_size(mut self, kmer_size: usize) -> Self {
        if kmer_size < 2 {
            println!("Warning: k-mer size must be at least 2. Using k=2 instead.");
            self.kmer_size = 2;
        } else {
            self.kmer_size = kmer_size;
            // Initialize the k-mer table if k > 2
            if kmer_size > 2 {
                self.kmer_table = Some(KmerTable::new(kmer_size));
            }
        }
        self
    }
    
    /// Set whether to use GPU acceleration for grammar construction
    /// Stores the provided GpuContext.
    pub fn with_gpu(mut self, gpu_context: Option<&GpuContext>) -> Self {
        // Clone the context if provided, otherwise set to None
        self.gpu_context = gpu_context.cloned();
        self
    }

    /// Initializes the builder with the input sequence.
    /// Converts the raw byte sequence into Terminal Symbols.
    fn initialize_sequence(&mut self, initial_sequence: &[EncodedBase], source_grammar_id: usize) {
        self.sequence = initial_sequence
            .iter()
            .enumerate()
            .map(|(i, &base)| {
                // Initial symbols are always on the '+' strand? Or does this need context?
                // Assume '+' for now.
                // Use 'i' (original position) as the symbol's ID for terminals.
                Symbol::terminal(i, base, Direction::Forward, Some(source_grammar_id), Some(i))
            })
            .collect();
        println!("Initialized sequence with {} symbols.", self.sequence.len());
        println!("DEBUG: First 20 base values after initialization:");
        for (_i, sym) in self.sequence.iter().take(20).enumerate() {
            if let SymbolType::Terminal(base) = sym.symbol_type {
                print!("{} ", base.0);
            }
        }
        println!("");
        use std::io::Write;
        std::io::stdout().flush().unwrap();
    }

    /// Populates the pattern table (digram or k-mer) based on the current sequence state.
    /// Uses parallel processing for large sequences.
    fn rebuild_pattern_table(&mut self) -> Result<()> {
        let start = Instant::now();
        
        if self.kmer_size == 2 {
            // Traditional digram processing
            self.digram_table = DigramTable::new();
            self.digram_table.populate(&self.sequence, self.reverse_aware);
        } else {
            // K-mer processing
            if self.kmer_table.is_none() {
                self.kmer_table = Some(KmerTable::new(self.kmer_size));
            }
            
            let kmer_table = self.kmer_table.as_mut().unwrap();
            kmer_table.clear();
            kmer_table.populate(&self.sequence, self.reverse_aware)?;
        }
        
        self.metrics.digram_table_time += start.elapsed();
        Ok(())
    }

    /// Performs a single step of the grammar inference algorithm.
    /// Returns true if a rule was created, false if no more patterns are frequent enough.
    fn step(&mut self) -> Result<bool> {
        self.metrics.step_count += 1;
        if self.metrics.step_count % 100 == 0 {
             println!("  [CPU DEBUG] Step: {}, SeqLen: {}, Rules: {}", 
                 self.metrics.step_count, self.sequence.len(), self.rules.len());
        }

        self.rebuild_pattern_table()?;

        let pattern_found = if self.kmer_size == 2 {
            // --- Digram Processing --- //
            println!("DEBUG: DigramTable contains {} unique digrams", self.digram_table.len());
            if !self.digram_table.is_empty() {
                let top_digrams = self.digram_table.get_top_digrams(5);
                println!("DEBUG: Top digrams:");
                for (i, (key, occurrences)) in top_digrams.iter().enumerate() {
                    println!("  {}. Key: {:?}, Count: {}", i + 1, key, occurrences.len());
                    if i == 0 && !occurrences.is_empty() {
                        println!("     First occurrence position: {}", occurrences[0].0);
                        println!("     Min usage required: {}", self.min_rule_usage);
                    }
                }
            }

            if let Some((key_tuple, count, occurrences)) = self.digram_table.find_most_frequent_digram() {
                println!("DEBUG: Found digram with key {:?}, count: {}", key_tuple, count);
                if count >= self.min_rule_usage {
                    println!("DEBUG: Creating rule for digram with count {}", count);

                    // Reconstruct the first digram instance to create the rule
                    // This requires getting the symbols from the sequence using the first position
                    let first_pos = occurrences[0].0;
                    if first_pos + 1 < self.sequence.len() {
                        let s1 = self.sequence[first_pos].clone();
                        let s2 = self.sequence[first_pos + 1].clone();
                        let digram_symbols = vec![s1, s2];

                        let new_rule = Rule {
                            id: self.next_rule_id,
                            symbols: digram_symbols,
                            usage_count: count,
                            positions: occurrences.iter().map(|(pos, _source)| *pos).collect(),
                            depth: None,
                            assembly_index: None, // Initialize assembly_index
                        };

                        self.rules.insert(self.next_rule_id, new_rule);

                        let start_replacement = Instant::now();
                        
                        let digram_occurrences_for_replacement: Vec<(usize, (Symbol, Symbol))> = occurrences.iter()
                            .filter_map(|(pos, _source)| {
                                if *pos + 1 < self.sequence.len() {
                                    Some((*pos, (self.sequence[*pos].clone(), self.sequence[*pos + 1].clone())))
                                } else {
                                    None
                                }
                             })
                            .collect();
                        self.replace_digram_occurrences(self.next_rule_id, &digram_occurrences_for_replacement)?;
                        self.metrics.replacement_time += start_replacement.elapsed();
                        self.metrics.replacement_count += 1;

                        self.next_rule_id += 1;
                        true // Pattern processed
                    } else {
                         println!("Warning: Most frequent digram at position {} is too close to sequence end.", first_pos);
                         false // Could not process this pattern
                    }
                } else {
                    println!("DEBUG: Digram count {} is below minimum usage {}", count, self.min_rule_usage);
                    false // No frequent pattern found
                }
            } else {
                println!("DEBUG: No frequent digrams found");
                false // No pattern found
            }

        } else {
            // --- K-mer Processing (k > 2) --- //
            let kmer_table = match &self.kmer_table {
                 Some(table) => table,
                 None => {
                    println!("Warning: K-mer size > 2 but KmerTable is not initialized.");
                    return Ok(false);
                 }
            };

            if let Some((kmer_key, count, occurrences)) = kmer_table.find_most_frequent(self.min_rule_usage) {
                println!("DEBUG: Found k-mer with key {:?}, count: {}", kmer_key, count);
                if count >= self.min_rule_usage {
                     println!("DEBUG: Creating rule for k-mer with count {}", count);
                     let representative_kmer = &occurrences[0].1;

                     let new_rule = Rule {
                         id: self.next_rule_id,
                         symbols: representative_kmer.clone(),
                         usage_count: count,
                         positions: occurrences.iter().map(|(pos, _)| *pos).collect(),
                         depth: None,
                         assembly_index: None, // Initialize assembly_index
                     };
                     self.rules.insert(self.next_rule_id, new_rule);

                     let start = Instant::now();
                     self.replace_kmer_occurrences(self.next_rule_id, &occurrences)?;
                     self.metrics.replacement_time += start.elapsed();
                     self.metrics.replacement_count += 1;

                     self.next_rule_id += 1;
                     true // Pattern processed
                } else {
                    println!("DEBUG: K-mer count {} is below minimum usage {}", count, self.min_rule_usage);
                    false // No frequent pattern found
                }
            } else {
                println!("DEBUG: No frequent k-mers found");
                false // No pattern found
            }
        };

        if pattern_found {
            // Check for rule eviction if necessary after processing a pattern
            if let Some(max_count) = self.max_rule_count {
                if self.rules.len() > max_count {
                    let evict_start = Instant::now();
                    self.evict_rules(max_count)?;
                    self.metrics.eviction_time += evict_start.elapsed();
                }
            }
            Ok(true)
        } else {
            Ok(false)
        }
    }
    
    /// Replaces a k-mer with a rule at all occurrences.
    /// Returns the number of replacements made.
    fn replace_kmer_occurrences(
        &mut self,
        rule_id: usize,
        occurrences: &[(usize, Vec<Symbol>)]
    ) -> Result<()> {
        if occurrences.is_empty() {
            return Ok(());
        }
        
        // Get the k-mer size
        let k = self.kmer_table.as_ref().map_or(2, |table| table.k());
        
        // Create a list of positions to replace
        let positions: Vec<usize> = occurrences.iter().map(|(pos, _)| *pos).collect();
        
        // Create a new sequence with the replacements
        let mut new_sequence = Vec::with_capacity(self.sequence.len());
        let mut i = 0;
        
        while i < self.sequence.len() {
            if positions.contains(&i) {
                // Add a reference to the new rule
                new_sequence.push(Symbol::non_terminal(i, rule_id, Direction::Forward));
                i += k; // Skip ahead by k symbols
            } else {
                // Keep the original symbol
                new_sequence.push(self.sequence[i]);
                i += 1;
            }
        }
        
        // Update the sequence
        self.sequence = new_sequence;
        
        Ok(())
    }
    
    /// Enables streaming mode for processing large sequences in chunks.
    pub fn enable_streaming_mode(mut self) -> Self {
        self.stream_mode = true;
        self
    }

    /// Process a chunk of the sequence incrementally
    pub fn process_sequence_chunk(&mut self, chunk: &[EncodedBase], source_grammar_id: usize, chunk_offset: usize) -> Result<()> {
        let chunk_start_time = Instant::now();
        self.chunk_count += 1;
        self.total_bases_processed += chunk.len();
        
        if self.sequence.is_empty() && chunk_offset == 0 { // Ensure it's truly the first chunk of a sequence
            // First chunk - initialize the sequence
            self.initialize_sequence(chunk, source_grammar_id);
        } else {
            // Add this chunk to the existing sequence
            let start_symbol_id = self.sequence.len(); // ID within the current builder's sequence context
            let chunk_symbols: Vec<Symbol> = chunk
                .iter()
                .enumerate()
                .map(|(i, &base)| {
                    // start_symbol_id + i is the ID in the context of this GrammarBuilder's sequence
                    // chunk_offset + i is the original position in the source file/chromosome
                    Symbol::terminal(start_symbol_id + i, base, Direction::Forward, Some(source_grammar_id), Some(chunk_offset + i))
                })
                .collect();
        
            self.sequence.extend(chunk_symbols);
        }
        
        // After adding the chunk, rebuild pattern table and process rules
        self.rebuild_pattern_table()?;
        
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
        while self.step()? && steps_performed < max_steps_per_chunk {
            steps_performed += 1;
        }
        
        // In streaming mode, evict rules if needed to control memory usage
        if self.stream_mode {
            if let Some(max_count) = self.max_rule_count {
                if self.rules.len() > max_count {
                    let evict_start = Instant::now();
                    self.evict_rules(max_count)?;
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
    
    /// Finalizes the grammar after all chunks have been processed.
    pub fn finalize_grammar(&mut self) -> Result<()> {
        // For streaming mode, finalize the construction
        if !self.stream_mode {
            println!("Not in streaming mode, nothing to finalize.");
            return Ok(());
        }
        
        println!("Finalizing grammar after {} chunks, {} total bases processed.", 
                 self.chunk_count, self.total_bases_processed);
                 
        // Continue processing until no more patterns meet the threshold
        let mut steps_performed = 0;
        while self.step()? {
            steps_performed += 1;
            
            // Report progress periodically
            if steps_performed % 100 == 0 {
                println!("  Finalization progress: {} additional steps", steps_performed);
            }
            
            // Periodic rule eviction during finalization, if enabled
            if let Some(max_count) = self.max_rule_count {
                if steps_performed % 50 == 0 && self.rules.len() > max_count {
                    let evict_start = Instant::now();
                    self.evict_rules(max_count)?;
                    self.metrics.eviction_time += evict_start.elapsed();
                }
            }
        }
        
        // Final inlining of single-use rules
        let inline_start = Instant::now();
        self.inline_single_use_rules();
        self.metrics.inlining_time += inline_start.elapsed();
        
        println!("Grammar finalization complete. Performed {} additional steps.", steps_performed);
        println!("Final stats: {} symbols, {} rules", self.sequence.len(), self.rules.len());
        
        Ok(())
    }

    /// Builds a grammar from the given sequence.
    /// This is the main entry point for grammar construction.
    pub fn build_grammar(&mut self, initial_sequence: &[EncodedBase], source_grammar_id: usize) -> Result<()> {
        // Use GPU acceleration if context is available
        if self.gpu_context.is_some() {
            match self.build_grammar_with_gpu(initial_sequence, source_grammar_id) {
                Ok(_) => return Ok(()),
                Err(e) => {
                    println!("GPU processing failed: {}. Falling back to CPU.", e);
                    // Fall through to CPU implementation
                }
            }
        }
        
        // CPU implementation follows
        let start = Instant::now();
        
        // Initialize the sequence if not in streaming mode
        if !self.stream_mode {
            self.initialize_sequence(initial_sequence, source_grammar_id);
            // Debug: print first 10 symbols after initialization
            println!("DEBUG: First 10 symbols after initialization: {:#?}", &self.sequence.iter().take(10).collect::<Vec<_>>());
        } else {
            // In streaming mode, we process the chunk directly
            // build_grammar is called with a single sequence, so chunk_offset is 0
            return self.process_sequence_chunk(initial_sequence, source_grammar_id, 0);
        }
        
        // Build the digram table
        let table_start = Instant::now();
        self.rebuild_pattern_table()?;
        // Debug: print top 5 digrams after building digram table
        let top_digrams = self.digram_table.get_top_digrams(5);
        println!("DEBUG: Top 5 digrams after table build:");
        for (i, (key, occurrences)) in top_digrams.iter().enumerate() {
            println!("  {}. Key: {:?}, Count: {}", i+1, key, occurrences.len());
        }
        self.metrics.digram_table_time += table_start.elapsed();
        
        // Iteratively replace digrams
        let mut steps = 0;
        while self.step()? {
            steps += 1;
            if steps % 1000 == 0 {
                println!("Step {}, Rule count: {}", steps, self.rules.len());
            }
            
            // Check for rule count limit and evict if necessary
            if let Some(max_rules) = self.max_rule_count {
                if self.rules.len() > max_rules {
                    let evict_start = Instant::now();
                    self.evict_least_used_rules(max_rules)?;
                    self.metrics.eviction_time += evict_start.elapsed();
                }
            }
        }
        
        // Inline rules that are used only once, as they don't contribute to compression
        let inline_start = Instant::now();
        self.inline_single_use_rules();
        self.metrics.inlining_time += inline_start.elapsed();
        
        // Calculate rule depths for hierarchy analysis
        self.calculate_rule_depths();
        
        // Calculate Assembly Indices for all rules
        match calculate_rule_assembly_indices(&mut self.rules) {
            Ok(_) => log::debug!("Assembly indices calculated successfully."),
            Err(e) => {
                log::warn!("Failed to calculate assembly indices: {}. Proceeding without them.", e);
                // Decide if this should be a hard error or just a warning.
                // For now, logging as warning and proceeding.
            }
        }

        // Print stats
        self.metrics.step_count = steps;
        println!("Grammar construction completed in {:?}", start.elapsed());
        println!("Rules: {}, Steps: {}", self.rules.len(), steps);
        log::info!("Grammar construction complete. Total rules: {}, Max depth: {}", self.rules.len(), self.get_max_rule_depth());
        Ok(())
    }
    
    /// Process the grammar using GPU acceleration, including digram finding and suffix array construction.
    pub fn build_grammar_with_gpu(&mut self, initial_sequence: &[EncodedBase], source_grammar_id: usize) -> Result<()> {
        if self.gpu_context.is_none() {
            log::warn!("GPU context not provided, falling back to CPU implementation.");
            return self.build_grammar(initial_sequence, source_grammar_id); // Pass source_grammar_id to fallback
        }
        log::info!("Starting grammar construction with GPU acceleration...");
        // Clone the GpuContext to avoid borrow issues
        let gpu_context = self.gpu_context.as_ref().unwrap().clone(); 

        self.initialize_sequence(initial_sequence, source_grammar_id);
        let mut current_rule_id_counter = 0; 

        loop {
            self.metrics.step_count += 1;
            if self.metrics.step_count % 10 == 0 {
                log::info!(
                    "  [GPU Step: {}] SeqLen: {}, Rules: {}, TotalBases: {}",
                    self.metrics.step_count,
                    self.sequence.len(),
                    self.rules.len(),
                    self.total_bases_processed
                );
            }

            let mut gpu_sequence = GpuSequence::from_symbols(&self.sequence)?;
            // Pass the cloned gpu_context
            gpu_sequence.upload_to_gpu(&gpu_context)?;

            match gpu_sequence.find_most_frequent_digram_enhanced(self.min_rule_usage, self.reverse_aware, &gpu_context)? {
                Some((_key_hash, occurrences_gpu)) => { 
                    if occurrences_gpu.is_empty() || occurrences_gpu.len() < self.min_rule_usage {
                        log::info!("No frequent digrams found by GPU or count below threshold. Finalizing grammar.");
                        break;
                    }

                    let first_occurrence_symbols = &occurrences_gpu[0].1;
                    let new_rule = Rule {
                        id: current_rule_id_counter, 
                        symbols: vec![first_occurrence_symbols.0.clone(), first_occurrence_symbols.1.clone()],
                        usage_count: occurrences_gpu.len(),
                        positions: occurrences_gpu.iter().map(|(pos, _syms)| *pos).collect(),
                        depth: None, 
                        assembly_index: None, // Initialize assembly_index
                    };
                    self.rules.insert(current_rule_id_counter, new_rule);

                    let replacement_start = Instant::now();
                    self.replace_digram_occurrences(current_rule_id_counter, &occurrences_gpu)?;
                    self.metrics.replacement_time += replacement_start.elapsed();
                    self.metrics.replacement_count += 1;

                    current_rule_id_counter += 1;
                }
                None => {
                    log::info!("No frequent digrams found by GPU. Finalizing grammar.");
                    break;
                }
            }
            if let Some(max_rules) = self.max_rule_count {
                if self.rules.len() > max_rules {
                    self.evict_rules(max_rules)?;
                }
            }
        }

        self.finalize_grammar()?;

        // Calculate rule depths for hierarchy analysis
        self.calculate_rule_depths();

        // Calculate Assembly Indices for all rules
        match calculate_rule_assembly_indices(&mut self.rules) {
            Ok(_) => log::debug!("Assembly indices calculated successfully after GPU build."),
            Err(e) => {
                log::warn!("Failed to calculate assembly indices after GPU build: {}. Proceeding without them.", e);
            }
        }

        log::info!("GPU-accelerated grammar construction complete. Total rules: {}, Max depth: {}", self.rules.len(), self.get_max_rule_depth());
        Ok(())
    }

    /// Calculate the hierarchical depth of each rule in the grammar using recursion and memoization.
    fn calculate_rule_depths(&mut self) {
        let mut depths = HashMap::new();
        let mut max_depth = 0;
        let rule_ids: Vec<usize> = self.rules.keys().cloned().collect(); // Collect keys first

        for rule_id in rule_ids {
            let mut visited_stack = HashSet::new(); // Reset cycle detection for each top-level call
            match self._calculate_rule_depth_recursive(rule_id, &mut depths, &mut visited_stack) {
                Ok(depth) => {
                    if let Some(rule) = self.rules.get_mut(&rule_id) {
                        rule.depth = Some(depth);
                        max_depth = max_depth.max(depth);
                    } else {
                        println!("Warning: Rule {} not found while setting depth.", rule_id);
                    }
                }
                Err(e) => {
                    println!("Warning: {} while calculating depth for rule {}. Assigning depth 0.", e, rule_id);
                    if let Some(rule) = self.rules.get_mut(&rule_id) {
                        rule.depth = Some(0); // Assign a default depth on cycle
                    }
                }
            }
        }
        // Store calculated depths in self.rule_depths as well, if needed elsewhere
        self.rule_depths = self.rules.iter()
            .filter_map(|(id, rule)| rule.depth.map(|d| (*id, d)))
            .collect();
            
        // Note: The Grammar struct's max_depth should be set when returning the grammar.
        // This function primarily focuses on setting the depth within each Rule struct.
    }

    /// Recursive helper function to calculate depth for a single rule with memoization and cycle detection.
    fn _calculate_rule_depth_recursive(
        &self,
        rule_id: usize,
        depth_cache: &mut HashMap<usize, Result<usize, String>>, // Cache stores Result to handle cycles
        visited_stack: &mut HashSet<usize> // Track nodes in current recursion path
    ) -> Result<usize, String> { // Return String error for simpler handling

        // Check cache first
        if let Some(cached_result) = depth_cache.get(&rule_id) {
            return cached_result.clone().map_err(|e| e.clone());
        }

        // Check for cycle
        if !visited_stack.insert(rule_id) {
            let error_msg = format!("Cycle detected involving rule {}", rule_id);
            depth_cache.insert(rule_id, Err(error_msg.clone())); // Cache cycle error
            return Err(error_msg);
        }

        let rule = match self.rules.get(&rule_id) {
            Some(r) => r,
            None => {
                 visited_stack.remove(&rule_id); // Backtrack
                 let error_msg = format!("Rule {} not found during depth calculation", rule_id);
                 // Don't cache error for non-existent rule?
                 return Err(error_msg);
            }
        };

        let mut max_sub_depth = 0;
        // Clone symbols to avoid borrowing self while calling recursively
        let rule_symbols = rule.symbols.clone(); 
        
        for symbol in &rule_symbols {
            if let SymbolType::NonTerminal(sub_rule_id) = symbol.symbol_type {
                match self._calculate_rule_depth_recursive(sub_rule_id, depth_cache, visited_stack) {
                    Ok(depth) => max_sub_depth = max_sub_depth.max(depth),
                    Err(e) => {
                        // Cycle detected deeper down, propagate the error
                        visited_stack.remove(&rule_id); // Backtrack
                        depth_cache.insert(rule_id, Err(e.clone())); // Cache the error
                        return Err(e);
                    }
                }
            }
        }

        visited_stack.remove(&rule_id); // Backtrack: remove from current path
        let current_depth = max_sub_depth + 1;
        depth_cache.insert(rule_id, Ok(current_depth)); // Cache successful result
        Ok(current_depth)
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
    pub fn process_bytes_chunk(&mut self, chunk: &[u8], source_grammar_id: usize, chunk_offset: usize) -> Result<()> {
        // Convert raw bytes to EncodedBase
        let encoded_chunk: Vec<EncodedBase> = chunk.iter()
            .filter_map(|&b| EncodedBase::from_base(b))
            .collect();
        
        // Use the regular process_sequence_chunk method
        self.process_sequence_chunk(&encoded_chunk, source_grammar_id, chunk_offset)
    }

    /// Evicts rules to keep the number of rules below a maximum.
    /// This method will inline the least frequently used rules.
    fn evict_rules(&mut self, max_count: usize) -> Result<()> {
        self.evict_least_used_rules(max_count)
    }
    
    /// Evicts the least used rules to keep the number of rules below the specified maximum.
    fn evict_least_used_rules(&mut self, max_count: usize) -> Result<()> {
        if self.rules.len() <= max_count {
            return Ok(()); // Nothing to do
        }
        
        let start = Instant::now();
        
        // Find rules to evict (rules with the lowest usage count)
        let mut rule_usages: Vec<(usize, usize)> = self.rules.iter()
            .map(|(id, rule)| (*id, rule.usage_count))
            .collect();
        
        // Sort by usage count (ascending)
        rule_usages.sort_by_key(|(_, count)| *count);
        
        // Evict rules until we're under the limit
        let num_to_evict = self.rules.len() - max_count;
        for i in 0..num_to_evict {
            let (rule_id, _) = rule_usages[i];
            
            // Inline rule at all occurrences
            // For simplicity, we'll just rebuild the sequence with this rule replaced
            self.inline_rule(rule_id)?;
            
            // Remove the rule
            self.rules.remove(&rule_id);
        }
        
        self.metrics.eviction_time += start.elapsed();
        println!("Evicted {} rules (keeping at most {})", num_to_evict, max_count);
        
        // Rebuild the pattern table after inlining and removal to reflect changes
        // self.rebuild_pattern_table()?;

        Ok(())
    }
    
    /// Inlines a specific rule at all its occurrences.
    fn inline_rule(&mut self, rule_id_to_inline: usize) -> Result<()> {
        let start_time = Instant::now();

        let rule_to_inline = match self.rules.get(&rule_id_to_inline) {
            Some(r) => r.clone(),
            None => {
                // Rule not found, nothing to inline.
                // This can happen if the rule was already inlined or removed.
                return Ok(());
            }
        };

        // --- Step 1: Identify inlining sites and collect digrams to be removed ---
        let mut sites_to_inline: Vec<usize> = Vec::new();
        for (i, sym) in self.sequence.iter().enumerate() {
            if let SymbolType::NonTerminal(id) = sym.symbol_type {
                if id == rule_id_to_inline {
                    sites_to_inline.push(i);
                }
            }
        }

        if sites_to_inline.is_empty() {
            // No occurrences of the rule in the current sequence to inline.
            return Ok(());
        }

        let mut digrams_to_remove: HashSet<(DigramKeyTuple, usize, DigramSource)> = HashSet::new();
        for &idx in &sites_to_inline {
            // Digram (S_{idx-1}, S_{idx}) where S_{idx} is the NT being inlined.
            if idx > 0 {
                // Ensure sequence[idx-1] is valid before accessing.
                if let Some(s1_ref) = self.sequence.get(idx - 1) {
                    let s2_val = self.sequence[idx]; // This is the NT being inlined.
                    let key = DigramTable::canonical_key((s1_ref, &s2_val), self.reverse_aware);
                    digrams_to_remove.insert((key, idx - 1, DigramSource::Original));
                }
            }
            // Digram (S_{idx}, S_{idx+1}) where S_{idx} is the NT being inlined.
            if idx + 1 < self.sequence.len() {
                let s1_val = self.sequence[idx]; // This is the NT being inlined.
                if let Some(s2_ref) = self.sequence.get(idx + 1) {
                    let key = DigramTable::canonical_key((&s1_val, s2_ref), self.reverse_aware);
                    digrams_to_remove.insert((key, idx, DigramSource::Original));
                }
            }
        }

        // --- Step 2: Perform sequence inlining ---
        let mut new_sequence = Vec::with_capacity(self.sequence.len());
        let mut current_new_idx = 0;
        // Store (start_new_sequence_idx, end_new_sequence_idx) for each block of inlined symbols.
        let mut inlined_block_boundaries_in_new_seq: Vec<(usize, usize)> = Vec::new();

        let mut last_original_idx_processed = 0;
        for &inline_site_orig_idx in &sites_to_inline {
            // Copy symbols from original sequence before the current inline site.
            if inline_site_orig_idx > last_original_idx_processed {
                for i in last_original_idx_processed..inline_site_orig_idx {
                    let mut sym_to_copy = self.sequence[i];
                    sym_to_copy.id = current_new_idx;
                    new_sequence.push(sym_to_copy);
                    current_new_idx += 1;
                }
            }

            // Inline the rule symbols.
            let block_start_new_idx = current_new_idx;
            if !rule_to_inline.symbols.is_empty() {
                for rule_sym in &rule_to_inline.symbols {
                    let mut sym_to_insert = *rule_sym;
                    sym_to_insert.id = current_new_idx;
                    new_sequence.push(sym_to_insert);
                    current_new_idx += 1;
                }
                inlined_block_boundaries_in_new_seq.push((block_start_new_idx, current_new_idx));
            }
            // If rule_to_inline.symbols is empty, the NT is effectively deleted.
            // No block is added to inlined_block_boundaries_in_new_seq.

            last_original_idx_processed = inline_site_orig_idx + 1;
        }

        // Copy any remaining symbols from the original sequence.
        if last_original_idx_processed < self.sequence.len() {
            for i in last_original_idx_processed..self.sequence.len() {
                let mut sym_to_copy = self.sequence[i];
                sym_to_copy.id = current_new_idx;
                new_sequence.push(sym_to_copy);
                current_new_idx += 1;
            }
        }
        
        self.sequence = new_sequence;

        // --- Step 3: Update DigramTable ---
        // Remove old digrams that were destroyed.
        for (key, pos, source) in digrams_to_remove {
            self.digram_table.remove_occurrence(&key, pos, source);
        }

        // Add new digrams formed by the inlining.
        // This includes digrams internal to the inlined blocks and boundary digrams.
        let mut scan_points = HashSet::new(); // Points around which to scan for new digrams

        for &(block_start, block_end) in &inlined_block_boundaries_in_new_seq {
            // Add points for scanning internal to the block and its boundaries
            if block_start > 0 {
                scan_points.insert(block_start - 1); // For digram ending at block_start - 1, starting block_start - 1
            }
            for i in block_start..block_end { // For digrams starting from block_start up to block_end - 1
                scan_points.insert(i);
            }
        }
        
        // Add digrams that might have been formed if an NT was simply deleted (empty rule inlined)
        // This involves checking positions around the original `sites_to_inline`.
        // This part is tricky due to coordinate changes. A simpler approach is to scan affected regions.
        // The current `scan_points` approach is more targeted.

        for &scan_pos in &scan_points {
            if scan_pos + 1 < self.sequence.len() {
                 // Ensure symbols being passed to add_digram have their `id` field (position) correct.
                 // self.sequence[scan_pos].id should be scan_pos, and self.sequence[scan_pos+1].id should be scan_pos+1
                 // This was handled during new_sequence construction.
                self.digram_table.add_digram(
                    scan_pos, // Position of the first symbol of the digram
                    self.sequence[scan_pos],
                    self.sequence[scan_pos+1],
                    self.reverse_aware
                );
            }
        }
        
        // If a non-terminal was replaced by an empty rule, a new digram might form
        // by joining S_{idx-1} and S_{idx+1} from the original sequence.
        // This needs careful handling of indices in the new sequence.
        // The current `scan_points` derived from `inlined_block_boundaries_in_new_seq`
        // might miss this if the rule was empty.

        // A more robust way for adding new digrams after any modification:
        // Determine a range [min_affected_idx, max_affected_idx] in the *new* sequence.
        // Then iterate from max(0, min_affected_idx - 1) up to min(len-2, max_affected_idx)
        // and add all digrams (s_i, s_{i+1}).
        
        // For now, the scan_points approach will be used.
        // It covers digrams within inserted blocks and those formed with one adjacent original symbol.

        self.metrics.inlining_time += start_time.elapsed();
        Ok(())
    }
    
    /// Inlines all rules that are used only once.
    fn inline_single_use_rules(&mut self) {
        let start = Instant::now();
        
        // Find rules used only once
        let single_use_rules: Vec<usize> = self.rules.iter()
            .filter(|(_, rule)| rule.usage_count == 1)
            .map(|(id, _)| *id)
            .collect();
        
        // Inline each rule
        for rule_id in single_use_rules {
            if let Err(e) = self.inline_rule(rule_id) {
                println!("Warning: Failed to inline rule {}: {}", rule_id, e);
                continue;
            }
            
            // Remove the rule
            self.rules.remove(&rule_id);
        }
        
        self.metrics.inlining_time += start.elapsed();
    }

    fn replace_digram_occurrences(
        &mut self,
        rule_id: usize,
        occurrences: &[(usize, (Symbol, Symbol))]
    ) -> Result<()> {
        if occurrences.is_empty() || self.sequence.is_empty() {
            return Ok(());
        }
        let mut new_sequence = Vec::with_capacity(self.sequence.len());
        let mut current_pos = 0;
        let mut sorted_occurrences = occurrences.to_vec();
        sorted_occurrences.sort_by_key(|(pos, _)| *pos);

        let mut occurrence_iter = sorted_occurrences.iter().peekable();

        while current_pos < self.sequence.len() {
            if let Some((occ_pos, _)) = occurrence_iter.peek() {
                if current_pos == *occ_pos {
                    // Replace digram with new rule symbol
                    new_sequence.push(Symbol::non_terminal(
                        current_pos, // Original position for the new symbol
                        rule_id,
                        Direction::Forward, // Assuming Forward for now
                    ));
                    // Advance past the two symbols that formed the digram
                    current_pos += 2;
                    occurrence_iter.next(); // Consume this occurrence

                    // Skip any overlapping occurrences
                    while let Some((next_occ_pos, _)) = occurrence_iter.peek() {
                        if *next_occ_pos < current_pos {
                            occurrence_iter.next(); // Consume overlapping occurrence
                        } else {
                            break;
                        }
                    }
                    continue;
                }
            }
            // Copy symbol if not part of a replaced digram
            new_sequence.push(self.sequence[current_pos].clone());
            current_pos += 1;
        }
        self.sequence = new_sequence;
        Ok(())
    }
}

#[derive(Debug, Clone)]
pub struct RuleCandidate {
    pub symbols: Vec<Symbol>,
    pub key_tuple: Option<DigramKeyTuple>,
    pub count: usize,
    pub positions: Vec<usize>,
    pub non_terminal_id: Option<usize>,
}

pub struct RuleBuilder {
    pub grammar: Grammar,
    pub table: DigramTable,
    pub kmer_table: KmerTable,
    pub next_rule_id: usize,
    pub min_count: usize,
}

impl RuleBuilder {
    pub fn new(sequence: Vec<Symbol>, min_count: usize, reverse_aware_digrams: bool, kmer_size: Option<usize>) -> Self {
        let initial_rules = HashMap::new();
        let grammar = Grammar {
            sequence: sequence.clone(),
            rules: initial_rules,
            max_depth: 0,
            origins: HashMap::new(), // Initialize origins
        };
        let table = DigramTable::build(&sequence, reverse_aware_digrams);
        let kmer_table = if let Some(k) = kmer_size {
            KmerTable::build(&sequence, k, false)
        } else {
            KmerTable::new(2)
        };
        Self {
            grammar,
            table,
            kmer_table,
            next_rule_id: 0,
            min_count,
        }
    }

    pub fn select_candidate_pair(&self) -> Option<RuleCandidate> {
        // Prioritize DigramTable for now
        if let Some((digram_key, count, occurrences)) = self.table.find_most_frequent_digram() {
            if count >= self.min_count && !occurrences.is_empty() {
                let first_pos = occurrences[0].0;
                if first_pos + 1 < self.grammar.sequence.len() {
                    let s1 = self.grammar.sequence[first_pos].clone();
                    let s2 = self.grammar.sequence[first_pos + 1].clone();
                    let positions = occurrences.iter().map(|(p, _source)| *p).collect::<Vec<usize>>();
                    
                    return Some(RuleCandidate {
                        symbols: vec![s1, s2],
                        key_tuple: Some(digram_key),
                        count,
                        positions,
                        non_terminal_id: None,
                    });
                } else {
                    log::warn!("DigramTable: Frequent digram's first position {} is too close to end of sequence len {}.", first_pos, self.grammar.sequence.len());
                    return None; 
                }
            }
        }
        None
    }

    // ... other methods ...
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
        let seq = encode_seq(b"AAAAAA");
        builder.build_grammar(&seq, 0)?;
        let (_, rules) = builder.get_grammar();
        // There should be a rule of length 2 with A then A
        let has_aa_rule = rules.values().any(|r| r.symbols.len() == 2 &&
            matches!(r.symbols[0].symbol_type, SymbolType::Terminal(t) if t == EncodedBase(0)) &&
            matches!(r.symbols[1].symbol_type, SymbolType::Terminal(t) if t == EncodedBase(0)));
        assert!(has_aa_rule, "No AA rule found");
        Ok(())
    }

    #[test]
    fn test_overlapping_patterns() -> Result<()> {
        let mut builder = GrammarBuilder::new(2, false);
        let seq = encode_seq(b"ABCABCABC");
        builder.build_grammar(&seq, 0)?;
        let (_final_sequence, rules) = builder.get_grammar();
        // There should be at least one rule of length 2 with valid terminal content
        let found = rules.values().any(|r| r.symbols.len() == 2 &&
            matches!(r.symbols[0].symbol_type, SymbolType::Terminal(_)) &&
            matches!(r.symbols[1].symbol_type, SymbolType::Terminal(_)));
        assert!(found, "No valid 2-symbol rule found");
        Ok(())
    }

    #[test]
    fn test_nested_rules() -> Result<()> {
        let mut builder = GrammarBuilder::new(2, false);
        let seq = encode_seq(b"ACACACACACACACAC");
        builder.build_grammar(&seq, 0)?;
        let (_final_sequence, rules) = builder.get_grammar();
        // There should be at least one rule that contains a non-terminal (nested rule)
        let found = rules.values().any(|r| r.symbols.iter().any(|s| matches!(s.symbol_type, SymbolType::NonTerminal(_))));
        assert!(found, "No nested rule found");
        Ok(())
    }

    #[test]
    fn test_rule_utility() -> Result<()> {
        let mut builder = GrammarBuilder::new(2, false);
        let seq = encode_seq(b"AAAACCCC");
        builder.build_grammar(&seq, 0)?;
        let (_final_sequence, rules) = builder.get_grammar();
        // Find rules by content
        let aa_rule = rules.values().find(|r| r.symbols.len() == 2 &&
            matches!(r.symbols[0].symbol_type, SymbolType::Terminal(t) if t == EncodedBase(0)) &&
            matches!(r.symbols[1].symbol_type, SymbolType::Terminal(t) if t == EncodedBase(0)));
        let cc_rule = rules.values().find(|r| r.symbols.len() == 2 &&
            matches!(r.symbols[0].symbol_type, SymbolType::Terminal(t) if t == EncodedBase(1)) &&
            matches!(r.symbols[1].symbol_type, SymbolType::Terminal(t) if t == EncodedBase(1)));
        assert!(aa_rule.is_some(), "No AA rule found");
        assert!(cc_rule.is_some(), "No CC rule found");
        Ok(())
    }

    #[test]
    fn test_no_frequent_digrams() -> Result<()> {
        let mut builder = GrammarBuilder::new(2, false);
        let seq = encode_seq(b"ACGT");
        builder.build_grammar(&seq, 0)?;
        let (final_sequence, rules) = builder.get_grammar();
        assert!(rules.is_empty());
        assert_eq!(final_sequence.len(), 4);
        for sym in final_sequence.iter() {
            assert!(matches!(sym.symbol_type, SymbolType::Terminal(_)));
        }
        Ok(())
    }

    #[test]
    fn test_empty_sequence() -> Result<()> {
        let mut builder = GrammarBuilder::new(2, false);
        let seq = encode_seq(b"");
        builder.build_grammar(&seq, 0)?;
        let (final_sequence, rules) = builder.get_grammar();
        assert!(rules.is_empty());
        assert!(final_sequence.is_empty());
        Ok(())
    }

    #[test]
    fn test_single_base_sequence() -> Result<()> {
        let mut builder = GrammarBuilder::new(2, false);
        let seq = encode_seq(b"A");
        builder.build_grammar(&seq, 0)?;
        let (final_sequence, rules) = builder.get_grammar();
        assert!(rules.is_empty());
        assert_eq!(final_sequence.len(), 1);
        assert!(matches!(final_sequence[0].symbol_type, SymbolType::Terminal(t) if t == EncodedBase(0)));
        Ok(())
    }

    #[test]
    fn test_reverse_complement_aware() -> Result<()> {
        let mut builder_aware = GrammarBuilder::new(2, true);
        let seq = encode_seq(b"ACGTACGT");
        builder_aware.build_grammar(&seq, 0)?;
        let (_final_sequence_aware, rules_aware) = builder_aware.get_grammar();
        assert!(!rules_aware.is_empty(), "Reverse aware should find repeating patterns");
        Ok(())
    }

    // Add more tests:
    // - Edge cases: single base repeats (AAAA), alternating bases (ABAB)
    // - Minimum rule usage > 2
    // - Sequences that generate cycles (needs careful handling or cycle detection)
    // - Very long sequences (if performance testing is desired here)
} 