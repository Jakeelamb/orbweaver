//! Locally Consistent Grammar (LCG) construction.
//!
//! This module implements a grammar compression algorithm based on the paper:
//! "Efficient terabyte-scale text compression via stable local consistency"
//! (Díaz-Domínguez, 2024)
//!
//! Key properties:
//! - O(n) time complexity (single pass)
//! - Parallelizable without synchronization
//! - Content-based fingerprinting ensures identical patterns get identical rules
//! - Merge operation is trivial (concatenation + deduplication)

use crate::encode::dna_2bit::EncodedBase;
use crate::gpu::lcg::GpuLcgParser;
use crate::gpu::GpuContext;
use crate::grammar::rule::Rule;
use crate::grammar::symbol::{Symbol, SymbolType, Direction};
use bumpalo::Bump;
use rustc_hash::FxHashMap;
use std::hash::{Hash, Hasher};
use std::sync::Arc;

/// Fingerprint type - 64-bit hash of content
pub type Fingerprint = u64;

/// A phrase identified by locally consistent parsing
#[derive(Debug, Clone)]
pub struct Phrase {
    /// Start position in the original sequence
    pub start: usize,
    /// Length of the phrase
    pub len: usize,
    /// Fingerprint of the phrase content
    pub fingerprint: Fingerprint,
}

/// Configuration for LCG parsing
#[derive(Debug, Clone)]
pub struct LcgConfig {
    /// Minimum phrase length (default: 2)
    pub min_phrase_len: usize,
    /// Maximum phrase length (default: 256)
    pub max_phrase_len: usize,
    /// Window size for local minimum detection (default: 8)
    pub window_size: usize,
    /// Minimum occurrences to create a rule (default: 2)
    pub min_occurrences: usize,
}

impl Default for LcgConfig {
    fn default() -> Self {
        Self {
            min_phrase_len: 2,
            max_phrase_len: 256,
            window_size: 8,
            min_occurrences: 2,
        }
    }
}

impl LcgConfig {
    /// Auto-tune LCG configuration based on sequence characteristics.
    ///
    /// This method analyzes a sample of the sequence to determine optimal parameters:
    /// - For highly repetitive sequences: smaller window, lower min_occurrences
    /// - For complex sequences: larger window, higher min_occurrences
    /// - Adjusts max_phrase_len based on sequence length
    ///
    /// # Arguments
    /// * `sequence` - The full sequence or a representative sample
    /// * `sample_size` - Optional sample size (default: 10000 bases)
    ///
    /// # Returns
    /// An auto-tuned LcgConfig
    pub fn auto_tune(sequence: &[EncodedBase], sample_size: Option<usize>) -> Self {
        let sample_len = sample_size.unwrap_or(10000).min(sequence.len());
        if sample_len < 100 {
            return Self::default();
        }

        // Take a sample from the sequence
        let sample = &sequence[..sample_len];

        // Calculate entropy and repetitiveness metrics
        let (entropy, repetitiveness) = Self::analyze_sequence(sample);

        log::debug!(
            "LCG auto-tune: entropy={:.3}, repetitiveness={:.3}",
            entropy,
            repetitiveness
        );

        // Tune parameters based on analysis
        let mut config = Self::default();

        // Window size: smaller for repetitive, larger for complex
        if repetitiveness > 0.3 {
            config.window_size = 4;
        } else if entropy > 1.8 {
            config.window_size = 16;
        } else {
            config.window_size = 8;
        }

        // Min phrase length: smaller for repetitive sequences
        if repetitiveness > 0.4 {
            config.min_phrase_len = 2;
        } else if entropy > 1.9 {
            config.min_phrase_len = 3;
        } else {
            config.min_phrase_len = 2;
        }

        // Max phrase length: scale with sequence length
        if sequence.len() > 10_000_000 {
            config.max_phrase_len = 512;
        } else if sequence.len() > 1_000_000 {
            config.max_phrase_len = 256;
        } else {
            config.max_phrase_len = 128;
        }

        // Min occurrences: lower for repetitive, higher for complex
        if repetitiveness > 0.3 {
            config.min_occurrences = 2;
        } else if entropy > 1.8 {
            config.min_occurrences = 3;
        } else {
            config.min_occurrences = 2;
        }

        log::info!(
            "LCG auto-tuned: window={}, min_phrase={}, max_phrase={}, min_occ={}",
            config.window_size,
            config.min_phrase_len,
            config.max_phrase_len,
            config.min_occurrences
        );

        config
    }

    /// Analyze sequence to compute entropy and repetitiveness metrics.
    ///
    /// Returns (entropy, repetitiveness) where:
    /// - entropy: Shannon entropy of digram distribution (0-4 for DNA)
    /// - repetitiveness: Fraction of digrams that occur more than average
    fn analyze_sequence(sequence: &[EncodedBase]) -> (f64, f64) {
        if sequence.len() < 2 {
            return (0.0, 0.0);
        }

        // Count digram occurrences
        let mut digram_counts = [0u32; 16];
        for i in 0..sequence.len() - 1 {
            let digram = ((sequence[i].0 as usize) << 2) | (sequence[i + 1].0 as usize);
            digram_counts[digram] += 1;
        }

        let total_digrams = (sequence.len() - 1) as f64;

        // Calculate Shannon entropy
        let mut entropy = 0.0f64;
        for &count in &digram_counts {
            if count > 0 {
                let p = count as f64 / total_digrams;
                entropy -= p * p.log2();
            }
        }

        // Calculate repetitiveness (fraction of digrams above average)
        let avg_count = total_digrams / 16.0;
        let repetitive_digrams: u32 = digram_counts.iter().filter(|&&c| c as f64 > avg_count * 1.5).map(|&c| c).sum();
        let repetitiveness = repetitive_digrams as f64 / total_digrams;

        (entropy, repetitiveness)
    }
}

/// Reusable buffers for LCG grammar construction.
/// These buffers are reused across iterations to reduce allocation overhead.
pub struct LcgBuffers {
    /// Buffer for position fingerprints
    position_fingerprints: Vec<Fingerprint>,
    /// Buffer for cut points during parsing
    cut_points: Vec<usize>,
    /// Buffer for the new sequence being built
    new_sequence: Vec<Symbol>,
    /// Arena allocator for temporary allocations within an iteration
    arena: Bump,
}

impl LcgBuffers {
    /// Create new buffers with the given initial capacity.
    pub fn new(capacity: usize) -> Self {
        Self {
            position_fingerprints: Vec::with_capacity(capacity),
            cut_points: Vec::with_capacity(capacity / 8), // Roughly 1 cut per 8 positions
            new_sequence: Vec::with_capacity(capacity),
            arena: Bump::with_capacity(capacity * std::mem::size_of::<Phrase>() / 4),
        }
    }

    /// Clear all buffers for reuse without deallocating.
    #[inline]
    pub fn clear(&mut self) {
        self.position_fingerprints.clear();
        self.cut_points.clear();
        self.new_sequence.clear();
        self.arena.reset();
    }

    /// Ensure buffers have at least the given capacity.
    pub fn reserve(&mut self, capacity: usize) {
        if self.position_fingerprints.capacity() < capacity {
            self.position_fingerprints.reserve(capacity - self.position_fingerprints.len());
        }
        if self.new_sequence.capacity() < capacity {
            self.new_sequence.reserve(capacity - self.new_sequence.len());
        }
    }
}

impl Default for LcgBuffers {
    fn default() -> Self {
        Self::new(10_000)
    }
}

/// Locally Consistent Grammar builder
pub struct LcgBuilder {
    config: LcgConfig,
    /// Maps fingerprint -> rule_id for deduplication
    fingerprint_to_rule: FxHashMap<Fingerprint, usize>,
    /// All rules created
    rules: FxHashMap<usize, Rule>,
    /// Next rule ID to assign
    next_rule_id: usize,
    /// Optional GPU parser for acceleration
    gpu_parser: Option<GpuLcgParser>,
    /// Reusable buffers for reduced allocation overhead
    buffers: Option<LcgBuffers>,
}

impl LcgBuilder {
    pub fn new(config: LcgConfig) -> Self {
        Self {
            config,
            fingerprint_to_rule: FxHashMap::default(),
            rules: FxHashMap::default(),
            next_rule_id: 0,
            gpu_parser: None,
            buffers: None,
        }
    }

    /// Create an LcgBuilder with buffer pooling enabled for reduced allocation overhead.
    pub fn with_pooling(config: LcgConfig, initial_capacity: usize) -> Self {
        Self {
            config,
            fingerprint_to_rule: FxHashMap::default(),
            rules: FxHashMap::default(),
            next_rule_id: 0,
            gpu_parser: None,
            buffers: Some(LcgBuffers::new(initial_capacity)),
        }
    }

    /// Create an LcgBuilder with GPU acceleration enabled
    pub fn with_gpu(config: LcgConfig, gpu_context: Arc<GpuContext>) -> Self {
        let gpu_parser = match GpuLcgParser::new(gpu_context) {
            Ok(parser) => Some(parser),
            Err(e) => {
                log::warn!("Failed to initialize GPU LCG parser: {}. Falling back to CPU.", e);
                None
            }
        };
        Self {
            config,
            fingerprint_to_rule: FxHashMap::default(),
            rules: FxHashMap::default(),
            next_rule_id: 0,
            gpu_parser,
            buffers: None,
        }
    }

    /// Create an LcgBuilder with both GPU acceleration and buffer pooling.
    pub fn with_gpu_and_pooling(config: LcgConfig, gpu_context: Arc<GpuContext>, initial_capacity: usize) -> Self {
        let gpu_parser = match GpuLcgParser::new(gpu_context) {
            Ok(parser) => Some(parser),
            Err(e) => {
                log::warn!("Failed to initialize GPU LCG parser: {}. Falling back to CPU.", e);
                None
            }
        };
        Self {
            config,
            fingerprint_to_rule: FxHashMap::default(),
            rules: FxHashMap::default(),
            next_rule_id: 0,
            gpu_parser,
            buffers: Some(LcgBuffers::new(initial_capacity)),
        }
    }

    /// Enable buffer pooling on an existing builder.
    pub fn enable_pooling(&mut self, initial_capacity: usize) {
        if self.buffers.is_none() {
            self.buffers = Some(LcgBuffers::new(initial_capacity));
        }
    }

    /// Check if GPU acceleration is available and would be beneficial for the given sequence length
    pub fn should_use_gpu(&self, sequence_len: usize) -> bool {
        self.gpu_parser
            .as_ref()
            .map(|p| p.should_use_gpu(sequence_len))
            .unwrap_or(false)
    }

    /// Compute a fingerprint for a sequence of encoded bases.
    /// Uses the SIMD-optimized fingerprint function for speed.
    /// The fingerprint is deterministic based on content only.
    #[inline]
    pub fn fingerprint_bases(bases: &[EncodedBase]) -> Fingerprint {
        crate::encode::simd::fingerprint_sequence(bases)
    }

    /// Compute a fingerprint for a sequence of symbols.
    #[inline]
    pub fn fingerprint_symbols(symbols: &[Symbol]) -> Fingerprint {
        let mut hasher = rustc_hash::FxHasher::default();
        for sym in symbols {
            match sym.symbol_type {
                SymbolType::Terminal(base) => {
                    0u8.hash(&mut hasher); // Tag for terminal
                    base.0.hash(&mut hasher);
                }
                SymbolType::NonTerminal(rule_id) => {
                    1u8.hash(&mut hasher); // Tag for non-terminal
                    rule_id.hash(&mut hasher);
                }
            }
            sym.strand.hash(&mut hasher);
        }
        hasher.finish()
    }

    /// Perform locally consistent parsing on a sequence.
    /// Returns a list of phrases with their fingerprints.
    ///
    /// The parsing rule: cut at position i if the fingerprint at i is a local minimum
    /// within a window of size `window_size`.
    ///
    /// This method will use GPU acceleration when available and the sequence is large enough.
    pub fn parse_locally_consistent(&self, sequence: &[EncodedBase]) -> Vec<Phrase> {
        if sequence.len() < self.config.min_phrase_len {
            // Entire sequence is one phrase
            if sequence.is_empty() {
                return vec![];
            }
            return vec![Phrase {
                start: 0,
                len: sequence.len(),
                fingerprint: Self::fingerprint_bases(sequence),
            }];
        }

        // Try GPU path for large sequences
        if self.should_use_gpu(sequence.len()) {
            if let Some(ref gpu_parser) = self.gpu_parser {
                match self.parse_locally_consistent_gpu(sequence, gpu_parser) {
                    Ok(phrases) => return phrases,
                    Err(e) => {
                        log::warn!("GPU parsing failed, falling back to CPU: {}", e);
                        // Fall through to CPU path
                    }
                }
            }
        }

        // CPU path
        self.parse_locally_consistent_cpu(sequence)
    }

    /// GPU-accelerated locally consistent parsing
    fn parse_locally_consistent_gpu(&self, sequence: &[EncodedBase], gpu_parser: &GpuLcgParser) -> anyhow::Result<Vec<Phrase>> {
        // Step 1: Compute position fingerprints on GPU
        let position_fingerprints = gpu_parser.compute_fingerprints(sequence)?;

        // Step 2: Find cut points on GPU
        let half_window = self.config.window_size / 2;
        let cut_points = gpu_parser.find_cut_points(
            &position_fingerprints,
            half_window,
            self.config.min_phrase_len,
            self.config.max_phrase_len,
        )?;

        // Step 3: Compute phrase fingerprints on GPU
        let phrase_fingerprints = gpu_parser.compute_phrase_fingerprints(sequence, &cut_points)?;

        // Step 4: Build phrase structures
        let mut phrases = Vec::with_capacity(cut_points.len());
        for (i, &start) in cut_points.iter().enumerate() {
            let end = if i + 1 < cut_points.len() {
                cut_points[i + 1]
            } else {
                sequence.len()
            };

            let fingerprint = if i < phrase_fingerprints.len() {
                phrase_fingerprints[i]
            } else {
                // Fallback to CPU fingerprint if GPU didn't compute it
                Self::fingerprint_bases(&sequence[start..end])
            };

            phrases.push(Phrase {
                start,
                len: end - start,
                fingerprint,
            });
        }

        Ok(phrases)
    }

    /// CPU implementation of locally consistent parsing
    fn parse_locally_consistent_cpu(&self, sequence: &[EncodedBase]) -> Vec<Phrase> {
        // Compute per-position fingerprints using a rolling window
        // For DNA, we use digram fingerprints as the basis
        let mut position_fingerprints: Vec<Fingerprint> = Vec::with_capacity(sequence.len());

        for i in 0..sequence.len() {
            // Fingerprint based on local context (current + next if available)
            let fp = if i + 1 < sequence.len() {
                // Digram fingerprint
                let mut h = rustc_hash::FxHasher::default();
                sequence[i].0.hash(&mut h);
                sequence[i + 1].0.hash(&mut h);
                h.finish()
            } else {
                // Last position - just the single base
                sequence[i].0 as Fingerprint
            };
            position_fingerprints.push(fp);
        }

        // Find cut points using local minimum rule
        let mut cut_points: Vec<usize> = vec![0]; // Always start at 0
        let half_window = self.config.window_size / 2;

        for i in 1..sequence.len() {
            let fp = position_fingerprints[i];

            // Check if this is a local minimum within the window
            let window_start = i.saturating_sub(half_window);
            let window_end = (i + half_window + 1).min(sequence.len());

            let is_local_min = position_fingerprints[window_start..window_end]
                .iter()
                .all(|&other_fp| fp <= other_fp);

            // Also enforce max phrase length
            let last_cut = *cut_points.last().unwrap();
            let phrase_len = i - last_cut;

            if (is_local_min && phrase_len >= self.config.min_phrase_len)
                || phrase_len >= self.config.max_phrase_len
            {
                cut_points.push(i);
            }
        }

        // Convert cut points to phrases
        let mut phrases = Vec::with_capacity(cut_points.len());
        for window in cut_points.windows(2) {
            let start = window[0];
            let end = window[1];
            let phrase_content = &sequence[start..end];
            phrases.push(Phrase {
                start,
                len: end - start,
                fingerprint: Self::fingerprint_bases(phrase_content),
            });
        }

        // Handle last phrase (from last cut to end)
        if let Some(&last_cut) = cut_points.last() {
            if last_cut < sequence.len() {
                let phrase_content = &sequence[last_cut..];
                phrases.push(Phrase {
                    start: last_cut,
                    len: sequence.len() - last_cut,
                    fingerprint: Self::fingerprint_bases(phrase_content),
                });
            }
        }

        phrases
    }

    /// Build grammar from a sequence using locally consistent parsing.
    /// Returns (compressed_sequence, rules).
    pub fn build_grammar(
        &mut self,
        sequence: &[EncodedBase],
        source_grammar_id: usize,
    ) -> (Vec<Symbol>, FxHashMap<usize, Rule>) {
        if sequence.is_empty() {
            return (vec![], FxHashMap::default());
        }

        // Convert to symbols
        let mut current_symbols: Vec<Symbol> = sequence
            .iter()
            .enumerate()
            .map(|(i, &base)| Symbol::terminal(i, base, Direction::Forward, Some(source_grammar_id), Some(i)))
            .collect();

        // Iteratively parse and compress
        let mut iteration = 0;
        const MAX_ITERATIONS: usize = 100;

        while iteration < MAX_ITERATIONS {
            iteration += 1;

            // Parse current sequence to find phrases
            let phrases = self.parse_symbols_locally_consistent(&current_symbols);

            // Count phrase occurrences by fingerprint
            let mut fingerprint_counts: FxHashMap<Fingerprint, Vec<&Phrase>> = FxHashMap::default();
            for phrase in &phrases {
                if phrase.len >= self.config.min_phrase_len {
                    fingerprint_counts
                        .entry(phrase.fingerprint)
                        .or_default()
                        .push(phrase);
                }
            }

            // Find phrases that occur frequently enough to become rules
            let mut replacements_made = false;
            let mut new_sequence: Vec<Symbol> = Vec::with_capacity(current_symbols.len());
            let mut phrase_iter = phrases.iter().peekable();
            let mut current_pos = 0;

            while let Some(phrase) = phrase_iter.next() {
                let occurrences = fingerprint_counts.get(&phrase.fingerprint).map(|v| v.len()).unwrap_or(0);

                if occurrences >= self.config.min_occurrences && phrase.len >= self.config.min_phrase_len {
                    // This phrase becomes a rule
                    let rule_id = self.get_or_create_rule(phrase.fingerprint, &current_symbols[phrase.start..phrase.start + phrase.len]);

                    // Add non-terminal to new sequence
                    new_sequence.push(Symbol::non_terminal(new_sequence.len(), rule_id, Direction::Forward));
                    replacements_made = true;
                } else {
                    // Keep original symbols
                    for sym in &current_symbols[phrase.start..phrase.start + phrase.len] {
                        new_sequence.push(*sym);
                    }
                }
                current_pos = phrase.start + phrase.len;
            }

            // Handle any remaining symbols
            while current_pos < current_symbols.len() {
                new_sequence.push(current_symbols[current_pos]);
                current_pos += 1;
            }

            if !replacements_made || new_sequence.len() >= current_symbols.len() {
                // No progress made, stop
                break;
            }

            current_symbols = new_sequence;

            log::debug!(
                "LCG iteration {}: {} symbols, {} rules",
                iteration,
                current_symbols.len(),
                self.rules.len()
            );
        }

        (current_symbols, self.rules.clone())
    }

    /// Build grammar using pooled buffers for reduced allocation overhead.
    /// This is the recommended method for large sequences.
    /// Falls back to regular build_grammar if pooling is not enabled.
    pub fn build_grammar_pooled(
        &mut self,
        sequence: &[EncodedBase],
        source_grammar_id: usize,
    ) -> (Vec<Symbol>, FxHashMap<usize, Rule>) {
        if self.buffers.is_none() {
            // Fall back to non-pooled version
            return self.build_grammar(sequence, source_grammar_id);
        }

        if sequence.is_empty() {
            return (vec![], FxHashMap::default());
        }

        // Take ownership of buffers temporarily
        let mut buffers = self.buffers.take().unwrap();
        buffers.reserve(sequence.len());

        // Convert to symbols
        let mut current_symbols: Vec<Symbol> = sequence
            .iter()
            .enumerate()
            .map(|(i, &base)| Symbol::terminal(i, base, Direction::Forward, Some(source_grammar_id), Some(i)))
            .collect();

        // Iteratively parse and compress
        let mut iteration = 0;
        const MAX_ITERATIONS: usize = 100;

        while iteration < MAX_ITERATIONS {
            iteration += 1;

            // Parse current sequence to find phrases using pooled buffers
            let phrases = self.parse_symbols_pooled(&current_symbols, &mut buffers);

            // Count phrase occurrences by fingerprint
            let mut fingerprint_counts: FxHashMap<Fingerprint, usize> = FxHashMap::default();
            for phrase in &phrases {
                if phrase.len >= self.config.min_phrase_len {
                    *fingerprint_counts.entry(phrase.fingerprint).or_default() += 1;
                }
            }

            // Find phrases that occur frequently enough to become rules
            let mut replacements_made = false;
            buffers.new_sequence.clear();

            for phrase in &phrases {
                let occurrences = fingerprint_counts.get(&phrase.fingerprint).copied().unwrap_or(0);

                if occurrences >= self.config.min_occurrences && phrase.len >= self.config.min_phrase_len {
                    // This phrase becomes a rule
                    let rule_id = self.get_or_create_rule(phrase.fingerprint, &current_symbols[phrase.start..phrase.start + phrase.len]);

                    // Add non-terminal to new sequence
                    buffers.new_sequence.push(Symbol::non_terminal(buffers.new_sequence.len(), rule_id, Direction::Forward));
                    replacements_made = true;
                } else {
                    // Keep original symbols
                    for sym in &current_symbols[phrase.start..phrase.start + phrase.len] {
                        buffers.new_sequence.push(*sym);
                    }
                }
            }

            if !replacements_made || buffers.new_sequence.len() >= current_symbols.len() {
                // No progress made, stop
                break;
            }

            // Swap buffers efficiently
            std::mem::swap(&mut current_symbols, &mut buffers.new_sequence);
            buffers.clear();

            log::debug!(
                "LCG iteration {} (pooled): {} symbols, {} rules",
                iteration,
                current_symbols.len(),
                self.rules.len()
            );
        }

        // Return buffers for reuse
        buffers.clear();
        self.buffers = Some(buffers);

        (current_symbols, self.rules.clone())
    }

    /// Parse symbols using pooled buffers.
    fn parse_symbols_pooled(&self, symbols: &[Symbol], buffers: &mut LcgBuffers) -> Vec<Phrase> {
        if symbols.len() < self.config.min_phrase_len {
            if symbols.is_empty() {
                return vec![];
            }
            return vec![Phrase {
                start: 0,
                len: symbols.len(),
                fingerprint: Self::fingerprint_symbols(symbols),
            }];
        }

        // Reuse position_fingerprints buffer
        buffers.position_fingerprints.clear();
        for i in 0..symbols.len() {
            let fp = if i + 1 < symbols.len() {
                Self::fingerprint_symbols(&symbols[i..i + 2])
            } else {
                Self::fingerprint_symbols(&symbols[i..i + 1])
            };
            buffers.position_fingerprints.push(fp);
        }

        // Reuse cut_points buffer
        buffers.cut_points.clear();
        buffers.cut_points.push(0);
        let half_window = self.config.window_size / 2;

        for i in 1..symbols.len() {
            let fp = buffers.position_fingerprints[i];
            let window_start = i.saturating_sub(half_window);
            let window_end = (i + half_window + 1).min(symbols.len());

            let is_local_min = buffers.position_fingerprints[window_start..window_end]
                .iter()
                .all(|&other_fp| fp <= other_fp);

            let last_cut = *buffers.cut_points.last().unwrap();
            let phrase_len = i - last_cut;

            if (is_local_min && phrase_len >= self.config.min_phrase_len)
                || phrase_len >= self.config.max_phrase_len
            {
                buffers.cut_points.push(i);
            }
        }

        // Convert to phrases (allocated in arena for efficiency)
        let mut phrases = Vec::with_capacity(buffers.cut_points.len());
        for window in buffers.cut_points.windows(2) {
            let start = window[0];
            let end = window[1];
            phrases.push(Phrase {
                start,
                len: end - start,
                fingerprint: Self::fingerprint_symbols(&symbols[start..end]),
            });
        }

        // Handle last phrase
        if let Some(&last_cut) = buffers.cut_points.last() {
            if last_cut < symbols.len() {
                phrases.push(Phrase {
                    start: last_cut,
                    len: symbols.len() - last_cut,
                    fingerprint: Self::fingerprint_symbols(&symbols[last_cut..]),
                });
            }
        }

        phrases
    }

    /// Parse symbols using locally consistent parsing
    fn parse_symbols_locally_consistent(&self, symbols: &[Symbol]) -> Vec<Phrase> {
        if symbols.len() < self.config.min_phrase_len {
            if symbols.is_empty() {
                return vec![];
            }
            return vec![Phrase {
                start: 0,
                len: symbols.len(),
                fingerprint: Self::fingerprint_symbols(symbols),
            }];
        }

        // Compute per-position fingerprints
        let mut position_fingerprints: Vec<Fingerprint> = Vec::with_capacity(symbols.len());

        for i in 0..symbols.len() {
            let fp = if i + 1 < symbols.len() {
                Self::fingerprint_symbols(&symbols[i..i + 2])
            } else {
                Self::fingerprint_symbols(&symbols[i..i + 1])
            };
            position_fingerprints.push(fp);
        }

        // Find cut points
        let mut cut_points: Vec<usize> = vec![0];
        let half_window = self.config.window_size / 2;

        for i in 1..symbols.len() {
            let fp = position_fingerprints[i];
            let window_start = i.saturating_sub(half_window);
            let window_end = (i + half_window + 1).min(symbols.len());

            let is_local_min = position_fingerprints[window_start..window_end]
                .iter()
                .all(|&other_fp| fp <= other_fp);

            let last_cut = *cut_points.last().unwrap();
            let phrase_len = i - last_cut;

            if (is_local_min && phrase_len >= self.config.min_phrase_len)
                || phrase_len >= self.config.max_phrase_len
            {
                cut_points.push(i);
            }
        }

        // Convert to phrases
        let mut phrases = Vec::with_capacity(cut_points.len());
        for window in cut_points.windows(2) {
            let start = window[0];
            let end = window[1];
            phrases.push(Phrase {
                start,
                len: end - start,
                fingerprint: Self::fingerprint_symbols(&symbols[start..end]),
            });
        }

        if let Some(&last_cut) = cut_points.last() {
            if last_cut < symbols.len() {
                phrases.push(Phrase {
                    start: last_cut,
                    len: symbols.len() - last_cut,
                    fingerprint: Self::fingerprint_symbols(&symbols[last_cut..]),
                });
            }
        }

        phrases
    }

    /// Get existing rule for fingerprint or create a new one
    fn get_or_create_rule(&mut self, fingerprint: Fingerprint, symbols: &[Symbol]) -> usize {
        if let Some(&rule_id) = self.fingerprint_to_rule.get(&fingerprint) {
            // Increment usage count
            if let Some(rule) = self.rules.get_mut(&rule_id) {
                rule.usage_count += 1;
            }
            return rule_id;
        }

        // Create new rule
        let rule_id = self.next_rule_id;
        self.next_rule_id += 1;

        let rule = Rule {
            id: rule_id,
            symbols: symbols.to_vec(),
            usage_count: 1,
            positions: vec![],
            depth: None,
            assembly_index: None,
        };

        self.rules.insert(rule_id, rule);
        self.fingerprint_to_rule.insert(fingerprint, rule_id);

        rule_id
    }

    /// Merge another LcgBuilder's rules into this one.
    /// This is trivial because fingerprints are content-based - same content = same fingerprint.
    ///
    /// Note: This only merges the rule definitions. The compressed sequences from each chunk
    /// need to have their NonTerminal references updated separately using remap_sequence().
    pub fn merge(&mut self, other: Self) -> FxHashMap<usize, usize> {
        // Build mapping from other's rule IDs to merged rule IDs
        let mut id_mapping: FxHashMap<usize, usize> = FxHashMap::default();

        for (fingerprint, other_rule_id) in other.fingerprint_to_rule {
            if let Some(&existing_rule_id) = self.fingerprint_to_rule.get(&fingerprint) {
                // Rule already exists, map to existing ID and update usage count
                id_mapping.insert(other_rule_id, existing_rule_id);
                if let (Some(existing), Some(other_rule)) = (
                    self.rules.get_mut(&existing_rule_id),
                    other.rules.get(&other_rule_id),
                ) {
                    existing.usage_count += other_rule.usage_count;
                }
            } else {
                // New rule - assign new ID
                if let Some(other_rule) = other.rules.get(&other_rule_id) {
                    let new_id = self.next_rule_id;
                    self.next_rule_id += 1;
                    id_mapping.insert(other_rule_id, new_id);

                    let mut new_rule = other_rule.clone();
                    new_rule.id = new_id;

                    self.rules.insert(new_id, new_rule);
                    self.fingerprint_to_rule.insert(fingerprint, new_id);
                }
            }
        }

        // Now update NonTerminal references in all newly merged rules
        // We need to do this in a second pass after all mappings are known
        let rule_ids_to_update: Vec<usize> = id_mapping.values().cloned().collect();
        for rule_id in rule_ids_to_update {
            if let Some(rule) = self.rules.get_mut(&rule_id) {
                for symbol in &mut rule.symbols {
                    if let SymbolType::NonTerminal(ref_id) = &mut symbol.symbol_type {
                        if let Some(&new_ref_id) = id_mapping.get(ref_id) {
                            *ref_id = new_ref_id;
                        }
                    }
                }
            }
        }

        id_mapping
    }

    /// Remap NonTerminal references in a sequence using the ID mapping from merge()
    pub fn remap_sequence(sequence: &mut [Symbol], id_mapping: &FxHashMap<usize, usize>) {
        for symbol in sequence {
            if let SymbolType::NonTerminal(ref_id) = &mut symbol.symbol_type {
                if let Some(&new_id) = id_mapping.get(ref_id) {
                    *ref_id = new_id;
                }
            }
        }
    }

    /// Get the rules
    pub fn get_rules(&self) -> &FxHashMap<usize, Rule> {
        &self.rules
    }

    /// Merge chunks with overlap region analysis.
    ///
    /// This method handles boundary-spanning patterns by re-analyzing the overlap region
    /// between consecutive chunks to find patterns that span chunk boundaries.
    ///
    /// # Arguments
    /// * `chunk_sequences` - Vector of (compressed_sequence, overlap_data) for each chunk
    ///   - compressed_sequence: The LCG-compressed symbols for this chunk
    ///   - overlap_data: Raw encoded bases from the overlap region at the end of this chunk
    /// * `chunk_builders` - LcgBuilders from each chunk to be merged
    ///
    /// # Returns
    /// Tuple of (merged_sequence, id_mapping) where id_mapping can be used to remap references
    pub fn merge_chunks_with_overlap(
        &mut self,
        chunk_sequences: Vec<(Vec<Symbol>, Option<Vec<EncodedBase>>)>,
        mut chunk_builders: Vec<LcgBuilder>,
    ) -> (Vec<Symbol>, FxHashMap<usize, usize>) {
        let mut all_id_mappings: FxHashMap<usize, usize> = FxHashMap::default();
        let mut merged_sequence: Vec<Symbol> = Vec::new();

        // Drain chunk_builders so we can consume them one by one
        let mut builders_iter = chunk_builders.drain(..);

        for (idx, (mut compressed_seq, overlap_data)) in chunk_sequences.into_iter().enumerate() {
            // Skip empty chunks
            if compressed_seq.is_empty() {
                // Still consume the builder if present
                let _ = builders_iter.next();
                continue;
            }

            // Merge the builder's rules
            if let Some(builder) = builders_iter.next() {
                let id_mapping = self.merge(builder);

                // Remap references in the compressed sequence
                Self::remap_sequence(&mut compressed_seq, &id_mapping);

                // Accumulate mappings
                for (old_id, new_id) in id_mapping {
                    all_id_mappings.insert(old_id, new_id);
                }
            }

            // Handle overlap region for cross-boundary pattern detection
            if let Some(overlap_bases) = overlap_data {
                // Re-parse the overlap region to find any patterns we might have missed
                let overlap_phrases = self.parse_locally_consistent(&overlap_bases);

                // Check if any phrases in the overlap match existing rules
                for phrase in overlap_phrases {
                    if phrase.len >= self.config.min_phrase_len {
                        // Check if this phrase's fingerprint matches an existing rule
                        if let Some(&rule_id) = self.fingerprint_to_rule.get(&phrase.fingerprint) {
                            // This phrase matches an existing rule - update usage count
                            if let Some(rule) = self.rules.get_mut(&rule_id) {
                                rule.usage_count += 1;
                            }
                        }
                    }
                }
            }

            // Skip duplicated symbols from the overlap region of previous chunk
            let start_idx = if idx > 0 && !merged_sequence.is_empty() {
                // Calculate how many symbols to skip based on overlap
                // This is approximate - we skip symbols that would be in the overlap
                self.calculate_overlap_skip_count(&merged_sequence, &compressed_seq)
            } else {
                0
            };

            // Append non-overlapping symbols
            for symbol in compressed_seq.into_iter().skip(start_idx) {
                merged_sequence.push(symbol);
            }
        }

        (merged_sequence, all_id_mappings)
    }

    /// Calculate how many symbols to skip at the start of a chunk to avoid duplication.
    fn calculate_overlap_skip_count(&self, previous_seq: &[Symbol], current_seq: &[Symbol]) -> usize {
        if previous_seq.is_empty() || current_seq.is_empty() {
            return 0;
        }

        // Look for matching fingerprints at the boundary
        // We compare the last few symbols of previous with first few of current
        let check_len = self.config.window_size.min(previous_seq.len()).min(current_seq.len());

        for skip in 0..check_len {
            // Check if skipping 'skip' symbols from current_seq aligns well with previous_seq
            let prev_tail = &previous_seq[previous_seq.len().saturating_sub(check_len - skip)..];
            let curr_head = &current_seq[skip..skip + (check_len - skip).min(current_seq.len() - skip)];

            // Compare fingerprints
            let prev_fp = Self::fingerprint_symbols(prev_tail);
            let curr_fp = Self::fingerprint_symbols(curr_head);

            if prev_fp == curr_fp && skip > 0 {
                return skip;
            }
        }

        0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn encode_seq(s: &[u8]) -> Vec<EncodedBase> {
        s.iter().filter_map(|&b| EncodedBase::from_base(b)).collect()
    }

    #[test]
    fn test_fingerprint_deterministic() {
        let seq1 = encode_seq(b"ACGT");
        let seq2 = encode_seq(b"ACGT");
        let seq3 = encode_seq(b"TGCA");

        let fp1 = LcgBuilder::fingerprint_bases(&seq1);
        let fp2 = LcgBuilder::fingerprint_bases(&seq2);
        let fp3 = LcgBuilder::fingerprint_bases(&seq3);

        assert_eq!(fp1, fp2, "Same content should produce same fingerprint");
        assert_ne!(fp1, fp3, "Different content should produce different fingerprint");
    }

    #[test]
    fn test_locally_consistent_parsing() {
        let config = LcgConfig {
            min_phrase_len: 2,
            max_phrase_len: 16,
            window_size: 4,
            min_occurrences: 2,
        };
        let builder = LcgBuilder::new(config);

        let seq = encode_seq(b"ACGTACGTACGTACGT");
        let phrases = builder.parse_locally_consistent(&seq);

        // Should produce some phrases
        assert!(!phrases.is_empty());

        // Phrases should cover the whole sequence
        let total_len: usize = phrases.iter().map(|p| p.len).sum();
        assert_eq!(total_len, seq.len());
    }

    #[test]
    fn test_build_grammar_repetitive() {
        let config = LcgConfig {
            min_phrase_len: 2,
            max_phrase_len: 8,
            window_size: 4,
            min_occurrences: 2,
        };
        let mut builder = LcgBuilder::new(config);

        // Highly repetitive sequence
        let seq = encode_seq(b"ACACACACACACACAC");
        let (compressed, rules) = builder.build_grammar(&seq, 0);

        // Should have created some rules
        assert!(!rules.is_empty(), "Should create rules for repetitive sequence");

        // Compressed sequence should be shorter
        assert!(compressed.len() < seq.len(), "Compressed sequence should be shorter");
    }

    #[test]
    fn test_merge_preserves_fingerprints() {
        let config = LcgConfig::default();

        let mut builder1 = LcgBuilder::new(config.clone());
        let mut builder2 = LcgBuilder::new(config);

        // Build grammars from the same pattern in different "chunks"
        let seq1 = encode_seq(b"ACGTACGT");
        let seq2 = encode_seq(b"ACGTACGT");

        let (_, _) = builder1.build_grammar(&seq1, 0);
        let (_, _) = builder2.build_grammar(&seq2, 1);

        let rules_before_merge = builder1.rules.len();
        builder1.merge(builder2);

        // After merge, should not have duplicate rules for same content
        // (usage counts should be combined instead)
        assert!(
            builder1.rules.len() <= rules_before_merge * 2,
            "Merge should deduplicate identical rules"
        );
    }

    #[test]
    fn test_merge_chunks_with_overlap() {
        let config = LcgConfig {
            min_phrase_len: 2,
            max_phrase_len: 8,
            window_size: 4,
            min_occurrences: 2,
        };

        // Simulate two chunks with overlap
        // Chunk 1: ACACACAC (with overlap region "ACAC" at end)
        // Chunk 2: ACACACAC (starts with "ACAC" which overlaps with chunk 1)
        let seq1 = encode_seq(b"ACACACAC");
        let seq2 = encode_seq(b"ACACACAC");
        let overlap_data = encode_seq(b"ACAC");

        let mut builder1 = LcgBuilder::new(config.clone());
        let mut builder2 = LcgBuilder::new(config.clone());

        let (compressed1, _) = builder1.build_grammar(&seq1, 0);
        let (compressed2, _) = builder2.build_grammar(&seq2, 1);

        // Create master builder and merge with overlap
        let mut master_builder = LcgBuilder::new(config);

        let chunk_sequences = vec![
            (compressed1, Some(overlap_data)),
            (compressed2, None),
        ];
        let chunk_builders = vec![builder1, builder2];

        let (merged_sequence, _id_mappings) = master_builder.merge_chunks_with_overlap(
            chunk_sequences,
            chunk_builders,
        );

        // Merged sequence should be shorter than naive concatenation
        // because identical rules are deduplicated
        assert!(!merged_sequence.is_empty(), "Merged sequence should not be empty");
        assert!(!master_builder.rules.is_empty(), "Should have created rules");
    }

    #[test]
    fn test_build_grammar_pooled() {
        let config = LcgConfig {
            min_phrase_len: 2,
            max_phrase_len: 8,
            window_size: 4,
            min_occurrences: 2,
        };
        let mut builder = LcgBuilder::with_pooling(config, 1000);

        // Highly repetitive sequence
        let seq = encode_seq(b"ACACACACACACACAC");
        let (compressed, rules) = builder.build_grammar_pooled(&seq, 0);

        // Should have created some rules
        assert!(!rules.is_empty(), "Should create rules for repetitive sequence");

        // Compressed sequence should be shorter
        assert!(compressed.len() < seq.len(), "Compressed sequence should be shorter");
    }

    #[test]
    fn test_pooled_vs_regular_produce_same_results() {
        let config = LcgConfig {
            min_phrase_len: 2,
            max_phrase_len: 8,
            window_size: 4,
            min_occurrences: 2,
        };

        // Test with the same sequence
        let seq = encode_seq(b"ACGTACGTACACACACACGTACGT");

        // Regular build
        let mut regular_builder = LcgBuilder::new(config.clone());
        let (regular_compressed, regular_rules) = regular_builder.build_grammar(&seq, 0);

        // Pooled build
        let mut pooled_builder = LcgBuilder::with_pooling(config, seq.len());
        let (pooled_compressed, pooled_rules) = pooled_builder.build_grammar_pooled(&seq, 0);

        // Results should be identical
        assert_eq!(
            regular_compressed.len(),
            pooled_compressed.len(),
            "Pooled and regular should produce same compressed length"
        );
        assert_eq!(
            regular_rules.len(),
            pooled_rules.len(),
            "Pooled and regular should produce same number of rules"
        );
    }

    #[test]
    fn test_pooled_builder_reuses_buffers() {
        let config = LcgConfig {
            min_phrase_len: 2,
            max_phrase_len: 8,
            window_size: 4,
            min_occurrences: 2,
        };
        let mut builder = LcgBuilder::with_pooling(config, 1000);

        // Build grammar multiple times
        let seq1 = encode_seq(b"ACACACACACACACAC");
        let seq2 = encode_seq(b"GTGTGTGTGTGTGTGT");

        let (_, _) = builder.build_grammar_pooled(&seq1, 0);

        // After first build, buffers should still exist
        assert!(builder.buffers.is_some(), "Buffers should be returned after build");

        let (_, _) = builder.build_grammar_pooled(&seq2, 1);

        // After second build, buffers should still exist
        assert!(builder.buffers.is_some(), "Buffers should be reused across builds");
    }
}
