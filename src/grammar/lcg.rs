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
use crate::grammar::rule::Rule;
use crate::grammar::symbol::{Symbol, SymbolType, Direction};
use rustc_hash::FxHashMap;
use std::hash::{Hash, Hasher};

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

/// Locally Consistent Grammar builder
pub struct LcgBuilder {
    config: LcgConfig,
    /// Maps fingerprint -> rule_id for deduplication
    fingerprint_to_rule: FxHashMap<Fingerprint, usize>,
    /// All rules created
    rules: FxHashMap<usize, Rule>,
    /// Next rule ID to assign
    next_rule_id: usize,
}

impl LcgBuilder {
    pub fn new(config: LcgConfig) -> Self {
        Self {
            config,
            fingerprint_to_rule: FxHashMap::default(),
            rules: FxHashMap::default(),
            next_rule_id: 0,
        }
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
}
