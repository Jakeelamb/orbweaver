use crate::grammar::builder::GrammarBuilder;
use crate::grammar::rule::Rule;
use crate::grammar::symbol::{Symbol, SymbolType, Direction};
use crate::grammar::engine::Grammar;
use crate::utils::export::expand_rule_to_string;
use anyhow::{Context, Result};
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

/// Writes the generated grammar to a GFAv1 file.
///
/// Args:
///     grammar_builder: The GrammarBuilder instance after build_grammar() has run.
///     output_path: The path to the output GFA file.
pub fn write_grammar_gfa(grammar_builder: &GrammarBuilder, output_path: &Path) -> Result<()> {
    println!("Writing grammar to GFAv1: {}", output_path.display());

    let (_final_sequence, rules) = grammar_builder.get_grammar();

    let file = File::create(output_path)
        .with_context(|| format!("Failed to create GFA output file: {}", output_path.display()))?;
    let mut writer = BufWriter::new(file);

    // --- Write GFA Header --- 
    writeln!(writer, "H\tVN:Z:1.0")
        .context("Failed to write GFA header")?;

    // --- Write S Lines (Segments for each Rule) --- 
    // This requires decoding the rule's symbols back into a base sequence.
    for rule in rules.values() {
        let segment_id = rule.id;
        let sequence = expand_rule_to_string(rule, &rules, Direction::Forward);
        // S <sid> <sequence> [slen] [tags...]
        // Using rule ID as segment ID. Sequence needs reconstruction.
        // slen tag is optional but recommended if sequence is long.
        writeln!(writer, "S\t{}\t{}\tLN:i:{}", segment_id, sequence, sequence.len())
             .with_context(|| format!("Failed to write S line for rule {}", rule.id))?;
    }

    // --- Write L Lines (Links representing rule structure R -> S1 S2) --- 
    for rule in rules.values() {
        let rule_id = rule.id;
        if rule.symbols.len() == 2 { // Assuming all rules replace digrams initially
            let sym1 = rule.symbols[0];
            let sym2 = rule.symbols[1];
            
            // Get segment IDs and orientations for the link
            let (seg1_id, seg1_orient) = get_segment_id_and_orientation(&sym1);
            let (seg2_id, seg2_orient) = get_segment_id_and_orientation(&sym2);

            // L <sid1> <orient1> <sid2> <orient2> <overlap> [tags...]
            // Link 1: Connects start of rule to first symbol
            // This isn't quite right for GFA representing grammar rules.
            // GFA usually links *ends* of segments.
            // How to represent R -> S1 S2?
            // Maybe a node for R, nodes for S1, S2 and edges?
            // Or just represent the final sequence structure? 

            // Let's rethink GFA representation. 
            // Maybe S lines for *terminals* and rules?
            // Or maybe S lines only for rules, and L lines show the R -> S1 S2 structure?
            // L <RuleID> <RuleOrient> <Symbol1ID> <Symbol1Orient> 0M (Overlap 0M)
            // L <Symbol1ID> <Symbol1Orient> <Symbol2ID> <Symbol2Orient> 0M (If S1 is a rule)
            
            // Simpler first approach: L line indicates S1 is followed by S2 within Rule R.
            // GFA viewers might not interpret this as intended for grammar structure.
            writeln!(writer, "L\t{}\t{}\t{}\t{}\t0M", seg1_id, seg1_orient, seg2_id, seg2_orient)
                 .with_context(|| format!("Failed to write L line for rule {}", rule_id))?;

        } else {
            // Handle rules with != 2 symbols if they exist later?
            eprintln!("Warning: Rule {} does not have exactly 2 symbols, skipping L line generation.", rule.id);
        }
    }

    // Alternative: Write L lines based on the *final sequence* connections?
    // For symbol[i] and symbol[i+1] in final_sequence:
    //   (id1, or1) = get_segment_id_and_orientation(symbol[i])
    //   (id2, or2) = get_segment_id_and_orientation(symbol[i+1])
    //   write L id1 or1 id2 or2 0M

    println!("Successfully wrote grammar to GFA (basic representation).");
    Ok(())
}

/// Determine the ID and orientation for a symbol in GFA
fn get_segment_id_and_orientation(symbol: &Symbol) -> (String, char) {
    match symbol.symbol_type {
        SymbolType::Terminal(base) => (
            format!("T_{}", base.to_char()), // Segment ID for terminal, e.g., T_A
            symbol.strand.to_char()       // Orientation
        ),
        SymbolType::NonTerminal(rule_id) => (
            format!("R{}", rule_id),        // Segment ID for non-terminal
            symbol.strand.to_char()       // Orientation
        )
    }
}

pub fn write_gfa_top_n_longest_rules(grammar: &Grammar, output_path: &Path, n: usize) -> Result<()> {
    if n == 0 {
        // Optionally, create an empty GFA or just print a message and return.
        // For now, just print and return.
        println!("Number of rules to keep (n) is 0. No GFA file will be generated for top N rules.");
        let file = File::create(output_path)
            .with_context(|| format!("Failed to create empty GFA output file: {}", output_path.display()))?;
        let mut writer = BufWriter::new(file);
        writeln!(writer, "H\tVN:Z:1.0").context("Failed to write GFA header for empty file")?;
        return Ok(());
    }

    println!(
        "Generating GFA for top {} longest rules (by expanded length) to: {}",
        n,
        output_path.display()
    );

    // 1. Calculate expanded lengths and identify top N rules
    let mut rule_lengths: Vec<(usize, usize, &Rule)> = Vec::new(); // (rule_id, length, rule_ref)
    for (id, rule) in &grammar.rules {
        // Use Direction::Forward for the canonical representation of the rule's sequence
        let expanded_seq = expand_rule_to_string(rule, &grammar.rules, Direction::Forward);
        rule_lengths.push((*id, expanded_seq.len(), rule));
    }

    rule_lengths.sort_by(|a, b| b.1.cmp(&a.1)); // Sort by length descending

    let top_n_rules_info: Vec<&Rule> = rule_lengths.iter().take(n).map(|&(_, _, rule_ref)| rule_ref).collect();

    if top_n_rules_info.is_empty() && !grammar.rules.is_empty() {
        println!("Top N rules list is empty, but grammar has rules. Check n ({}).", n);
    } else if grammar.rules.is_empty() {
        println!("No rules found in the grammar. GFA file will be effectively empty (header only).");
    }


    // 2. Collect all segments (rules and terminals) needed for S-lines
    let mut needed_rule_segments: HashMap<usize, &Rule> = HashMap::new();
    let mut needed_terminal_chars: HashSet<char> = HashSet::new(); // Store base char for terminals

    for top_rule in &top_n_rules_info {
        needed_rule_segments.insert(top_rule.id, top_rule); // The top rule itself
        for symbol_in_def in &top_rule.symbols {
            match symbol_in_def.symbol_type {
                SymbolType::Terminal(base) => {
                    needed_terminal_chars.insert(base.to_char());
                }
                SymbolType::NonTerminal(sub_rule_id) => {
                    if let Some(sub_rule) = grammar.rules.get(&sub_rule_id) {
                        needed_rule_segments.insert(sub_rule_id, sub_rule); // Direct sub-rules
                    } else {
                        // This case should ideally not happen if grammar is consistent
                        eprintln!("Warning: Sub-rule R{} referenced by R{} not found in grammar.rules.", sub_rule_id, top_rule.id);
                    }
                }
            }
        }
    }

    // 3. Write GFA file
    let file = File::create(output_path)
        .with_context(|| format!("Failed to create GFA output file: {}", output_path.display()))?;
    let mut writer = BufWriter::new(file);

    // GFA Header
    writeln!(writer, "H\tVN:Z:1.0").context("Failed to write GFA header")?;

    // S-Lines for Rules
    for (id, rule) in &needed_rule_segments {
        let expanded_seq = expand_rule_to_string(rule, &grammar.rules, Direction::Forward);
        writeln!(writer, "S\tR{}\t{}\tLN:i:{}", id, expanded_seq, expanded_seq.len())
            .with_context(|| format!("Failed to write S-line for rule R{}", id))?;
    }

    // S-Lines for Terminals
    for base_char in &needed_terminal_chars {
        let terminal_seg_id = format!("T_{}", base_char);
        // Sequence for terminal S-line is just the base itself
        writeln!(writer, "S\t{}\t{}\tLN:i:1", terminal_seg_id, base_char)
            .with_context(|| format!("Failed to write S-line for terminal T_{}", base_char))?;
    }

    // L-Lines (defining the internal structure of the top N rules)
    for top_rule in &top_n_rules_info {
        // Rules are R -> S1 S2 ... Sn. Create links S1->S2, S2->S3 ...
        if top_rule.symbols.len() >= 2 {
            for i in 0..(top_rule.symbols.len() - 1) {
                let sym1 = top_rule.symbols[i];
                let sym2 = top_rule.symbols[i+1];

                // Use the existing helper to get standardized segment names and orientations
                let (seg1_id_str, seg1_orient_char) = get_segment_id_and_orientation(&sym1);
                let (seg2_id_str, seg2_orient_char) = get_segment_id_and_orientation(&sym2);

                writeln!(
                    writer,
                    "L\t{}\t{}\t{}\t{}\t0M", // 0M for adjacency
                    seg1_id_str, seg1_orient_char, seg2_id_str, seg2_orient_char
                ).with_context(|| format!("Failed to write L-line for internal structure of rule R{}", top_rule.id))?;
            }
        }
    }

    println!(
        "Successfully wrote GFA for top {} longest rules to: {}",
        top_n_rules_info.len(), // Actual number written, might be less than n if grammar is small
        output_path.display()
    );
    Ok(())
} 