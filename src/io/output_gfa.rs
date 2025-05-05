use crate::grammar::builder::GrammarBuilder;
use crate::grammar::rule::Rule;
use crate::grammar::symbol::{Symbol, SymbolType};
use anyhow::{Context, Result};
use std::collections::HashMap;
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
        let sequence = reconstruct_rule_sequence(rule.id, rules)?;
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
            let (seg1_id, seg1_orient) = get_segment_id_and_orientation(sym1);
            let (seg2_id, seg2_orient) = get_segment_id_and_orientation(sym2);

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

/// Helper to get a segment ID and orientation string for a symbol.
/// Terminals could be represented as segments themselves (e.g., 'T_A', 'T_C'), 
/// or we could skip links involving terminals if only rules are segments.
/// For now, let terminals be represented by their base char.
fn get_segment_id_and_orientation(symbol: Symbol) -> (String, char) {
    match symbol.symbol_type {
        SymbolType::Terminal(base) => (format!("T_{}", base as char), symbol.strand), // e.g., T_A, T_C
        SymbolType::NonTerminal(rule_id) => (rule_id.to_string(), symbol.strand),
    }
}

/// Recursively reconstructs the base sequence of a rule.
/// Needs access to the full rule map.
fn reconstruct_rule_sequence(rule_id: usize, all_rules: &HashMap<usize, Rule>) -> Result<String> {
    let rule = all_rules.get(&rule_id)
        .ok_or_else(|| anyhow::anyhow!("Rule ID {} not found during sequence reconstruction", rule_id))?;
    
    let mut sequence = String::new();
    let mut visited = std::collections::HashSet::new(); // Detect cycles
    visited.insert(rule_id);

    fn append_symbol_sequence(
        symbol: Symbol,
        all_rules: &HashMap<usize, Rule>,
        current_sequence: &mut String,
        visited_stack: &mut std::collections::HashSet<usize>,
    ) -> Result<()> {
        match symbol.symbol_type {
            SymbolType::Terminal(base) => {
                current_sequence.push(base as char);
            }
            SymbolType::NonTerminal(sub_rule_id) => {
                if !visited_stack.insert(sub_rule_id) {
                    // Cycle detected
                    return Err(anyhow::anyhow!(
                        "Cycle detected during sequence reconstruction involving rule {}",
                        sub_rule_id
                    ));
                }
                let sub_rule = all_rules.get(&sub_rule_id).ok_or_else(|| {
                    anyhow::anyhow!(
                        "Rule ID {} (dependency of {}) not found during sequence reconstruction",
                        sub_rule_id, "current_rule_id" // Need context here
                    )
                })?;
                for &sub_symbol in &sub_rule.symbols {
                    append_symbol_sequence(sub_symbol, all_rules, current_sequence, visited_stack)?;
                }
                visited_stack.remove(&sub_rule_id); // Remove after successfully processing children
            }
        }
        Ok(())
    }

    for &symbol in &rule.symbols {
        // Pass a mutable reference to visited set for cycle detection within the recursive call.
        append_symbol_sequence(symbol, all_rules, &mut sequence, &mut visited)?;
    }

    // TODO: Handle strand properly during reconstruction?
    // If a symbol has strand '-', should its sequence be reverse complemented?
    // This adds significant complexity.

    Ok(sequence)
} 