use crate::grammar::builder::GrammarBuilder;
use crate::grammar::symbol::{Symbol, SymbolType};
use anyhow::{Context, Result};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

/// Generates a unique node ID for DOT graph representation.
fn dot_node_id(symbol: Symbol) -> String {
    match symbol.symbol_type {
        // Include instance ID and strand in terminal nodes for uniqueness?
        // Maybe just base and strand is sufficient if we don't draw multiple instances.
        SymbolType::Terminal(base) => format!("T_{}{}", base as char, symbol.strand),
        SymbolType::NonTerminal(rule_id) => format!("R{}{}", rule_id, symbol.strand),
    }
}

/// Generates a label for a DOT graph node.
fn dot_node_label(symbol: Symbol) -> String {
    match symbol.symbol_type {
        SymbolType::Terminal(base) => format!("{}{}", base as char, symbol.strand),
        SymbolType::NonTerminal(rule_id) => format!("R{}{}", rule_id, symbol.strand),
    }
}

/// Writes the generated grammar structure to a DOT file for visualization.
///
/// Args:
///     grammar_builder: The GrammarBuilder instance after build_grammar() has run.
///     output_path: The path to the output DOT file.
pub fn write_grammar_dot(grammar_builder: &GrammarBuilder, output_path: &Path) -> Result<()> {
    println!("Writing grammar to DOT: {}", output_path.display());

    let (_final_sequence, rules) = grammar_builder.get_grammar();

    let file = File::create(output_path)
        .with_context(|| format!("Failed to create DOT output file: {}", output_path.display()))?;
    let mut writer = BufWriter::new(file);

    writeln!(writer, "digraph Grammar {{")
        .context("Failed to write DOT header")?;
    writeln!(writer, "  rankdir=LR;") // Left-to-right layout often better for sequences
        .context("Failed to write DOT rankdir")?;
    writeln!(writer, "  node [shape=box, style=filled, fillcolor=lightgrey];")
         .context("Failed to write DOT node style")?;

    // Keep track of declared nodes to avoid duplicates
    let mut declared_nodes = std::collections::HashSet::new();

    // Define Rule Nodes
    for rule in rules.values() {
        let rule_node_id = format!("R{}{}", rule.id, '+'); // Represent rule node as always '+'?
        let rule_label = format!("R{} (Usage={})", rule.id, rule.usage_count);
        if declared_nodes.insert(rule_node_id.clone()) {
             writeln!(writer, "  {} [label=\"{}\", shape=box, fillcolor=lightblue];", rule_node_id, rule_label)
                .with_context(|| format!("Failed to write DOT node for rule {}", rule.id))?;
        }
    }

    // Define Edges and potentially Terminal Nodes based on Rules
    for rule in rules.values() {
        let rule_node_id = format!("R{}{}", rule.id, '+'); // Match node definition
        
        for symbol in &rule.symbols {
            let symbol_node_id = dot_node_id(*symbol);
            let symbol_label = dot_node_label(*symbol);

            // Declare terminal nodes on demand
            if let SymbolType::Terminal(_) = symbol.symbol_type {
                 if declared_nodes.insert(symbol_node_id.clone()) {
                     writeln!(writer, "  {} [label=\"{}\", shape=ellipse, fillcolor=white];", symbol_node_id, symbol_label)
                         .context("Failed to write DOT node for terminal")?;
                 }
            } else {
                 // Ensure non-terminal nodes used in rules are declared if somehow missed
                 // (Should have been declared in the first loop, but for safety)
                 if !declared_nodes.contains(&symbol_node_id) {
                      let nt_rule_id = match symbol.symbol_type { SymbolType::NonTerminal(id) => id, _ => unreachable!() };
                      let nt_label = format!("R{} {}", nt_rule_id, symbol.strand);
                       if declared_nodes.insert(symbol_node_id.clone()) {
                            writeln!(writer, "  {} [label=\"{}\", shape=box, fillcolor=lightblue];", symbol_node_id, nt_label)
                                .context("Failed to write DOT node for non-terminal symbol")?;
                       }
                 }
            }

            // Write edge from rule to symbol
            writeln!(writer, "  {} -> {};", rule_node_id, symbol_node_id)
                .with_context(|| format!("Failed to write DOT edge for rule {}", rule.id))?;
        }
    }

    // Optionally, add edges representing the final sequence?
    // This could make the graph very large and linear.

    writeln!(writer, "}}")
        .context("Failed to write DOT footer")?;

    println!("Successfully wrote grammar to DOT.");
    Ok(())
} 