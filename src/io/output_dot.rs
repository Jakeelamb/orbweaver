use crate::grammar::builder::GrammarBuilder;
use crate::grammar::symbol::{Symbol, SymbolType};
use crate::grammar::engine::Grammar;
use anyhow::{Context, Result};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

/// Options for DOT graph visualization
pub struct DotOptions {
    /// Whether to include terminals in the graph
    pub include_terminals: bool,
    /// Whether to include rule usage counts
    pub include_usage_counts: bool,
    /// Whether to color nodes by rule depth
    pub color_by_depth: bool,
}

/// Generates a unique node ID for DOT graph representation.
fn dot_node_id(symbol: Symbol) -> String {
    match symbol.symbol_type {
        // Include instance ID and strand in terminal nodes for uniqueness?
        // Maybe just base and strand is sufficient if we don't draw multiple instances.
        SymbolType::Terminal(base) => format!("T_{}{}", base.to_char(), symbol.strand),
        SymbolType::NonTerminal(rule_id) => format!("R{}{}", rule_id, symbol.strand),
    }
}

/// Generates a label for a DOT graph node.
fn dot_node_label(symbol: Symbol) -> String {
    match symbol.symbol_type {
        SymbolType::Terminal(base) => format!("{}{}", base.to_char(), symbol.strand),
        SymbolType::NonTerminal(rule_id) => format!("R{}{}", rule_id, symbol.strand),
    }
}

/// Writes the generated grammar structure to a DOT file for visualization.
///
/// Args:
///     grammar_builder: The GrammarBuilder instance after build_grammar() has run.
///     output_path: The path to the output DOT file.
pub fn write_grammar_dot(grammar_builder: &GrammarBuilder, output_path: &Path, options: &DotOptions) -> Result<()> {
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

/// Write a grammar to a DOT file for visualization
pub fn write_grammar_dot_from_grammar(path: &Path, grammar: &Grammar, options: &DotOptions) -> Result<()> {
    let mut file = File::create(path)
        .context(format!("Failed to create DOT file: {}", path.display()))?;
    
    // Write DOT header
    writeln!(file, "digraph Grammar {{")?;
    writeln!(file, "  rankdir=LR;")?;
    writeln!(file, "  node [shape=box, style=filled, fillcolor=lightblue];")?;
    
    // Write nodes for each rule
    for (rule_id, rule) in &grammar.rules {
        let rule_depth = rule.depth.unwrap_or(0);
        let color = if options.color_by_depth {
            // Color by depth: deeper rules are darker
            let hue = 0.6; // Blue
            let saturation = 0.8;
            let lightness = (1.0 - (rule_depth as f32 * 0.1).min(0.8)).max(0.2);
            format!("\"#{:02x}{:02x}{:02x}\"", 
                (hue * 255.0) as u8, 
                (saturation * 255.0) as u8, 
                (lightness * 255.0) as u8)
        } else {
            "lightblue".to_string()
        };
        
        let label = if options.include_usage_counts {
            format!("R{} (used: {})", rule_id, rule.usage_count)
        } else {
            format!("R{}", rule_id)
        };
        
        writeln!(file, "  R{} [label=\"{}\", fillcolor={}];", rule_id, label, color)?;
    }
    
    // Write edges for rule references
    for (rule_id, rule) in &grammar.rules {
        let mut pos = 0;
        for symbol in &rule.symbols {
            match symbol.symbol_type {
                SymbolType::Terminal(base) if options.include_terminals => {
                    let base_char = match base.0 {
                        0 => 'A',
                        1 => 'C',
                        2 => 'G',
                        3 => 'T',
                        _ => 'N',
                    };
                    
                    let terminal_id = format!("R{}_T{}", rule_id, pos);
                    let strand_char = symbol.strand.to_char();
                    
                    writeln!(file, "  {} [label=\"{}{}\", shape=ellipse, fillcolor=lightgreen];", 
                        terminal_id, base_char, strand_char)?;
                    writeln!(file, "  R{} -> {};", rule_id, terminal_id)?;
                },
                SymbolType::NonTerminal(ref_rule_id) => {
                    let strand_char = symbol.strand.to_char();
                    writeln!(file, "  R{} -> R{} [label=\"{}\"];", rule_id, ref_rule_id, strand_char)?;
                },
                _ => {}
            }
            pos += 1;
        }
    }
    
    // Write DOT footer
    writeln!(file, "}}")?;
    
    Ok(())
} 