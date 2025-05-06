use anyhow::{Context, Result};
use crate::grammar::engine::Grammar;
use crate::grammar::symbol::SymbolType;
use std::io::{BufWriter, Write};
use std::path::Path;
use std::fs::File;

/// Options for DOT graph visualization
#[derive(Debug, Clone)]
pub struct DotOptions {
    /// Whether to include terminal nodes in the graph
    pub include_terminals: bool,
    /// Whether to include rule usage counts
    pub include_usage_counts: bool,
    /// Whether to color-code nodes by rule depth
    pub color_by_depth: bool,
}

impl Default for DotOptions {
    fn default() -> Self {
        Self {
            include_terminals: true,
            include_usage_counts: true,
            color_by_depth: true,
        }
    }
}

/// Generate a DOT representation of the grammar suitable for Graphviz
pub fn grammar_to_dot(grammar: &Grammar, options: &DotOptions) -> Result<String> {
    let mut output = String::new();
    
    // Start the digraph
    output.push_str("digraph Grammar {\n");
    output.push_str("  graph [rankdir=LR];\n");
    output.push_str("  node [shape=box, style=filled, fontname=\"Arial\"];\n");
    
    // Generate a subgraph for the main sequence
    output.push_str("  subgraph cluster_sequence {\n");
    output.push_str("    label=\"Main Sequence\";\n");
    output.push_str("    style=filled;\n");
    output.push_str("    color=lightblue;\n");
    
    // Create nodes for sequence symbols
    for (i, symbol) in grammar.sequence.iter().enumerate() {
        match symbol.symbol_type {
            SymbolType::Terminal(base) => {
                if options.include_terminals {
                    let color = match base.0 {
                        0 => "lightgreen", // A
                        1 => "lightblue",  // C
                        2 => "orange",     // G
                        3 => "pink",       // T
                        _ => "gray",       // Others
                    };
                    
                    output.push_str(&format!(
                        "    seq_{} [label=\"{}\", color={}];\n",
                        i,
                        base.to_char(),
                        color
                    ));
                }
            }
            SymbolType::NonTerminal(rule_id) => {
                let label = if options.include_usage_counts {
                    if let Some(rule) = grammar.rules.get(&rule_id) {
                        format!("R{}({})", rule_id, rule.usage_count)
                    } else {
                        format!("R{}", rule_id)
                    }
                } else {
                    format!("R{}", rule_id)
                };
                
                let color = if options.color_by_depth {
                    // Color by rule depth
                    let depth = get_rule_depth(rule_id, &grammar.rules);
                    get_depth_color(depth, grammar.max_depth)
                } else {
                    "lightyellow"
                };
                
                output.push_str(&format!(
                    "    seq_{} [label=\"{}\", color={}];\n",
                    i,
                    label,
                    color
                ));
            }
        }
    }
    
    // Add edges between sequence symbols
    for i in 0..(grammar.sequence.len().saturating_sub(1)) {
        output.push_str(&format!("    seq_{} -> seq_{};\n", i, i + 1));
    }
    
    output.push_str("  }\n\n");
    
    // Create subgraphs for each rule
    for (&rule_id, rule) in &grammar.rules {
        output.push_str(&format!("  subgraph cluster_rule_{} {{\n", rule_id));
        output.push_str(&format!("    label=\"Rule {}\";\n", rule_id));
        output.push_str("    style=filled;\n");
        output.push_str("    color=lightgrey;\n");
        
        // Create nodes for rule symbols
        for (i, symbol) in rule.symbols.iter().enumerate() {
            match symbol.symbol_type {
                SymbolType::Terminal(base) => {
                    if options.include_terminals {
                        let color = match base.0 {
                            0 => "lightgreen", // A
                            1 => "lightblue",  // C
                            2 => "orange",     // G
                            3 => "pink",       // T
                            _ => "gray",       // Others
                        };
                        
                        output.push_str(&format!(
                            "    rule_{}_{} [label=\"{}\", color={}];\n",
                            rule_id,
                            i,
                            base.to_char(),
                            color
                        ));
                    }
                }
                SymbolType::NonTerminal(child_rule_id) => {
                    let label = if options.include_usage_counts {
                        if let Some(child_rule) = grammar.rules.get(&child_rule_id) {
                            format!("R{}({})", child_rule_id, child_rule.usage_count)
                        } else {
                            format!("R{}", child_rule_id)
                        }
                    } else {
                        format!("R{}", child_rule_id)
                    };
                    
                    let color = if options.color_by_depth {
                        // Color by rule depth
                        let depth = get_rule_depth(child_rule_id, &grammar.rules);
                        get_depth_color(depth, grammar.max_depth)
                    } else {
                        "lightyellow"
                    };
                    
                    output.push_str(&format!(
                        "    rule_{}_{} [label=\"{}\", color={}];\n",
                        rule_id,
                        i,
                        label,
                        color
                    ));
                }
            }
        }
        
        // Add edges between rule symbols
        for i in 0..(rule.symbols.len().saturating_sub(1)) {
            output.push_str(&format!(
                "    rule_{}_{} -> rule_{}_{};\n",
                rule_id, i, rule_id, i + 1
            ));
        }
        
        output.push_str("  }\n\n");
    }
    
    // Add edges between rules (for non-terminals)
    for (&rule_id, rule) in &grammar.rules {
        for (i, symbol) in rule.symbols.iter().enumerate() {
            match symbol.symbol_type {
                SymbolType::NonTerminal(child_rule_id) => {
                    if grammar.rules.contains_key(&child_rule_id) {
                        output.push_str(&format!(
                            "  rule_{}_{} -> rule_{}_0 [style=dashed, color=blue];\n",
                            rule_id, i, child_rule_id
                        ));
                    }
                }
                _ => {}
            }
        }
    }
    
    // Close the digraph
    output.push_str("}\n");
    
    Ok(output)
}

/// Helper function to get a rule's depth
fn get_rule_depth(rule_id: usize, rules: &std::collections::HashMap<usize, crate::grammar::rule::Rule>) -> usize {
    fn recurse(
        id: usize,
        rules: &std::collections::HashMap<usize, crate::grammar::rule::Rule>,
        visited: &mut std::collections::HashSet<usize>,
    ) -> usize {
        if visited.contains(&id) {
            return 0; // Prevent infinite recursion
        }
        
        visited.insert(id);
        
        if let Some(rule) = rules.get(&id) {
            let mut max_child_depth = 0;
            
            for symbol in &rule.symbols {
                match symbol.symbol_type {
                    SymbolType::NonTerminal(child_id) => {
                        let child_depth = recurse(child_id, rules, visited);
                        max_child_depth = max_child_depth.max(child_depth);
                    }
                    _ => {}
                }
            }
            
            // This rule's depth is 1 + the max depth of its children
            1 + max_child_depth
        } else {
            0
        }
    }
    
    let mut visited = std::collections::HashSet::new();
    recurse(rule_id, rules, &mut visited)
}

/// Get color based on depth
fn get_depth_color(depth: usize, max_depth: usize) -> &'static str {
    if max_depth == 0 {
        return "lightyellow";
    }
    
    // Normalize depth to a value between 0 and 1
    let normalized = depth as f64 / max_depth as f64;
    
    // Green (low depth) to Red (high depth)
    if normalized < 0.2 {
        "lightgreen"
    } else if normalized < 0.4 {
        "palegreen"
    } else if normalized < 0.6 {
        "khaki"
    } else if normalized < 0.8 {
        "lightsalmon"
    } else {
        "lightcoral"
    }
}

/// Write the DOT representation to a file
pub fn write_grammar_dot<P: AsRef<Path>>(
    path: P,
    grammar: &Grammar,
    options: &DotOptions,
) -> Result<()> {
    let dot_content = grammar_to_dot(grammar, options)?;
    let file = File::create(path).context("Failed to create DOT output file")?;
    let mut writer = BufWriter::new(file);
    writer.write_all(dot_content.as_bytes())?;
    Ok(())
}

/// Generate a Graphical Fragment Assembly (GFA) format representation
pub fn grammar_to_gfa(grammar: &Grammar) -> Result<String> {
    let mut output = String::new();
    
    // Header
    output.push_str("H\tVN:Z:1.0\n");
    
    // Map to store sequences for each segment
    let mut segments = std::collections::HashMap::new();
    
    // Generate segments for terminals
    for (i, symbol) in grammar.sequence.iter().enumerate() {
        if let SymbolType::Terminal(base) = symbol.symbol_type {
            let segment_id = format!("T{}", i);
            let segment_seq = format!("{}", base.to_char());
            segments.insert(segment_id, segment_seq);
        }
    }
    
    // Generate segments for rules
    for (&rule_id, rule) in &grammar.rules {
        let segment_id = format!("R{}", rule_id);
        let mut segment_seq = String::new();
        
        for symbol in &rule.symbols {
            if let SymbolType::Terminal(base) = symbol.symbol_type {
                segment_seq.push(base.to_char());
            } else {
                // For non-terminals in GFA, we'll represent them as 'N'
                segment_seq.push('N');
            }
        }
        
        segments.insert(segment_id, segment_seq);
    }
    
    // Write segments (S lines)
    for (id, seq) in &segments {
        output.push_str(&format!("S\t{}\t{}\n", id, seq));
    }
    
    // Write links for the main sequence (L lines)
    for i in 0..(grammar.sequence.len().saturating_sub(1)) {
        let source = match grammar.sequence[i].symbol_type {
            SymbolType::Terminal(_) => format!("T{}", i),
            SymbolType::NonTerminal(rule_id) => format!("R{}", rule_id),
        };
        
        let target = match grammar.sequence[i+1].symbol_type {
            SymbolType::Terminal(_) => format!("T{}", i+1),
            SymbolType::NonTerminal(rule_id) => format!("R{}", rule_id),
        };
        
        output.push_str(&format!("L\t{}\t+\t{}\t+\t0M\n", source, target));
    }
    
    // Write links for rules (L lines)
    for (&rule_id, rule) in &grammar.rules {
        for i in 0..(rule.symbols.len().saturating_sub(1)) {
            let source = match rule.symbols[i].symbol_type {
                SymbolType::Terminal(_) => continue, // Skip terminals in rules
                SymbolType::NonTerminal(source_rule_id) => format!("R{}", source_rule_id),
            };
            
            let target = match rule.symbols[i+1].symbol_type {
                SymbolType::Terminal(_) => continue, // Skip terminals in rules
                SymbolType::NonTerminal(target_rule_id) => format!("R{}", target_rule_id),
            };
            
            output.push_str(&format!("L\t{}\t+\t{}\t+\t0M\n", source, target));
        }
    }
    
    Ok(output)
}

/// Write the GFA representation to a file
pub fn write_grammar_gfa<P: AsRef<Path>>(
    path: P,
    grammar: &Grammar,
) -> Result<()> {
    let gfa_content = grammar_to_gfa(grammar)?;
    let file = File::create(path).context("Failed to create GFA output file")?;
    let mut writer = BufWriter::new(file);
    writer.write_all(gfa_content.as_bytes())?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::grammar::engine::Grammar;
    use crate::grammar::rule::Rule;
    use crate::grammar::symbol::{Symbol, Direction};
    use crate::encode::dna_2bit::EncodedBase;
    use crate::utils::export::format_symbol;
    use std::collections::HashMap;

    fn create_test_grammar() -> Grammar {
        // Create a simple grammar:
        // Rule 0: A+ C+
        // Rule 1: G+ T+
        // Final sequence: R0+ R1+ R0-
        
        let mut rules = HashMap::new();
        
        // Rule 0: A+ C+
        let rule0 = Rule {
            id: 0,
            symbols: vec![
                Symbol::terminal(0, EncodedBase(0), Direction::Forward), // A+
                Symbol::terminal(1, EncodedBase(1), Direction::Forward), // C+
            ],
            usage_count: 2,
            depth: Some(1),
            positions: vec![],
        };
        
        // Rule 1: G+ T+
        let rule1 = Rule {
            id: 1,
            symbols: vec![
                Symbol::terminal(2, EncodedBase(2), Direction::Forward), // G+
                Symbol::terminal(3, EncodedBase(3), Direction::Forward), // T+
            ],
            usage_count: 1,
            depth: Some(1),
            positions: vec![],
        };
        
        rules.insert(0, rule0);
        rules.insert(1, rule1);
        
        // Final sequence: R0+ R1+ R0-
        let sequence = vec![
            Symbol::non_terminal(10, 0, Direction::Forward), // R0+
            Symbol::non_terminal(11, 1, Direction::Forward), // R1+
            Symbol::non_terminal(12, 0, Direction::Reverse), // R0-
        ];
        
        Grammar {
            sequence,
            rules,
            max_depth: 1,
        }
    }
    
    #[test]
    fn test_grammar_to_dot() {
        let grammar = create_test_grammar();
        let options = DotOptions::default();
        
        let dot = grammar_to_dot(&grammar, &options).unwrap();
        
        // Check that the DOT output contains essential elements
        assert!(dot.contains("digraph Grammar"));
        assert!(dot.contains("R0")); // Rule 0 node
        assert!(dot.contains("R1")); // Rule 1 node
        // Add more checks if necessary, e.g., links
    }
    
    #[test]
    fn test_grammar_to_gfa() {
        let grammar = create_test_grammar();
        
        let gfa = grammar_to_gfa(&grammar).unwrap();
        
        // Check that the GFA output contains essential elements
        assert!(gfa.contains("VN:Z:1.0"));
        assert!(gfa.contains("S\tR0\tNN"), "Missing segment for R0"); // R0 has NN based on current impl
        assert!(gfa.contains("S\tR1\tNN"), "Missing segment for R1"); // R1 has NN based on current impl
        assert!(gfa.contains("S\tT_A\tA"), "Missing segment for Terminal A");
        assert!(gfa.contains("S\tT_C\tC"), "Missing segment for Terminal C");
        // ... check others

        // Check L lines (links)
        assert!(gfa.contains("L\tR0\t+\tR1\t+\t0M"), "Missing link R0->R1 in sequence");
        assert!(gfa.contains("L\tR1\t+\tR0\t-\t0M"), "Missing link R1->R0- in sequence");
    }
    
    #[test]
    fn test_format_symbol_test() { // Renamed from test_format_symbol to avoid conflict
        let grammar = create_test_grammar();
        
        // Test terminal symbols
        let term1 = Symbol::terminal(0, EncodedBase(0), Direction::Forward);
        assert_eq!(format_symbol(&term1, &grammar), "A+");
        
        let term2 = Symbol::terminal(1, EncodedBase(2), Direction::Reverse);
        assert_eq!(format_symbol(&term2, &grammar), "G-");
        
        // Test non-terminal symbols
        let nonterm1 = Symbol::non_terminal(2, 0, Direction::Forward);
        assert_eq!(format_symbol(&nonterm1, &grammar), "R0+");
        
        let nonterm2 = Symbol::non_terminal(3, 1, Direction::Reverse);
        assert_eq!(format_symbol(&nonterm2, &grammar), "R1-");
    }
} 