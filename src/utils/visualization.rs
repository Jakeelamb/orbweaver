use anyhow::{Context, Result};
use crate::grammar::engine::Grammar;
use crate::grammar::symbol::{SymbolType, Direction};
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
    
    // Write segments (S lines) for each rule
    for (rule_id, rule) in &grammar.rules {
        // Expand the rule to its full DNA sequence
        let expanded_sequence = crate::utils::export::expand_rule_to_string(rule, &grammar.rules, Direction::Forward);
        output.push_str(&format!("S\t{}\t{}\tLN:i:{}\n", 
            rule_id, // Use rule ID directly as segment ID
            expanded_sequence, 
            expanded_sequence.len()
        ));
    }

    // Write segments for terminal symbols used in the final sequence?
    // This might make the graph too complex. Let's skip for now.

    // Write links (L lines) representing the final sequence structure
    for i in 0..(grammar.sequence.len().saturating_sub(1)) {
        let sym1 = &grammar.sequence[i];
        let sym2 = &grammar.sequence[i + 1];

        let seg1_id = match sym1.symbol_type {
            SymbolType::Terminal(b) => format!("T_{}{}", b.to_char(), sym1.strand.to_char()), // Placeholder for terminal segment
            SymbolType::NonTerminal(id) => id.to_string(),
        };
        let seg1_orient = sym1.strand.to_char();

        let seg2_id = match sym2.symbol_type {
            SymbolType::Terminal(b) => format!("T_{}{}", b.to_char(), sym2.strand.to_char()), // Placeholder for terminal segment
            SymbolType::NonTerminal(id) => id.to_string(),
        };
        let seg2_orient = sym2.strand.to_char();

        // If a symbol is terminal, we currently don't have an S line for it.
        // We only create S lines for rules. This linking strategy might be incomplete.
        // A better approach might involve Walk lines (W) or Containment lines (C),
        // but sticking to basic S/L for now.
        // We'll only write links between rules (NonTerminals).
        if let (SymbolType::NonTerminal(_), SymbolType::NonTerminal(_)) = 
            (&sym1.symbol_type, &sym2.symbol_type) {
            output.push_str(&format!("L\t{}\t{}\t{}\t{}\t0M\n", 
                seg1_id, seg1_orient, seg2_id, seg2_orient
            ));
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
    use crate::grammar::rule::Rule;
    use crate::grammar::symbol::{Symbol, SymbolType, Direction};
    use crate::encode::dna_2bit::EncodedBase;
    use std::collections::HashMap;

    fn create_test_grammar() -> Grammar {
        let mut rules = HashMap::new();
        let rule0 = Rule {
            id: 0,
            symbols: vec![
                Symbol::terminal(0, EncodedBase(0), Direction::Forward),
                Symbol::terminal(1, EncodedBase(1), Direction::Forward),
            ],
            usage_count: 4,
            positions: vec![0, 2, 4, 6],
            depth: Some(1),
        };
        let rule1 = Rule {
            id: 1,
            symbols: vec![
                Symbol::non_terminal(2, 0, Direction::Forward),
                Symbol::terminal(3, EncodedBase(2), Direction::Forward),
            ],
            usage_count: 2,
            positions: vec![1, 5],
            depth: Some(2),
        };
        rules.insert(0, rule0);
        rules.insert(1, rule1);
        Grammar {
            sequence: vec![
                Symbol::non_terminal(10, 0, Direction::Forward),
                Symbol::non_terminal(11, 1, Direction::Forward),
            ],
            rules,
            max_depth: 2,
            origins: HashMap::new(),
        }
    }

    #[test]
    fn test_grammar_to_dot() {
        let grammar = create_test_grammar();
        let options = DotOptions {
            include_terminals: true,
            include_usage_counts: true,
            color_by_depth: true,
        };
        let dot = grammar_to_dot(&grammar, &options).unwrap();
        assert!(dot.contains("digraph Grammar"));
        assert!(dot.contains("R0"));
        assert!(dot.contains("R1"));
    }

    #[test]
    fn test_grammar_to_gfa() {
        let grammar = create_test_grammar();
        let gfa = grammar_to_gfa(&grammar).unwrap();
        println!("Generated GFA:\n{}", gfa); // Print for debugging
        assert!(gfa.contains("H\tVN:Z:1.0"));
        // Check for S lines using rule IDs directly
        assert!(gfa.contains("S\t0\tAC\tLN:i:2"), "GFA missing S line for rule 0"); 
        assert!(gfa.contains("S\t1\tACG\tLN:i:3"), "GFA missing S line for rule 1");
        // Check for the L line connecting R0 and R1 in the final sequence
        assert!(gfa.contains("L\t0\t+\t1\t+\t0M"), "GFA missing L line for sequence R0 -> R1");
    }

    #[test]
    fn test_format_symbol_test() {
        let s = Symbol::terminal(0, EncodedBase(0), Direction::Forward);
        let mut rules = HashMap::new();
        let rule = Rule {
            id: 0,
            symbols: vec![s],
            usage_count: 1,
            positions: vec![0],
            depth: Some(1),
        };
        rules.insert(0, rule);
        let g = Grammar {
            sequence: vec![s],
            rules,
            max_depth: 1,
            origins: HashMap::new(),
        };
        let formatted = format_symbol(&s, &g);
        assert!(formatted.contains("A"));
    }
} 