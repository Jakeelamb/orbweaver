use anyhow::{Context, Result};
use crate::grammar::engine::Grammar;
use crate::grammar::symbol::{Symbol, SymbolType, Direction};
use crate::encode::dna_2bit::EncodedBase;
use std::io::{BufWriter, Write};
use std::path::Path;
use std::fs::File;
// use std::process::Command; // No longer needed if convert_dot_to_graphml is gone
// use crate::utils::predictable_random_color; // Commented out due to unresolved import
use quick_xml::events::{BytesEnd, BytesStart, Event};
use quick_xml::Writer as XmlWriter; // Renamed to avoid conflict with std::io::Write

/// Options for DOT graph visualization
#[derive(Debug, Clone)]
pub struct DotOptions {
    /// Whether to include terminal nodes in the graph
    pub include_terminals: bool,
    /// Maximum depth of rules to display in graph visualizations.
    pub max_depth: Option<usize>,
    /// Skip rules above this depth in graph visualizations.
    pub skip_rules_above_depth: Option<usize>,
    /// Use a transparent background for graph visualizations (PNG/SVG).
    pub transparent_background: bool,
    /// Use dark mode styling for graph visualizations.
    pub dark_mode: bool,
    /// Whether to include rule usage counts
    pub include_usage_counts: bool,
    /// Whether to color-code nodes by rule depth
    pub color_by_depth: bool,
    /// Graphviz layout engine (e.g., dot, sfdp).
    pub engine: String,
}

impl Default for DotOptions {
    fn default() -> Self {
        Self {
            include_terminals: true,
            max_depth: None,
            skip_rules_above_depth: None,
            transparent_background: false,
            dark_mode: false,
            include_usage_counts: true,
            color_by_depth: true,
            engine: "sfdp".to_string(),
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

/// Generate a GraphML representation of the grammar directly.
pub fn grammar_to_graphml_direct(grammar: &Grammar, options: &DotOptions) -> Result<String> {
    let mut buffer = Vec::new();
    let mut xml_writer = XmlWriter::new_with_indent(&mut buffer, b' ', 2);

    xml_writer.write_event(Event::Start(BytesStart::new("graphml")
        .with_attributes(vec![
            ("xmlns".as_bytes(), "http://graphml.graphdrawing.org/xmlns".as_bytes()),
            ("xmlns:xsi".as_bytes(), "http://www.w3.org/2001/XMLSchema-instance".as_bytes()),
            ("xsi:schemaLocation".as_bytes(), "http://graphml.graphdrawing.org/xmlns http://graphml.graphdrawing.org/xmlns/1.0/graphml.xsd".as_bytes())
        ])))?;

    // Define keys for node attributes
    let node_keys = vec![
        ("d_label", "label", "string"),      // e.g., R1, A+
        ("d_type", "type", "string"),       // "rule", "terminal"
        ("d_usage", "usage_count", "int"), // For rules
        ("d_depth", "depth", "int"),       // For rules
    ];
    for (id, name, type_) in node_keys {
        xml_writer.write_event(Event::Start(BytesStart::new("key")
            .with_attributes(vec![("id".as_bytes(), id.as_bytes()), ("for".as_bytes(), "node".as_bytes()), ("attr.name".as_bytes(), name.as_bytes()), ("attr.type".as_bytes(), type_.as_bytes())])))?;
        xml_writer.write_event(Event::End(BytesEnd::new("key")))?;
    }

    // Define keys for edge attributes (optional, e.g., relationship type)
    // let edge_keys = vec![("e_type", "relation", "string")]; // Example: "expands_to", "sequence_link"
    // for (id, name, type_) in edge_keys {
    //     xml_writer.write_event(Event::Start(BytesStart::new("key").with_attributes(vec![...])))?;
    //     xml_writer.write_event(Event::End(BytesEnd::new("key")))?;
    // }

    xml_writer.write_event(Event::Start(BytesStart::new("graph").with_attributes(vec![("id".as_bytes(), "G".as_bytes()), ("edgedefault".as_bytes(), "directed".as_bytes())])))?;

    // --- Create Nodes and Edges for Rules ---
    let mut defined_terminals = std::collections::HashSet::new();

    for (&rule_id, rule) in &grammar.rules {
        // Filter by depth if options specify
        if let Some(max_d) = options.max_depth {
            if rule.depth.unwrap_or(0) > max_d {
                continue;
            }
        }
        if let Some(min_d) = options.skip_rules_above_depth { // This logic is inverted, skip_rules_ABOVE_depth
            if rule.depth.unwrap_or(0) > min_d { // So if rule_depth > skip_depth, we skip. Correct.
                continue;
            }
        }

        let rule_node_id = format!("R{}", rule_id);
        xml_writer.write_event(Event::Start(BytesStart::new("node").with_attributes(vec![("id".as_bytes(), rule_node_id.as_bytes())])))?;
        xml_writer.write_event(Event::Start(BytesStart::new("data").with_attributes(vec![("key".as_bytes(), "d_label".as_bytes())])))?;
        xml_writer.write_event(Event::Text(quick_xml::events::BytesText::new(&rule_node_id)))?;
        xml_writer.write_event(Event::End(BytesEnd::new("data")))?;
        xml_writer.write_event(Event::Start(BytesStart::new("data").with_attributes(vec![("key".as_bytes(), "d_type".as_bytes())])))?;
        xml_writer.write_event(Event::Text(quick_xml::events::BytesText::new("rule")))?;
        xml_writer.write_event(Event::End(BytesEnd::new("data")))?;
        if options.include_usage_counts {
            xml_writer.write_event(Event::Start(BytesStart::new("data").with_attributes(vec![("key".as_bytes(), "d_usage".as_bytes())])))?;
            xml_writer.write_event(Event::Text(quick_xml::events::BytesText::new(&rule.usage_count.to_string())))?;
            xml_writer.write_event(Event::End(BytesEnd::new("data")))?;
        }
        if options.color_by_depth { // Using this option to decide if depth data is written
            xml_writer.write_event(Event::Start(BytesStart::new("data").with_attributes(vec![("key".as_bytes(), "d_depth".as_bytes())])))?;
            xml_writer.write_event(Event::Text(quick_xml::events::BytesText::new(&rule.depth.unwrap_or(0).to_string())))?;
            xml_writer.write_event(Event::End(BytesEnd::new("data")))?;
        }
        xml_writer.write_event(Event::End(BytesEnd::new("node")))?;

        // For each symbol in the rule, create an edge and potentially a terminal node
        for (symbol_idx, symbol) in rule.symbols.iter().enumerate() {
            match symbol.symbol_type {
                SymbolType::Terminal(base) => {
                    if options.include_terminals {
                        let terminal_label = format!("{}{}", base.to_char(), symbol.strand.to_char());
                        let terminal_node_id = format!("T_{}_{}{}", base.to_char(), symbol.strand.to_char(), symbol.id); // Ensure unique enough ID, include original symbol.id
                        
                        if !defined_terminals.contains(&terminal_node_id) {
                            xml_writer.write_event(Event::Start(BytesStart::new("node").with_attributes(vec![("id".as_bytes(), terminal_node_id.as_bytes())])))?;
                            xml_writer.write_event(Event::Start(BytesStart::new("data").with_attributes(vec![("key".as_bytes(), "d_label".as_bytes())])))?;
                            xml_writer.write_event(Event::Text(quick_xml::events::BytesText::new(&terminal_label)))?;
                            xml_writer.write_event(Event::End(BytesEnd::new("data")))?;
                            xml_writer.write_event(Event::Start(BytesStart::new("data").with_attributes(vec![("key".as_bytes(), "d_type".as_bytes())])))?;
                            xml_writer.write_event(Event::Text(quick_xml::events::BytesText::new("terminal")))?;
                            xml_writer.write_event(Event::End(BytesEnd::new("data")))?;
                            xml_writer.write_event(Event::End(BytesEnd::new("node")))?;
                            defined_terminals.insert(terminal_node_id.clone());
                        }
                        // Edge from rule to terminal
                        let edge_id = format!("E_R{}_to_T{}_{}", rule_id, terminal_node_id, symbol_idx);
                        xml_writer.write_event(Event::Start(BytesStart::new("edge")
                            .with_attributes(vec![
                                ("id".as_bytes(), edge_id.as_bytes()),
                                ("source".as_bytes(), rule_node_id.as_bytes()),
                                ("target".as_bytes(), terminal_node_id.as_bytes())
                            ])))?;
                        xml_writer.write_event(Event::End(BytesEnd::new("edge")))?;
                    }
                }
                SymbolType::NonTerminal(child_rule_id) => {
                    let child_rule_node_id = format!("R{}", child_rule_id);
                    // Edge from parent rule to child rule symbol instance
                    let edge_id = format!("E_R{}_to_R{}_{}", rule_id, child_rule_id, symbol_idx);
                     xml_writer.write_event(Event::Start(BytesStart::new("edge")
                        .with_attributes(vec![
                            ("id".as_bytes(), edge_id.as_bytes()),
                            ("source".as_bytes(), rule_node_id.as_bytes()),
                            ("target".as_bytes(), child_rule_node_id.as_bytes())
                        ])))?;
                    xml_writer.write_event(Event::End(BytesEnd::new("edge")))?;
                }
            }
        }
    }

    // --- Create Nodes and Edges for the Main Sequence ---
    if !grammar.sequence.is_empty() {
        for i in 0..grammar.sequence.len() {
            let symbol = &grammar.sequence[i];
            let current_symbol_node_id: String;

            match symbol.symbol_type {
                SymbolType::Terminal(base) => {
                    if options.include_terminals {
                        let terminal_label = format!("{}{}", base.to_char(), symbol.strand.to_char());
                        // ID needs to be consistent with how terminals from rules are made, using symbol.id for uniqueness.
                        let node_id = format!("T_{}_{}{}", base.to_char(), symbol.strand.to_char(), symbol.id);
                        current_symbol_node_id = node_id.clone();

                        if !defined_terminals.contains(&node_id) {
                            xml_writer.write_event(Event::Start(BytesStart::new("node").with_attributes(vec![("id".as_bytes(), node_id.as_bytes())])))?;
                            xml_writer.write_event(Event::Start(BytesStart::new("data").with_attributes(vec![("key".as_bytes(), "d_label".as_bytes())])))?;
                            xml_writer.write_event(Event::Text(quick_xml::events::BytesText::new(&terminal_label)))?;
                            xml_writer.write_event(Event::End(BytesEnd::new("data")))?;
                            xml_writer.write_event(Event::Start(BytesStart::new("data").with_attributes(vec![("key".as_bytes(), "d_type".as_bytes())])))?;
                            xml_writer.write_event(Event::Text(quick_xml::events::BytesText::new("terminal")))?;
                            xml_writer.write_event(Event::End(BytesEnd::new("data")))?;
                            xml_writer.write_event(Event::End(BytesEnd::new("node")))?;
                            defined_terminals.insert(node_id.clone());
                        }
                    } else {
                        // If terminals are not included, skip creating edges for them in sequence
                        continue; 
                    }
                }
                SymbolType::NonTerminal(rule_id) => {
                    // Rule nodes are expected to be defined already by the rule processing loop above.
                    // Filter by depth if options specify (consistency with rule node generation)
                    if let Some(rule) = grammar.rules.get(&rule_id) {
                        if let Some(max_d) = options.max_depth {
                            if rule.depth.unwrap_or(0) > max_d {
                                continue; // Skip this symbol in sequence if its rule is too deep
                            }
                        }
                        if let Some(min_d) = options.skip_rules_above_depth {
                            if rule.depth.unwrap_or(0) > min_d {
                                continue; // Skip this symbol if its rule is shallower than skip_depth
                            }
                        }
                    }
                    current_symbol_node_id = format!("R{}", rule_id);
                }
            }

            // Add edge from the previous sequence symbol to the current one
            if i > 0 {
                let prev_symbol = &grammar.sequence[i-1];
                let prev_symbol_node_id: String;
                
                // Determine previous symbol's node ID, applying same filtering as current
                match prev_symbol.symbol_type {
                    SymbolType::Terminal(base) => {
                        if !options.include_terminals { continue; }
                        prev_symbol_node_id = format!("T_{}_{}{}", base.to_char(), prev_symbol.strand.to_char(), prev_symbol.id);
                    }
                    SymbolType::NonTerminal(rule_id) => {
                        if let Some(rule) = grammar.rules.get(&rule_id) {
                            if let Some(max_d) = options.max_depth {
                                if rule.depth.unwrap_or(0) > max_d { continue; }
                            }
                            if let Some(min_d) = options.skip_rules_above_depth {
                                if rule.depth.unwrap_or(0) > min_d { continue; }
                            }
                        }
                        prev_symbol_node_id = format!("R{}", rule_id);
                    }
                }

                // Ensure current_symbol_node_id was actually set (not skipped by continue)
                // This check might be tricky if current_symbol_node_id was skipped due to its own filtering.
                // For simplicity, if either prev or current was filtered out by options, this edge won't be drawn.
                // A more robust way would be to build a list of valid sequence node IDs first.
                // However, given the current loop structure, we rely on both being present from their respective iterations.

                let edge_id = format!("S_edge_{}_to_{}", i-1, i);
                xml_writer.write_event(Event::Start(BytesStart::new("edge")
                    .with_attributes(vec![
                        ("id".as_bytes(), edge_id.as_bytes()),
                        ("source".as_bytes(), prev_symbol_node_id.as_bytes()),
                        ("target".as_bytes(), current_symbol_node_id.as_bytes())
                    ])))?;
                xml_writer.write_event(Event::End(BytesEnd::new("edge")))?;
            }
        }
    }

    xml_writer.write_event(Event::End(BytesEnd::new("graph")))?;
    xml_writer.write_event(Event::End(BytesEnd::new("graphml")))?;

    Ok(String::from_utf8(buffer)?)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::grammar::rule::Rule;
    // use crate::grammar::symbol::{Symbol, SymbolType, Direction}; // Already imported via super::*
    use crate::encode::dna_2bit::EncodedBase;
    use std::collections::HashMap;

    fn create_test_grammar() -> Grammar { // Ensure this test helper is robust
        let mut rules = HashMap::new();
        let rule0 = Rule {
            id: 0,
            symbols: vec![
                Symbol { id: 0, symbol_type: SymbolType::Terminal(EncodedBase::from_base(b'A').unwrap()), strand: Direction::Forward },
                Symbol { id: 1, symbol_type: SymbolType::Terminal(EncodedBase::from_base(b'C').unwrap()), strand: Direction::Forward },
            ],
            usage_count: 4,
            positions: Vec::new(), // Added default
            depth: Some(1),
        };
        let rule1 = Rule {
            id: 1,
            symbols: vec![
                Symbol { id: 2, symbol_type: SymbolType::NonTerminal(0), strand: Direction::Forward }, 
                Symbol { id: 3, symbol_type: SymbolType::Terminal(EncodedBase::from_base(b'G').unwrap()), strand: Direction::Forward },
            ],
            usage_count: 2,
            positions: Vec::new(), // Added default
            depth: Some(2),
        };
        rules.insert(0, rule0);
        rules.insert(1, rule1);
        Grammar {
            sequence: vec![
                Symbol { id: 10, symbol_type: SymbolType::NonTerminal(0), strand: Direction::Forward },
                Symbol { id: 11, symbol_type: SymbolType::NonTerminal(1), strand: Direction::Forward },
            ],
            rules,
            max_depth: 2, 
            origins: HashMap::new(), // Added default
        }
    }

    #[test]
    fn test_grammar_to_dot() {
        let grammar = create_test_grammar();
        let options = DotOptions {
            include_terminals: true,
            include_usage_counts: true,
            color_by_depth: true,
            ..Default::default()
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
        let g = create_test_grammar(); // Use the helper to get a Grammar instance
        let base_a = EncodedBase::from_base(b'A').expect("Failed to encode base A");
        let s = Symbol {
            id: 0, // Arbitrary ID for the test symbol itself
            symbol_type: SymbolType::Terminal(base_a),
            strand: Direction::Forward,
        };
        assert_eq!(crate::utils::export::format_symbol(&s, &g), "A+"); 
    }
} 