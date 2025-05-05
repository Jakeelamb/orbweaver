use crate::grammar::builder::GrammarBuilder;
use crate::grammar::rule::Rule;
use crate::grammar::symbol::Symbol;
use anyhow::{Context, Result};
use serde::Serialize;
use std::collections::HashMap;
use std::fs::File;
use std::io::BufWriter;
use std::path::Path;

/// Structure representing the final grammar for JSON serialization.
#[derive(Serialize)]
struct JsonGrammar<'a> {
    final_sequence: &'a Vec<Symbol>,
    rules: &'a HashMap<usize, Rule>,
}

/// Writes the generated grammar (final sequence and rules) to a JSON file.
///
/// Args:
///     grammar_builder: The GrammarBuilder instance after build_grammar() has run.
///     output_path: The path to the output JSON file.
pub fn write_grammar_json(grammar_builder: &GrammarBuilder, output_path: &Path) -> Result<()> {
    println!("Writing grammar to JSON: {}", output_path.display());

    let (final_sequence, rules) = grammar_builder.get_grammar();

    let json_grammar = JsonGrammar {
        final_sequence,
        rules,
    };

    let file = File::create(output_path)
        .with_context(|| format!("Failed to create JSON output file: {}", output_path.display()))?;
    let writer = BufWriter::new(file);

    // Use serde_json::to_writer_pretty for readable output
    serde_json::to_writer_pretty(writer, &json_grammar)
        .with_context(|| format!("Failed to serialize grammar to JSON: {}", output_path.display()))?;

    println!("Successfully wrote grammar to JSON.");
    Ok(())
} 