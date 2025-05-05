use crate::grammar::builder::GrammarBuilder;
use crate::grammar::symbol::{Symbol, SymbolType};
use anyhow::{Context, Result};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

/// Formats a Symbol for text output.
fn format_symbol(symbol: Symbol) -> String {
    match symbol.symbol_type {
        SymbolType::Terminal(base) => format!("{}{}", base as char, symbol.strand),
        SymbolType::NonTerminal(rule_id) => format!("R{}{}", rule_id, symbol.strand),
    }
}

/// Writes the generated grammar in a human-readable text format.
///
/// Args:
///     grammar_builder: The GrammarBuilder instance after build_grammar() has run.
///     output_path: The path to the output text file.
pub fn write_grammar_text(grammar_builder: &GrammarBuilder, output_path: &Path) -> Result<()> {
    println!("Writing grammar to Text: {}", output_path.display());

    let (final_sequence, rules) = grammar_builder.get_grammar();

    let file = File::create(output_path)
        .with_context(|| format!("Failed to create text output file: {}", output_path.display()))?;
    let mut writer = BufWriter::new(file);

    // --- Write Final Sequence --- 
    writeln!(writer, "== Final Sequence ({} symbols) ==", final_sequence.len())
        .context("Failed to write text header")?;
    let seq_string: Vec<String> = final_sequence.iter().map(|&s| format_symbol(s)).collect();
    // Write sequence, potentially wrapping lines for readability
    let max_line_len = 80;
    for chunk in seq_string.join(" ").as_bytes().chunks(max_line_len) {
         writeln!(writer, "{}", std::str::from_utf8(chunk).unwrap_or("ERROR"))
             .context("Failed to write final sequence chunk")?;
    }
    writeln!(writer).context("Failed to write newline")?;

    // --- Write Rules --- 
    writeln!(writer, "== Rules ({} total) ==", rules.len())
        .context("Failed to write rules header")?;
    
    // Sort rules by ID for consistent output
    let mut sorted_rules: Vec<_> = rules.values().collect();
    sorted_rules.sort_by_key(|r| r.id);

    for rule in sorted_rules {
        let symbols_str: Vec<String> = rule.symbols.iter().map(|&s| format_symbol(s)).collect();
        writeln!(
            writer,
            "R{} [Usage={}] -> {}",
            rule.id,
            rule.usage_count,
            symbols_str.join(" ")
        )
        .with_context(|| format!("Failed to write rule R{}", rule.id))?;
    }

    println!("Successfully wrote grammar to Text.");
    Ok(())
} 