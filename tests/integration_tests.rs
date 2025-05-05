// Integration tests for the Orbweaver application

use anyhow::Result;
use assert_cmd::prelude::*; // Add methods on commands
use predicates::prelude::*; // Used for writing assertions
use std::fs;
use std::path::PathBuf;
use std::process::Command; // Run programs
use tempfile::TempDir;
use std::path::Path;

// Helper function to get the path to the compiled binary
fn get_orbweaver_bin() -> std::path::PathBuf {
    // Get the output directory where cargo puts the compiled binary
    let manifest_dir = std::env::var("CARGO_MANIFEST_DIR").unwrap();
    let target_dir = std::path::Path::new(&manifest_dir).join("target").join("debug");
    
    // Construct the path to the binary
    let orbweaver_path = target_dir.join("orbweaver");
    assert!(orbweaver_path.exists(), "orbweaver binary not found at: {:?}", orbweaver_path);
    orbweaver_path
}

// Helper to create a temporary test FASTA file
fn create_test_fasta(dir: &tempfile::TempDir, filename: &str, content: &str) -> Result<PathBuf> {
    let path = dir.path().join(filename);
    fs::write(&path, content)?;
    Ok(path)
}

#[test]
fn test_tiny_fasta_to_json() {
    let temp_dir = TempDir::new().unwrap();
    let output_path = temp_dir.path().join("output.json");
    
    // Create a simple FASTA file
    let fasta_path = temp_dir.path().join("test.fasta");
    fs::write(&fasta_path, ">test\nACGTACGTACGT\n").unwrap();
    
    // Run the application on the FASTA file
    let status = Command::new(get_orbweaver_bin())
        .args(&["-i", fasta_path.to_str().unwrap(), "-j", output_path.to_str().unwrap()])
        .status()
        .unwrap();
    
    assert!(status.success());
    assert!(Path::new(&output_path).exists());
    
    // Basic JSON validation
    let json_content = fs::read_to_string(&output_path).unwrap();
    assert!(json_content.starts_with("{"), "JSON should start with opening brace");
    assert!(json_content.contains("\"final_sequence\":"), "JSON should contain final_sequence key");
    assert!(json_content.contains("\"rules\":"), "JSON should contain rules key");
}

#[test]
fn test_comprehensive_functionality() {
    let temp_dir = TempDir::new().unwrap();
    
    // Create input FASTA with a sequence containing repeating patterns
    let fasta_path = temp_dir.path().join("test_repeats.fasta");
    fs::write(&fasta_path, ">repeats\nACGTACGTATATATGCGCGCGCACACACAC\n").unwrap();
    
    // Output file paths
    let json_path = temp_dir.path().join("grammar.json");
    let text_path = temp_dir.path().join("grammar.txt");
    let dot_path = temp_dir.path().join("graph.dot");
    let gfa_path = temp_dir.path().join("grammar.gfa");
    let fasta_blocks_path = temp_dir.path().join("blocks.fasta");
    
    // Run the application with multiple output formats and analysis
    let status = Command::new(get_orbweaver_bin())
        .args(&[
            "-i", fasta_path.to_str().unwrap(),
            "-j", json_path.to_str().unwrap(),
            "--output-text", text_path.to_str().unwrap(),
            "--visualize", dot_path.to_str().unwrap(),
            "--output-gfa", gfa_path.to_str().unwrap(),
            "--export-blocks", fasta_blocks_path.to_str().unwrap(),
            "--kmer-size", "3",
            "--min-rule-usage", "2",
            "--stats"
        ])
        .status()
        .unwrap();
    
    // Verify the command succeeded
    assert!(status.success());
    
    // Verify all output files exist
    assert!(Path::new(&json_path).exists(), "JSON output file not created");
    assert!(Path::new(&text_path).exists(), "Text output file not created");
    assert!(Path::new(&dot_path).exists(), "DOT output file not created");
    assert!(Path::new(&gfa_path).exists(), "GFA output file not created");
    assert!(Path::new(&fasta_blocks_path).exists(), "FASTA blocks output file not created");
    
    // Read and check contents of some output files
    let json_content = fs::read_to_string(&json_path).unwrap();
    assert!(json_content.contains("rules"), "JSON output doesn't contain rules");
    
    let text_content = fs::read_to_string(&text_path).unwrap();
    assert!(text_content.contains("Rule"), "Text output doesn't contain rules");
    
    let dot_content = fs::read_to_string(&dot_path).unwrap();
    assert!(dot_content.contains("digraph"), "DOT file doesn't have correct format");
    
    let gfa_content = fs::read_to_string(&gfa_path).unwrap();
    assert!(gfa_content.contains("S\t"), "GFA file doesn't contain segment lines");
    
    let fasta_blocks_content = fs::read_to_string(&fasta_blocks_path).unwrap();
    assert!(fasta_blocks_content.contains(">"), "FASTA blocks file doesn't contain sequences");
}

// TODO: Add more integration tests:
// - Test with different arguments (--skip-ns, --reverse-aware, kmer size)
// - Test GFA, Text, FASTA outputs
// - Test error handling (e.g., non-existent input file) 