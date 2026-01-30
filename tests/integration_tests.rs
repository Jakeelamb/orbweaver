// Integration tests for the Orbweaver application

use std::fs;
use std::process::Command;
use tempfile::TempDir;

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

#[test]
fn test_tiny_fasta_to_json() {
    let temp_dir = TempDir::new().unwrap();

    // Create a simple FASTA file
    let fasta_path = temp_dir.path().join("test.fasta");
    fs::write(&fasta_path, ">test\nACGTACGTACGT\n").unwrap();

    // Run the application on the FASTA file from the temp directory
    // so the run output directory is created there
    let status = Command::new(get_orbweaver_bin())
        .args(&["-i", fasta_path.to_str().unwrap()])
        .current_dir(temp_dir.path())
        .status()
        .unwrap();

    assert!(status.success());

    // Find the run output directory
    let run_dirs: Vec<_> = fs::read_dir(temp_dir.path())
        .unwrap()
        .filter_map(|e| e.ok())
        .filter(|e| e.path().is_dir() && e.file_name().to_str().map_or(false, |s| !s.starts_with('.')))
        .filter(|e| e.file_name().to_str().map_or(false, |s| s != "test.fasta"))
        .collect();

    assert!(!run_dirs.is_empty(), "No run output directory created");
    let run_dir = run_dirs[0].path();

    // Check for grammar.json in the run directory
    let json_path = run_dir.join("grammar.json");
    assert!(json_path.exists(), "grammar.json not found in run directory");

    // Basic JSON validation
    let json_content = fs::read_to_string(&json_path).unwrap();
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

    // Run the application with analysis flags
    // The application creates output files in a run-specific directory
    let status = Command::new(get_orbweaver_bin())
        .args(&[
            "-i", fasta_path.to_str().unwrap(),
            "--kmer-size", "3",
            "--min-rule-usage", "2",
            "--stats"
        ])
        .current_dir(temp_dir.path())
        .status()
        .unwrap();

    // Verify the command succeeded
    assert!(status.success());

    // Find the run output directory (it's generated with a timestamp/UUID)
    let run_dirs: Vec<_> = fs::read_dir(temp_dir.path())
        .unwrap()
        .filter_map(|e| e.ok())
        .filter(|e| e.path().is_dir() && e.file_name().to_str().map_or(false, |s| !s.starts_with('.')))
        .collect();

    assert!(!run_dirs.is_empty(), "No run output directory created");
    let run_dir = run_dirs[0].path();

    // Check for expected files in the run directory
    let grammar_json_path = run_dir.join("grammar.json");
    let grammar_gfa_path = run_dir.join("grammar.gfa");
    let grammar_dot_path = run_dir.join("grammar.dot");

    assert!(grammar_json_path.exists(), "JSON output file not created in run directory");
    assert!(grammar_gfa_path.exists(), "GFA output file not created in run directory");
    assert!(grammar_dot_path.exists(), "DOT output file not created in run directory");

    // Read and check contents of output files
    let json_content = fs::read_to_string(&grammar_json_path).unwrap();
    assert!(json_content.contains("rules"), "JSON output doesn't contain rules");

    let dot_content = fs::read_to_string(&grammar_dot_path).unwrap();
    assert!(dot_content.contains("digraph"), "DOT file doesn't have correct format");

    let gfa_content = fs::read_to_string(&grammar_gfa_path).unwrap();
    assert!(gfa_content.contains("S\t"), "GFA file doesn't contain segment lines");
}

// TODO: Add more integration tests:
// - Test with different arguments (--skip-ns, --reverse-aware, kmer size)
// - Test GFA, Text, FASTA outputs
// - Test error handling (e.g., non-existent input file) 