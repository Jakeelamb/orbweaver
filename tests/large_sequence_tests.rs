// tests/large_sequence_tests.rs

use anyhow::Result;
use assert_cmd::Command;
use assert_cmd::assert::OutputAssertExt;
use predicates::prelude::*;
use std::fs;
use std::io::Write;
use tempfile::tempdir;
use std::path::PathBuf;

// Helper function to run orbweaver with given arguments on a fasta file
// Returns the command output and the path to the output JSON file.
fn run_orbweaver_on_fasta(fasta_content: &str, args: &[&str]) -> Result<(std::process::Output, PathBuf)> {
    let dir = tempdir()?;
    let input_path = dir.path().join("large_test.fasta");
    let output_json_file_path = dir.path().join("large_test_output.json");

    let mut file = fs::File::create(&input_path)?;
    writeln!(file, "{}", fasta_content)?;

    let mut cmd = Command::cargo_bin("orbweaver")?;
    cmd.arg("-i")
        .arg(&input_path)
        .arg("-j")
        .arg(&output_json_file_path); // Use the actual path for output
    
    for arg in args {
        cmd.arg(arg);
    }
    
    let output = cmd.ok()?; // Use ok() to allow checking status and output
    Ok((output, output_json_file_path))
}

// Helper function to generate a repetitive sequence
fn generate_repetitive_sequence(pattern: &str, repetitions: usize, header: &str) -> String {
    let mut sequence = String::new();
    sequence.push_str(&format!(">{}\n", header));
    for _ in 0..repetitions {
        sequence.push_str(pattern);
    }
    sequence.push('\n');
    sequence
}


#[test]
fn test_10mb_highly_repetitive_sequence_standard_mode() -> Result<()> {
    // Approx 1MB pattern, repeated 10 times = 10MB
    let pattern_unit = "ACGTN"; // 5 bytes
    let one_mb_repetitions = (1024 * 1024) / pattern_unit.len(); // ~200k repetitions for 1MB
    
    let mut pattern_1mb = String::new();
    for _ in 0..one_mb_repetitions {
        pattern_1mb.push_str(pattern_unit);
    }
    
    let fasta_content = generate_repetitive_sequence(&pattern_1mb, 10, "highly_repetitive_10mb");
    assert!(fasta_content.len() >= 10 * 1024 * 1024, "Generated FASTA content should be at least 10MB");

    println!("Generated 10MB repetitive sequence. Running Orbweaver...");
    let start_time = std::time::Instant::now();

    let (output, output_json_path) = run_orbweaver_on_fasta(&fasta_content, &["--stats"])?;
    
    let duration = start_time.elapsed();
    println!("Orbweaver processing took: {:?}", duration);

    output.assert().success().stdout(predicate::str::contains("Grammar construction completed"));
    
    // Check that output json was created and is not empty
    assert!(output_json_path.exists(), "Output JSON file should exist");
    let metadata = fs::metadata(&output_json_path)?;
    assert!(metadata.len() > 0, "Output JSON file should not be empty");
    
    // Further assertions could be:
    // - Check JSON output for expected number of rules (should be few for highly repetitive)
    // - Check compression ratio
    // - Check processing time against a baseline (if established)
    
    Ok(())
}

#[test]
fn test_10mb_highly_repetitive_sequence_streaming_mode() -> Result<()> {
    let pattern_unit = "GATTACA"; // 7 bytes
    let one_mb_repetitions = (1024 * 1024) / pattern_unit.len(); 
    
    let mut pattern_1mb = String::new();
    for _ in 0..one_mb_repetitions {
        pattern_1mb.push_str(pattern_unit);
    }

    let fasta_content = generate_repetitive_sequence(&pattern_1mb, 10, "highly_repetitive_10mb_stream");
    assert!(fasta_content.len() >= 10 * 1024 * 1024);

    println!("Generated 10MB repetitive sequence for streaming test. Running Orbweaver...");
    let start_time = std::time::Instant::now();
    
    // Use a specific chunk size for streaming, e.g., 1MB
    let (output, output_json_path) = run_orbweaver_on_fasta(&fasta_content, &["--streaming", "--chunk-size-streaming", "1048576", "--stats"])?;
    
    let duration = start_time.elapsed();
    println!("Orbweaver streaming processing took: {:?}", duration);

    output.assert().success().stdout(predicate::str::contains("Grammar construction completed"));
    
    assert!(output_json_path.exists(), "Output JSON file should exist for streaming test");
    let metadata = fs::metadata(&output_json_path)?;
    assert!(metadata.len() > 0, "Output JSON file should not be empty for streaming test");

    Ok(())
}

// TODO: Add test for a less repetitive 10MB sequence (e.g., from a generator script or real data subset)
// TODO: Add test for 10MB with GPU enabled and disabled, comparing time or results if applicable. 