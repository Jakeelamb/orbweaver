use anyhow::{Context, Result};
use clap::{Parser, Subcommand};
use orbweaver::encode::dna_2bit::EncodedBase;
use orbweaver::fasta::reader::FastaStream;
use orbweaver::grammar::engine::{Grammar, Sequitur};
use orbweaver::parallel::engine::{parallel_sequitur, ParallelMetrics};
use orbweaver::parallel::chunking::ChunkingConfig;
use orbweaver::utils::export::{export_grammar, calculate_compression_stats, print_compression_stats};
use orbweaver::utils::io::OutputFormat;
use orbweaver::utils::stats::SequenceStats;
use std::path::{Path, PathBuf};
use std::time::Instant;
use log::info;

#[derive(Parser)]
#[command(
    author = "Orbweaver Team",
    version = "0.1.0",
    about = "Grammar-based genomic sequence analysis",
    long_about = "Orbweaver is a tool for analyzing genomic sequences using context-free grammars."
)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Build a grammar from a FASTA file
    Build {
        /// Input FASTA file
        #[arg(short, long)]
        input: PathBuf,

        /// Output JSON file for the grammar
        #[arg(short, long)]
        output: PathBuf,

        /// Output text file for the grammar (human-readable)
        #[arg(long)]
        output_text: Option<PathBuf>,

        /// Output GFA file for the grammar (graph format)
        #[arg(long)]
        output_gfa: Option<PathBuf>,
        
        /// Export rule sequences as FASTA
        #[arg(long)]
        export_blocks: Option<PathBuf>,

        /// Minimum number of times a digram must appear to create a rule
        #[arg(short, long, default_value = "100")]
        min_rule_usage: usize,

        /// Consider reverse complements as equivalent
        #[arg(long, default_value = "true")]
        reverse_aware: bool,

        /// Use parallel processing
        #[arg(short, long, default_value = "false")]
        parallel: bool,

        /// Number of threads to use (default: number of CPU cores)
        #[arg(short, long)]
        threads: Option<usize>,

        /// Chunk size for parallel processing
        #[arg(long, default_value = "100000")]
        chunk_size: usize,

        /// Overlap size between chunks
        #[arg(long, default_value = "1000")]
        chunk_overlap: usize,
        
        /// Print statistics about the grammar. (Enabled by default)
        /// Use --no-stats to disable.
        #[arg(long, default_value_t = true)]
        stats: bool,
        
        /// Skip N bases in the input sequence
        #[arg(long, default_value = "true")]
        skip_ns: bool,
        
        /// Maximum k-mer size for statistics
        #[arg(short, long, default_value = "21")]
        kmer_size: usize,
    },
    
    /// Display information about a grammar file
    Info {
        /// Input grammar file (JSON)
        #[arg(short, long)]
        input: PathBuf,
        
        /// Calculate and display sequence statistics
        #[arg(long, default_value = "false")]
        sequence_stats: bool,
        
        /// Show compression statistics
        #[arg(long, default_value = "true")]
        compression: bool,
    },
    
    /// Convert a grammar file from one format to another
    Convert {
        /// Input grammar file
        #[arg(short, long)]
        input: PathBuf,
        
        /// Output file
        #[arg(short, long)]
        output: PathBuf,
        
        /// Input format (json, text, gfa, dot, fasta)
        #[arg(long)]
        input_format: Option<String>,
        
        /// Output format (json, text, gfa, dot, fasta)
        #[arg(long)]
        output_format: Option<String>,
    },
}

/// Main entry point for the CLI
pub fn main() -> Result<()> {
    // Initialize logging
    env_logger::init();
    
    // Parse command line arguments
    let cli = Cli::parse();
    
    match cli.command {
        Commands::Build {
            input,
            output,
            output_text,
            output_gfa,
            export_blocks,
            min_rule_usage,
            reverse_aware,
            parallel,
            threads,
            chunk_size,
            chunk_overlap,
            stats,
            skip_ns,
            kmer_size,
        } => {
            build_grammar(
                &input,
                &output,
                output_text.as_deref(),
                output_gfa.as_deref(),
                export_blocks.as_deref(),
                min_rule_usage,
                reverse_aware,
                parallel,
                threads,
                chunk_size,
                chunk_overlap,
                stats,
                skip_ns,
                kmer_size,
            )
        }
        Commands::Info { 
            input,
            sequence_stats,
            compression,
        } => {
            display_grammar_info(&input, sequence_stats, compression)
        }
        Commands::Convert {
            input,
            output,
            input_format,
            output_format,
        } => {
            convert_grammar(&input, &output, input_format.as_deref(), output_format.as_deref())
        }
    }
}

/// Build a grammar from a FASTA file and write it to output files
fn build_grammar(
    input_path: &Path,
    output_path: &Path,
    output_text_path: Option<&Path>,
    output_gfa_path: Option<&Path>,
    export_blocks_path: Option<&Path>,
    min_rule_usage: usize,
    reverse_aware: bool,
    parallel: bool,
    threads: Option<usize>,
    chunk_size: usize,
    chunk_overlap: usize,
    stats: bool,
    skip_ns: bool,
    kmer_size: usize,
) -> Result<()> {
    println!("Building grammar from {}", input_path.display());
    println!("  Min rule usage: {}", min_rule_usage);
    println!("  Reverse aware: {}", reverse_aware);
    println!("  Skip N bases: {}", skip_ns);
    
    let start = Instant::now();
    
    // Create a FASTA stream
    let mut stream = FastaStream::new(input_path)
        .with_context(|| format!("Failed to open FASTA file: {}", input_path.display()))?;
    
    if skip_ns {
        stream = stream.skip_ns();
    }
    
    // Collect the sequence
    println!("Reading sequence...");
    let read_start = Instant::now();
    let bases: Vec<EncodedBase> = stream.collect();
    let read_duration = read_start.elapsed();
    println!("Read {} bases in {:.2?}", bases.len(), read_duration);
    
    // Calculate sequence statistics if requested
    if stats {
        println!("\nSequence Statistics:");
        let seq_stats = SequenceStats::from_encoded_bases(&bases, Some(kmer_size));
        println!("{}", seq_stats);
    }
    
    // Process the sequence
    let grammar = if parallel {
        let num_threads = threads.unwrap_or_else(|| {
            std::thread::available_parallelism().map(|p| p.get()).unwrap_or(1)
        });
        
        println!("Using parallel processing with {} threads", num_threads);
        println!("  Chunk size: {}", chunk_size);
        println!("  Chunk overlap: {}", chunk_overlap);
        
        let chunking_config = ChunkingConfig {
            chunk_size,
            overlap_size: chunk_overlap,
            min_rule_usage,
            reverse_aware,
            num_threads,
            show_progress: true,
            adaptive_chunking: false,
            max_memory_per_chunk: None,
        };
        
        let (grammar, metrics) = parallel_sequitur(&bases, chunking_config)?;
        
        print_parallel_metrics(&metrics);
        
        grammar
    } else {
        println!("Using serial processing");
        
        let mut sequitur = Sequitur::new(min_rule_usage, reverse_aware);
        let process_start = Instant::now();
        let result = sequitur.build_grammar(&bases)?;
        let process_duration = process_start.elapsed();
        println!("Grammar built in {:.2?}", process_duration);
        
        result
    };
    
    // Export the grammar in various formats
    println!("Exporting grammar...");
    
    // JSON output is always required
    export_grammar(&grammar, output_path, OutputFormat::Json)?;
    println!("  JSON grammar written to {}", output_path.display());
    
    // Optional text output
    if let Some(text_path) = output_text_path {
        export_grammar(&grammar, text_path, OutputFormat::Text)?;
        println!("  Text grammar written to {}", text_path.display());
    }
    
    // Optional GFA output
    if let Some(gfa_path) = output_gfa_path {
        export_grammar(&grammar, gfa_path, OutputFormat::Gfa)?;
        println!("  GFA graph written to {}", gfa_path.display());
    }
    
    // Optional FASTA output for rule sequences
    if let Some(fasta_path) = export_blocks_path {
        export_grammar(&grammar, fasta_path, OutputFormat::Fasta)?;
        println!("  Rule sequences exported to {}", fasta_path.display());
    }
    
    // Print statistics
    if stats {
        println!("\nGrammar Statistics:");
        let compression_stats = calculate_compression_stats(&grammar);
        print_compression_stats(&compression_stats)?;
    } else {
        // Print minimal stats
        let elapsed = start.elapsed();
        println!("Finished in {:.2?}", elapsed);
        println!("Final sequence length: {}", grammar.sequence.len());
        println!("Number of rules: {}", grammar.rules.len());
        println!("Maximum rule depth: {}", grammar.max_depth);
    }
    
    Ok(())
}

/// Display information about a grammar file
fn display_grammar_info(
    grammar_path: &Path, 
    sequence_stats: bool,
    compression: bool,
) -> Result<()> {
    println!("Grammar information for {}", grammar_path.display());
    
    // Read the grammar file
    let file = std::fs::File::open(grammar_path)
        .with_context(|| format!("Failed to open grammar file: {}", grammar_path.display()))?;
    
    // Parse the JSON
    let grammar: Grammar = serde_json::from_reader(file)
        .with_context(|| format!("Failed to parse grammar file: {}", grammar_path.display()))?;
    
    // Display basic information
    println!("Final sequence length: {}", grammar.sequence.len());
    println!("Number of rules: {}", grammar.rules.len());
    println!("Maximum rule depth: {}", grammar.max_depth);
    
    // Calculate and display compression statistics if requested
    if compression {
        println!("\nCompression Statistics:");
        let compression_stats = calculate_compression_stats(&grammar);
        print_compression_stats(&compression_stats)?;
    }
    
    // Calculate and display sequence statistics if requested
    if sequence_stats {
        println!("\nSequence Statistics:");
        let expanded = orbweaver::utils::export::expand_sequence_to_string(&grammar.sequence, &grammar);
        let stats = SequenceStats::from_string(&expanded, None);
        println!("{}", stats);
    }
    
    Ok(())
}

/// Convert a grammar file from one format to another
fn convert_grammar(
    input_path: &Path,
    output_path: &Path,
    input_format_str: Option<&str>,
    output_format_str: Option<&str>,
) -> Result<()> {
    // Determine input format
    let input_format = match input_format_str {
        Some(format_str) => OutputFormat::from_str(format_str)
            .with_context(|| format!("Invalid input format: {}", format_str))?,
        None => orbweaver::utils::io::format_from_extension(input_path)
            .with_context(|| format!("Could not determine format from extension: {}", input_path.display()))?,
    };
    
    // Determine output format
    let output_format = match output_format_str {
        Some(format_str) => OutputFormat::from_str(format_str)
            .with_context(|| format!("Invalid output format: {}", format_str))?,
        None => orbweaver::utils::io::format_from_extension(output_path)
            .with_context(|| format!("Could not determine format from extension: {}", output_path.display()))?,
    };
    
    println!("Converting grammar from {} to {}", 
             input_path.display(), output_path.display());
    println!("  Input format: {:?}", input_format);
    println!("  Output format: {:?}", output_format);
    
    // Currently only support JSON as input format
    if input_format != OutputFormat::Json {
        anyhow::bail!("Currently only JSON is supported as input format");
    }
    
    // Read the grammar
    let file = std::fs::File::open(input_path)
        .with_context(|| format!("Failed to open grammar file: {}", input_path.display()))?;
    
    let grammar: Grammar = serde_json::from_reader(file)
        .with_context(|| format!("Failed to parse grammar file: {}", input_path.display()))?;
    
    // Export to the desired format
    export_grammar(&grammar, output_path, output_format)?;
    
    println!("Conversion successful!");
    
    Ok(())
}

/// Print performance metrics for parallel processing
fn print_parallel_metrics(metrics: &ParallelMetrics) {
    println!("Parallel processing metrics:");
    println!("  Number of chunks: {}", metrics.chunk_count);
    println!("  Chunk processing time: {:.2?}", metrics.chunk_processing_time);
    println!("  Rule deduplication time: {:.2?}", metrics.deduplication_time);
    println!("  Merge time: {:.2?}", metrics.merge_time);
    println!("  Total time: {:.2?}", metrics.total_time);
    println!("  Rules before merge: {}", metrics.rules_before_merge);
    println!("  Rules after merge: {}", metrics.rules_after_merge);
    println!("  Rule reduction: {:.2}%", 
        100.0 * (metrics.rules_before_merge - metrics.rules_after_merge) as f64 / metrics.rules_before_merge as f64);
} 