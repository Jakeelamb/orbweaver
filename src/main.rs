use anyhow::{Context, Result, bail};
use clap::{Parser, command};
use std::path::PathBuf;
use std::time::Instant;
use orbweaver::fasta::reader::read_fasta_sequences;
use orbweaver::encode::dna_2bit::{EncodedBase};
use orbweaver::grammar::builder::GrammarBuilder;
use orbweaver::grammar::engine::{Grammar, Sequitur};
use orbweaver::io::output_json::write_grammar_json;
use orbweaver::analysis::stats::calculate_and_print_stats;
use orbweaver::io::output_dot::DotOptions;
use orbweaver::parallel::chunking::ChunkingConfig;
use orbweaver::parallel::engine::parallel_sequitur;
use orbweaver::encode::bitvec;
use rayon::prelude::*;
use orbweaver::utils;

// Use the configuration module if needed, or args directly
// use orbweaver::config::Config; 

// Use other library modules as needed later
// use orbweaver::io; // Example
// use orbweaver::fasta;
// use orbweaver::grammar;
// use orbweaver::analysis; // If motifs analysis is re-integrated later

// Configure jemalloc for better memory management if feature is enabled
#[cfg(feature = "jemalloc")]
#[global_allocator]
static ALLOC: tikv_jemallocator::Jemalloc = tikv_jemallocator::Jemalloc;

/// Orbweaver: A grammar-based approach to genomic sequence analysis.
/// 
/// Orbweaver processes DNA sequences from FASTA files to build context-free grammars
/// by identifying repeating patterns. It can output the grammar in various formats and
/// provide statistics about the compression and structure.
#[derive(Parser, Debug)]
#[command(author = "Orbweaver Team", version, about, long_about = None)]
#[command(help_template = "\
{before-help}{name} {version}
{author-with-newline}{about-with-newline}
{usage-heading} {usage}

{all-args}{after-help}
")]
struct OrbweaverArgs {
    /// Input FASTA file path (.fa, .fasta, .fna).
    /// 
    /// Path to the file containing DNA sequences in FASTA format.
    /// Currently processes the first sequence in multi-sequence files.
    #[clap(short, long, value_parser, required = true)]
    input: PathBuf,

    // --- Output Format Options ---
    
    /// Output JSON file path for the grammar.
    /// 
    /// Writes a detailed JSON representation of the grammar including all
    /// rules and the final compressed sequence.
    #[clap(short = 'j', long, value_parser)]
    output_json: Option<PathBuf>,

    /// Output human-readable text representation of the grammar.
    /// 
    /// Writes a simplified text format showing rules and the final sequence
    /// in a more readable format than JSON.
    #[clap(long, value_parser)]
    output_text: Option<PathBuf>,

    /// Output GFAv1 representation of the grammar.
    /// 
    /// Exports the grammar as a graph in GFA (Graphical Fragment Assembly) format,
    /// which can be visualized with tools like Bandage.
    #[clap(long, value_parser)]
    output_gfa: Option<PathBuf>,

    /// Generate a .dot file for visualizing the grammar.
    /// 
    /// Creates a DOT file for visualization with Graphviz tools like dot, neato, etc.
    /// Example usage: dot -Tpng grammar.dot -o grammar.png
    #[clap(long, value_parser)]
    visualize: Option<PathBuf>,

    /// Export grammar rules as sequences in FASTA format.
    /// 
    /// Writes each grammar rule as a separate FASTA record, where each record
    /// contains the fully expanded DNA sequence for that rule.
    #[clap(long, value_parser)]
    export_blocks: Option<PathBuf>,

    /// Print statistics about the generated grammar.
    /// 
    /// Displays metrics about the grammar: rule count, depth, compression ratio, etc.
    #[clap(long, action = clap::ArgAction::SetTrue)]
    stats: bool,

    // --- Grammar Construction Options ---

    /// K-mer size for processing.
    /// 
    /// Used in some analysis algorithms. The grammar builder currently uses
    /// digrams (k=2) internally regardless of this setting.
    #[clap(short, long, value_parser, default_value_t = 21)]
    kmer_size: usize,

    /// Minimum usage count for a rule to be kept.
    /// 
    /// A digram must appear at least this many times to be replaced by a rule.
    /// Higher values lead to fewer rules with more usage each.
    #[clap(long, value_parser, default_value_t = 2)]
    min_rule_usage: usize,

    /// Maximum number of rules allowed (triggers eviction).
    /// 
    /// When the number of rules exceeds this limit, less frequently used rules
    /// are evicted and inlined. If not specified, no limit is enforced.
    #[clap(long, value_parser)]
    max_rule_count: Option<usize>,

    /// Perform reverse complement canonicalization.
    /// 
    /// When enabled, a digram and its reverse complement are treated as equivalent.
    /// For example, "AC" and "GT" (reverse complement) are considered the same pattern.
    #[clap(long, default_value_t = true, action = clap::ArgAction::Set)]
    reverse_aware: bool,

    // --- Input Processing Options ---

    /// Skip N bases in FASTA input.
    /// 
    /// When enabled, 'N' or 'n' characters in the input are skipped,
    /// as they typically represent unknown bases.
    #[clap(long, default_value_t = true, action = clap::ArgAction::Set)]
    skip_ns: bool,

    /// Process genome in chunks of this size.
    /// 
    /// Divides large sequences into smaller chunks for processing.
    /// Useful for very large genomes that exceed available memory.
    #[clap(long, value_parser)]
    chunk_size: Option<usize>,

    /// Overlap between chunks when chunking is enabled.
    /// 
    /// Specifies the number of bases that overlap between adjacent chunks.
    /// Helps ensure patterns spanning chunk boundaries are detected.
    /// Note: Only relevant when chunking is enabled.
    #[clap(long, value_parser, default_value_t = 1000)]
    chunk_overlap: usize,
    
    /// Number of threads to use for parallel processing.
    /// 
    /// Sets the number of worker threads for parallelized operations.
    /// Default: number of logical CPU cores.
    #[clap(long, value_parser)]
    threads: Option<usize>,

    /// Enable performance profiling.
    /// 
    /// Generates a CPU profile and flame graph of the application execution.
    /// Results are saved to the 'profile' directory.
    #[clap(long, action = clap::ArgAction::SetTrue)]
    profile: bool,

    /// Enable 2-bit encoding to reduce memory usage.
    /// 
    /// Use 2-bit encoding for DNA sequences (A=00, C=01, G=10, T=11),
    /// reducing memory usage by up to 75% compared to ASCII storage.
    #[clap(long, default_value_t = true, action = clap::ArgAction::Set)]
    use_encoding: bool,
    
    /// Enable streaming mode for low memory usage.
    /// 
    /// Process the FASTA file in a streaming fashion, without loading
    /// the entire sequence into memory at once. Recommended for large genomes.
    #[clap(long, action = clap::ArgAction::SetTrue)]
    streaming: bool,
}

fn main() -> Result<()> {
    let args = OrbweaverArgs::parse();
    
    // Start profiling if requested
    #[cfg(feature = "profiling")]
    let guard = if args.profile {
        println!("Starting profiling");
        let guard = pprof::ProfilerGuard::new(100)?;
        Some(guard)
    } else {
        None
    };

    println!("Orbweaver: Processing file {}", args.input.display());
    
    let start = Instant::now();

    // Determine execution mode: parallel chunking or standard
    let grammar = if let Some(chunk_size) = args.chunk_size {
        // Process in parallel using chunking
        process_chunked_mode(&args, chunk_size)?
    } else {
        // Process in standard mode (load entire sequence)
        process_standard_mode(&args)?
    };
    
    println!("Grammar construction completed in {:?}", start.elapsed());
    
    // Write outputs in requested formats
    if let Some(ref path) = args.output_json {
        write_grammar_json(path, &grammar.sequence, &grammar.rules)?;
        println!("Wrote JSON grammar to {}", path.display());
    }
    
    if let Some(ref path) = args.output_text {
        orbweaver::utils::export::export_grammar_text(&grammar, &mut utils::io::open_file_for_writing(path)?)?;
        println!("Wrote text grammar to {}", path.display());
    }
    
    if let Some(ref path) = args.output_gfa {
        orbweaver::utils::export::export_grammar_gfa(&grammar, &mut utils::io::open_file_for_writing(path)?)?;
        println!("Wrote GFA to {}", path.display());
    }
    
    if let Some(ref path) = args.visualize {
        let options = DotOptions {
            include_terminals: true,
            include_usage_counts: true,
            color_by_depth: true,
        };
        
        orbweaver::utils::export::export_grammar_dot(&grammar, &mut utils::io::open_file_for_writing(path)?)?;
        println!("Wrote DOT visualization to {}", path.display());
    }
    
    if let Some(ref path) = args.export_blocks {
        orbweaver::utils::export::export_grammar_fasta(&grammar, &mut utils::io::open_file_for_writing(path)?)?;
        println!("Wrote FASTA blocks to {}", path.display());
    }
    
    // Print statistics if requested
    if args.stats {
        calculate_and_print_stats(&grammar)?;
    }
    
    #[cfg(feature = "profiling")]
    if let Some(guard) = guard {
        if let Ok(report) = guard.report().build() {
            let file = File::create("flamegraph.svg").unwrap();
            report.flamegraph(file).unwrap();
            println!("Wrote profiling data to flamegraph.svg");
        }
    }
    
    Ok(())
}

fn process_standard_mode(args: &OrbweaverArgs) -> Result<Grammar> {
    println!("Using standard mode (loading entire sequence)");
    
    // Read the FASTA file
    let sequences = read_fasta_sequences(&args.input, args.skip_ns)
        .context("Failed to read FASTA file")?;
    
    if sequences.is_empty() {
        bail!("No sequences found in input file");
    }
    
    let (_, bases) = &sequences[0];
    println!("Processing sequence: {} bases", bases.len());
    
    // Process with 2-bit encoding if enabled
    if args.use_encoding {
        let encoded_bases = encode_dna(bases);
        println!("Using 2-bit encoding (converted {} bytes to {} encoded bases)", 
            bases.len(), encoded_bases.len());
        
        // Build grammar using sequitur algorithm
        let mut grammar_builder = GrammarBuilder::new(args.min_rule_usage, args.reverse_aware);
        
        // Set rule count limit if specified
        if let Some(max_rules) = args.max_rule_count {
            grammar_builder = grammar_builder.with_max_rules(max_rules);
        }
        
        grammar_builder.build_grammar(&encoded_bases)?;
        
        // Get the resulting grammar
        let (sequence, rules) = grammar_builder.get_grammar();
        let max_depth = grammar_builder.get_max_rule_depth();
        
        // Calculate memory savings from 2-bit encoding
        let (saving_pct, description) = bitvec::estimate_memory_savings(bases.len(), encoded_bases.len());
        println!("Memory savings from 2-bit encoding: {:.1}% ({})", saving_pct, description);
        
        Ok(Grammar {
            sequence: sequence.clone(),
            rules: rules.clone(),
            max_depth,
        })
    } else {
        // Process without encoding (less memory efficient)
        println!("WARNING: Processing without 2-bit encoding (higher memory usage)");
        // Convert Vec<u8> to Vec<EncodedBase> for Sequitur
        let encoded_bases_fallback = encode_dna(bases);
        let mut sequitur = Sequitur::new(args.min_rule_usage, args.reverse_aware);
        sequitur.build_grammar(&encoded_bases_fallback)
    }
}

fn process_chunked_mode(args: &OrbweaverArgs, chunk_size: usize) -> Result<Grammar> {
    println!("Using parallel chunking mode with chunk size: {}", chunk_size);
    
    // Configure the chunking parameters
    let config = ChunkingConfig {
        chunk_size,
        overlap_size: args.chunk_overlap,
        min_rule_usage: args.min_rule_usage,
        reverse_aware: args.reverse_aware,
        num_threads: args.threads.unwrap_or_else(num_cpus::get),
        show_progress: true,
    };
    
    // Read and process the sequence
    if args.use_encoding {
        // Read the sequence with encoding
        let sequences = read_fasta_sequences(&args.input, args.skip_ns)?;
        if sequences.is_empty() {
            bail!("No sequences found in input file");
        }
        
        let (_, bases) = &sequences[0];
        println!("Processing sequence: {} bases", bases.len());
        
        // Convert to EncodedBase for more efficient processing
        let encoded_bases = encode_dna(bases);
        
        // Run the parallel sequitur algorithm
        let (grammar, metrics) = parallel_sequitur(&encoded_bases, config)?;
        
        // Display execution metrics
        println!("Parallel execution metrics:");
        println!("  Chunk count: {}", metrics.chunk_count);
        println!("  Chunk processing time: {:?}", metrics.chunk_processing_time);
        println!("  Merge time: {:?}", metrics.merge_time);
        println!("  Deduplication time: {:?}", metrics.deduplication_time);
        println!("  Total time: {:?}", metrics.total_time);
        println!("  Rule count (before merge): {}", metrics.rules_before_merge);
        println!("  Rule count (after merge): {}", metrics.rules_after_merge);
        
        // Calculate memory savings from 2-bit encoding
        let (saving_pct, description) = bitvec::estimate_memory_savings(bases.len(), encoded_bases.len());
        println!("Memory savings from 2-bit encoding: {:.1}% ({})", saving_pct, description);
        
        Ok(grammar)
    } else {
        bail!("Chunked mode requires use_encoding to be enabled for memory efficiency.")
    }
}

// Add this helper function for encoding DNA
fn encode_dna(bases: &[u8]) -> Vec<EncodedBase> {
    bases.iter()
        .filter_map(|&b| EncodedBase::from_base(b))
        .collect()
}
