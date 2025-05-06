use anyhow::{Context, Result, bail};
use clap::{Parser, command};
use std::path::PathBuf;
use std::time::Instant;
use orbweaver::fasta::reader::{read_fasta_sequences, read_fasta_sequences_async};
use orbweaver::encode::dna_2bit::{EncodedBase};
use orbweaver::grammar::builder::GrammarBuilder;
use orbweaver::grammar::engine::{Grammar};
use orbweaver::io::output_json::write_grammar_json;
use orbweaver::analysis::stats::calculate_and_print_stats;
use orbweaver::utils::visualization::DotOptions;
use orbweaver::parallel::chunking::ChunkingConfig;
use orbweaver::parallel::engine::{parallel_sequitur, merge_grammars};
use orbweaver::encode::bitvec;
use rayon::prelude::*;
use orbweaver::utils;
use sysinfo::{System};
use std::sync::Arc;

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
#[derive(Parser, Debug, Clone)]
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
    /// Processes all sequences in multi-sequence files by default.
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
    /// If set to 0, will determine size dynamically.
    #[clap(long, value_parser, default_value_t = 0)]
    chunk_size: usize,

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
    
    /// Enable adaptive chunk sizing.
    /// 
    /// Dynamically adjust chunk sizes based on sequence complexity
    /// and available memory. Requires chunked mode.
    #[clap(long, action = clap::ArgAction::SetTrue)]
    adaptive_chunking: bool,
    
    /// Process only specific sequences by index (0-based).
    /// 
    /// For multi-sequence FASTA files, process only the sequences
    /// with the specified indices. Default: process all sequences.
    #[clap(long, value_parser, value_delimiter = ',')]
    sequence_indices: Option<Vec<usize>>,
    
    /// Maximum memory usage per chunk (in MB).
    /// 
    /// Limits the memory usage of each chunk when adaptive chunking
    /// is enabled. Helps prevent out-of-memory errors.
    #[clap(long, value_parser)]
    max_memory_per_chunk_mb: Option<usize>,
}

#[tokio::main]
async fn main() -> Result<()> {
    let args = OrbweaverArgs::parse();
    
    // Start profiling if requested
    #[cfg(feature = "profiling")]
    let guard = if args.profile {
        println!("Starting profiling");
        let guard = pprof::ProfilerGuardBuilder::default().frequency(100).build()?;
        Some(guard)
    } else {
        None
    };

    println!("Orbweaver: Processing file {}", args.input.display());
    
    let start = Instant::now();

    // Determine execution mode: parallel chunking or standard
    let grammar = if args.chunk_size > 0 || args.adaptive_chunking {
        let args_clone = args.clone();
        tokio::task::spawn_blocking(move || process_chunked_mode(&args_clone))
            .await?
            .context("Chunked processing failed")?
    } else {
        process_standard_mode(&args).await?
    };
    
    println!("Grammar construction completed in {:?}", start.elapsed());
    
    // Write outputs in requested formats
    generate_outputs(&grammar, &args)?;
    
    // Stop profiling if it was started
    #[cfg(feature = "profiling")]
    if let Some(guard) = guard {
        if let Ok(report) = guard.report().build() {
            if let Some(parent) = std::path::Path::new("profile").parent() {
                 std::fs::create_dir_all(parent)?;
             }
             let file = File::create("profile/flamegraph.svg")
                 .context("Failed to create flamegraph file")?;
            report.flamegraph(file).context("Failed to write flamegraph")?;
            println!("Wrote profiling data to profile/flamegraph.svg");
        } else {
            eprintln!("Failed to build profiling report.");
        }
    }
    
    Ok(())
}

/// Generates all requested output files from the final grammar.
fn generate_outputs(grammar: &Grammar, args: &OrbweaverArgs) -> Result<()> {
    // Output the grammar in various formats based on CLI args
    if let Some(json_path) = &args.output_json {
        println!("Writing grammar to JSON file: {}", json_path.display());
        write_grammar_json(json_path, &grammar.sequence, &grammar.rules)?;
    }
    
    if let Some(text_path) = &args.output_text {
        println!("Writing grammar to text file: {}", text_path.display());
        utils::export::write_grammar_text(text_path, &grammar.sequence, &grammar.rules)?;
    }
    
    if let Some(gfa_path) = &args.output_gfa {
        println!("Writing grammar to GFA file: {}", gfa_path.display());
        utils::visualization::write_grammar_gfa(gfa_path, grammar)?;
    }
    
    if let Some(dot_path) = &args.visualize {
        let options = DotOptions {
            include_terminals: true,
            include_usage_counts: true,
            color_by_depth: true,
        };
        println!("Writing grammar visualization to DOT file: {}", dot_path.display());
        utils::visualization::write_grammar_dot(dot_path, grammar, &options)?;
    }
    
    if let Some(fasta_path) = &args.export_blocks {
        println!("Exporting grammar rules to FASTA file: {}", fasta_path.display());
        utils::export::export_grammar(grammar, fasta_path, utils::io::OutputFormat::Fasta)?;
    }
    
    if args.stats {
        // Calculate and print statistics
        calculate_and_print_stats(grammar)?;
        
        // Calculate compression stats
        let stats = utils::export::calculate_compression_stats(grammar);
        utils::export::print_compression_stats(&stats)?;
        
        // Memory usage estimation with 2-bit encoding
        if args.use_encoding {
            let original_size = stats.original_size;
            let (savings_pct, savings_desc) = bitvec::estimate_memory_savings(original_size, original_size);
            println!("\n2-Bit Encoding Memory Savings:");
            println!("  Reduced memory usage by approximately {:.1}% {}", savings_pct, savings_desc);
            println!("  Original: {} bytes, With 2-bit encoding: ~{} bytes", 
                original_size, (original_size + 3) / 4);
        }
    }
    
    Ok(())
}

async fn process_standard_mode(args: &OrbweaverArgs) -> Result<Grammar> {
    println!("Using standard mode (loading entire sequence asynchronously)");
    
    let all_sequences = read_fasta_sequences_async(&args.input, args.skip_ns).await
        .context("Failed to read FASTA file asynchronously")?;
    
    if all_sequences.is_empty() {
        bail!("No sequences found in input file");
    }
    
    println!("Read {} sequences.", all_sequences.len());
    
    // Filter sequences based on user-provided indices (if any)
    let selected_sequences = if let Some(ref indices) = args.sequence_indices {
        println!("Processing selected sequences only: {:?}", indices);
        
        let mut filtered = Vec::new();
        for &idx in indices {
            if idx < all_sequences.len() {
                filtered.push(all_sequences[idx].clone());
            } else {
                println!("Warning: Sequence index {} is out of range (max: {}), skipping", 
                         idx, all_sequences.len() - 1);
            }
        }
        
        if filtered.is_empty() {
            bail!("No valid sequences selected for processing");
        }
        
        filtered
    } else {
        all_sequences
    };
    
    println!("Processing {} sequences in total", selected_sequences.len());
    
    let args_clone = Arc::new(args.clone());
    
    // Process each sequence in parallel using Rayon
    let mut grammars = vec![];
    
    for (seq_id, (record_id, bases)) in selected_sequences.iter().enumerate() {
        println!("Processing sequence {} ({}): {} bases", seq_id, record_id, bases.len());
        
        let args_clone = Arc::clone(&args_clone);
        let bases_clone = bases.clone();
        
        // Process each sequence in a separate blocking task
        let grammar_result = tokio::task::spawn_blocking(move || -> Result<Grammar, anyhow::Error> {
            if args_clone.use_encoding {
                let encoded_bases = encode_dna(&bases_clone);
                println!("  Using 2-bit encoding ({} bases -> {} encoded bases)", 
                         bases_clone.len(), encoded_bases.len());
                
                let mut grammar_builder = GrammarBuilder::new(args_clone.min_rule_usage, args_clone.reverse_aware);
                if let Some(max_rules) = args_clone.max_rule_count {
                    grammar_builder = grammar_builder.with_max_rules(max_rules);
                }
                grammar_builder.build_grammar(&encoded_bases)?;
                let (sequence, rules) = grammar_builder.get_grammar();
                let max_depth = grammar_builder.get_max_rule_depth();
                
                Ok(Grammar {
                    sequence: sequence.clone(),
                    rules: rules.clone(),
                    max_depth,
                })
            } else {
                println!("  WARNING: Processing without 2-bit encoding.");
                let encoded_bases_fallback = encode_dna(&bases_clone);
                let mut grammar_builder = GrammarBuilder::new(args_clone.min_rule_usage, args_clone.reverse_aware);
                grammar_builder.build_grammar(&encoded_bases_fallback)?;
                let (sequence, rules) = grammar_builder.get_grammar();
                let max_depth = grammar_builder.get_max_rule_depth();
                Ok(Grammar { sequence: sequence.clone(), rules: rules.clone(), max_depth })
            }
        }).await?
          .with_context(|| format!("Failed to build grammar for sequence {}", record_id))?;
        
        grammars.push(grammar_result);
    }
    
    println!("Merging {} individual grammars...", grammars.len());
    
    // Merge all grammars
    let dummy_config = ChunkingConfig::default();
    let total_len_estimate = selected_sequences.iter().map(|(_, b)| b.len()).sum();
    let (merged_grammar, merge_metrics) = merge_grammars(grammars, &dummy_config, total_len_estimate)?;
    
    println!("Merging complete. Final grammar has {} rules.", merged_grammar.rules.len());
    println!("Total merging time: {:?}", merge_metrics.merge_time);
    
    Ok(merged_grammar)
}

fn process_chunked_mode(args: &OrbweaverArgs) -> Result<Grammar> {
    // Determine chunk size: use CLI arg, calculate dynamically, or use adaptive sizing
    let effective_chunk_size = if args.chunk_size > 0 {
        println!("Using provided chunk size: {}", args.chunk_size);
        args.chunk_size
    } else {
        // Dynamic chunk sizing based on available memory
        let mut sys = System::new_all();
        sys.refresh_memory(); // Refresh memory information
        let available_memory = sys.available_memory() as usize; // In bytes
        
        // Heuristic: Use 1/4 of available memory, with min/max bounds
        let calculated_size = (available_memory / 4) as usize;
        let min_chunk_size = 1_000_000; // 1MB
        let max_chunk_size = 500_000_000; // 500MB
        
        let dynamic_size = calculated_size.clamp(min_chunk_size, max_chunk_size);
        println!("Dynamically calculated chunk size based on available memory ({:.2} GB): {} bytes", 
                 available_memory as f64 / 1_000_000_000.0,
                 dynamic_size);
        dynamic_size
    };

    println!("Using parallel chunking mode with effective chunk size: {}", effective_chunk_size);
    
    // Configure the chunking parameters
    let max_memory_per_chunk = args.max_memory_per_chunk_mb.map(|mb| mb * 1024 * 1024);
    
    let config = ChunkingConfig {
        chunk_size: effective_chunk_size,
        overlap_size: args.chunk_overlap,
        min_rule_usage: args.min_rule_usage,
        reverse_aware: args.reverse_aware,
        num_threads: args.threads.unwrap_or_else(num_cpus::get),
        show_progress: true,
        adaptive_chunking: args.adaptive_chunking,
        max_memory_per_chunk,
    };
    
    if !args.use_encoding {
        bail!("Chunked mode requires use_encoding to be enabled for memory efficiency.");
    }
    
    // Read all sequences from the FASTA file
    let all_sequences = read_fasta_sequences(&args.input, args.skip_ns)
        .context("Failed to read FASTA sequence(s) for chunked mode")?;
    
    if all_sequences.is_empty() {
        bail!("No sequences found in input file for chunking");
    }
    
    // Filter sequences based on user-provided indices (if any)
    let selected_sequences = if let Some(ref indices) = args.sequence_indices {
        println!("Processing selected sequences only: {:?}", indices);
        
        let mut filtered = Vec::new();
        for &idx in indices {
            if idx < all_sequences.len() {
                filtered.push(all_sequences[idx].clone());
            } else {
                println!("Warning: Sequence index {} is out of range (max: {}), skipping", 
                         idx, all_sequences.len() - 1);
            }
        }
        
        if filtered.is_empty() {
            bail!("No valid sequences selected for processing");
        }
        
        filtered
    } else {
        all_sequences
    };
    
    println!("Processing {} sequences in chunked mode", selected_sequences.len());
    
    // Process each sequence
    let mut grammars = Vec::new();
    let mut total_bases = 0;
    
    for (seq_idx, (record_id, bases)) in selected_sequences.iter().enumerate() {
        println!("Processing sequence {} ({}): {} bases", seq_idx, record_id, bases.len());
        total_bases += bases.len();
        
        // Convert to EncodedBase for parallel processing
        let encoded_bases = encode_dna(bases);
        
        // Run the parallel sequitur algorithm for this sequence
        let (grammar, metrics) = parallel_sequitur(&encoded_bases, config.clone())?;
        
        // Display execution metrics for this sequence
        println!("Parallel execution metrics for sequence {}:", seq_idx);
        println!("  Chunk count: {}", metrics.chunk_count);
        println!("  Chunk processing time: {:?}", metrics.chunk_processing_time);
        println!("  Merge time: {:?}", metrics.merge_time);
        println!("  Deduplication time: {:?}", metrics.deduplication_time);
        println!("  Total time: {:?}", metrics.total_time);
        println!("  Rule count (before merge): {}", metrics.rules_before_merge);
        println!("  Rule count (after merge): {}", metrics.rules_after_merge);
        
        grammars.push(grammar);
    }
    
    // Calculate memory savings from 2-bit encoding
    let (saving_pct, description) = bitvec::estimate_memory_savings(total_bases, total_bases / 4);
    println!("Memory savings from 2-bit encoding: {:.1}% ({})", saving_pct, description);
    
    // If we processed multiple sequences, merge their grammars
    if grammars.len() > 1 {
        println!("Merging {} sequence grammars...", grammars.len());
        let merge_start = Instant::now();
        let (merged_grammar, merge_metrics) = merge_grammars(grammars, &config, total_bases)?;
        
        println!("Multi-sequence merging complete in {:?}", merge_start.elapsed());
        println!("Final merged grammar has {} rules", merged_grammar.rules.len());
        
        Ok(merged_grammar)
    } else if grammars.len() == 1 {
        // Return the single grammar directly
        Ok(grammars.remove(0))
    } else {
        bail!("No grammars were produced")
    }
}

fn encode_dna(bases: &[u8]) -> Vec<EncodedBase> {
    bases.iter()
        .filter_map(|&b| EncodedBase::from_base(b))
        .collect()
}
