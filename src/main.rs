use anyhow::{Context, Result};
use clap::{Parser, command};
use std::path::PathBuf;
use std::time::Instant;
use std::fs::File;
use std::io::Write;
use orbweaver::fasta::reader::{read_fasta_sequences, InMemoryFastaReader, ChunkedSequenceIterator};
use orbweaver::fasta::encoder::{encode_dna_2bit, decode_dna_2bit};
use orbweaver::grammar::builder::GrammarBuilder;
use orbweaver::io::output_json::write_grammar_json;
use orbweaver::io::output_gfa::write_grammar_gfa;
use orbweaver::io::output_text::write_grammar_text;
use orbweaver::io::output_fasta::export_rules_as_fasta;
use orbweaver::analysis::stats::calculate_and_print_stats;
use orbweaver::io::output_dot::write_grammar_dot;

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
    #[clap(short, long, value_parser, default_value_t = 2)]
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
    
    /// Use 2-bit encoding to reduce memory usage.
    /// 
    /// When enabled, the input DNA sequence is encoded using 2 bits per base
    /// before being processed, significantly reducing memory consumption.
    /// This is especially useful for large genomes.
    #[clap(long, default_value_t = true, action = clap::ArgAction::Set)]
    use_encoding: bool,

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
}

fn main() -> Result<()> {
    // Initialize logger for better error and warning messages
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info")).init();
    
    let args = OrbweaverArgs::parse();
    
    // Set up thread pool for parallel processing
    if let Some(thread_count) = args.threads {
        rayon::ThreadPoolBuilder::new()
            .num_threads(thread_count)
            .build_global()
            .context("Failed to configure thread pool")?;
        println!("Using {} threads for parallel operations", thread_count);
    } else {
        println!("Using default thread count: {} threads", rayon::current_num_threads());
    }

    #[cfg(feature = "profiling")]
    let _guard = if args.profile {
        println!("Profiling enabled, generating flamegraph");
        Some(profiling::start_profiling("orbweaver")?)
    } else {
        None
    };

    let start_time = Instant::now();

    println!("Starting Orbweaver...");
    println!("Input file: {}", args.input.display());
    if let Some(out) = &args.output_json {
        println!("Output JSON: {}", out.display());
    }
    println!("K-mer size: {}", args.kmer_size);
    println!("Skip Ns: {}", args.skip_ns);
    println!("Reverse Aware: {}", args.reverse_aware);
    println!("Using 2-bit encoding: {}", args.use_encoding);
    if let Some(cs) = args.chunk_size {
        println!("Chunk Size: {}", cs);
        println!("Chunk Overlap: {}", args.chunk_overlap);
    }
    if let Some(max) = args.max_rule_count {
        println!("Max rule count: {}", max);
    }

    // --- Input Validation ---
    if !args.input.exists() {
        anyhow::bail!("Input file not found: {}", args.input.display());
    }
    if !args.input.is_file() {
        anyhow::bail!("Input path is not a file: {}", args.input.display());
    }

    if args.kmer_size == 0 {
        anyhow::bail!("K-mer size must be greater than 0.");
    }
    
    if args.kmer_size > 255 {
        anyhow::bail!("K-mer size must not exceed 255.");
    }
    
    if let Some(chunk_s) = args.chunk_size {
         if chunk_s <= args.chunk_overlap {
            anyhow::bail!("Chunk overlap ({}) must be smaller than chunk size ({}).", 
                         args.chunk_overlap, chunk_s);
        }
    }
    
    if let Some(max_rules) = args.max_rule_count {
        if max_rules > 0 && max_rules < args.min_rule_usage {
             // Use eprintln for warnings
             eprintln!("Warning: --max-rule-count ({}) is less than --min-rule-usage ({}). This might lead to unexpected eviction behavior.", max_rules, args.min_rule_usage);
        }
        if max_rules == 0 {
            anyhow::bail!("--max-rule-count cannot be zero.");
        }
    }
    if args.min_rule_usage == 0 {
        anyhow::bail!("--min-rule-usage must be at least 1.");
    }

    println!("Configuration validated.");
    println!("----------------------------------------");

    // --- Core Logic --- 
    
    let read_start = Instant::now();

    // 1. Read FASTA sequence(s)
    let sequences = read_fasta_sequences(&args.input, args.skip_ns)
        .context("Failed to read FASTA file")?;
    
    if sequences.is_empty() {
        println!("No sequences found in the input file after processing.");
        return Ok(());
    }

    let read_time = read_start.elapsed();
    println!("FASTA reading completed in {:.2?}", read_time);

    // For now, process only the first sequence.
    let (first_seq_id, first_seq_data) = &sequences[0];
    println!(
        "Processing first sequence: '{}' ({} bases)",
        first_seq_id,
        first_seq_data.len()
    );
    if sequences.len() > 1 {
         println!("Warning: Input FASTA has multiple sequences. Only the first one ('{}') will be processed.", first_seq_id);
    }

    let initial_len = first_seq_data.len(); // Store initial length for stats

    // Create the grammar builder with the specified parameters
    let grammar_start = Instant::now();
    let mut grammar_builder = if let Some(max_rules) = args.max_rule_count {
        GrammarBuilder::new(args.min_rule_usage, args.reverse_aware)
            .with_max_rules(max_rules)
    } else {
        GrammarBuilder::new(args.min_rule_usage, args.reverse_aware)
    };

    // 2. Process the sequence based on encoding and chunking options
    if let Some(chunk_size) = args.chunk_size {
        println!("Processing sequence in chunks of size {} with {} bases overlap", 
                 chunk_size, args.chunk_overlap);
                 
        // Create a reader for chunked processing
        let mut reader = InMemoryFastaReader::new(vec![(first_seq_id.clone(), first_seq_data.clone())]);
        let chunk_iterator = ChunkedSequenceIterator::new(&mut reader, chunk_size, args.chunk_overlap, 0);
        
        let mut chunks_processed = 0;
        let mut total_bases_processed = 0;
        
        // Process each chunk
        for chunk in chunk_iterator {
            let chunk_start = Instant::now();
            chunks_processed += 1;
            total_bases_processed += chunk.data.len();
            
            println!("Processing chunk {}: bases {} to {} ({} bases)", 
                     chunks_processed, chunk.start_pos, chunk.end_pos, chunk.data.len());
            
            let sequence_to_process = if args.use_encoding {
                match encode_dna_2bit(&chunk.data) {
                    Ok(encoded) => {
                        let decoded = decode_dna_2bit(&encoded, chunk.data.len());
                        decoded
                    },
                    Err(e) => {
                        println!("Warning: Failed to encode chunk: {}. Using raw data.", e);
                        chunk.data
                    }
                }
            } else {
                chunk.data
            };
            
            // Process this chunk with the same grammar builder (continuing from previous chunks)
            grammar_builder.build_grammar(&sequence_to_process)
                .with_context(|| format!("Failed processing chunk {}", chunks_processed))?;
            
            let chunk_time = chunk_start.elapsed();
            println!("Completed chunk {}: current grammar has {} rules (took {:.2?})", 
                     chunks_processed, grammar_builder.get_grammar().1.len(), chunk_time);
                     
            if chunk.is_last {
                println!("Reached end of sequence after {} chunks ({} bases processed)", 
                         chunks_processed, total_bases_processed);
            }
        }
    } else {
        // Process the entire sequence at once
        let encode_start = Instant::now();
        let sequence_to_process = if args.use_encoding {
            println!("Encoding sequence using 2-bit representation to reduce memory usage...");
            match encode_dna_2bit(first_seq_data) {
                Ok(encoded) => {
                    println!("  Original size: {} bytes, Encoded size: {} bytes ({}% reduction)",
                        first_seq_data.len(),
                        encoded.len(),
                        (100.0 * (first_seq_data.len() - encoded.len()) as f64 / first_seq_data.len() as f64) as u32
                    );
                    // Decode back for processing
                    let decoded = decode_dna_2bit(&encoded, first_seq_data.len());
                    decoded
                },
                Err(e) => {
                    println!("Warning: Failed to encode sequence: {}. Falling back to raw processing.", e);
                    first_seq_data.clone()
                }
            }
        } else {
            first_seq_data.clone()
        };
        println!("Encoding completed in {:.2?}", encode_start.elapsed());

        // Build the grammar from the full sequence
        let build_start = Instant::now();
        grammar_builder.build_grammar(&sequence_to_process)
            .context("Failed during grammar construction")?;
        println!("Grammar construction completed in {:.2?}", build_start.elapsed());
    }

    let grammar_time = grammar_start.elapsed();
    // 4. Get final grammar (sequence and rules)
    let (final_sequence, rules) = grammar_builder.get_grammar();
    println!("----------------------------------------");
    println!("Grammar construction complete in {:.2?}.", grammar_time);
    println!("  Final sequence length: {}", final_sequence.len());
    println!("  Number of rules generated: {}", rules.len());
    
    // Print rule depth information if available
    println!("  Maximum rule depth: {}", grammar_builder.get_max_rule_depth());
    println!("  Average rule depth: {:.2}", grammar_builder.get_avg_rule_depth());

    // 5. Write Output Files
    let output_start = Instant::now();
    if let Some(json_path) = &args.output_json {
        let json_start = Instant::now();
        write_grammar_json(&grammar_builder, json_path)
            .context("Failed to write grammar to JSON")?;
        println!("JSON output completed in {:.2?}", json_start.elapsed());
    }
    if let Some(gfa_path) = &args.output_gfa {
        let gfa_start = Instant::now();
        write_grammar_gfa(&grammar_builder, gfa_path)
            .context("Failed to write grammar to GFA")?;
        println!("GFA output completed in {:.2?}", gfa_start.elapsed());
    }
    if let Some(text_path) = &args.output_text {
        let text_start = Instant::now();
        write_grammar_text(&grammar_builder, text_path)
            .context("Failed to write grammar to Text")?;
        println!("Text output completed in {:.2?}", text_start.elapsed());
    }
    if let Some(fasta_path) = &args.export_blocks {
         let fasta_start = Instant::now();
         export_rules_as_fasta(&grammar_builder, fasta_path)
            .context("Failed to export rules to FASTA")?;
         println!("FASTA output completed in {:.2?}", fasta_start.elapsed());
    }
    println!("All outputs completed in {:.2?}", output_start.elapsed());

    // 6. Calculate and Print Stats
    if args.stats {
        let stats_start = Instant::now();
        calculate_and_print_stats(&grammar_builder, initial_len)
            .context("Failed to calculate statistics")?;
        println!("Stats calculation completed in {:.2?}", stats_start.elapsed());
    }

    // 7. Write Visualization
    if let Some(dot_path) = &args.visualize {
        let dot_start = Instant::now();
        write_grammar_dot(&grammar_builder, dot_path)
            .context("Failed to write grammar to DOT")?;
        println!("DOT visualization completed in {:.2?}", dot_start.elapsed());
    }

    #[cfg(feature = "profiling")]
    if args.profile {
        profiling::finish_profiling(_guard, "orbweaver")?;
    }

    println!("----------------------------------------");
    println!("Orbweaver finished in {:.2?}.", start_time.elapsed());
    
    // Write a simple performance summary file 
    let summary_path = args.input.with_file_name(format!(
        "{}_perf_summary.txt", 
        args.input.file_stem().unwrap_or_default().to_string_lossy()
    ));
    
    let mut summary_file = File::create(&summary_path)?;
    writeln!(summary_file, "Orbweaver Performance Summary")?;
    writeln!(summary_file, "------------------------")?;
    writeln!(summary_file, "Input file: {}", args.input.display())?;
    writeln!(summary_file, "Sequence length: {} bases", initial_len)?;
    writeln!(summary_file, "Threads used: {}", rayon::current_num_threads())?;
    writeln!(summary_file, "Total runtime: {:.2?}", start_time.elapsed())?;
    writeln!(summary_file, "FASTA reading: {:.2?} ({:.1}%)", 
        read_time, 
        100.0 * read_time.as_secs_f64() / start_time.elapsed().as_secs_f64()
    )?;
    writeln!(summary_file, "Grammar construction: {:.2?} ({:.1}%)", 
        grammar_time,
        100.0 * grammar_time.as_secs_f64() / start_time.elapsed().as_secs_f64()
    )?;
    writeln!(summary_file, "------------------------")?;
    writeln!(summary_file, "Rules generated: {}", rules.len())?;
    writeln!(summary_file, "Final sequence length: {}", final_sequence.len())?;
    writeln!(summary_file, "Compression ratio: {:.4}", final_sequence.len() as f64 / initial_len as f64)?;
    
    println!("Performance summary written to: {}", summary_path.display());
    
    Ok(())
}
