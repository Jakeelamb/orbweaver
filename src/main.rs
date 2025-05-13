use anyhow::{Context, Result, bail};
use clap::{Parser, command};
use std::path::PathBuf;
use std::time::Instant;
use std::fs;
use chrono::{Utc, DateTime};
use serde::{Serialize, Deserialize};
use orbweaver::fasta::reader::{read_sequences_from_multiple_files, read_sequences_from_multiple_files_async};
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
use sysinfo::{System, SystemExt};
use std::sync::Arc;
use orbweaver::gpu::GpuContext;
use std::collections::HashMap;
use std::fs::File;
use bincode;

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
#[derive(Parser, Debug, Clone, Serialize, Deserialize)]
#[command(author = "Orbweaver Team", version, about, long_about = None)]
#[command(help_template = "\
{before-help}{name} {version}
{author-with-newline}{about-with-newline}
{usage-heading} {usage}

{all-args}{after-help}
")]
struct OrbweaverArgs {
    /// Input FASTA file paths (.fa, .fasta, .fna), comma-separated.
    ///
    /// Paths to files containing DNA sequences in FASTA format.
    /// Processes all sequences in multi-sequence files by default.
    #[clap(short, long, value_parser, required = true, value_delimiter = ',')]
    input_files: Vec<PathBuf>,

    // --- Workflow and Output Organization ---
    /// Base directory for all outputs.
    ///
    /// All run-specific outputs will be placed in subdirectories under this path.
    #[clap(short = 'o', long, value_parser, required = true)]
    output_dir: PathBuf,

    /// Species identifier (e.g., "homo_sapiens", "drosophila_melanogaster").
    ///
    /// Used to organize outputs by species.
    #[clap(long, value_parser, required = true)]
    species_id: String,

    /// Assembly identifier (e.g., "GRCh38.p13", "dm6").
    ///
    /// Used to organize outputs by a specific genome assembly version.
    #[clap(long, value_parser, required = true)]
    assembly_id: String,

    /// Custom run identifier (optional).
    ///
    /// If not provided, a timestamp-based ID will be generated.
    /// Used to distinguish different runs on the same assembly.
    #[clap(long, value_parser)]
    run_id: Option<String>,

    // --- Resumption Options ---
    /// Resume a previous run.
    ///
    /// If this flag is set, --resume-run-dir must also be provided.
    /// Most other arguments will be ignored and loaded from the metadata of the run to be resumed.
    #[clap(long, action = clap::ArgAction::SetTrue)]
    resume: bool,

    /// Path to the specific run directory to resume.
    ///
    /// Only used if --resume is active. This should be the full path to a previous run's output directory
    /// (e.g., output_base/species/assembly/run_id_timestamp/).
    #[clap(long, value_parser, requires = "resume")]
    resume_run_dir: Option<PathBuf>,

    /// Force re-run if a completed run with the same ID is found.
    ///
    /// If set, will overwrite existing completed run data and re-process.
    /// Also applies when --resume is used on a completed run.
    #[clap(long, action = clap::ArgAction::SetTrue)]
    force_rerun: bool,

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

    /// Output tabular summary of repeats.
    ///
    /// Writes a text file summarizing identified repeats, sorted by size or frequency.
    #[clap(long, value_parser)]
    output_repeats: Option<PathBuf>,

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
    
    /// Chunk size for FastaStream when streaming mode is enabled (bytes).
    /// Default: 1MB (1024 * 1024).
    #[clap(long, value_parser)]
    chunk_size_streaming: Option<usize>,
    
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

    /// Disable GPU acceleration (CPU fallback).
    ///
    /// GPU acceleration is used by default. This flag forces CPU-only processing.
    #[clap(long, action = clap::ArgAction::SetTrue)]
    no_gpu: bool,
}

#[derive(Serialize, Deserialize, Debug, Clone)]
struct RunMetadata {
    args: OrbweaverArgs, // Storing a clone of the args used for this run
    tool_version: String,
    status: String, // e.g., "starting", "completed", "failed", "checkpointed"
    start_time: DateTime<Utc>,
    end_time: Option<DateTime<Utc>>,
    run_id: String, // The specific ID for this run
    processed_sequence_checkpoints: HashMap<String, PathBuf>, // Maps sequence ID to its Grammar checkpoint file
    // Potentially add last_checkpoint_path: Option<PathBuf> later for other modes
}

#[tokio::main]
async fn main() -> Result<()> {
    let mut args = OrbweaverArgs::parse();
    let current_tool_version = command!().get_version().unwrap_or("unknown").to_string();

    let run_specific_output_dir: PathBuf;
    let run_id_for_metadata: String;
    let effective_args: OrbweaverArgs; // Args that will be used for the run
    let mut current_metadata: Option<RunMetadata> = None; // Declare at higher scope

    if args.resume {
        println!("Attempting to resume run...");
        let resume_dir = args.resume_run_dir.clone().context("--resume flag requires --resume-run-dir to be specified.")?;
        
        if !resume_dir.is_dir() {
            bail!("Resume directory not found or is not a directory: {:?}", resume_dir);
        }
        run_specific_output_dir = resume_dir;

        let metadata_path = run_specific_output_dir.join("run_metadata.json");
        if !metadata_path.exists() {
            bail!("Metadata file (run_metadata.json) not found in resume directory: {:?}", metadata_path);
        }

        let content = fs::read_to_string(&metadata_path)
            .with_context(|| format!("Failed to read metadata file from: {:?}", metadata_path))?;
        let mut loaded_metadata = serde_json::from_str::<RunMetadata>(&content)
            .with_context(|| format!("Failed to parse metadata from: {:?}", metadata_path))?;

        println!("Resuming run_id: {} from species: {}, assembly: {}.", 
                 loaded_metadata.run_id, loaded_metadata.args.species_id, loaded_metadata.args.assembly_id);
        println!("Original run started at: {}. Original tool version: {}", loaded_metadata.start_time, loaded_metadata.tool_version);

        if loaded_metadata.tool_version != current_tool_version {
            println!(
                "Warning: Resuming with tool version {} but original run used version {}.",
                current_tool_version,
                loaded_metadata.tool_version
            );
        }

        if loaded_metadata.status == "completed" {
            if !args.force_rerun {
                println!("Run {} was already completed successfully on {}. Output at: {:?}", 
                         loaded_metadata.run_id, 
                         loaded_metadata.end_time.map_or_else(|| "N/A".to_string(), |t| t.to_rfc3339()), 
                         run_specific_output_dir);
                println!("Use --force-rerun to re-process this completed run.");
                return Ok(());
            } else {
                println!("Warning: --force-rerun specified. Re-processing completed run {}.", loaded_metadata.run_id);
                // Proceed to treat as a normal resumption, metadata will be updated to starting next.
            }
        }

        effective_args = loaded_metadata.args.clone();
        // However, the output directory related fields in the *original* args object
        // (output_dir, species_id, assembly_id, run_id) are now defined by the resume_run_dir context.
        // We need to ensure these are correctly set for any logic that might use them directly from `args`
        // rather than `effective_args` or `run_specific_output_dir`.
        // For clarity, let's update the original `args` to reflect the resumed context.
        // This might be important if any part of the code *only* looks at the initially parsed `args`.
        args.output_dir = run_specific_output_dir.parent().unwrap().parent().unwrap().parent().unwrap().to_path_buf();
        args.species_id = run_specific_output_dir.parent().unwrap().parent().unwrap().file_name().unwrap().to_string_lossy().into_owned();
        args.assembly_id = run_specific_output_dir.parent().unwrap().file_name().unwrap().to_string_lossy().into_owned();
        args.run_id = Some(run_specific_output_dir.file_name().unwrap().to_string_lossy().into_owned());

        run_id_for_metadata = loaded_metadata.run_id.clone(); // Use the run_id from the metadata
        // Update metadata for this new attempt
        loaded_metadata.status = "starting".to_string();
        loaded_metadata.start_time = Utc::now(); // New start time for this resume attempt
        loaded_metadata.end_time = None;
        loaded_metadata.tool_version = current_tool_version.clone(); // Update to current tool version

        fs::write(
            &metadata_path,
            serde_json::to_string_pretty(&loaded_metadata)?
        ).with_context(|| format!("Failed to write updated metadata to {:?}", metadata_path))?;
        println!("Updated metadata for resumed run: status - starting");
        current_metadata = Some(loaded_metadata); // For use later in success/failure update

    } else {
        // This is a new run (not resuming)
        if args.output_dir == PathBuf::from("") || args.species_id.is_empty() || args.assembly_id.is_empty() {
            // This check is a bit redundant due to `required=true` on these args in clap,
            // but good for explicit safety if those were ever removed.
            bail!("--output-dir, --species-id, and --assembly-id are required for a new run.");
        }
        run_id_for_metadata = args.run_id.clone().unwrap_or_else(|| Utc::now().format("%Y%m%d_%H%M%S").to_string());
        run_specific_output_dir = args.output_dir.join(&args.species_id).join(&args.assembly_id).join(&run_id_for_metadata);
        effective_args = args.clone(); // For a new run, effective_args are the parsed args

        if !run_specific_output_dir.exists() {
            fs::create_dir_all(&run_specific_output_dir)
                .with_context(|| format!("Failed to create output directory: {:?}", run_specific_output_dir))?;
            println!("Created output directory: {:?}", run_specific_output_dir);
        } else {
            println!("Using existing output directory for new run: {:?}", run_specific_output_dir);
        }
        
        let metadata_path = run_specific_output_dir.join("run_metadata.json");
        let new_run_metadata = RunMetadata {
            args: effective_args.clone(),
            tool_version: current_tool_version.clone(),
            status: "starting".to_string(),
            start_time: Utc::now(),
            end_time: None,
            run_id: run_id_for_metadata.clone(),
            processed_sequence_checkpoints: HashMap::new(),
        };

        if metadata_path.exists() {
            if args.force_rerun {
                println!("Warning: --force-rerun specified. Overwriting existing metadata/output for new run {}.", run_id_for_metadata);
            } else {
                // Attempt to read existing metadata to check its status if not forcing rerun
                let existing_meta_read_attempt = fs::read_to_string(&metadata_path);
                match existing_meta_read_attempt {
                    Ok(content) => {
                        match serde_json::from_str::<RunMetadata>(&content) {
                            Ok(existing_meta) => {
                                if existing_meta.status == "completed" {
                                    println!("A run with ID {} already completed successfully in this directory.", existing_meta.run_id);
                                    println!("Output at: {:?}", run_specific_output_dir);
                                    println!("Use --force-rerun to overwrite and re-process, or specify a different --run-id.");
                                    return Ok(());
                                } else {
                                    println!("Warning: Found existing metadata for an incomplete run ({}). Will overwrite and start fresh for new run ID {}.", existing_meta.status, run_id_for_metadata);
                                }
                            }
                            Err(e) => { // serde_json::Error
                                println!("Warning: Could not parse existing metadata file at {:?}: {}. Overwriting.", metadata_path, e);
                            }
                        }
                    }
                    Err(e) => { // std::io::Error
                        println!("Warning: Could not read existing metadata file at {:?}: {}. Overwriting.", metadata_path, e);
                    }
                }
            }
        }

        fs::write(
            &metadata_path,
            serde_json::to_string_pretty(&new_run_metadata)?
        ).with_context(|| format!("Failed to write initial metadata to {:?}", metadata_path))?;
        println!("Wrote initial metadata: status - starting");
        current_metadata = Some(new_run_metadata);
    }

    // The rest of the `main` function will now use `effective_args` for all processing decisions
    // and `run_specific_output_dir` for all output locations.
    // The `current_metadata` variable (which is Option<RunMetadata>) needs to be passed or updated for final status.

    // Re-assign `args` to `effective_args` to simplify the rest of the function if it refers to `args`.
    // This is a bit of a shadow, but makes diffs smaller for the rest of the code.
    // Alternatively, pass `effective_args` down to all functions.
    args = effective_args;

    // Update output paths in `args` to be relative to the `run_specific_output_dir`
    // This needs to happen for both new and resumed runs, using the `args` that will actually be processed.
    if let Some(path_opt) = args.output_json.as_mut() {
        let file_name = path_opt.file_name().unwrap_or_else(|| std::ffi::OsStr::new("grammar.json"));
        *path_opt = run_specific_output_dir.join(file_name);
    }
    if let Some(path_opt) = args.output_text.as_mut() {
        let file_name = path_opt.file_name().unwrap_or_else(|| std::ffi::OsStr::new("grammar.txt"));
        *path_opt = run_specific_output_dir.join(file_name);
    }
    if let Some(path_opt) = args.output_gfa.as_mut() {
        let file_name = path_opt.file_name().unwrap_or_else(|| std::ffi::OsStr::new("grammar.gfa"));
        *path_opt = run_specific_output_dir.join(file_name);
    }
    if let Some(path_opt) = args.visualize.as_mut() {
        let file_name = path_opt.file_name().unwrap_or_else(|| std::ffi::OsStr::new("grammar.dot"));
        *path_opt = run_specific_output_dir.join(file_name);
    }
    if let Some(path_opt) = args.export_blocks.as_mut() {
        let file_name = path_opt.file_name().unwrap_or_else(|| std::ffi::OsStr::new("blocks.fasta"));
        *path_opt = run_specific_output_dir.join(file_name);
    }
    if let Some(path_opt) = args.output_repeats.as_mut() {
        let file_name = path_opt.file_name().unwrap_or_else(|| std::ffi::OsStr::new("repeats.txt"));
        *path_opt = run_specific_output_dir.join(file_name);
    }

    // The `current_metadata` variable now holds the metadata struct that needs to be updated upon completion or failure.
    // The original metadata handling section after this point needs to be adapted.
    // Specifically, instead of creating a new `metadata` variable, we'll use `current_metadata.unwrap()` (after ensuring it's Some).
    
    // --- This section replaces the previous metadata handling logic ---
    let mut metadata_to_update = current_metadata.expect("Metadata should have been initialized for new or resumed run");
    let metadata_path_for_updates = run_specific_output_dir.join("run_metadata.json"); // ensure correct path

    let mut use_gpu = !args.no_gpu;
    let gpu_context: Option<GpuContext> = if use_gpu {
        match GpuContext::new() {
            Ok(ctx) => {
                println!("GPU context initialized successfully.");
                Some(ctx)
            }
            Err(e) => {
                println!("Warning: Failed to initialize GPU context: {}. Falling back to CPU.", e);
                use_gpu = false; // Force CPU mode
                None
            }
        }
    } else {
        None
    };

    // --- Main Processing Logic --- 
    let result: Result<Grammar> = {
        // Conditional profiling setup
        #[cfg(feature = "profiling")]
        let _guard = if args.profile {
            println!("Starting profiling");
            // Adjust pprof output directory if desired, e.g., to run_specific_output_dir.join("profile")
            // For now, it uses the default "profile" directory in CWD.
            let guard = pprof::ProfilerGuardBuilder::default()
                .frequency(100)
                // .blocklist(&["libc", "libgcc", "pthread", "vdso"])
                .build()?;
            Some(guard)
        } else {
            None
        };

        // Configure number of threads for Rayon
        if let Some(threads) = args.threads {
            rayon::ThreadPoolBuilder::new().num_threads(threads).build_global()?;
        }

        // Memory usage logging (consider moving to a more central place or a utility function)
        let mut sys = System::new_all();
        sys.refresh_all();
        println!("Total memory: {} MB", sys.total_memory() / 1024 / 1024);
        println!("Used memory: {} MB", sys.used_memory() / 1024 / 1024);

        let processing_start_time = Instant::now();

        let grammar_result_tuple: Result<(Grammar, HashMap<String, PathBuf>)>;
        let grammar_result_single: Result<Grammar>;

        let final_grammar: Grammar;
        let mut temp_checkpoints: Option<HashMap<String, PathBuf>> = None;

        if args.chunk_size > 0 || args.adaptive_chunking {
            println!("Processing in chunked mode.");
            let chunked_result_tuple = process_chunked_mode(&args, &run_specific_output_dir, &metadata_to_update.processed_sequence_checkpoints);
            match chunked_result_tuple {
                Ok((g, cp)) => {
                    final_grammar = g;
                    temp_checkpoints = Some(cp);
                }
                Err(e) => return Err(e), 
            }
        } else if args.streaming {
            println!("Processing in streaming mode.");
            // TODO: Implement checkpointing for streaming mode with GrammarBuilder
            // For now, it will behave like standard mode for checkpointing demonstration
            grammar_result_tuple = process_standard_mode(&args, use_gpu, gpu_context.as_ref(), &run_specific_output_dir, &metadata_to_update.processed_sequence_checkpoints).await;
            match grammar_result_tuple {
                Ok((g, cp)) => {
                    final_grammar = g;
                    temp_checkpoints = Some(cp);
                }
                Err(e) => return Err(e),
            }
        } else {
            println!("Processing in standard (in-memory) mode.");
            grammar_result_tuple = process_standard_mode(&args, use_gpu, gpu_context.as_ref(), &run_specific_output_dir, &metadata_to_update.processed_sequence_checkpoints).await;
            match grammar_result_tuple {
                Ok((g, cp)) => {
                    final_grammar = g;
                    temp_checkpoints = Some(cp);
                }
                Err(e) => return Err(e),
            }
        };

        // Update metadata with any new checkpoints from standard/streaming (mocked) mode
        if let Some(checkpoints) = temp_checkpoints {
            metadata_to_update.processed_sequence_checkpoints = checkpoints;
            // Optionally, save metadata here if you want to reflect checkpoints immediately 
            // even if the full run doesn't complete. For now, saved at very end.
        }

        // At this point, final_grammar is populated. We return it as Ok.
        Ok(final_grammar) // This makes the block an expression that yields Result<Grammar>
    };

    // --- Handle Results and Finalize Metadata ---
    match result {
        Ok(grammar) => {
            println!("Grammar construction successful.");
            // Use `args` (which is `effective_args`) for stats and outputs
            if args.stats {
                calculate_and_print_stats(&grammar)?;
            }
            // Pass `args` (which is `effective_args`) to generate_outputs
            generate_outputs(&grammar, &args)?;
            
            metadata_to_update.status = "completed".to_string();
            metadata_to_update.end_time = Some(Utc::now());
            fs::write(
                &metadata_path_for_updates, // use the correctly scoped path
                serde_json::to_string_pretty(&metadata_to_update)?
            ).with_context(|| format!("Failed to write final metadata to {:?}", metadata_path_for_updates))?;
            println!("Run {} completed. Metadata updated.", metadata_to_update.run_id);
        }
        Err(e) => {
            eprintln!("Error during grammar construction: {:?}", e);
            metadata_to_update.status = "failed".to_string();
            metadata_to_update.end_time = Some(Utc::now());
            fs::write(
                &metadata_path_for_updates, // use the correctly scoped path
                serde_json::to_string_pretty(&metadata_to_update)?
            ).with_context(|| format!("Failed to write error metadata to {:?}", metadata_path_for_updates))?;
            println!("Run {} failed. Metadata updated.", metadata_to_update.run_id);
            return Err(e); // Propagate the error
        }
    }

    // Stop profiling if it was started
    #[cfg(feature = "profiling")]
    if let Some(guard) = _guard { // This assumes _guard is in scope
        println!("Stopping profiling");
        if let Ok(report) = guard.report().build() {
            let profile_dir = run_specific_output_dir.join("profile");
            fs::create_dir_all(&profile_dir).context("Failed to create profile directory")?;
            // Use `run_id_for_metadata` for the flamegraph name to ensure consistency
            let flamegraph_path = profile_dir.join(format!("flamegraph_{}.svg", run_id_for_metadata));
            let file = fs::File::create(&flamegraph_path)
                .with_context(|| format!("Failed to create flamegraph file: {:?}", flamegraph_path))?;
            report.flamegraph(file)
                .with_context(|| format!("Failed to write flamegraph to {:?}", flamegraph_path))?;
            println!("Flamegraph saved to {:?}", flamegraph_path);
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
    
    // Add output for repeats summary
    if let Some(repeats_path) = &args.output_repeats {
        println!("Writing repeat summary to: {}", repeats_path.display());
        utils::export::write_repeat_summary(repeats_path, grammar)?;
    }
    
    Ok(())
}

async fn process_standard_mode(
    args: &OrbweaverArgs, 
    use_gpu: bool, 
    gpu_context: Option<&GpuContext>,
    run_specific_output_dir: &PathBuf, 
    initial_checkpoints: &HashMap<String, PathBuf>
) -> Result<(Grammar, HashMap<String, PathBuf>)> { 
    println!("Using standard mode (loading entire sequence asynchronously)");

    let checkpoints_dir = run_specific_output_dir.join("checkpoints");
    if !checkpoints_dir.exists() {
        fs::create_dir_all(&checkpoints_dir)
            .with_context(|| format!("Failed to create checkpoints directory: {:?}", checkpoints_dir))?;
    }

    let mut current_run_checkpoints = initial_checkpoints.clone();
    
    let all_sequences = read_sequences_from_multiple_files_async(&args.input_files, args.skip_ns).await
        .context("Failed to read FASTA files asynchronously")?;
    
    if all_sequences.is_empty() {
        bail!("No sequences found in input files");
    }
    
    println!("Read {} sequences.", all_sequences.len());
    
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
    
    let total_len_estimate = selected_sequences.iter().map(|(_, b)| b.len()).sum();

    let mut grammars = vec![];
    let args_arc = Arc::new(args.clone());

    for (seq_idx, (record_id_owned, bases_owned)) in selected_sequences.into_iter().enumerate() {
        let seq_identifier = format!("seq_{}_{}", seq_idx, record_id_owned.replace(|c: char| !c.is_alphanumeric(), "_"));
        let checkpoint_file_name = format!("{}.bincode", seq_identifier);
        let checkpoint_path = checkpoints_dir.join(&checkpoint_file_name);

        if let Some(existing_checkpoint_path) = current_run_checkpoints.get(&seq_identifier) {
            if existing_checkpoint_path.exists() {
                println!("Loading grammar for sequence {} (ID: {}) from checkpoint {:?}", seq_idx, record_id_owned, existing_checkpoint_path);
                match File::open(existing_checkpoint_path) {
                    Ok(file) => {
                        match bincode::deserialize_from(file) {
                            Ok(grammar) => {
                                grammars.push(grammar);
                                continue; // Successfully loaded, move to next sequence
                            }
                            Err(e) => {
                                println!("Warning: Failed to deserialize checkpoint {:?}: {}. Recomputing.", existing_checkpoint_path, e);
                            }
                        }
                    }
                    Err(e) => {
                        println!("Warning: Failed to open checkpoint file {:?}: {}. Recomputing.", existing_checkpoint_path, e);
                    }
                }
            } else {
                println!("Warning: Checkpoint file for sequence {} not found at {:?}, though listed in metadata. Recomputing.", seq_identifier, existing_checkpoint_path);
            }
        }
        
        // If not loaded from checkpoint, process the sequence
        println!("Processing sequence {} (ID: {}): {} bases", seq_idx, record_id_owned, bases_owned.len());
        
        let task_args = Arc::clone(&args_arc);
        let task_gpu_context = gpu_context.cloned();
        
        let grammar_result: Grammar = tokio::task::spawn_blocking(move || -> Result<Grammar, anyhow::Error> {
            if task_args.use_encoding {
                let encoded_bases = encode_dna(&bases_owned); 
                println!("  Using 2-bit encoding ({} bases -> {} encoded bases)", 
                         bases_owned.len(), encoded_bases.len());
                
                let mut grammar_builder = GrammarBuilder::new(task_args.min_rule_usage, task_args.reverse_aware);
                if let Some(max_rules) = task_args.max_rule_count {
                    grammar_builder = grammar_builder.with_max_rules(max_rules);
                }

                let use_gpu_override = false; 

                if use_gpu_override && task_gpu_context.is_some() {
                    grammar_builder = grammar_builder.with_gpu(task_gpu_context.as_ref());
                    grammar_builder.build_grammar_with_gpu(&encoded_bases)?;
                } else {
                    if use_gpu && task_gpu_context.is_none() {
                        println!("Warning: GPU mode selected, but GPU context is unavailable for sequence {}. Falling back to CPU.", record_id_owned);
                    } else if use_gpu_override {
                        println!("INFO: GPU path forced off for debugging.");
                    }
                    grammar_builder.build_grammar(&encoded_bases)?;
                }
                
                let (sequence, rules) = grammar_builder.get_grammar();
                let max_depth = grammar_builder.get_max_rule_depth();
                
                Ok(Grammar {
                    sequence: sequence.clone(),
                    rules: rules.clone(),
                    max_depth,
                    origins: HashMap::new(), 
                })
            } else {
                println!("  WARNING: Processing without 2-bit encoding.");
                let encoded_bases_fallback = encode_dna(&bases_owned); 
                let mut grammar_builder = GrammarBuilder::new(task_args.min_rule_usage, task_args.reverse_aware);

                if use_gpu && task_gpu_context.is_some() {
                    grammar_builder = grammar_builder.with_gpu(task_gpu_context.as_ref());
                    grammar_builder.build_grammar_with_gpu(&encoded_bases_fallback)?;
                } else {
                    if use_gpu && task_gpu_context.is_none() {
                         println!("Warning: GPU mode selected (no encoding), but GPU context is unavailable for sequence {}. Falling back to CPU.", record_id_owned);
                    }
                    grammar_builder.build_grammar(&encoded_bases_fallback)?;
                }
                let (sequence, rules) = grammar_builder.get_grammar();
                let max_depth = grammar_builder.get_max_rule_depth();
                Ok(Grammar {
                    sequence: sequence.clone(), 
                    rules: rules.clone(), 
                    max_depth, 
                    origins: HashMap::new(), 
                })
            }
        }).await?
          .with_context(|| format!("Failed to build grammar for sequence index {}", seq_idx))?;
        
        // Save checkpoint for the newly built grammar
        match File::create(&checkpoint_path) {
            Ok(file_writer) => {
                if let Err(e) = bincode::serialize_into(file_writer, &grammar_result) {
                    println!("Warning: Failed to save checkpoint for sequence {}: {}. Proceeding without checkpoint for this sequence.", seq_identifier, e);
                } else {
                    println!("Saved checkpoint for sequence {} to {:?}", seq_identifier, checkpoint_path);
                    current_run_checkpoints.insert(seq_identifier.clone(), checkpoint_path.clone());
                }
            }
            Err(e) => {
                println!("Warning: Failed to create checkpoint file {:?}: {}. Proceeding without checkpoint for this sequence.", checkpoint_path, e);
            }
        }
        grammars.push(grammar_result);
    }
    
    println!("Merging {} individual grammars...", grammars.len());
    
    let dummy_config = ChunkingConfig::default(); // Assuming this is okay for now
    let (merged_grammar, _merge_metrics) = merge_grammars(grammars, &dummy_config, total_len_estimate)?;
    
    println!("Merging complete. Final grammar has {} rules.", merged_grammar.rules.len());
    // println!("Total merging time: {:?}", merge_metrics.merge_time);
    
    Ok((merged_grammar, current_run_checkpoints))
}

fn process_chunked_mode(
    args: &OrbweaverArgs, 
    run_specific_output_dir: &PathBuf, 
    initial_checkpoints: &HashMap<String, PathBuf>
) -> Result<(Grammar, HashMap<String, PathBuf>)> {
    let use_gpu = !args.no_gpu;
    let mut actual_use_gpu = use_gpu;
    let gpu_context_chunked: Option<GpuContext> = if use_gpu {
        match GpuContext::new() {
            Ok(ctx) => {
                println!("GPU context initialized successfully for chunked mode.");
                Some(ctx)
            }
            Err(e) => {
                println!("Warning: Failed to initialize GPU context for chunked mode: {}. Falling back to CPU.", e);
                actual_use_gpu = false;
                None
            }
        }
    } else {
        None
    };

    let checkpoints_dir = run_specific_output_dir.join("checkpoints");
    if !checkpoints_dir.exists() {
        fs::create_dir_all(&checkpoints_dir)
            .with_context(|| format!("Failed to create checkpoints directory: {:?}", checkpoints_dir))?;
    }
    let mut current_run_checkpoints = initial_checkpoints.clone();

    let effective_chunk_size = if args.chunk_size > 0 {
        println!("Using provided chunk size: {}", args.chunk_size);
        args.chunk_size
    } else {
        // Dynamic chunk sizing based on available memory
        let mut sys = System::new();
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
        use_gpu: actual_use_gpu, // Use the possibly updated flag
    };
    
    if !args.use_encoding {
        bail!("Chunked mode requires use_encoding to be enabled for memory efficiency.");
    }
    
    // Read all sequences from the multiple FASTA files
    let all_sequences = read_sequences_from_multiple_files(&args.input_files, args.skip_ns)
        .context("Failed to read FASTA sequence(s) for chunked mode")?;
    
    if all_sequences.is_empty() {
        bail!("No sequences found in input files for chunking");
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
    
    // Calculate total length estimate *before* consuming selected_sequences
    let total_len_estimate = selected_sequences.iter().map(|(_, b)| b.len()).sum();

    // Process each sequence
    let mut grammars = Vec::new();
    let mut total_bases = 0;
    
    for (seq_idx, (record_id, bases)) in selected_sequences.iter().enumerate() {
        let seq_identifier = format!("seq_{}_{}", seq_idx, record_id.replace(|c: char| !c.is_alphanumeric(), "_"));
        let checkpoint_file_name = format!("{}.bincode", seq_identifier);
        let checkpoint_path = checkpoints_dir.join(&checkpoint_file_name);

        if let Some(existing_checkpoint_path) = current_run_checkpoints.get(&seq_identifier) {
            if existing_checkpoint_path.exists() {
                println!("Loading grammar for chunked sequence {} (ID: {}) from checkpoint {:?}", seq_idx, record_id, existing_checkpoint_path);
                match File::open(existing_checkpoint_path) {
                    Ok(file) => {
                        match bincode::deserialize_from(file) {
                            Ok(grammar) => {
                                grammars.push(grammar);
                                total_bases += bases.len(); // Still count towards total for merge estimate
                                continue; 
                            }
                            Err(e) => {
                                println!("Warning: Failed to deserialize checkpoint {:?}: {}. Recomputing.", existing_checkpoint_path, e);
                            }
                        }
                    }
                    Err(e) => {
                        println!("Warning: Failed to open checkpoint file {:?}: {}. Recomputing.", existing_checkpoint_path, e);
                    }
                }
            } else {
                println!("Warning: Checkpoint file for chunked sequence {} not found at {:?}. Recomputing.", seq_identifier, existing_checkpoint_path);
            }
        }

        println!("Processing sequence {} ({}): {} bases in chunked mode", seq_idx, record_id, bases.len());
        total_bases += bases.len();
        
        let encoded_bases = encode_dna(bases);
        let (grammar_result, metrics) = parallel_sequitur(&encoded_bases, config.clone(), gpu_context_chunked.as_ref())?;
        
        println!("Parallel execution metrics for sequence {}:", seq_idx);
        println!("  Chunk count: {}", metrics.chunk_count);
        println!("  Chunk processing time: {:?}", metrics.chunk_processing_time);
        println!("  Merge time: {:?}", metrics.merge_time);
        println!("  Deduplication time: {:?}", metrics.deduplication_time);
        println!("  Total time: {:?}", metrics.total_time);
        println!("  Rule count (before merge): {}", metrics.rules_before_merge);
        println!("  Rule count (after merge): {}", metrics.rules_after_merge);
        
        match File::create(&checkpoint_path) {
            Ok(file_writer) => {
                if let Err(e) = bincode::serialize_into(file_writer, &grammar_result) {
                    println!("Warning: Failed to save checkpoint for chunked sequence {}: {}. Proceeding without checkpoint.", seq_identifier, e);
                } else {
                    println!("Saved checkpoint for chunked sequence {} to {:?}", seq_identifier, checkpoint_path);
                    current_run_checkpoints.insert(seq_identifier.clone(), checkpoint_path.clone());
                }
            }
            Err(e) => {
                println!("Warning: Failed to create checkpoint file {:?}: {}. Proceeding without checkpoint.", checkpoint_path, e);
            }
        }
        grammars.push(grammar_result);
    }
    
    if grammars.len() > 1 {
        println!("Merging {} sequence grammars (chunked mode)...", grammars.len());
        let merge_start = Instant::now();
        let (merged_grammar, _merge_metrics) = merge_grammars(grammars, &config, total_len_estimate)?;
        
        println!("Multi-sequence merging (chunked mode) complete in {:?}", merge_start.elapsed());
        println!("Final merged grammar (chunked mode) has {} rules", merged_grammar.rules.len());
        
        Ok((merged_grammar, current_run_checkpoints))
    } else if grammars.len() == 1 {
        Ok((grammars.remove(0), current_run_checkpoints))
    } else {
        if selected_sequences.is_empty() {
             bail!("No sequences selected for processing in chunked mode.")
        } else {
            bail!("No grammars were produced in chunked mode from the provided sequences.")
        }
    }
}

fn encode_dna(bases: &[u8]) -> Vec<EncodedBase> {
    bases.iter()
        .filter_map(|&b| EncodedBase::from_base(b))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use clap::CommandFactory;

    #[test]
    fn test_default_streaming_enabled() {
        let cmd = OrbweaverArgs::command();
        let matches = cmd.get_matches_from(vec!["orbweaver", "-i", "input.fa"]);
        let streaming_arg = matches.get_flag("streaming");
        assert!(streaming_arg, "Streaming should be enabled by default if not otherwise specified.");
    }
}