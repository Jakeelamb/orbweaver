#![feature(iter_next_chunk)]
use anyhow::{Context, Result, bail};
use clap::{Parser, command};
use std::path::PathBuf;
use std::fs;
use chrono::{Utc, DateTime};
use serde::{Serialize, Deserialize};
use orbweaver::fasta::reader::{read_sequences_from_multiple_files_async};
use orbweaver::encode::dna_2bit::{EncodedBase};
use orbweaver::grammar::builder::GrammarBuilder;
use orbweaver::grammar::engine::{Grammar};
use orbweaver::analysis::stats::calculate_and_print_stats;
use orbweaver::utils::visualization::DotOptions;
use orbweaver::parallel::chunking::ChunkingConfig;
use orbweaver::parallel::engine::{merge_grammars};
use std::sync::Arc;
use orbweaver::gpu::GpuContext;
use std::collections::HashMap;
use std::fs::File;
use bincode;
use log::{info};
use orbweaver::utils::export::{write_grammar_text, write_grammar_json as io_write_grammar_json, write_grammar_graphml};
use orbweaver::utils::visualization::{write_grammar_gfa, write_grammar_dot};
use orbweaver::io::output_fasta::write_grammar_fasta;
use std::fs::write as fs_write;
use orbweaver::utils::profiling; // Added for profiling RAII guard
use uuid;

// Use the configuration module if needed, or args directly
// use orbweaver::config::Config; 

// Use other library modules as needed later
// use orbweaver::io; // Example
// use orbweaver::fasta;
// use orbweaver::grammar;
// use orbweaver::analysis; // If motifs analysis is re-integrated later

// Configure jemalloc for better memory management if feature is enabled
#[cfg(feature = "profiling")]
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
    #[clap(long, value_parser, default_value_t = 100)]
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

    /// Graph-related arguments
    #[clap(long, value_parser, default_value_t = String::from("sfdp"))]
    graph_engine: String,

    /// Include terminal symbols in graph visualizations.
    #[clap(long, value_parser, default_value_t = true)]
    graph_include_terminals: bool,

    /// Maximum depth of rules to display in graph visualizations.
    #[clap(long, value_parser)]
    graph_max_depth: Option<usize>,

    /// Skip rules above this depth in graph visualizations (useful for very deep grammars).
    #[clap(long, value_parser)]
    graph_skip_rules_above_depth: Option<usize>,

    /// Use a transparent background for graph visualizations (PNG/SVG).
    #[clap(long, value_parser, default_value_t = false)]
    graph_transparent_background: bool,

    /// Use dark mode styling for graph visualizations.
    #[clap(long, value_parser, default_value_t = false)]
    graph_dark_mode: bool,

    /// Show usage counts on nodes in graph visualizations.
    #[clap(long, value_parser, default_value_t = true)]
    graph_show_usage_counts: bool,

    /// Color nodes by depth and show depth information in graph visualizations.
    #[clap(long, value_parser, default_value_t = true)]
    pub graph_show_depth: bool,
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

impl RunMetadata {
    fn save(&self, path: &PathBuf) -> Result<()> {
        let json_string = serde_json::to_string_pretty(self)
            .context("Failed to serialize RunMetadata to JSON")?;
        fs::write(path, json_string)
            .context(format!("Failed to write RunMetadata to file: {}", path.display()))
    }

    fn load(path: &PathBuf) -> Result<Self> {
        let json_string = fs::read_to_string(path)
            .context(format!("Failed to read RunMetadata from file: {}", path.display()))?;
        serde_json::from_str(&json_string)
            .context(format!("Failed to deserialize RunMetadata from JSON: {}", path.display()))
    }
}

#[tokio::main]
async fn main() -> Result<()> {
    // Initialize logging
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info")).init();

    let mut args = OrbweaverArgs::parse();
    let cli_force_rerun = args.force_rerun;

    // --- Profiling Setup ---
    let _profiler_session_guard = if args.profile {
        #[cfg(feature = "profiling")]
        {
            match profiling::start_profiling("orbweaver_run") {
                Ok(guard) => Some(guard),
                Err(e) => {
                    eprintln!("Failed to start profiling: {}. Continuing without profiling.", e);
                    None
                }
            }
        }
        #[cfg(not(feature = "profiling"))]
        {
            println!("Profiling requested via --profile flag, but Orbweaver was not compiled with the 'profiling' feature.");
            println!("To enable profiling, recompile with: cargo build --features profiling");
            match profiling::start_profiling("orbweaver_run_disabled") {
                 Ok(guard) => Some(guard),
                 Err(_) => None, 
            }
        }
    } else {
        None
    };

    // --- Run ID and Output Directory Setup ---
    let _run_id_str = args.run_id.clone().unwrap_or_else(|| {
        uuid::Uuid::new_v4().to_string()
    });

    let first_input_file = args.input_files.get(0).context("No input files provided.")?;
    let _input_file_stem = first_input_file.file_stem()
        .context("Could not extract file stem from input file.")?
        .to_string_lossy();

    let mut metadata_to_save;
    let mut initial_checkpoints: HashMap<String, PathBuf> = HashMap::new();
    let run_specific_output_dir: PathBuf;

    let base_output_dir = PathBuf::from("."); 
    let determined_run_id = args.run_id.clone().unwrap_or_else(|| Utc::now().format("%Y%m%d_%H%M%S").to_string());

    // --- Resumption Logic ---
    if args.resume {
        // Clone resume_run_dir to take ownership *before* args is potentially reassigned.
        if let Some(owned_resume_dir) = args.resume_run_dir.clone() {
            let metadata_path = owned_resume_dir.join("run_metadata.json");
            if metadata_path.exists() {
                println!("Resuming run from: {}", owned_resume_dir.display());
                match RunMetadata::load(&metadata_path) {
                    Ok(loaded_metadata) => {
                        if loaded_metadata.status == "completed" && !cli_force_rerun {
                            println!("Run {} was already completed. Use --force-rerun to re-process.", loaded_metadata.run_id);
                            return Ok(());
                        }
                        
                        let mut original_loaded_args = loaded_metadata.args.clone();
                        
                        if args.output_json.is_some() { original_loaded_args.output_json = args.output_json.clone(); }
                        if args.output_text.is_some() { original_loaded_args.output_text = args.output_text.clone(); }
                        if args.output_gfa.is_some() { original_loaded_args.output_gfa = args.output_gfa.clone(); }

                        args = original_loaded_args; 
                        
                        args.resume = true; 
                        // Use the owned_resume_dir here, which is from the *original* args before overwrite.
                        args.resume_run_dir = Some(owned_resume_dir.clone()); 
                        args.force_rerun = cli_force_rerun; 
                        
                        args.run_id = Some(loaded_metadata.run_id.clone()); 
                        run_specific_output_dir = owned_resume_dir.clone(); // Output to the resumed run's directory.

                        metadata_to_save = loaded_metadata;
                        metadata_to_save.status = if args.force_rerun { "starting (forced rerun on resumed)".to_string() } else { "resuming".to_string() };
                        metadata_to_save.start_time = Utc::now(); 
                        metadata_to_save.end_time = None;
                        if !args.force_rerun { 
                           initial_checkpoints = metadata_to_save.processed_sequence_checkpoints.clone();
                        } else {
                            metadata_to_save.processed_sequence_checkpoints = HashMap::new(); 
                        }

                        println!("Successfully loaded and configured for resuming run: {}", metadata_to_save.run_id);
                        println!("Effective arguments for this resumed run: {:#?}", args);
                    }
                    Err(e) => {
                        bail!("Failed to load metadata from {}: {}. Cannot resume.", metadata_path.display(), e);
                    }
                }
            } else {
                bail!("Metadata file run_metadata.json not found in resume directory: {}. Cannot resume.", owned_resume_dir.display());
            }
        } else {
            bail!("--resume flag was set, but --resume-run-dir was not provided.");
        }
    } else {
        // This is a new run or a non-resumed forced re-run of an existing ID
        run_specific_output_dir = base_output_dir.join(&determined_run_id);
        if run_specific_output_dir.exists() && !cli_force_rerun { 
            println!(
                "Output directory {} already exists. Use --force-rerun to overwrite or a different --run-id.", 
                run_specific_output_dir.display()
            );
            return Ok(());
        }
        if run_specific_output_dir.exists() && cli_force_rerun { 
            fs::remove_dir_all(&run_specific_output_dir)
                 .with_context(|| format!("Failed to remove existing output directory for forced rerun: {}", run_specific_output_dir.display()))?;
        }
        fs::create_dir_all(&run_specific_output_dir)
            .with_context(|| format!("Failed to create run output directory: {}", run_specific_output_dir.display()))?;

        metadata_to_save = RunMetadata {
            args: args.clone(), 
            tool_version: command!().get_version().unwrap_or("unknown").to_string(),
            status: "starting".to_string(),
            start_time: Utc::now(),
            end_time: None,
            run_id: determined_run_id.clone(), 
            processed_sequence_checkpoints: HashMap::new(),
        };
        args.force_rerun = cli_force_rerun;
        args.run_id = Some(determined_run_id); 
    }

    // Save initial metadata (either new or updated for resume)
    let metadata_path = run_specific_output_dir.join("run_metadata.json");
    metadata_to_save.save(&metadata_path)
        .with_context(|| format!("Failed to save initial run metadata to {}", metadata_path.display()))?;

    // Log the arguments being used for this run
    info!("Starting Orbweaver run with arguments: {:#?}", args);
    info!("Run ID: {}", metadata_to_save.run_id);
    info!("Output will be saved to: {}", run_specific_output_dir.display());

    // --- Main Processing Logic ---
    let result: Result<(Grammar, HashMap<String, PathBuf>, Vec<(String, usize)>)> = if args.streaming {
        process_streaming_mode(&args, None, &run_specific_output_dir, &initial_checkpoints).await
    } else if args.chunk_size > 0 || args.adaptive_chunking {
        // TODO: Pass GPU context to chunked mode if/when it supports it
        Ok(process_chunked_mode(&args, &run_specific_output_dir, &initial_checkpoints)?)
    } else {
        process_standard_mode(&args, false, None, &run_specific_output_dir, &initial_checkpoints).await
    };

    match result {
        Ok((final_grammar, final_checkpoints, chromosome_details)) => {
            info!("Grammar construction completed successfully.");
            generate_outputs(&final_grammar, &args, &run_specific_output_dir, &chromosome_details)?;
            
            // Update and save metadata as completed
            metadata_to_save.status = "completed".to_string();
            metadata_to_save.end_time = Some(Utc::now());
            metadata_to_save.processed_sequence_checkpoints = final_checkpoints;
            metadata_to_save.save(&metadata_path)
                .with_context(|| format!("Failed to save final run metadata to {}", metadata_path.display()))?;
            info!("Orbweaver run {} completed successfully.", metadata_to_save.run_id);

            // No explicit finish_profiling call needed, _profiler_session_guard will handle it on drop
            Ok(())
        }
        Err(e) => {
            eprintln!("Error during Orbweaver processing: {:?}", e);
            // Update and save metadata as failed
            metadata_to_save.status = "failed".to_string();
            metadata_to_save.end_time = Some(Utc::now());
            // Checkpoints might be partially updated, save what we have
            // If result was Err, final_checkpoints isn't available in this scope.
            // We should ideally get the checkpoints from the error or pass metadata_to_save into processing functions.
            // For now, just save metadata without updating checkpoints on error.
            if let Err(save_err) = metadata_to_save.save(&metadata_path) {
                 eprintln!("Additionally, failed to save error state metadata to {}: {}", metadata_path.display(), save_err);
            }
            // No explicit finish_profiling call needed here either
            Err(e) // Propagate the original error
        }
    }
}

/// Generates all requested output files from the final grammar.
fn generate_outputs(grammar: &Grammar, args: &OrbweaverArgs, run_specific_output_dir: &PathBuf, chromosome_details: &[(String, usize)]) -> Result<()> {
    info!("Generating outputs in directory: {:?}", run_specific_output_dir);

    let dot_options = DotOptions {
        engine: args.graph_engine.clone(),
        include_terminals: args.graph_include_terminals,
        max_depth: args.graph_max_depth,
        skip_rules_above_depth: args.graph_skip_rules_above_depth,
        transparent_background: args.graph_transparent_background,
        dark_mode: args.graph_dark_mode,
        include_usage_counts: args.graph_show_usage_counts, // Fixed field name
        color_by_depth: args.graph_show_depth,             // Fixed field name
    };

    // Always generate GraphML output
    let graphml_path = run_specific_output_dir.join("grammar.graphml");
    info!("Writing grammar to GraphML: {:?}", graphml_path);
    write_grammar_graphml(&graphml_path, grammar, &dot_options) // Corrected arguments
        .with_context(|| format!("Failed to write GraphML to {:?}", graphml_path))?;

    // Always generate FASTA export of grammar rules (blocks)
    let fasta_export_path = run_specific_output_dir.join("rules.fasta");
    info!("Exporting grammar blocks to FASTA: {:?}", fasta_export_path);
    write_grammar_fasta(&fasta_export_path, grammar)
        .with_context(|| format!("Failed to export rules as FASTA to {:?}", fasta_export_path))?;

    // Always generate DOT file output
    let dot_path = run_specific_output_dir.join("grammar.dot");
    info!("Writing grammar to DOT: {:?}", dot_path);
    write_grammar_dot(&dot_path, grammar, &dot_options) // Use write_grammar_dot directly
        .with_context(|| format!("Failed to write DOT file to {:?}", dot_path))?;

    // JSON output: Always generate, using specified path or defaulting to grammar.json
    let json_output_path = match args.output_json {
        Some(ref user_path) => run_specific_output_dir.join(user_path),
        None => run_specific_output_dir.join("grammar.json"),
    };
    info!("Writing grammar to JSON: {:?}", json_output_path);
    io_write_grammar_json(&json_output_path, &grammar.sequence, &grammar.rules)
        .with_context(|| format!("Failed to write JSON grammar to {:?}", json_output_path))?;

    // Text output (conditional based on user argument)
    if let Some(ref user_path) = args.output_text {
        let final_path = run_specific_output_dir.join(user_path);
        info!("Writing grammar to text: {:?}", final_path);
        write_grammar_text(&final_path, &grammar.sequence, &grammar.rules) // Corrected arguments
            .with_context(|| format!("Failed to write text grammar to {:?}", final_path))?;
    }

    // Always generate GFA output
    let gfa_path = match args.output_gfa {
        Some(ref user_path) => run_specific_output_dir.join(user_path),
        None => run_specific_output_dir.join("grammar.gfa"),
    };
    info!("Writing grammar to GFA: {:?}", gfa_path);
    write_grammar_gfa(&gfa_path, grammar)
        .with_context(|| format!("Failed to write GFA to {:?}", gfa_path))?;

    if let Some(ref user_path) = args.output_repeats {
        let final_path = run_specific_output_dir.join(user_path);
        info!("Writing repeats summary to: {:?}", final_path);
        // Placeholder for actual repeat writing logic
        // grammar.write_repeats_summary(&final_path)?;
        fs_write(&final_path, "Repeat summary placeholder.") // Replace with actual call
            .with_context(|| format!("Failed to write repeats summary to {:?}", final_path))?;
    }

    if args.stats {
        info!("Calculating and printing statistics...");
        let stats_path = run_specific_output_dir.join("grammar_stats.txt");
        calculate_and_print_stats(grammar, &stats_path, args.stats, chromosome_details)
            .with_context(|| format!("Failed to calculate or write stats to {:?}", stats_path))?;
    }

    info!("All requested outputs generated.");
    Ok(())
}

async fn process_standard_mode(
    args: &OrbweaverArgs, 
    use_gpu: bool, 
    gpu_context: Option<&GpuContext>,
    run_specific_output_dir: &PathBuf,
    initial_checkpoints: &HashMap<String, PathBuf>
) -> Result<(Grammar, HashMap<String, PathBuf>, Vec<(String, usize)>)> { 
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

    let chromosome_details: Vec<(String, usize)> = selected_sequences.iter()
        .map(|(name, bases)| (name.clone(), bases.len()))
        .collect();

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
                    grammar_builder.build_grammar_with_gpu(&encoded_bases, seq_idx)?;
                } else {
                    if use_gpu && task_gpu_context.is_none() {
                        println!("Warning: GPU mode selected, but GPU context is unavailable for sequence {}. Falling back to CPU.", record_id_owned);
                    } else if use_gpu_override {
                        println!("INFO: GPU path forced off for debugging.");
                    }
                    grammar_builder.build_grammar(&encoded_bases, seq_idx)?;
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
                    grammar_builder.build_grammar_with_gpu(&encoded_bases_fallback, seq_idx)?;
                } else {
                    if use_gpu && task_gpu_context.is_none() {
                         println!("Warning: GPU mode selected (no encoding), but GPU context is unavailable for sequence {}. Falling back to CPU.", record_id_owned);
                    }
                    grammar_builder.build_grammar(&encoded_bases_fallback, seq_idx)?;
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
    
    Ok((merged_grammar, current_run_checkpoints, chromosome_details))
}

async fn process_streaming_mode(
    args: &OrbweaverArgs,
    gpu_context: Option<&GpuContext>,
    run_specific_output_dir: &PathBuf,
    initial_checkpoints: &HashMap<String, PathBuf>
) -> Result<(Grammar, HashMap<String, PathBuf>, Vec<(String, usize)>)> {
    println!("Using streaming mode.");

    let checkpoints_dir = run_specific_output_dir.join("checkpoints");
    if !checkpoints_dir.exists() {
        fs::create_dir_all(&checkpoints_dir)
            .with_context(|| format!("Failed to create checkpoints directory: {:?}", checkpoints_dir))?;
    }
    let mut current_run_checkpoints = initial_checkpoints.clone();

    let mut grammars: Vec<Grammar> = Vec::new();
    let mut total_bases_processed_estimate: usize = 0;

    for (file_idx, input_file_path) in args.input_files.iter().enumerate() {
        let file_identifier = format!(
            "file_{}_{}",
            file_idx,
            input_file_path.file_stem().unwrap_or_default().to_string_lossy()
        );
        let checkpoint_file_name = format!("{}_stream.grammar", file_identifier);
        let checkpoint_path = checkpoints_dir.join(&checkpoint_file_name);

        if let Some(existing_checkpoint_path) = current_run_checkpoints.get(&file_identifier) {
            if existing_checkpoint_path.exists() && !args.force_rerun {
                println!("Loading grammar for file {} ({}) from checkpoint {:?}", file_idx, input_file_path.display(), existing_checkpoint_path);
                match File::open(existing_checkpoint_path) {
                    Ok(file) => match bincode::deserialize_from(file) {
                        Ok(grammar) => { grammars.push(grammar); continue; }
                        Err(e) => println!("Warning: Failed to deserialize checkpoint {:?}: {}. Recomputing.", existing_checkpoint_path, e),
                    },
                    Err(e) => println!("Warning: Failed to open checkpoint file {:?}: {}. Recomputing.", existing_checkpoint_path, e),
                }
            }
        }

        println!("Processing file {} ({}) in streaming mode...", file_idx, input_file_path.display());

        let mut grammar_builder = GrammarBuilder::new(args.min_rule_usage, args.reverse_aware).enable_streaming_mode();
        if let Some(max_rules) = args.max_rule_count {
            grammar_builder = grammar_builder.with_max_rules(max_rules);
        }
        if !args.no_gpu && gpu_context.is_some(){
            grammar_builder = grammar_builder.with_gpu(gpu_context);
        }

        let file = File::open(input_file_path)
            .with_context(|| format!("Failed to open FASTA file: {}", input_file_path.display()))?;
        let fasta_reader = bio::io::fasta::Reader::new(std::io::BufReader::new(file));
        
        let mut record_processed_bases = 0;

        for result in fasta_reader.records() {
            let record = result.with_context(|| format!("Failed to read FASTA record from {}", input_file_path.display()))?;
            println!("Processing record: {} from file {}", record.id(), input_file_path.display());
            
            let seq_data = record.seq();
            let mut current_pos = 0;
            let stream_chunk_size = args.chunk_size_streaming.unwrap_or(1024 * 1024); // Default to 1MB if not set

            while current_pos < seq_data.len() {
                let end = std::cmp::min(current_pos + stream_chunk_size, seq_data.len());
                let chunk = &seq_data[current_pos..end];
                
                let encoded_chunk_bases: Vec<EncodedBase>;
                if args.skip_ns {
                    encoded_chunk_bases = chunk.iter().filter(|&&b| b != b'N' && b != b'n').filter_map(|&b| EncodedBase::from_base(b)).collect();
                } else {
                    encoded_chunk_bases = chunk.iter().filter_map(|&b| EncodedBase::from_base(b)).collect();
                }

                if !encoded_chunk_bases.is_empty() {
                    grammar_builder.process_sequence_chunk(&encoded_chunk_bases, file_idx, record_processed_bases)
                        .with_context(|| format!("Failed to process sequence chunk for record {} from file {}", record.id(), input_file_path.display()))?;
                    record_processed_bases += encoded_chunk_bases.len();
                }
                current_pos = end;
            }
            println!("Finished record: {}. Processed approx. {} bases for this record.", record.id(), record_processed_bases);
            total_bases_processed_estimate += record_processed_bases;
            record_processed_bases = 0; // Reset for next record in the same file
        }
        
        grammar_builder.finalize_grammar()
            .with_context(|| format!("Failed to finalize grammar for file {}", input_file_path.display()))?;
        
        let (sequence_ref, rules_ref) = grammar_builder.get_grammar();
        let max_depth = grammar_builder.get_max_rule_depth();
        let file_grammar = Grammar { 
            sequence: sequence_ref.clone(),
            rules: rules_ref.clone(),
            max_depth, 
            origins: HashMap::new() 
        };

        // Save checkpoint logic (assuming this part is mostly okay for now)
        match File::create(&checkpoint_path) {
            Ok(file_writer) => {
                if let Err(e) = bincode::serialize_into(file_writer, &file_grammar) {
                    println!("Warning: Failed to save checkpoint for file {}: {}. Proceeding without checkpoint.", file_identifier, e);
                } else {
                    println!("Saved checkpoint for file {} to {:?}", file_identifier, checkpoint_path);
                    current_run_checkpoints.insert(file_identifier.clone(), checkpoint_path.clone());
                }
            }
            Err(e) => {
                 println!("Warning: Failed to create checkpoint file {:?}: {}. Proceeding without checkpoint.", checkpoint_path, e);
            }
        }
        grammars.push(file_grammar);
    }

    if grammars.is_empty() {
        bail!("No grammars were built. Input files might be empty or failed to process.");
    }

    if grammars.len() == 1 {
        Ok((grammars.remove(0), current_run_checkpoints, Vec::new())) // Return empty vec for chromosome_details
    } else {
        println!("Merging {} grammars from streaming mode...", grammars.len());
        let dummy_config = ChunkingConfig::default();
        merge_grammars(grammars, &dummy_config, total_bases_processed_estimate)
             .map(|(g, _metrics)| (g, current_run_checkpoints, Vec::new())) // Return empty vec for chromosome_details
            .context("Failed to merge grammars from streaming mode")
    }
}

fn process_chunked_mode(
    args: &OrbweaverArgs, 
    run_specific_output_dir: &PathBuf,
    initial_checkpoints: &HashMap<String, PathBuf>
) -> Result<(Grammar, HashMap<String, PathBuf>, Vec<(String, usize)>)> {
    if args.input_files.len() > 1 {
         return Err(anyhow::anyhow!("Chunked mode currently supports only a single input FASTA file. Please provide one file or use standard/streaming mode for multiple files."));
    }
    let input_path = args.input_files.first().ok_or_else(|| anyhow::anyhow!("No input file provided for chunked mode."))?;
    println!("Processing (chunked) input file: {}", input_path.display());

    let seq_id = input_path.file_stem().unwrap_or_default().to_string_lossy().to_string();
    let checkpoint_filename = format!("{}_chunked_final.grammar", seq_id);
    let _checkpoint_path = run_specific_output_dir.join(&checkpoint_filename);

    if let Some(resumed_grammar_path) = initial_checkpoints.get(&format!("{}_final", seq_id)) { // Check for final checkpoint
        if resumed_grammar_path.exists() && !args.force_rerun {
            println!("Chunked mode: Resuming final grammar for {} from checkpoint: {}", seq_id, resumed_grammar_path.display());
            // TODO: Load the fully merged grammar.
            // let grammar = ... load from resumed_grammar_path ...
            // return Ok((grammar, initial_checkpoints.clone()));
            println!("Skipping {} due to existing final checkpoint (full chunked resume not yet implemented here).", seq_id);
             return Err(anyhow::anyhow!("Chunked mode final checkpoint exists, resume not fully implemented.")); // Placeholder
        }
    }
    
    let sequences = orbweaver::fasta::reader::read_fasta_sequences(input_path, args.skip_ns)?;
    if sequences.is_empty() {
        return Err(anyhow::anyhow!("No sequences found in FASTA file: {}", input_path.display()));
    }
    
    // For chunked mode, typically operate on the first sequence if multiple are present,
    // or concatenate them (concatenation needs careful thought for biological meaning).
    // Here, we'll use the first sequence.
    let (_header, sequence_bytes) = &sequences[0];
    
    let encoded_sequence = if args.use_encoding {
        encode_dna(sequence_bytes)
    } else {
        // This path is less common as GrammarBuilder expects EncodedBase.
        // If use_encoding is false, we'd need a different builder path or convert here.
        // For now, assume use_encoding=true is typical for builder.
        // If not, `GrammarBuilder::process_bytes_chunk` would be relevant.
        eprintln!("Warning: Chunked mode typically uses encoded sequences. Behavior with raw bytes might be limited.");
        // Convert to Vec<EncodedBase> anyway, perhaps with a lossy conversion for non-ACGT.
        sequence_bytes.iter().filter_map(|b| EncodedBase::from_base(*b)).collect()
    };

    let num_threads = args.threads.unwrap_or_else(num_cpus::get);
    let chunking_config = ChunkingConfig {
        chunk_size: args.chunk_size, // If 0, parallel::engine will use adaptive/default
        overlap_size: args.chunk_overlap,
        min_rule_usage: args.min_rule_usage,
        reverse_aware: args.reverse_aware,
        num_threads,
        show_progress: true, // Or make this an arg
        adaptive_chunking: args.adaptive_chunking,
        max_memory_per_chunk: args.max_memory_per_chunk_mb.map(|mb| mb * 1024 * 1024),
        use_gpu: !args.no_gpu, // Assuming ChunkingConfig takes this
    };

    // `parallel_sequitur` would need to be adapted to save/load chunk-level checkpoints
    // and a final merged checkpoint, using `run_specific_output_dir`.
    // It would also need to accept the `initial_checkpoints` and return `current_checkpoints`.
    // For now, this is a simplified call.
    let (grammar, _parallel_metrics) = orbweaver::parallel::engine::parallel_sequitur(
        &encoded_sequence, 
        chunking_config,
        None // Pass GpuContext here if available and `parallel_sequitur` supports it
    )?;
    
    // Save final merged grammar checkpoint
    // utils::export::write_grammar_json(&checkpoint_path, &grammar.sequence, &grammar.rules)?;
    let final_checkpoints = initial_checkpoints.clone();
    // final_checkpoints.insert(format!("{}_final", seq_id), checkpoint_path);
    // println!("Saved final chunked grammar checkpoint for {} to {}", seq_id, checkpoint_filename);

    Ok((grammar, final_checkpoints, Vec::new())) // Return empty vec for chromosome_details
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