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
use orbweaver::analysis::stats::calculate_and_print_stats;
use orbweaver::utils::visualization::DotOptions;
use orbweaver::parallel::chunking::ChunkingConfig;
use orbweaver::parallel::engine::{parallel_sequitur, merge_grammars};
use rayon::prelude::*;
use orbweaver::utils;
use sysinfo::{System, SystemExt};
use std::sync::Arc;
use orbweaver::gpu::GpuContext;
use std::collections::HashMap;
use std::fs::File;
use bincode;
use log::{info, warn};
use orbweaver::utils::export::{write_grammar_text, write_grammar_json as io_write_grammar_json, write_grammar_graphml};
use orbweaver::utils::visualization::{write_grammar_dot, write_grammar_gfa};
use std::io::BufWriter;

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

    /// Output GraphML representation of the grammar.
    /// 
    /// Exports the grammar as a graph in GraphML format,
    /// suitable for import into various graph analysis tools.
    #[clap(long, value_parser)]
    output_graphml: Option<PathBuf>,

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
    #[clap(long, value_parser, default_value_t = 10)]
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
    graph_show_depth: bool,
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
    // --- Argument Parsing ---
    let mut current_args = OrbweaverArgs::parse();

    // --- Derive species_id from input file ---
    // This needs to happen early, but run_specific_output_dir for *this* run uses it.
    // For resume, the resume_run_dir is explicit.
    let derived_species_id_for_new_run = if let Some(first_input_file) = current_args.input_files.first() {
        first_input_file.file_stem()
            .and_then(|stem| stem.to_str())
            .map(|s| s.to_string())
            .unwrap_or_else(|| {
                eprintln!("Warning: Could not derive species ID from input file {:?}. Using 'unknown_species'.", first_input_file);
                "unknown_species".to_string()
            })
    } else {
        eprintln!("Error: No input files provided.");
        return Err(anyhow::anyhow!("No input files provided to derive species ID for new run."));
    };
    // We print this later, after we know if we are resuming or not.

    // --- Load original arguments if resuming ---
    let original_args_for_resume: OrbweaverArgs = if current_args.resume {
        if let Some(resume_dir_path) = &current_args.resume_run_dir { // resume_dir_path is the user-provided absolute path
            let metadata_path = resume_dir_path.join("run_metadata.json");
            if metadata_path.exists() {
                let metadata_str = fs::read_to_string(&metadata_path)
                    .context(format!("Failed to read metadata from resumed run: {}", metadata_path.display()))?;
                let loaded_metadata: RunMetadata = serde_json::from_str(&metadata_str)
                    .context(format!("Failed to parse metadata from resumed run: {}", metadata_path.display()))?;
                loaded_metadata.args
            } else {
                // Construct the path string for the warning message here, as metadata_path is out of scope
                let non_existent_metadata_path_display = resume_dir_path.join("run_metadata.json").display().to_string();
                println!("Warning: Could not find run_metadata.json in resume directory: {}. Defaulting behavior for outputs may apply if not specified in current command.", non_existent_metadata_path_display);
                current_args.clone()
            }
        } else {
            // This case should be prevented by clap's `requires` attribute.
            println!("Warning: --resume was set but --resume-run-dir was not. This is unexpected.");
            current_args.clone()
        }
    } else {
        current_args.clone()
    };

    // --- Determine run-specific output directory for *this* execution ---
    // If resuming, actual_run_id should ideally come from the resumed run's metadata to ensure consistency,
    // unless --force-rerun or some other overriding condition is met.
    // For simplicity now, if current_args.run_id is given, it's used, otherwise a new timestamp.
    // This means resuming a run and giving a *new* run_id will create a *new* output dir for the resumed data.
    let actual_run_id = current_args.run_id.clone().unwrap_or_else(|| Utc::now().format("%Y%m%d_%H%M%S").to_string());
    
    // The species_id for the output path is based on the *current* command's input files, even if resuming.
    // This seems logical as the resumed run might be on a *different* (but compatible) input.
    println!("Derived species ID for this run's output: {}", derived_species_id_for_new_run);
    let base_output_dir = PathBuf::from(".").join(&derived_species_id_for_new_run);
    let run_specific_output_dir = base_output_dir.join(&actual_run_id);

    if !run_specific_output_dir.exists() {
        fs::create_dir_all(&run_specific_output_dir)
            .context(format!("Failed to create run-specific output directory: {}", run_specific_output_dir.display()))?;
    }
    println!("Run-specific output directory: {}", run_specific_output_dir.display());

    // --- Default Output File Paths ---
    // These defaults are applied if the user doesn't specify a path,
    // and intelligently handles resumes (respects original run's choices for specific file outputs).

    // Default .json output
    if current_args.output_json.is_none() {
        let mut skip_defaulting = false;
        if current_args.resume {
            // Check if the original run *explicitly specified* an output_json path.
            // If original_args_for_resume.output_json was Some, it means the user set it.
            // We don't want to override that with a new default if they are resuming and *didn't* specify one now.
            // However, if original_args_for_resume.output_json was None, it means it wasn't set in the original run.
            // In that case, if current_args.output_json is *also* None, it's okay to apply the new default.
            if original_args_for_resume.output_json.is_some() && current_args.output_json.is_none() {
                 // User specified it in original, not in current resume command: respect original.
                 // Effectively, we "copy" the original path if it existed and none is given now.
                 // But if we are *here*, current_args.output_json IS None.
                 // The question is if the *original* had a specific one, or if it also defaulted.
                 // The current logic: if current is None, default it, UNLESS (resume AND original_had_one)
                 // This needs to be: if current_args.output_json is None, default it,
                 // UNLESS (current_args.resume AND original_args_for_resume.output_json was Some(...))
                 // The current skip_defaulting logic seems okay.
            }
        }
        if !skip_defaulting { // skip_defaulting is effectively true if (resume AND original_args_for_resume.output_json.is_some())
            current_args.output_json = Some(run_specific_output_dir.join("grammar.json"));
        }
    }

    // Default .fasta (blocks) output
    if current_args.export_blocks.is_none() {
        let mut skip_defaulting = false;
        if current_args.resume {
            if original_args_for_resume.export_blocks.is_some() {
                skip_defaulting = true;
            }
        }
        if !skip_defaulting {
            current_args.export_blocks = Some(run_specific_output_dir.join("rules.fasta"));
        }
    }

    // Default .tsv (repeats summary) output
    if current_args.output_repeats.is_none() {
        let mut skip_defaulting = false;
        if current_args.resume {
            if original_args_for_resume.output_repeats.is_some() {
                skip_defaulting = true;
            }
        }
        if !skip_defaulting {
            current_args.output_repeats = Some(run_specific_output_dir.join("repeats_summary.tsv"));
        }
    }

    // Default .dot output
    if current_args.visualize.is_none() {
        let mut skip_defaulting = false;
        if current_args.resume {
            if original_args_for_resume.visualize.is_some() {
                 skip_defaulting = true;
            }
        }
        if !skip_defaulting {
            current_args.visualize = Some(run_specific_output_dir.join("grammar.dot"));
        }
    }

    // Default .graphml output
    if current_args.output_graphml.is_none() {
        let mut skip_defaulting = false;
        if current_args.resume {
            if original_args_for_resume.output_graphml.is_some() {
                 skip_defaulting = true;
            }
        }
        if !skip_defaulting {
            current_args.output_graphml = Some(run_specific_output_dir.join("grammar.graphml"));
        }
    }

    // --- Save initial/updated metadata for this run attempt ---
    let run_start_time = Utc::now();
    let mut metadata_to_save = RunMetadata {
        args: current_args.clone(), // This will now store args without species_id, assembly_id, output_dir
        tool_version: env!("CARGO_PKG_VERSION").to_string(),
        status: "starting".to_string(), // Will be updated to "completed" or "failed"
        start_time: run_start_time,
        end_time: None,
        run_id: actual_run_id.clone(), // actual_run_id is already determined
        processed_sequence_checkpoints: HashMap::new(), // Use initial checkpoints from potential resume
    };

    // The rest of the `main` function will now use `current_args` for all processing decisions
    // and `run_specific_output_dir` for all output locations.
    // The `metadata_to_save` variable now holds the metadata struct that needs to be updated upon completion or failure.
    
    // --- Main Processing Logic --- 
    let result: Result<Grammar> = {
        // Conditional profiling setup
        #[cfg(feature = "profiling")]
        let _guard = if current_args.profile {
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
        if let Some(threads) = current_args.threads {
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

        if current_args.chunk_size > 0 || current_args.adaptive_chunking {
            println!("Processing in chunked mode.");
            let chunked_result_tuple = process_chunked_mode(&current_args, &run_specific_output_dir, &metadata_to_save.processed_sequence_checkpoints);
            match chunked_result_tuple {
                Ok((g, cp)) => {
                    final_grammar = g;
                    temp_checkpoints = Some(cp);
                }
                Err(e) => return Err(e), 
            }
        } else if current_args.streaming {
            println!("Processing in streaming mode.");
            // TODO: Implement checkpointing for streaming mode with GrammarBuilder
            // For now, it will behave like standard mode for checkpointing demonstration
            grammar_result_tuple = process_standard_mode(&current_args, !current_args.no_gpu, None, &run_specific_output_dir, &metadata_to_save.processed_sequence_checkpoints).await;
            match grammar_result_tuple {
                Ok((g, cp)) => {
                    final_grammar = g;
                    temp_checkpoints = Some(cp);
                }
                Err(e) => return Err(e),
            }
        } else {
            println!("Processing in standard (in-memory) mode.");
            grammar_result_tuple = process_standard_mode(&current_args, !current_args.no_gpu, None, &run_specific_output_dir, &metadata_to_save.processed_sequence_checkpoints).await;
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
            metadata_to_save.processed_sequence_checkpoints = checkpoints;
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
            // Use `current_args` (which is `original_args_for_resume`) for stats and outputs
            if current_args.stats {
                calculate_and_print_stats(&grammar)?;
            }
            // Pass `current_args` (which is `original_args_for_resume`) to generate_outputs
            generate_outputs(&grammar, &current_args)?;
            
            metadata_to_save.status = "completed".to_string();
            metadata_to_save.end_time = Some(Utc::now());
            let metadata_path = run_specific_output_dir.join("run_metadata.json");
            fs::write(
                &metadata_path,
                serde_json::to_string_pretty(&metadata_to_save)?
            ).with_context(|| format!("Failed to write final metadata to {:?}", metadata_path))?;
            println!("Run {} completed. Metadata updated.", metadata_to_save.run_id);
        }
        Err(e) => {
            eprintln!("Error during grammar construction: {:?}", e);
            metadata_to_save.status = "failed".to_string();
            metadata_to_save.end_time = Some(Utc::now());
            let metadata_path = run_specific_output_dir.join("run_metadata.json");
            fs::write(
                &metadata_path,
                serde_json::to_string_pretty(&metadata_to_save)?
            ).with_context(|| format!("Failed to write error metadata to {:?}", metadata_path))?;
            println!("Run {} failed. Metadata updated.", metadata_to_save.run_id);
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
            // Use `actual_run_id` for the flamegraph name to ensure consistency
            let flamegraph_path = profile_dir.join(format!("flamegraph_{}.svg", actual_run_id));
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
    info!("Generating requested output files...");

    if let Some(path) = &args.output_json {
        info!("Writing grammar to JSON file: {}", path.display());
        io_write_grammar_json(path, &grammar.sequence, &grammar.rules)
            .with_context(|| format!("Failed to create JSON grammar file at {:?}", path))?;
    }

    if let Some(path) = &args.output_text {
        info!("Writing grammar to text file: {}", path.display());
        write_grammar_text(path, &grammar.sequence, &grammar.rules)
            .with_context(|| format!("Failed to create text grammar file at {:?}", path))?;
    }

    if let Some(path) = &args.output_gfa {
        info!("Writing grammar to GFA file: {}", path.display());
        write_grammar_gfa(path, grammar)
            .with_context(|| format!("Failed to create GFA file at {:?}", path))?;
    }

    // Common DotOptions for GraphML and DOT, derived from args.graph_*
    let common_dot_options = DotOptions {
        include_terminals: args.graph_include_terminals,
        max_depth: args.graph_max_depth,
        skip_rules_above_depth: args.graph_skip_rules_above_depth,
        transparent_background: args.graph_transparent_background,
        dark_mode: args.graph_dark_mode,
        include_usage_counts: args.graph_show_usage_counts,
        color_by_depth: args.graph_show_depth,
        engine: args.graph_engine.clone(),
    };

    if let Some(graphml_path) = &args.output_graphml {
        info!("Writing grammar to GraphML file: {}", graphml_path.display());
        write_grammar_graphml(graphml_path, grammar, &common_dot_options)
            .with_context(|| format!("Failed to create GraphML file at {:?}", graphml_path))?;
    }

    if let Some(dot_path) = &args.visualize {
        info!("Writing grammar visualization to DOT file: {}", dot_path.display());
        write_grammar_dot(dot_path, grammar, &common_dot_options)
            .with_context(|| format!("Failed to create DOT file at {:?}", dot_path))?;
    }

    if args.stats {
        info!("Calculating and printing statistics...");
        calculate_and_print_stats(grammar)?;
    }

    if let Some(path) = &args.export_blocks {
         info!("Exporting grammar blocks to FASTA file: {}", path.display());
         let file = File::create(path)
             .with_context(|| format!("Failed to create FASTA export file at {:?}", path))?;
         let mut writer = BufWriter::new(file);
         utils::export::export_grammar_fasta(grammar, &mut writer)
             .with_context(|| format!("Failed to export grammar blocks to FASTA at {:?}", path))?;
    }
    
    if let Some(path) = &args.output_repeats {
        info!("Writing repeats summary to: {}", path.display());
        warn!("Repeat summary output is not fully implemented yet for path: {}", path.display());
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