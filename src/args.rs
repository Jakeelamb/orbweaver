use clap::Parser;
use serde::{Serialize, Deserialize};
use std::path::PathBuf;

/// Orbweaver: A grammar-based approach to genomic sequence analysis.
/// 
/// Orbweaver processes DNA sequences from FASTA files to build context-free grammars
/// by identifying repeating patterns. It can output the grammar in various formats and
/// provide statistics about the compression and structure.
#[derive(Parser, Debug, Clone, Serialize, Deserialize, Default)]
#[command(author = "Orbweaver Team", version, about = "Orbweaver: Grammar-based genomic sequence analysis tool.\n\nOutputs are saved in ./<species_id>/<run_id>/ where <species_id> is derived from the input filename and <run_id> is a timestamp or user-provided.", long_about = None)]
#[command(help_template = "\n{before-help}{name} {version}\n{author-with-newline}{about-with-newline}\n{usage-heading} {usage}\n\n{all-args}{after-help}\n")]
pub struct OrbweaverArgs {
    /// Input FASTA file paths (.fa, .fasta, .fna), comma-separated.
    ///
    /// Paths to files containing DNA sequences in FASTA format.
    /// Processes all sequences in multi-sequence files by default.
    #[clap(short, long, value_parser, required = true, value_delimiter = ',')]
    pub input_files: Vec<PathBuf>,

    // --- Workflow and Output Organization ---
    /// Base directory for all outputs. (Will be derived if not provided if paths are relative, or used as absolute if given)
    #[clap(short = 'o', long, value_parser)]
    pub output_dir: Option<PathBuf>,

    /// Species identifier. (Will be derived from the first input file if not provided)
    #[clap(long, value_parser)]
    pub species_id: Option<String>,

    /// Assembly identifier. (Optional)
    #[clap(long, value_parser)]
    pub assembly_id: Option<String>,

    // --- Output Format Options ---
    
    /// Output JSON file path for the grammar.
    /// 
    /// Writes a detailed JSON representation of the grammar including all
    /// rules and the final compressed sequence.
    #[clap(short = 'j', long, value_parser)]
    pub output_json: Option<PathBuf>,

    /// Output human-readable text representation of the grammar.
    /// 
    /// Writes a simplified text format showing rules and the final sequence
    /// in a more readable format than JSON.
    #[clap(long, value_parser)]
    pub output_text: Option<PathBuf>,

    /// Output GFAv1 representation of the grammar.
    /// 
    /// Exports the grammar as a graph in GFA (Graphical Fragment Assembly) format,
    /// which can be visualized with tools like Bandage.
    #[clap(long, value_parser)]
    pub output_gfa: Option<PathBuf>,

    /// Generate a .dot file for visualizing the grammar.
    /// 
    /// Creates a DOT file for visualization with Graphviz tools like dot, neato, etc.
    /// Example usage: dot -Tpng grammar.dot -o grammar.png
    #[clap(long, value_parser)]
    pub visualize: Option<PathBuf>,

    /// Export grammar rules as sequences in FASTA format.
    /// 
    /// Writes each grammar rule as a separate FASTA record, where each record
    /// contains the fully expanded DNA sequence for that rule.
    #[clap(long, value_parser)]
    pub export_blocks: Option<PathBuf>,

    /// Output tabular summary of repeats.
    ///
    /// Writes a text file summarizing identified repeats, sorted by size or frequency.
    #[clap(long, value_parser)]
    pub output_repeats: Option<PathBuf>,

    /// Output GraphML representation of the grammar.
    /// 
    /// Exports the grammar as a graph in GraphML format,
    /// suitable for import into various graph analysis tools.
    #[clap(long, value_parser)]
    pub output_graphml: Option<PathBuf>,

    /// Print statistics about the generated grammar.
    /// 
    /// Displays metrics about the grammar: rule count, depth, compression ratio, etc.
    #[clap(long, action = clap::ArgAction::SetTrue)]
    pub stats: bool,

    // --- Grammar Construction Options ---

    /// K-mer size for processing.
    /// 
    /// Used in some analysis algorithms. The grammar builder currently uses
    /// digrams (k=2) internally regardless of this setting.
    #[clap(short, long, value_parser, default_value_t = 21)]
    pub kmer_size: usize,

    /// Minimum usage count for a rule to be kept.
    /// 
    /// A digram must appear at least this many times to be replaced by a rule.
    /// Higher values lead to fewer rules with more usage each.
    #[clap(long, value_parser, default_value_t = 10)]
    pub min_rule_usage: usize,

    /// Maximum number of rules allowed (triggers eviction).
    /// 
    /// When the number of rules exceeds this limit, less frequently used rules
    /// are evicted and inlined. If not specified, no limit is enforced.
    #[clap(long, value_parser)]
    pub max_rule_count: Option<usize>,

    /// Perform reverse complement canonicalization.
    /// 
    /// When enabled, a digram and its reverse complement are treated as equivalent.
    /// For example, "AC" and "GT" (reverse complement) are considered the same pattern.
    #[clap(long, default_value_t = true, action = clap::ArgAction::Set)]
    pub reverse_aware: bool,

    // --- Input Processing Options ---

    /// Skip N bases in FASTA input.
    /// 
    /// When enabled, 'N' or 'n' characters in the input are skipped,
    /// as they typically represent unknown bases.
    #[clap(long, default_value_t = true, action = clap::ArgAction::Set)]
    pub skip_ns: bool,

    /// Process genome in chunks of this size.
    /// 
    /// Divides large sequences into smaller chunks for processing.
    /// Useful for very large genomes that exceed available memory.
    /// If set to 0, will determine size dynamically.
    #[clap(long, value_parser, default_value_t = 0)]
    pub chunk_size: usize,

    /// Overlap between chunks when chunking is enabled.
    /// 
    /// Specifies the number of bases that overlap between adjacent chunks.
    /// Helps ensure patterns spanning chunk boundaries are detected.
    /// Note: Only relevant when chunking is enabled.
    #[clap(long, value_parser, default_value_t = 1000)]
    pub chunk_overlap: usize,
    
    /// Number of threads to use for parallel processing.
    /// 
    /// Sets the number of worker threads for parallelized operations.
    /// Default: number of logical CPU cores.
    #[clap(long, value_parser)]
    pub threads: Option<usize>,

    /// Enable performance profiling.
    /// 
    /// Generates a CPU profile and flame graph of the application execution.
    /// Results are saved to the 'profile' directory.
    #[clap(long, action = clap::ArgAction::SetTrue)]
    pub profile: bool,

    /// Enable 2-bit encoding to reduce memory usage.
    /// 
    /// Use 2-bit encoding for DNA sequences (A=00, C=01, G=10, T=11),
    /// reducing memory usage by up to 75% compared to ASCII storage.
    #[clap(long, default_value_t = true, action = clap::ArgAction::Set)]
    pub use_encoding: bool,
    
    /// Enable streaming mode for low memory usage.
    /// 
    /// Process the FASTA file in a streaming fashion, without loading
    /// the entire sequence into memory at once. Recommended for large genomes.
    #[clap(long, action = clap::ArgAction::SetTrue)]
    pub streaming: bool,
    
    /// Chunk size for FastaStream when streaming mode is enabled (bytes).
    /// Default: 1MB (1024 * 1024).
    #[clap(long, value_parser)]
    pub chunk_size_streaming: Option<usize>,
    
    /// Enable adaptive chunk sizing.
    /// 
    /// Dynamically adjust chunk sizes based on sequence complexity
    /// and available memory. Requires chunked mode.
    #[clap(long, action = clap::ArgAction::SetTrue)]
    pub adaptive_chunking: bool,
    
    /// Process only specific sequences by index (0-based).
    /// 
    /// For multi-sequence FASTA files, process only the sequences
    /// with the specified indices. Default: process all sequences.
    #[clap(long, value_parser, value_delimiter = ',')]
    pub sequence_indices: Option<Vec<usize>>,
    
    /// Maximum memory usage per chunk (in MB).
    /// 
    /// Limits the memory usage of each chunk when adaptive chunking
    /// is enabled. Helps prevent out-of-memory errors.
    #[clap(long, value_parser)]
    pub max_memory_per_chunk_mb: Option<usize>,

    /// Disable GPU acceleration (CPU fallback).
    ///
    /// GPU acceleration is used by default. This flag forces CPU-only processing.
    #[clap(long, action = clap::ArgAction::SetTrue)]
    pub no_gpu: bool,

    /// Output JSON in compact format (no pretty-printing).
    #[clap(long, action = clap::ArgAction::SetTrue)]
    pub compact_json: bool,

    // --- Resume/Checkpoint Options ---
    /// Unique identifier for a specific run. 
    /// If provided, Orbweaver will attempt to resume this run if it exists in the `checkpoint_dir`.
    /// If not provided or if the ID is new, a new run (with a generated UUID as ID if none was given) will be started.
    /// It is recommended to use descriptive IDs (e.g., `species_condition_kmer21`) for better run management.
    #[clap(long, value_parser)]
    pub run_id: Option<String>,

    /// Directory to store run manifests, checkpoints, and outputs.
    /// Each run will create a subdirectory within this path, named by its `run_id`.
    /// Defaults to `./orbweaver_runs` if not specified.
    #[clap(long, value_parser, default_value = "./orbweaver_runs")]
    pub checkpoint_dir: PathBuf,

    /// Configures periodic checkpointing. Format: "<N>s" for seconds or "<N>steps" for grammar construction steps.
    /// Examples: `3600s` (every hour), `1000steps` (every 1000 rule creations/replacements).
    /// If omitted but `checkpoint_dir` is specified, a final checkpoint is typically saved upon successful completion.
    #[clap(long, value_parser)]
    pub checkpoint_interval: Option<String>,

    /// Graph-related arguments
    #[clap(long, value_parser, default_value_t = String::from("sfdp"))]
    pub graph_engine: String,

    /// Include terminal symbols in graph visualizations.
    #[clap(long, value_parser, default_value_t = true)]
    pub graph_include_terminals: bool,

    /// Maximum depth of rules to display in graph visualizations.
    #[clap(long, value_parser)]
    pub graph_max_depth: Option<usize>,

    /// Skip rules above this depth in graph visualizations (useful for very deep grammars).
    #[clap(long, value_parser)]
    pub graph_skip_rules_above_depth: Option<usize>,

    /// Use a transparent background for graph visualizations (PNG/SVG).
    #[clap(long, value_parser, default_value_t = false)]
    pub graph_transparent_background: bool,

    /// Use dark mode styling for graph visualizations.
    #[clap(long, value_parser, default_value_t = false)]
    pub graph_dark_mode: bool,

    /// Show usage counts on nodes in graph visualizations.
    #[clap(long, value_parser, default_value_t = true)]
    pub graph_show_usage_counts: bool,

    /// Color nodes by depth and show depth information in graph visualizations.
    #[clap(long, value_parser, default_value_t = true)]
    pub graph_show_depth: bool,
} 