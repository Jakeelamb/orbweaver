use anyhow::{Context, Result};
use clap::{Parser, Subcommand};
use std::path::PathBuf;

// Use the modules declared in src/lib.rs
use orbweaver::ncbi;
use orbweaver::encoding;
use orbweaver::slp;
use orbweaver::motifs;

/// Orbweaver: A tool for processing genomic data, including NCBI fetching, encoding, and SLP construction.
#[derive(Parser, Debug)]
#[clap(author, version, about, long_about = None)]
struct Cli {
    #[clap(subcommand)]
    command: Commands,
}

#[derive(Subcommand, Debug)]
enum Commands {
    /// Fetch data from NCBI, process, and optionally download sample genomes.
    FetchNcbi {
        #[clap(long, help = "Force fetch data from NCBI and overwrite existing CSV files.")]
        generate: bool,

        #[clap(long, value_parser, default_value = "metazoan_chromosome_assemblies_with_lineage.csv", help = "Path to the main CSV file for loading/saving assembly data.")]
        csv_file: PathBuf,
    },

    /// Encode a FASTA file into VBQ format.
    EncodeVbq {
        #[clap(value_parser, help = "Input FASTA file path (.fa, .fasta, .fna)")]
        input: PathBuf,

        #[clap(value_parser, help = "Output VBQ file path (.vbq)")]
        output: PathBuf,

        #[clap(short, long, value_parser, default_value_t = 0, help = "Compression level (0-9)")]
        level: u8,
    },

    /// Encode a FASTA file into the custom Orbweaver binary format (.orb).
    EncodeOrb {
        #[clap(value_parser, help = "Input FASTA file path (.fa, .fasta, .fna)")]
        input: PathBuf,

        #[clap(value_parser, help = "Output ORB file path (.orb)")]
        output: PathBuf,
    },

    /// Build a Straight-Line Program (SLP) from a FASTA file (first sequence).
    BuildSlp {
        #[clap(short, long, value_parser, help = "Input FASTA file path.")]
        input: PathBuf,

        #[clap(short, long, value_parser, default_value_t = 1_000_000.0, help = "Normalization factor N for SLP stopping condition (freq < 1 / (len * N))")]
        n_factor: f64,

        // Placeholder for potential output file argument
        // #[clap(short, long, value_parser, help = "Output file for the generated SLP rules (optional).")]
        // output: Option<PathBuf>,
    },

    /// Find repeating substrings in a custom Orbweaver binary file (.orb).
    FindRepeats {
        #[clap(short, long, value_parser, help = "Input ORB file path (.orb).")]
        input: PathBuf,

        #[clap(short, long, value_parser, default_value_t = 10, help = "Minimum length of repeats to report.")]
        min_len: usize,
    },
}

fn main() -> Result<()> {
    let cli = Cli::parse();

    // Match on the subcommand
    match cli.command {
        Commands::FetchNcbi { generate, csv_file } => {
            println!("Running NCBI Fetch command...");
            ncbi::run_ncbi_fetch(generate, csv_file)
                .context("NCBI Fetch command failed")?
        }
        Commands::EncodeVbq { input, output, level } => {
            println!("Running Encode VBQ command...");
            encoding::run_fasta_to_vbq(input, output, level)
                .context("Encode VBQ command failed")?
        }
        Commands::EncodeOrb { input, output } => {
            println!("Running Encode ORB command...");
            encoding::run_fasta_to_orb_bin(input, output)
                .context("Encode ORB command failed")?
        }
        Commands::BuildSlp { input, n_factor } => {
            println!("Running Build SLP command...");
            slp::run_slp_build(input, n_factor)
                .context("Build SLP command failed")?
            // TODO: Add logic to handle optional SLP output file
        }
        Commands::FindRepeats { input, min_len } => {
            println!("Running Find Repeats command (on .orb file)...");

            // Calls the modified motifs::find_repeats which expects .orb format
            let repeats = motifs::find_repeats(&input, min_len)
                .context("Failed to find repeats")?;

            // Print the results (potentially large output)
            println!("\nFound {} unique repeats (length >= {}):", repeats.len(), min_len);
             let mut sorted_repeats: Vec<_> = repeats.into_iter().collect();
             sorted_repeats.sort_by(|a, b| b.1.cmp(&a.1).then_with(|| a.0.cmp(&b.0)));

            let max_to_print = 50; // Limit output
            for (repeat, count) in sorted_repeats.iter().take(max_to_print) {
                 println!("  - \"{}\" : {}", repeat, count);
            }
            if sorted_repeats.len() > max_to_print {
                println!("  ... (output truncated, {} repeats total)", sorted_repeats.len());
            }
        }
    }

    Ok(())
}
