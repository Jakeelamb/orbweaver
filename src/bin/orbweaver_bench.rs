use anyhow::{bail, Context, Result};
use chrono::Utc;
use clap::Parser;
use serde::Serialize;
use std::fs;
use std::mem::MaybeUninit;
use std::os::unix::process::ExitStatusExt;
use std::path::{Path, PathBuf};
use std::process::{Command, ExitStatus, Stdio};
use std::time::Instant;

#[derive(Parser, Debug)]
#[command(author, version, about = "Run the focused Orbweaver benchmark path and save metrics.")]
struct BenchArgs {
    #[arg(short, long)]
    input: PathBuf,

    #[arg(long, default_value = "target/release/orbweaver")]
    orbweaver_bin: PathBuf,

    #[arg(long, default_value = "benchmark_results")]
    output_dir: PathBuf,

    #[arg(long)]
    run_id: Option<String>,

    #[arg(long, default_value_t = 10_000_000)]
    mmap_chunk_size: usize,

    #[arg(long, default_value_t = 1000)]
    chunk_overlap: usize,

    #[arg(long, default_value_t = 1000)]
    streaming_output_flush_interval: usize,

    #[arg(long)]
    gpu: bool,
}

#[derive(Debug, Serialize)]
struct BenchmarkSummary {
    input: String,
    command: Vec<String>,
    exit_code: Option<i32>,
    success: bool,
    wall_time_ms: u128,
    peak_rss_kb: Option<u64>,
    run_dir: String,
    grammar_json: Artifact,
    rules_fasta: Artifact,
    motif_table: Artifact,
    graphml: GraphArtifact,
    dot: Artifact,
    gfa: Artifact,
    stats: ParsedStats,
}

#[derive(Debug, Serialize)]
struct Artifact {
    path: String,
    bytes: Option<u64>,
}

#[derive(Debug, Serialize)]
struct GraphArtifact {
    path: String,
    bytes: Option<u64>,
    nodes: Option<usize>,
    edges: Option<usize>,
}

#[derive(Debug, Default, Serialize)]
struct ParsedStats {
    number_of_rules: Option<usize>,
    compressed_sequence_symbols: Option<usize>,
    compressed_sequence_bp: Option<usize>,
    maximum_rule_depth: Option<usize>,
    total_symbols_in_rules: Option<usize>,
    motif_rows: Option<usize>,
}

fn main() -> Result<()> {
    let args = BenchArgs::parse();
    if !args.input.exists() {
        bail!("input does not exist: {}", args.input.display());
    }
    if !args.orbweaver_bin.exists() {
        bail!(
            "orbweaver binary does not exist: {}. Build it first with `cargo build --release`.",
            args.orbweaver_bin.display()
        );
    }
    let orbweaver_bin = args
        .orbweaver_bin
        .canonicalize()
        .with_context(|| format!("failed to resolve {}", args.orbweaver_bin.display()))?;

    fs::create_dir_all(&args.output_dir)
        .with_context(|| format!("failed to create {}", args.output_dir.display()))?;

    let run_id = args.run_id.clone().unwrap_or_else(|| {
        format!("bench_{}", Utc::now().format("%Y%m%d_%H%M%S"))
    });
    let run_dir = args.output_dir.join(&run_id);

    let mut command_args = vec![
        "--input-files".to_string(),
        args.input.display().to_string(),
        "--run-id".to_string(),
        run_id.clone(),
        "--mmap".to_string(),
        "--hierarchical-merge".to_string(),
        "--auto-tune-lcg".to_string(),
        "--streaming-output".to_string(),
        "--streaming-output-flush-interval".to_string(),
        args.streaming_output_flush_interval.to_string(),
        "--mmap-chunk-size".to_string(),
        args.mmap_chunk_size.to_string(),
        "--chunk-overlap".to_string(),
        args.chunk_overlap.to_string(),
        "--stats".to_string(),
    ];
    if args.gpu {
        command_args.push("--gpu".to_string());
    }

    let mut command_for_record = vec![orbweaver_bin.display().to_string()];
    command_for_record.extend(command_args.clone());

    let start = Instant::now();
    let (status, peak_rss_kb) = run_and_measure_peak_rss(&orbweaver_bin, &command_args, &args.output_dir)?;
    let wall_time_ms = start.elapsed().as_millis();
    let summary = build_summary(
        &args,
        command_for_record,
        status.code(),
        status.success(),
        wall_time_ms,
        peak_rss_kb,
        &run_dir,
    )?;
    write_summary(&run_dir, &summary)?;
    if !status.success() {
        bail!("orbweaver benchmark command failed with status {}", status);
    }
    Ok(())
}

fn build_summary(
    args: &BenchArgs,
    command: Vec<String>,
    exit_code: Option<i32>,
    success: bool,
    wall_time_ms: u128,
    peak_rss_kb: Option<u64>,
    run_dir: &Path,
) -> Result<BenchmarkSummary> {
    let grammar_stats = run_dir.join("grammar_stats.txt");
    let motif_table = run_dir.join("motif_table.tsv");
    let graphml = run_dir.join("grammar.graphml");
    let mut stats = parse_stats(&grammar_stats).unwrap_or_default();
    stats.motif_rows = count_tsv_rows(&motif_table);

    let (nodes, edges) = parse_graphml_size(&graphml);

    Ok(BenchmarkSummary {
        input: args.input.display().to_string(),
        command,
        exit_code,
        success,
        wall_time_ms,
        peak_rss_kb,
        run_dir: run_dir.display().to_string(),
        grammar_json: artifact(run_dir.join("grammar.json")),
        rules_fasta: artifact(run_dir.join("rules.fasta")),
        motif_table: artifact(motif_table),
        graphml: GraphArtifact {
            path: graphml.display().to_string(),
            bytes: file_len(&graphml),
            nodes,
            edges,
        },
        dot: artifact(run_dir.join("grammar.dot")),
        gfa: artifact(run_dir.join("grammar.gfa")),
        stats,
    })
}

fn write_summary(run_dir: &Path, summary: &BenchmarkSummary) -> Result<()> {
    fs::create_dir_all(run_dir)
        .with_context(|| format!("failed to create run dir {}", run_dir.display()))?;

    let json = serde_json::to_string_pretty(summary)?;
    fs::write(run_dir.join("benchmark_summary.json"), json)?;

    let tsv = format!(
        "input\trun_dir\tsuccess\texit_code\twall_time_ms\tpeak_rss_kb\tnumber_of_rules\tcompressed_sequence_symbols\tcompressed_sequence_bp\tmaximum_rule_depth\tmotif_rows\tgraphml_nodes\tgraphml_edges\n{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
        summary.input,
        summary.run_dir,
        summary.success,
        summary.exit_code.map(|v| v.to_string()).unwrap_or_else(|| "NA".to_string()),
        summary.wall_time_ms,
        summary.peak_rss_kb.map(|v| v.to_string()).unwrap_or_else(|| "NA".to_string()),
        opt_usize(summary.stats.number_of_rules),
        opt_usize(summary.stats.compressed_sequence_symbols),
        opt_usize(summary.stats.compressed_sequence_bp),
        opt_usize(summary.stats.maximum_rule_depth),
        opt_usize(summary.stats.motif_rows),
        opt_usize(summary.graphml.nodes),
        opt_usize(summary.graphml.edges),
    );
    fs::write(run_dir.join("benchmark_summary.tsv"), tsv)?;
    Ok(())
}

fn artifact(path: PathBuf) -> Artifact {
    Artifact {
        path: path.display().to_string(),
        bytes: file_len(&path),
    }
}

fn file_len(path: &Path) -> Option<u64> {
    fs::metadata(path).ok().map(|m| m.len())
}

fn opt_usize(value: Option<usize>) -> String {
    value.map(|v| v.to_string()).unwrap_or_else(|| "NA".to_string())
}

fn parse_stats(path: &Path) -> Result<ParsedStats> {
    let contents = fs::read_to_string(path)
        .with_context(|| format!("failed to read {}", path.display()))?;
    let mut stats = ParsedStats::default();
    for line in contents.lines() {
        if let Some(value) = parse_stat_line(line, "Number of Rules:") {
            stats.number_of_rules = Some(value);
        } else if let Some(value) = parse_stat_line(line, "Compressed Sequence Length (Symbols):") {
            stats.compressed_sequence_symbols = Some(value);
        } else if let Some(value) = parse_stat_line(line, "Compressed Sequence Length (Base Pairs):") {
            stats.compressed_sequence_bp = Some(value);
        } else if let Some(value) = parse_stat_line(line, "Maximum Rule Depth:") {
            stats.maximum_rule_depth = Some(value);
        } else if let Some(value) = parse_stat_line(line, "Total Symbols in Rule Definitions:") {
            stats.total_symbols_in_rules = Some(value);
        }
    }
    Ok(stats)
}

fn parse_stat_line(line: &str, key: &str) -> Option<usize> {
    line.trim()
        .strip_prefix(key)?
        .trim()
        .parse::<usize>()
        .ok()
}

fn count_tsv_rows(path: &Path) -> Option<usize> {
    let contents = fs::read_to_string(path).ok()?;
    Some(contents.lines().skip(1).filter(|line| !line.trim().is_empty()).count())
}

fn parse_graphml_size(path: &Path) -> (Option<usize>, Option<usize>) {
    let contents = match fs::read_to_string(path) {
        Ok(contents) => contents,
        Err(_) => return (None, None),
    };
    let nodes = contents.match_indices("<node ").count();
    let edges = contents.match_indices("<edge ").count();
    (Some(nodes), Some(edges))
}

fn run_and_measure_peak_rss(
    orbweaver_bin: &Path,
    command_args: &[String],
    output_dir: &Path,
) -> Result<(ExitStatus, Option<u64>)> {
    let child = Command::new(orbweaver_bin)
        .args(command_args)
        .current_dir(output_dir)
        .stdout(Stdio::inherit())
        .stderr(Stdio::inherit())
        .spawn()
        .with_context(|| format!("failed to launch {}", orbweaver_bin.display()))?;

    wait_for_child_with_usage(child.id())
}

fn wait_for_child_with_usage(pid: u32) -> Result<(ExitStatus, Option<u64>)> {
    let mut raw_status = 0;
    let mut usage = MaybeUninit::<RUsage>::zeroed();
    let waited = unsafe { wait4(pid as i32, &mut raw_status, 0, usage.as_mut_ptr()) };
    if waited < 0 {
        return Err(std::io::Error::last_os_error()).context("wait4 failed");
    }
    let usage = unsafe { usage.assume_init() };
    let peak_rss_kb = if usage.ru_maxrss > 0 {
        Some(usage.ru_maxrss as u64)
    } else {
        None
    };
    Ok((ExitStatus::from_raw(raw_status), peak_rss_kb))
}

#[repr(C)]
#[derive(Clone, Copy, Debug)]
struct TimeVal {
    tv_sec: i64,
    tv_usec: i64,
}

#[repr(C)]
#[derive(Clone, Copy, Debug)]
struct RUsage {
    ru_utime: TimeVal,
    ru_stime: TimeVal,
    ru_maxrss: i64,
    ru_ixrss: i64,
    ru_idrss: i64,
    ru_isrss: i64,
    ru_minflt: i64,
    ru_majflt: i64,
    ru_nswap: i64,
    ru_inblock: i64,
    ru_oublock: i64,
    ru_msgsnd: i64,
    ru_msgrcv: i64,
    ru_nsignals: i64,
    ru_nvcsw: i64,
    ru_nivcsw: i64,
}

extern "C" {
    fn wait4(pid: i32, status: *mut i32, options: i32, rusage: *mut RUsage) -> i32;
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_stat_lines() {
        assert_eq!(parse_stat_line("Number of Rules: 42", "Number of Rules:"), Some(42));
        assert_eq!(parse_stat_line("Number of Rules: nope", "Number of Rules:"), None);
    }
}
