use serde::{Serialize, Deserialize};
use std::path::PathBuf;
use crate::args::OrbweaverArgs; // For storing original arguments
use chrono::{DateTime, Utc}; // For timestamps

/// Represents the possible states of an Orbweaver run.
#[derive(Serialize, Deserialize, Debug, Clone, PartialEq, Eq)]
pub enum RunStatus {
    /// The run is currently being set up.
    Initializing,
    /// The run is actively processing.
    Running,
    /// The run has saved a checkpoint and can be resumed.
    Checkpointed,
    /// The run has finished all processing successfully.
    Completed,
    /// The run encountered an error and stopped prematurely.
    Failed,
}

/// Stores metadata and state for a single Orbweaver run.
/// This information is saved to `manifest.json` within the run's dedicated output directory.
#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct RunManifest {
    /// Unique identifier for this run.
    pub run_id: String,
    /// The original command-line arguments used to start this run.
    pub original_args: OrbweaverArgs,
    /// The current execution status of the run.
    pub status: RunStatus,
    /// Path to the last successfully saved checkpoint file, relative to the `output_dir_manifest`.
    pub last_checkpoint_file: Option<PathBuf>,
    /// Timestamp of when the last checkpoint was saved.
    pub last_checkpoint_time: Option<DateTime<Utc>>,
    /// Timestamp of when the run was initially started.
    pub run_start_time: DateTime<Utc>,
    /// Absolute path to this run's dedicated output directory where the manifest and all run artifacts are stored.
    pub output_dir_manifest: PathBuf, 
    // Add more fields as needed, e.g., total_chunks, processed_chunks for streaming/chunked mode progress
}

impl RunManifest {
    /// Creates a new `RunManifest` for a fresh run.
    pub fn new(run_id: String, args: OrbweaverArgs, run_output_dir: PathBuf) -> Self {
        RunManifest {
            run_id,
            original_args: args,
            status: RunStatus::Initializing,
            last_checkpoint_file: None,
            last_checkpoint_time: None,
            run_start_time: Utc::now(),
            output_dir_manifest: run_output_dir,
        }
    }

    /// Saves the manifest to the specified path as a JSON file.
    pub fn save(&self, manifest_path: &PathBuf) -> anyhow::Result<()> {
        use std::fs::File;
        use std::io::BufWriter;
        use anyhow::Context;

        let file = File::create(manifest_path)
            .with_context(|| format!("Failed to create manifest file at {:?}", manifest_path))?;
        let writer = BufWriter::new(file);
        serde_json::to_writer_pretty(writer, self)
            .with_context(|| format!("Failed to serialize manifest to {:?}", manifest_path))?;
        Ok(())
    }

    /// Loads a `RunManifest` from a JSON file at the specified path.
    pub fn load(manifest_path: &PathBuf) -> anyhow::Result<Self> {
        use std::fs::File;
        use std::io::BufReader;
        use anyhow::Context;

        let file = File::open(manifest_path)
            .with_context(|| format!("Failed to open manifest file at {:?}", manifest_path))?;
        let reader = BufReader::new(file);
        let manifest: Self = serde_json::from_reader(reader)
            .with_context(|| format!("Failed to deserialize manifest from {:?}", manifest_path))?;
        Ok(manifest)
    }
} 