// Declare the library modules
pub mod encoding;
pub mod motifs;
pub mod ncbi;
pub mod slp;

// Optional: Re-export key functions/structs if needed for easier access
// pub use ncbi::run_ncbi_fetch;
// pub use encoding::run_fasta_to_vbq;
// pub use slp::build_slp_for_sequence;

// Note: Ensure all necessary dependencies used by these modules
// are listed in Cargo.toml. 