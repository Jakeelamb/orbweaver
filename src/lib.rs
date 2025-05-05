// Declare the library modules
pub mod analysis;
pub mod fasta;
pub mod grammar;
pub mod io;
pub mod utils;
// pub mod tests; // Typically tests are not declared as a library module

// Optional: Re-export key functions/structs if needed for easier access
// pub use io::ncbi::run_ncbi_fetch;
// pub use io::encoding::run_fasta_to_vbq;
// pub use grammar::slp::build_slp_for_sequence;

// Note: Ensure all necessary dependencies used by these modules
// are listed in Cargo.toml. 