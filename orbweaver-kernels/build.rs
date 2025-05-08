fn main() {
    // No build steps are needed with embedded kernels
    println!("cargo:rerun-if-changed=src/lib.rs");
} 