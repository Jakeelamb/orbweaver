use serde::{Serialize, Deserialize};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum AlphabetType {
    Dna,        // A, C, G, T, N
    Protein,    // 20 standard amino acids, plus others
    Utf8,       // General UTF-8 characters
    Byte,       // Raw bytes (0-255)
}

impl Default for AlphabetType {
    fn default() -> Self {
        AlphabetType::Dna // Default to DNA for now
    }
}

impl AlphabetType {
    // We can add methods here later, e.g., to get max symbol value
    // or to validate symbols.
    pub fn from_str(s: &str) -> Result<Self, String> {
        match s.to_lowercase().as_str() {
            "dna" => Ok(AlphabetType::Dna),
            "protein" => Ok(AlphabetType::Protein),
            "utf8" | "utf-8" => Ok(AlphabetType::Utf8),
            "byte" | "bytes" => Ok(AlphabetType::Byte),
            _ => Err(format!("Unknown alphabet type: {}", s)),
        }
    }
}

// Implement std::str::FromStr to allow direct parsing from string by clap
impl std::str::FromStr for AlphabetType {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        AlphabetType::from_str(s)
    }
}
