use anyhow::{bail, Context, Result};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Read, Write};
use std::path::{Path, PathBuf};

/// Supported file formats for grammar output
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OutputFormat {
    /// JSON format
    Json,
    /// Text format
    Text,
    /// GFA (Graphical Fragment Assembly) format
    Gfa,
    /// DOT format for visualization
    Dot,
    /// FASTA format for rule sequences
    Fasta,
}

impl OutputFormat {
    /// Parse a string to get the output format
    pub fn from_str(s: &str) -> Result<Self> {
        match s.to_lowercase().as_str() {
            "json" => Ok(Self::Json),
            "text" => Ok(Self::Text),
            "gfa" => Ok(Self::Gfa),
            "dot" => Ok(Self::Dot),
            "fasta" => Ok(Self::Fasta),
            _ => bail!("Unsupported output format: {}", s),
        }
    }

    /// Get the file extension for this format
    pub fn extension(&self) -> &'static str {
        match self {
            Self::Json => "json",
            Self::Text => "txt",
            Self::Gfa => "gfa",
            Self::Dot => "dot",
            Self::Fasta => "fasta",
        }
    }
}

/// Get the file format from a file extension
pub fn format_from_extension(path: &Path) -> Result<OutputFormat> {
    let ext = path
        .extension()
        .and_then(|e| e.to_str())
        .map(|e| e.to_lowercase());

    match ext.as_deref() {
        Some("json") => Ok(OutputFormat::Json),
        Some("txt") | Some("text") => Ok(OutputFormat::Text),
        Some("gfa") => Ok(OutputFormat::Gfa),
        Some("dot") => Ok(OutputFormat::Dot),
        Some("fasta") | Some("fa") | Some("fna") => Ok(OutputFormat::Fasta),
        _ => bail!("Unknown or missing file extension for path: {}", path.display()),
    }
}

/// Open a file for reading
pub fn open_file_for_reading(path: &Path) -> Result<BufReader<File>> {
    let file = File::open(path)
        .with_context(|| format!("Failed to open file for reading: {}", path.display()))?;
    Ok(BufReader::new(file))
}

/// Open a file for writing
pub fn open_file_for_writing(path: &Path) -> Result<BufWriter<File>> {
    // Ensure parent directory exists
    if let Some(parent) = path.parent() {
        std::fs::create_dir_all(parent)
            .with_context(|| format!("Failed to create directory: {}", parent.display()))?;
    }
    
    let file = File::create(path)
        .with_context(|| format!("Failed to open file for writing: {}", path.display()))?;
    Ok(BufWriter::new(file))
}

/// Ensure a file has the correct extension
pub fn ensure_extension(path: &Path, format: OutputFormat) -> PathBuf {
    let extension = format.extension();
    
    let path_str = path.to_string_lossy();
    let path_without_ext = match path.extension() {
        Some(_) => {
            // Remove existing extension
            let pos = path_str.rfind('.').unwrap_or(path_str.len());
            path_str[..pos].to_string()
        }
        None => path_str.to_string(),
    };
    
    PathBuf::from(format!("{}.{}", path_without_ext, extension))
}

/// Read a file as a string
pub fn read_to_string(path: &Path) -> Result<String> {
    let mut file = File::open(path)
        .with_context(|| format!("Failed to open file: {}", path.display()))?;
    let mut content = String::new();
    file.read_to_string(&mut content)
        .with_context(|| format!("Failed to read file: {}", path.display()))?;
    Ok(content)
}

/// Write a string to a file
pub fn write_string_to_file(content: &str, path: &Path) -> Result<()> {
    let mut writer = open_file_for_writing(path)?;
    writer.write_all(content.as_bytes())
        .with_context(|| format!("Failed to write to file: {}", path.display()))?;
    writer.flush()
        .with_context(|| format!("Failed to flush file: {}", path.display()))?;
    Ok(())
}

/// Count lines in a file
pub fn count_lines(path: &Path) -> Result<usize> {
    let file = File::open(path)
        .with_context(|| format!("Failed to open file: {}", path.display()))?;
    let reader = BufReader::new(file);
    let count = reader.lines().count();
    Ok(count)
}

/// Check if a file exists
pub fn file_exists(path: &Path) -> bool {
    path.exists() && path.is_file()
}

/// Create a temporary file
pub fn create_temp_file(prefix: &str, extension: &str) -> Result<(PathBuf, BufWriter<File>)> {
    let temp_dir = std::env::temp_dir();
    let filename = format!("{}_{}.{}", prefix, uuid::Uuid::new_v4(), extension);
    let path = temp_dir.join(filename);
    
    let file = File::create(&path)
        .with_context(|| format!("Failed to create temporary file: {}", path.display()))?;
    
    Ok((path, BufWriter::new(file)))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::tempdir;
    
    #[test]
    fn test_output_format() {
        assert_eq!(OutputFormat::from_str("json").unwrap(), OutputFormat::Json);
        assert_eq!(OutputFormat::from_str("JSON").unwrap(), OutputFormat::Json);
        assert!(OutputFormat::from_str("unknown").is_err());
        
        assert_eq!(OutputFormat::Json.extension(), "json");
        assert_eq!(OutputFormat::Text.extension(), "txt");
    }
    
    #[test]
    fn test_format_from_extension() {
        let json_path = Path::new("test.json");
        assert_eq!(format_from_extension(json_path).unwrap(), OutputFormat::Json);
        
        let fasta_path = Path::new("test.fasta");
        assert_eq!(format_from_extension(fasta_path).unwrap(), OutputFormat::Fasta);
        
        let unknown_path = Path::new("test.unknown");
        assert!(format_from_extension(unknown_path).is_err());
    }
    
    #[test]
    fn test_ensure_extension() {
        let path = Path::new("test");
        let with_ext = ensure_extension(path, OutputFormat::Json);
        assert_eq!(with_ext, PathBuf::from("test.json"));
        
        let path_with_ext = Path::new("test.txt");
        let changed_ext = ensure_extension(path_with_ext, OutputFormat::Json);
        assert_eq!(changed_ext, PathBuf::from("test.json"));
    }
    
    #[test]
    fn test_read_write_string() -> Result<()> {
        let dir = tempdir()?;
        let file_path = dir.path().join("test.txt");
        
        let content = "Test content";
        write_string_to_file(content, &file_path)?;
        
        let read_content = read_to_string(&file_path)?;
        assert_eq!(read_content, content);
        
        Ok(())
    }
    
    #[test]
    fn test_count_lines() -> Result<()> {
        let dir = tempdir()?;
        let file_path = dir.path().join("test.txt");
        
        {
            let mut file = File::create(&file_path)?;
            writeln!(file, "Line 1")?;
            writeln!(file, "Line 2")?;
            writeln!(file, "Line 3")?;
        }
        
        let line_count = count_lines(&file_path)?;
        assert_eq!(line_count, 3);
        
        Ok(())
    }
} 