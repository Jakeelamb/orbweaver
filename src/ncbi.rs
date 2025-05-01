use anyhow::{Context, Result};
use csv::{ReaderBuilder, WriterBuilder};
use dotenv::dotenv;
use quick_xml::de::from_str;
use reqwest::blocking::Client;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::env;
use std::fs::{self, File};
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};
use std::thread;
use std::time::Duration;
use glob;

// --- Configuration & Constants ---
const EUTILS_BASE_URL: &str = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/";
const NCBI_RATE_LIMIT_DELAY: Duration = Duration::from_millis(200); // ~5 requests/sec with API key
const DEFAULT_MAX_IDS: usize = 10000;
const DEFAULT_CSV_FILE: &str = "metazoan_chromosome_assemblies_with_lineage.csv";
const GCF_CSV_FILE: &str = "gcf_assemblies_by_size.csv";
const SAMPLE_TSV_FILE: &str = "sample.tsv";
const GENOMES_DIR: &str = "genomes";
const TAXONOMY_RANKS: [&str; 7] = ["superkingdom", "kingdom", "phylum", "class", "order", "family", "genus"];

// --- Structs for Deserialization & Serialization ---
// (Keep these internal unless needed elsewhere)

#[derive(Deserialize, Debug)]
struct ESearchResult {
    count: String,
    idlist: Vec<String>,
}

#[derive(Deserialize, Debug)]
struct ESearchResponse {
    esearchresult: ESearchResult,
}

#[derive(Deserialize, Debug, Clone, Serialize)]
struct AssemblySummary {
    #[serde(rename = "assemblyaccession")]
    assembly_accession: String,
    #[serde(rename = "taxid")]
    taxid: String,
    #[serde(rename = "organism")]
    organism: String,
    #[serde(rename = "assemblyname")]
    assembly_name: String,
    #[serde(rename = "assemblystatus")]
    assembly_status: String,
    #[serde(rename = "meta")]
    meta: Option<String>,
    #[serde(rename = "contign50")]
    contig_n50: Option<serde_json::Value>,
    #[serde(rename = "scaffoldn50")]
    scaffold_n50: Option<serde_json::Value>,
}

#[derive(Deserialize, Debug)]
struct ESummaryResult {
    uids: Vec<String>,
    #[serde(flatten)]
    summaries: HashMap<String, AssemblySummary>,
}

#[derive(Deserialize, Debug)]
struct ESummaryResponse {
    result: ESummaryResult,
}

#[derive(Deserialize, Debug, Default, Clone, Serialize)]
struct TaxonomyLineage {
    superkingdom: Option<String>,
    kingdom: Option<String>,
    phylum: Option<String>,
    class: Option<String>,
    order: Option<String>,
    family: Option<String>,
    genus: Option<String>,
}

#[derive(Deserialize, Debug)]
struct TaxonInfo {
    #[serde(rename = "ScientificName")]
    scientific_name: String,
    #[serde(rename = "Rank")]
    rank: String,
}

#[derive(Deserialize, Debug)]
struct LineageTaxon {
    #[serde(rename = "TaxId")]
    tax_id: i64,
    #[serde(rename = "ScientificName")]
    scientific_name: String,
    #[serde(rename = "Rank")]
    rank: String,
}

#[derive(Deserialize, Debug)]
struct LineageEx {
    #[serde(rename = "Taxon", default)]
    taxons: Vec<LineageTaxon>,
}

#[derive(Deserialize, Debug)]
struct Taxon {
    #[serde(rename = "TaxId")]
    tax_id: i64,
    #[serde(rename = "ScientificName")]
    scientific_name: String,
    #[serde(rename = "Rank")]
    rank: String,
    #[serde(rename = "LineageEx")]
    lineage_ex: Option<LineageEx>,
}

#[derive(Deserialize, Debug)]
struct TaxaSet {
    #[serde(rename = "Taxon", default)]
    taxons: Vec<Taxon>,
}

#[derive(Debug, Clone, Serialize, Deserialize)] // Add Deserialize for reading existing CSV
pub struct AssemblyRecord { // Made public for potential use elsewhere, though maybe not needed
    // From AssemblySummary
    assembly_accession: String,
    taxid: String,
    organism: String,
    assembly_name: String,
    assembly_status: String,
    meta: Option<String>,
    contig_n50: Option<String>,
    scaffold_n50: Option<String>,

    // From TaxonomyLineage
    superkingdom: Option<String>,
    kingdom: Option<String>,
    phylum: Option<String>,
    class: Option<String>,
    order: Option<String>,
    family: Option<String>,
    genus: Option<String>,

    // Add assembly_size for sorting GCF records
    #[serde(skip_serializing_if = "Option::is_none")]
    #[serde(default)] // Needed for deserializing if field is missing
    assembly_size: Option<u64>,
}


// Helper to convert serde_json::Value (Number or String) to Option<String>
fn json_value_to_string(value: &Option<serde_json::Value>) -> Option<String> {
     value.as_ref().map(|v| {
        if v.is_string() {
            v.as_str().unwrap_or("").to_string()
        } else if v.is_number() {
            v.to_string()
        } else {
            String::new()
        }
    }).filter(|s| !s.is_empty())
}

// Helper to run external commands
pub fn run_command(program: &str, args: &[&str], cwd: Option<&Path>) -> Result<()> {
    println!("Running command: {} {} (in {:?})", program, args.join(" "), cwd.unwrap_or_else(|| Path::new(".")));
    let mut command = Command::new(program);
    command.args(args);
    if let Some(dir) = cwd {
        command.current_dir(dir);
    }
    command.stdout(Stdio::piped());
    command.stderr(Stdio::piped());
    let output = command.output().context(format!("Failed to execute command: {}", program))?;
    if !output.status.success() {
        let stdout = String::from_utf8_lossy(&output.stdout);
        let stderr = String::from_utf8_lossy(&output.stderr);
        anyhow::bail!(
            "Command '{}' failed with status {}\nStdout:\n{}\nStderr:\n{}",
            program, output.status, stdout, stderr
        );
    }
    println!("Command '{}' executed successfully.", program);
    Ok(())
}

// --- NCBI API Functions ---

pub fn get_metazoan_chromosome_assembly_ids(api_key: &str, max_ids: usize) -> Result<Vec<String>> {
    println!("Searching for Metazoan chromosome-level assemblies...");
    let search_term = r#"("Metazoa"[Organism]) AND "Chromosome"[Assembly Level]"#;
    let search_url = format!(
        "{}esearch.fcgi?db=assembly&term={}&retmax={}&retmode=json&api_key={}",
        EUTILS_BASE_URL, search_term, max_ids, api_key
    );

    let client = Client::new();
    let response = client.get(&search_url).send()
        .context(format!("Failed to send ESearch request to {}", search_url))?;

    if !response.status().is_success() {
        anyhow::bail!(
            "ESearch request failed with status: {}. Response: {}",
            response.status(),
            response.text().unwrap_or_else(|_| "<failed to read response body>".to_string())
        );
    }

    let search_data: ESearchResponse = response.json().context("Failed to parse ESearch JSON response")?;
    let assembly_ids = search_data.esearchresult.idlist;
    let count: usize = search_data.esearchresult.count.parse().unwrap_or(0);
    println!("Found {} assembly IDs matching criteria (fetched max {}).", count, assembly_ids.len());
    Ok(assembly_ids)
}


pub fn get_assembly_summaries(assembly_ids: &[String], api_key: &str, batch_size: usize) -> Result<Vec<AssemblySummary>> {
    println!("Fetching assembly summaries for {} IDs...", assembly_ids.len());
    let mut summaries = Vec::new();
    let client = Client::new();

    for batch in assembly_ids.chunks(batch_size) {
        let ids_str = batch.join(",");
        let summary_url = format!(
            "{}esummary.fcgi?db=assembly&id={}&retmode=json&api_key={}",
            EUTILS_BASE_URL, ids_str, api_key
        );

        thread::sleep(NCBI_RATE_LIMIT_DELAY); // Rate limit

        let response = client.get(&summary_url).send()
            .context(format!("Failed to send ESummary request to {}", summary_url))?;

        if !response.status().is_success() {
            eprintln!(
                "Warning: ESummary request for batch failed with status: {}. URL: {}",
                response.status(), summary_url
            );
             eprintln!("Response body: {}", response.text().unwrap_or_else(|_| "<failed to read body>".to_string()));
            continue; // Skip this batch on error
        }

         // Read the response text first for better debugging
        let response_text = response.text().context("Failed to read ESummary response text")?;

        // Attempt to parse the JSON
        match serde_json::from_str::<ESummaryResponse>(&response_text) {
             Ok(summary_data) => {
                // Extract summaries based on the UIDs returned in the response
                for uid in &summary_data.result.uids {
                    if let Some(summary) = summary_data.result.summaries.get(uid) {
                        summaries.push(summary.clone());
                    } else {
                         eprintln!("Warning: Summary for UID {} not found in response batch.", uid);
                    }
                 }
             }
             Err(e) => {
                 eprintln!("Warning: Failed to parse ESummary JSON response: {}. Skipping batch.", e);
                 eprintln!("Response text was: {}", response_text); // Log the problematic JSON
                 continue;
             }
         }
    }
    println!("Successfully fetched {} assembly summaries.", summaries.len());
    Ok(summaries)
}


pub fn get_taxonomy_lineages(taxids: &[i64], api_key: &str, batch_size: usize) -> Result<HashMap<String, TaxonomyLineage>> {
    println!("Fetching taxonomy lineages for {} taxids...", taxids.len());
    let mut lineages = HashMap::new();
    let client = Client::new();

    for batch in taxids.chunks(batch_size) {
        let ids_str = batch.iter().map(|id| id.to_string()).collect::<Vec<_>>().join(",");
        let lineage_url = format!(
            "{}efetch.fcgi?db=taxonomy&id={}&retmode=xml&api_key={}",
            EUTILS_BASE_URL, ids_str, api_key
        );

        thread::sleep(NCBI_RATE_LIMIT_DELAY); // Rate limit

        let response = client.get(&lineage_url).send()
             .context(format!("Failed to send EFetch (taxonomy) request to {}", lineage_url))?;

        if !response.status().is_success() {
            eprintln!(
                "Warning: EFetch (taxonomy) request failed with status: {}. URL: {}",
                response.status(), lineage_url
            );
            eprintln!("Response body: {}", response.text().unwrap_or_else(|_| "<failed to read body>".to_string()));

            continue; // Skip this batch on error
        }

        let response_text = response.text().context("Failed to read EFetch (taxonomy) response text")?;

        match from_str::<TaxaSet>(&response_text) {
            Ok(taxa_set) => {
                for taxon in taxa_set.taxons {
                    let mut lineage = TaxonomyLineage::default();
                    if let Some(lineage_ex) = taxon.lineage_ex {
                        for item in lineage_ex.taxons {
                            match item.rank.as_str() {
                                "superkingdom" => lineage.superkingdom = Some(item.scientific_name),
                                "kingdom" => lineage.kingdom = Some(item.scientific_name),
                                "phylum" => lineage.phylum = Some(item.scientific_name),
                                "class" => lineage.class = Some(item.scientific_name),
                                "order" => lineage.order = Some(item.scientific_name),
                                "family" => lineage.family = Some(item.scientific_name),
                                "genus" => lineage.genus = Some(item.scientific_name),
                                _ => {} // Ignore other ranks
                            }
                        }
                    }
                     // Use the tax_id from the main Taxon element as the key
                     lineages.insert(taxon.tax_id.to_string(), lineage);
                }
            }
            Err(e) => {
                eprintln!("Warning: Failed to parse EFetch (taxonomy) XML response: {}. Skipping batch.", e);
                eprintln!("Response text was: {}", response_text);
                continue;
            }
        }
    }
    println!("Successfully fetched {} taxonomy lineages.", lineages.len());
    Ok(lineages)
}


pub fn merge_data(summaries: Vec<AssemblySummary>, lineages: HashMap<String, TaxonomyLineage>) -> Vec<AssemblyRecord> {
    println!("Merging assembly summaries and taxonomy lineages...");
    summaries.into_iter().map(|summary| {
        let lineage = lineages.get(&summary.taxid).cloned().unwrap_or_default();
        let assembly_size = parse_total_length_from_meta(&summary.meta); // Calculate assembly size here
        AssemblyRecord {
            assembly_accession: summary.assembly_accession,
            taxid: summary.taxid,
            organism: summary.organism,
            assembly_name: summary.assembly_name,
            assembly_status: summary.assembly_status,
            meta: summary.meta,
            contig_n50: json_value_to_string(&summary.contig_n50),
            scaffold_n50: json_value_to_string(&summary.scaffold_n50),
            superkingdom: lineage.superkingdom,
            kingdom: lineage.kingdom,
            phylum: lineage.phylum,
            class: lineage.class,
            order: lineage.order,
            family: lineage.family,
            genus: lineage.genus,
            assembly_size, // Add the calculated size
        }
    }).collect()
}


pub fn save_records_to_csv<T: Serialize>(records: &[T], path: &Path) -> Result<()> {
    println!("Saving {} records to CSV: {}", records.len(), path.display());
    fs::create_dir_all(path.parent().unwrap_or_else(|| Path::new(".")))?; // Ensure directory exists
    let file = File::create(path).context(format!("Failed to create CSV file: {}", path.display()))?;
    let mut wtr = WriterBuilder::new().from_writer(file);
    for record in records {
        wtr.serialize(record)?;
    }
    wtr.flush()?;
    println!("Successfully saved CSV: {}", path.display());
    Ok(())
}

// Helper function to parse N50 values robustly
fn parse_n50(value_str: &Option<String>) -> Option<u64> {
    value_str.as_ref().and_then(|s| s.parse::<u64>().ok())
}

// Helper function to parse 'total-length' from meta string
fn parse_total_length_from_meta(meta: &Option<String>) -> Option<u64> {
    meta.as_ref().and_then(|m| {
        m.split(';')
            .find(|s| s.trim().starts_with("total-length="))
            .and_then(|s| s.split('=').nth(1))
            .and_then(|val| val.trim().parse::<u64>().ok())
    })
}


pub fn create_gcf_assemblies_by_size(records: &[AssemblyRecord]) -> Vec<AssemblyRecord> {
    println!("Filtering GCF assemblies and calculating sizes...");
    let mut gcf_records: Vec<AssemblyRecord> = records.iter()
        .filter(|r| r.assembly_accession.starts_with("GCF_"))
        .cloned() // Clone records that match
        .collect();

    println!("Found {} GCF assemblies.", gcf_records.len());

    // Sort by assembly size (descending), requires size calculation first
    gcf_records.sort_by(|a, b| {
        b.assembly_size.cmp(&a.assembly_size) // Sort descending by size
    });

    println!("Sorted GCF assemblies by size.");
    gcf_records
}


pub fn create_sample_tsv(gcf_records: &[AssemblyRecord], path: &Path) -> Result<()> {
     println!("Creating sample TSV file: {}", path.display());
    fs::create_dir_all(path.parent().unwrap_or_else(|| Path::new(".")))?; // Ensure directory exists

    let mut file = File::create(path)?;
    writeln!(file, "assembly_accession	organism_name	assembly_size")?; // Header

    let mut count = 0;
    for record in gcf_records.iter().take(10) { // Take top 10 largest
        if let (Some(size), Some(org)) = (record.assembly_size, Some(&record.organism)) {
             // Sanitize organism name for file system use (replace spaces, etc.)
             let sanitized_org_name = org.replace(|c: char| !c.is_alphanumeric() && c != '_', "_");
            writeln!(
                file,
                "{}	{}	{}",
                record.assembly_accession,
                sanitized_org_name, // Use sanitized name
                size
            )?;
             count += 1;
        }
    }

    println!("Saved top {} largest GCF assemblies to {}", count, path.display());
    Ok(())
}


pub fn download_sample_genome(sample_tsv_path: &Path, genomes_dir: &Path) -> Result<()> {
    println!("Downloading sample genome based on {}", sample_tsv_path.display());
    fs::create_dir_all(genomes_dir)?; // Ensure genomes directory exists

    // Read the sample TSV file
    let file = File::open(sample_tsv_path)
        .context(format!("Failed to open sample TSV file: {}", sample_tsv_path.display()))?;
    let mut rdr = ReaderBuilder::new().delimiter(b'\t').from_reader(file);

    let mut records = rdr.deserialize::<HashMap<String, String>>();

    // Find the record with the largest assembly size (should be the first record if sorted)
    if let Some(result) = records.next() {
        let record = result?;
        let accession = record.get("assembly_accession")
            .context("Missing 'assembly_accession' in sample TSV")?;
        // let organism = record.get("organism_name")
        //     .context("Missing 'organism_name' in sample TSV")?; // Not strictly needed for download URL

        println!("Selected largest assembly for download: {}", accession);

        // Construct the download URL (example for FASTA format)
        // Base URL: https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/xxx/xxx/xxx/GCF_xxxxxxxx.x_Assembly Name/
        // File: GCF_xxxxxxxx.x_Assembly Name_genomic.fna.gz

        // Need to query ESummary again for the full assembly name or FTP path if not stored
        // For now, let's *assume* a predictable pattern (which is often unreliable)
        // A more robust way involves querying esummary/elink for the FTP path.

        // Simplified approach: Use NCBI datasets tool (if available)
        let datasets_path = "datasets"; // Assume datasets CLI tool is in PATH
        let output_filename = format!("{}.zip", accession);
        let output_path = genomes_dir.join(&output_filename);
        let output_path_str = output_path.to_str().context("Invalid output path")?;

        println!("Attempting download using NCBI datasets tool...");
        run_command(
            datasets_path,
            &[
                "download",
                "genome",
                "accession",
                accession,
                "--include", "genome", // Download genome sequence
                "--filename", output_path_str,
            ],
            None // Run in current directory
        ).context(format!("Failed to download genome {} using datasets tool", accession))?;

        println!("Downloaded {} to {}", accession, output_path.display());

        // Unzip the downloaded file
        let unzip_output_dir = genomes_dir.join(accession); // Unzip into a directory named after the accession
        fs::create_dir_all(&unzip_output_dir)?;
        let unzip_output_dir_str = unzip_output_dir.to_str().context("Invalid unzip directory path")?;

        println!("Unzipping {} to {}", output_path.display(), unzip_output_dir.display());
        run_command(
            "unzip",
            &[
                "-o", // Overwrite files without prompting
                output_path_str,
                "-d", unzip_output_dir_str,
            ],
            None
        ).context(format!("Failed to unzip {}", output_path.display()))?;

        // --- Find the FASTA file ---
        // The unzipped structure is typically ncbi_dataset/data/GCF_xxxx/sequence.fna
        let fna_pattern = unzip_output_dir.join("ncbi_dataset/data").join(accession).join("*.fna");
        println!("Searching for FASTA file matching pattern: {}", fna_pattern.display());

        let fna_files: Vec<PathBuf> = glob::glob(fna_pattern.to_str().unwrap())?
                                        .filter_map(Result::ok)
                                        .collect();

        if let Some(fna_path) = fna_files.first() {
            println!("Found FASTA file: {}", fna_path.display());
            // Optionally move or rename the FASTA file here if needed
             // Example: Move to genomes/chrom_1.fasta (specific for this project?)
             let target_fasta_path = genomes_dir.join("chrom_1.fasta");
             println!("Copying {} to {}", fna_path.display(), target_fasta_path.display());
             fs::copy(fna_path, &target_fasta_path)
                 .context(format!("Failed to copy FASTA file to {}", target_fasta_path.display()))?;
             println!("Copied FASTA file for chromosome 1.");

        } else {
            eprintln!("Warning: Could not find the genomic FASTA file (.fna) after unzipping.");
             println!("Contents of {}:", unzip_output_dir.display());
             // Try listing contents for debugging
             match fs::read_dir(&unzip_output_dir) {
                 Ok(entries) => {
                     for entry in entries {
                         if let Ok(entry) = entry {
                             println!("  - {}", entry.path().display());
                         }
                     }
                 }
                 Err(e) => eprintln!("  Could not list directory contents: {}", e),
             }
        }


        // Clean up the zip file? Optional.
        // fs::remove_file(&output_path)?;

    } else {
        eprintln!("Warning: No records found in sample TSV file to download.");
    }

    Ok(())
}

// Main entry point for the NCBI fetch subcommand
pub fn run_ncbi_fetch(generate: bool, csv_path: PathBuf) -> Result<()> {
    dotenv().ok(); // Load .env file if present
    let api_key = env::var("NCBI_API_KEY").context("NCBI_API_KEY not found in environment or .env file")?;

    let records = if !csv_path.exists() || generate {
        println!("Generating new data from NCBI...");

        // 1. Fetch Assembly IDs
        let assembly_ids = get_metazoan_chromosome_assembly_ids(&api_key, DEFAULT_MAX_IDS)?;
        if assembly_ids.is_empty() {
            anyhow::bail!("No assembly IDs found. Exiting.");
        }

        // 2. Fetch Assembly Summaries
        let summaries = get_assembly_summaries(&assembly_ids, &api_key, 100)?; // Batch size 100 for summaries

        // 3. Extract TaxIDs
        let taxids: Vec<i64> = summaries.iter()
            .filter_map(|s| s.taxid.parse::<i64>().ok()) // Parse taxid to i64, filter out errors
            .collect();

        // 4. Fetch Taxonomy Lineages
        let lineages = get_taxonomy_lineages(&taxids, &api_key, 200)?; // Batch size 200 for taxonomy

        // 5. Merge Data
        let merged_records = merge_data(summaries, lineages);

        // 6. Save Combined Data
        save_records_to_csv(&merged_records, &csv_path)?;

        merged_records // Return the newly fetched records

    } else {
        println!("Loading existing data from CSV: {}", csv_path.display());
        let file = File::open(&csv_path)
            .context(format!("Failed to open existing CSV file: {}", csv_path.display()))?;
        let mut rdr = ReaderBuilder::new().from_reader(file);
        let records: Vec<AssemblyRecord> = rdr.deserialize().collect::<Result<_, _>>()?;
        println!("Loaded {} records.", records.len());
        records // Return the loaded records
    };

    // --- Post-processing Steps ---

     // Ensure the genomes directory exists
     fs::create_dir_all(GENOMES_DIR)?;
     let genomes_dir_path = PathBuf::from(GENOMES_DIR);

    // 7. Create GCF list sorted by size
    let gcf_records = create_gcf_assemblies_by_size(&records);
    let gcf_csv_path = genomes_dir_path.join(GCF_CSV_FILE);
    save_records_to_csv(&gcf_records, &gcf_csv_path)?;

    // 8. Create Sample TSV (Top 10 largest GCF)
    let sample_tsv_path = genomes_dir_path.join(SAMPLE_TSV_FILE);
    create_sample_tsv(&gcf_records, &sample_tsv_path)?;

    // 9. Download Sample Genome (Largest one from TSV)
    // Only download if generate flag is set OR the target fasta doesn't exist?
    // For now, let's tie it to the `generate` flag to avoid repeated downloads.
    if generate {
        download_sample_genome(&sample_tsv_path, &genomes_dir_path)?;
    } else {
        println!("Skipping sample genome download (use --generate to force download).");
    }


    println!("NCBI data processing complete.");
    Ok(())
} 