use anyhow::{Context, Result};
use log::{debug, info};
use neo4rs::{query, Graph, Row};
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::fs::File;
use std::io::Write;
use std::path::Path;

/// Summary statistics for a group of genomes
#[derive(Debug, Serialize, Deserialize)]
pub struct GroupStatistics {
    /// Group identifier
    pub group_id: String,
    
    /// Number of genomes in the group
    pub genome_count: usize,
    
    /// Average assembly index across all genomes in the group
    pub avg_assembly_index: f64,
    
    /// Standard deviation of assembly index across all genomes in the group
    pub std_dev_assembly_index: f64,
    
    /// Average graph depth across all chromosomes in the group
    pub avg_graph_depth: f64,
    
    /// Average motif length across all chromosomes in the group
    pub avg_motif_length: f64,
    
    /// Minimum genome size in the group (bp)
    pub min_genome_size: usize,
    
    /// Maximum genome size in the group (bp)
    pub max_genome_size: usize,
    
    /// Average genome size in the group (bp)
    pub avg_genome_size: f64,
}

/// Data point for correlation analysis
#[derive(Debug, Serialize, Deserialize)]
pub struct CorrelationPoint {
    /// Genome ID
    pub genome_id: String,
    
    /// Species name
    pub species_name: String,
    
    /// Related group
    pub related_group: String,
    
    /// Genome size in base pairs
    pub genome_size: usize,
    
    /// Assembly index
    pub assembly_index: f64,
}

/// Correlation analysis results
#[derive(Debug, Serialize, Deserialize)]
pub struct CorrelationAnalysis {
    /// Data points for the correlation
    pub data_points: Vec<CorrelationPoint>,
    
    /// Pearson correlation coefficient
    pub pearson_r: f64,
    
    /// Spearman correlation coefficient
    pub spearman_r: f64,
    
    /// P-value for the Pearson correlation
    pub pearson_p: f64,
    
    /// P-value for the Spearman correlation
    pub spearman_p: f64,
}

/// Performs various analyses on the graph database
pub struct Analyzer {
    /// The Neo4j graph connection
    graph: Graph,
    
    /// Output directory for analysis results
    output_dir: String,
}

impl Analyzer {
    /// Create a new analyzer
    pub async fn new(uri: &str, username: &str, password: &str, output_dir: &str) -> Result<Self> {
        info!("Connecting to Neo4j at {}", uri);
        
        let graph = Graph::new(uri, username, password)
            .await
            .context("Failed to connect to Neo4j database for analysis")?;
        
        Ok(Self {
            graph,
            output_dir: output_dir.to_string(),
        })
    }
    
    /// Calculate summary statistics for each related group
    pub async fn calculate_group_statistics(&self) -> Result<Vec<GroupStatistics>> {
        info!("Calculating group statistics");
        
        // Get list of distinct groups
        let groups_query = query("MATCH (g:Genome) RETURN DISTINCT g.related_group AS group_id");
        let groups = self.graph.run(groups_query).await?;
        
        let mut group_ids = Vec::new();
        for row in groups {
            let group_id: String = row.get("group_id").context("Failed to get group_id")?;
            group_ids.push(group_id);
        }
        
        debug!("Found {} distinct groups", group_ids.len());
        
        // Calculate statistics for each group
        let mut group_stats = Vec::new();
        
        for group_id in group_ids {
            // Get summary statistics for this group
            let stats_query = query(
                "MATCH (g:Genome {related_group: $group_id})
                OPTIONAL MATCH (c:Chromosome)-[:PART_OF]->(g)
                WITH 
                    g.related_group AS group_id,
                    COUNT(DISTINCT g) AS genome_count,
                    AVG(g.assembly_index_avg) AS avg_assembly_index,
                    STDEV(g.assembly_index_avg) AS std_dev_assembly_index,
                    AVG(c.graph_depth) AS avg_graph_depth,
                    AVG(c.avg_motif_length) AS avg_motif_length,
                    MIN(g.total_processed_size_bp) AS min_genome_size,
                    MAX(g.total_processed_size_bp) AS max_genome_size,
                    AVG(g.total_processed_size_bp) AS avg_genome_size
                RETURN
                    group_id, genome_count, avg_assembly_index, std_dev_assembly_index,
                    avg_graph_depth, avg_motif_length,
                    min_genome_size, max_genome_size, avg_genome_size"
            ).param("group_id", &group_id);
            
            let result = self.graph.run(stats_query).await?;
            
            if let Some(row) = result.into_iter().next() {
                let stats = GroupStatistics {
                    group_id: row.get("group_id").unwrap_or_else(|_| group_id.clone()),
                    genome_count: row.get::<i64>("genome_count").unwrap_or(0) as usize,
                    avg_assembly_index: row.get("avg_assembly_index").unwrap_or(0.0),
                    std_dev_assembly_index: row.get("std_dev_assembly_index").unwrap_or(0.0),
                    avg_graph_depth: row.get("avg_graph_depth").unwrap_or(0.0),
                    avg_motif_length: row.get("avg_motif_length").unwrap_or(0.0),
                    min_genome_size: row.get::<i64>("min_genome_size").unwrap_or(0) as usize,
                    max_genome_size: row.get::<i64>("max_genome_size").unwrap_or(0) as usize,
                    avg_genome_size: row.get("avg_genome_size").unwrap_or(0.0),
                };
                
                group_stats.push(stats);
            }
        }
        
        // Save results to JSON file
        let output_path = Path::new(&self.output_dir).join("group_statistics.json");
        let file = File::create(&output_path)
            .context(format!("Failed to create output file: {}", output_path.display()))?;
        serde_json::to_writer_pretty(file, &group_stats)
            .context("Failed to write group statistics to JSON file")?;
        
        info!("Group statistics written to {}", output_path.display());
        
        Ok(group_stats)
    }
    
    /// Perform correlation analysis of genome size vs. assembly index
    pub async fn analyze_size_vs_index_correlation(&self) -> Result<CorrelationAnalysis> {
        info!("Analyzing correlation between genome size and assembly index");
        
        // Get data points for correlation
        let query = query(
            "MATCH (g:Genome)
            RETURN
                g.genome_id AS genome_id,
                g.species_name AS species_name,
                g.related_group AS related_group,
                g.total_processed_size_bp AS genome_size,
                g.assembly_index_avg AS assembly_index
            ORDER BY g.genome_id"
        );
        
        let result = self.graph.run(query).await?;
        
        let mut data_points = Vec::new();
        let mut sizes = Vec::new();
        let mut indices = Vec::new();
        
        for row in result {
            let point = CorrelationPoint {
                genome_id: row.get("genome_id").context("Failed to get genome_id")?,
                species_name: row.get("species_name").context("Failed to get species_name")?,
                related_group: row.get("related_group").context("Failed to get related_group")?,
                genome_size: row.get::<i64>("genome_size").context("Failed to get genome_size")? as usize,
                assembly_index: row.get("assembly_index").context("Failed to get assembly_index")?,
            };
            
            sizes.push(point.genome_size as f64);
            indices.push(point.assembly_index);
            data_points.push(point);
        }
        
        // Calculate correlation coefficients
        // Note: In a real implementation, you'd use a statistics library to calculate these
        let pearson_r = calculate_pearson_correlation(&sizes, &indices);
        let spearman_r = calculate_spearman_correlation(&sizes, &indices);
        
        // Save results to JSON file
        let output_path = Path::new(&self.output_dir).join("correlation_analysis.json");
        
        let correlation = CorrelationAnalysis {
            data_points,
            pearson_r,
            spearman_r,
            pearson_p: 0.0, // Would require a proper statistical calculation
            spearman_p: 0.0, // Would require a proper statistical calculation
        };
        
        let file = File::create(&output_path)
            .context(format!("Failed to create output file: {}", output_path.display()))?;
        serde_json::to_writer_pretty(file, &correlation)
            .context("Failed to write correlation analysis to JSON file")?;
        
        info!("Correlation analysis written to {}", output_path.display());
        
        Ok(correlation)
    }
    
    /// Generate a CSV file with chromosome metrics for visualization
    pub async fn export_chromosome_metrics_csv(&self) -> Result<()> {
        info!("Exporting chromosome metrics to CSV");
        
        // Query for chromosome metrics
        let query = query(
            "MATCH (c:Chromosome)-[:PART_OF]->(g:Genome)
            RETURN
                g.genome_id AS genome_id,
                g.species_name AS species_name,
                g.related_group AS related_group,
                c.chromosome_id AS chromosome_id,
                c.processed_length_bp AS length,
                c.assembly_index AS assembly_index,
                c.graph_depth AS graph_depth,
                c.avg_motif_length AS avg_motif_length
            ORDER BY g.related_group, g.genome_id, c.chromosome_id"
        );
        
        let result = self.graph.run(query).await?;
        
        // Write to CSV
        let output_path = Path::new(&self.output_dir).join("chromosome_metrics.csv");
        let mut file = File::create(&output_path)
            .context(format!("Failed to create output file: {}", output_path.display()))?;
        
        // Write header
        writeln!(file, "genome_id,species_name,related_group,chromosome_id,length,assembly_index,graph_depth,avg_motif_length")?;
        
        // Write data rows
        for row in result {
            writeln!(
                file,
                "{},{},{},{},{},{},{},{}",
                row.get::<String>("genome_id").unwrap_or_default(),
                row.get::<String>("species_name").unwrap_or_default(),
                row.get::<String>("related_group").unwrap_or_default(),
                row.get::<String>("chromosome_id").unwrap_or_default(),
                row.get::<i64>("length").unwrap_or(0),
                row.get::<i64>("assembly_index").unwrap_or(0),
                row.get::<i64>("graph_depth").unwrap_or(0),
                row.get::<f64>("avg_motif_length").unwrap_or(0.0)
            )?;
        }
        
        info!("Chromosome metrics exported to {}", output_path.display());
        
        Ok(())
    }
}

/// Calculate Pearson correlation coefficient
fn calculate_pearson_correlation(x: &[f64], y: &[f64]) -> f64 {
    if x.len() != y.len() || x.is_empty() {
        return 0.0;
    }
    
    let n = x.len() as f64;
    let sum_x: f64 = x.iter().sum();
    let sum_y: f64 = y.iter().sum();
    let sum_xy: f64 = x.iter().zip(y.iter()).map(|(a, b)| a * b).sum();
    let sum_x_sq: f64 = x.iter().map(|a| a * a).sum();
    let sum_y_sq: f64 = y.iter().map(|a| a * a).sum();
    
    let numerator = n * sum_xy - sum_x * sum_y;
    let denominator = ((n * sum_x_sq - sum_x * sum_x) * (n * sum_y_sq - sum_y * sum_y)).sqrt();
    
    if denominator == 0.0 {
        0.0
    } else {
        numerator / denominator
    }
}

/// Calculate Spearman correlation coefficient
fn calculate_spearman_correlation(x: &[f64], y: &[f64]) -> f64 {
    if x.len() != y.len() || x.is_empty() {
        return 0.0;
    }
    
    // This is a simplistic implementation
    // A proper implementation would convert to ranks and then use Pearson
    
    // Create vectors of (value, index) pairs
    let mut x_indexed: Vec<(f64, usize)> = x.iter().copied().enumerate().map(|(i, v)| (v, i)).collect();
    let mut y_indexed: Vec<(f64, usize)> = y.iter().copied().enumerate().map(|(i, v)| (v, i)).collect();
    
    // Sort by value
    x_indexed.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
    y_indexed.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
    
    // Assign ranks
    let mut x_ranks = vec![0.0; x.len()];
    let mut y_ranks = vec![0.0; y.len()];
    
    for (rank, &(_, idx)) in x_indexed.iter().enumerate() {
        x_ranks[idx] = rank as f64 + 1.0;
    }
    
    for (rank, &(_, idx)) in y_indexed.iter().enumerate() {
        y_ranks[idx] = rank as f64 + 1.0;
    }
    
    // Calculate Pearson correlation on the ranks
    calculate_pearson_correlation(&x_ranks, &y_ranks)
} 