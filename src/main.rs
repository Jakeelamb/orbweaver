use anyhow::Result;
use clap::Parser;
use log::{debug, info, warn, LevelFilter};
use rayon::prelude::*;
use std::path::Path;

mod genome;
mod slp;
mod graph;
mod database;
mod analysis;
mod utils;

use crate::database::DatabaseConfig;
use crate::genome::{Chromosome, Genome};
use crate::graph::MotifGraph;
use crate::slp::SLPBuilder;
use crate::utils::ensure_dir_exists;

/// Orbweaver: A tool for computing assembly indices and building substring motif graphs
/// for genome assemblies from NCBI.
#[derive(Parser, Debug)]
#[clap(author, version, about)]
struct Args {
    /// Path to input FASTA file containing chromosomes/scaffolds
    #[clap(short, long)]
    input: String,

    /// Output directory for storing results
    #[clap(short, long, default_value = "output")]
    output: String,

    /// Neo4j connection URI
    #[clap(long, default_value = "bolt://localhost:7687")]
    neo4j_uri: String,

    /// Neo4j username
    #[clap(long, default_value = "neo4j")]
    neo4j_user: String,

    /// Neo4j password
    #[clap(long)]
    neo4j_password: String,

    /// Genome ID
    #[clap(long)]
    genome_id: String,

    /// Species name
    #[clap(long)]
    species_name: String,

    /// Related group (taxonomic group, clade, etc.)
    #[clap(long)]
    related_group: String,

    /// SLP threshold factor (higher means more aggressive compression)
    #[clap(long, default_value = "1000000")]
    threshold_factor: f64,

    /// Number of threads to use for parallel processing
    #[clap(short, long, default_value = "0")]
    threads: usize,

    /// Skip database operations
    #[clap(long)]
    skip_db: bool,

    /// Debug mode
    #[clap(short, long)]
    debug: bool,
}

#[tokio::main]
async fn main() -> Result<()> {
    // Parse command line arguments
    let args = Args::parse();

    // Setup logging
    let log_level = if args.debug {
        LevelFilter::Debug
    } else {
        LevelFilter::Info
    };
    
    env_logger::builder()
        .filter_level(log_level)
        .init();

    info!("Starting Orbweaver...");
    
    // Set up thread pool size if specified
    if args.threads > 0 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(args.threads)
            .build_global()?;
        info!("Using {} threads for parallel processing", args.threads);
    }

    // Create output directories
    let output_dir = Path::new(&args.output);
    ensure_dir_exists(output_dir)?;

    // Process genome to get chromosomes with VBQ files
    info!("Processing input genome from {}", args.input);
    let mut genome_processor = genome::GenomeProcessor::new(&args.input, &args.output)?;
    let mut genome = genome_processor.process()?;
    
    // Update genome with command line arguments
    genome.id = args.genome_id;
    genome.species_name = args.species_name;
    genome.related_group = args.related_group;
    
    info!("Processed genome: {} ({}) with {} chromosomes", 
        genome.species_name, genome.id, genome.chromosomes.len());

    // Process each chromosome to build SLPs and calculate metrics
    info!("Building SLPs and calculating assembly indices for each chromosome");
    
    // Use ParallelIterator to process chromosomes in parallel
    let chromosome_results: Vec<Result<(Chromosome, MotifGraph)>> = genome.chromosomes
        .par_iter()
        .map(|chromosome| {
            info!("Processing chromosome: {}", chromosome.id);
            
            // Build SLP
            let mut slp_builder = SLPBuilder::new()
                .with_threshold_factor(args.threshold_factor);
            
            let slp = slp_builder.build_from_vbq(chromosome)?;
            
            // Convert SLP to graph
            let motif_graph = MotifGraph::from_slp(&slp)?;
            
            // Create updated chromosome with metrics
            let mut updated_chromosome = chromosome.clone();
            updated_chromosome.assembly_index = Some(slp.assembly_index);
            updated_chromosome.graph_depth = Some(motif_graph.depth);
            updated_chromosome.avg_motif_length = Some(motif_graph.avg_motif_length);
            
            info!(
                "Chromosome {} processed: assembly index = {}, graph depth = {}, avg motif length = {:.2}",
                updated_chromosome.id,
                updated_chromosome.assembly_index.unwrap(),
                updated_chromosome.graph_depth.unwrap(),
                updated_chromosome.avg_motif_length.unwrap()
            );
            
            Ok((updated_chromosome, motif_graph))
        })
        .collect();
    
    // Process results and update genome
    let mut chromosome_graphs = Vec::new();
    let mut updated_chromosomes = Vec::new();
    let mut total_assembly_index = 0.0;
    let mut total_weight = 0.0;
    
    for result in chromosome_results {
        match result {
            Ok((chromosome, graph)) => {
                // Calculate weighted assembly index
                let weight = chromosome.processed_length as f64;
                total_assembly_index += chromosome.assembly_index.unwrap() as f64 * weight;
                total_weight += weight;
                
                updated_chromosomes.push(chromosome);
                chromosome_graphs.push(graph);
            }
            Err(e) => {
                warn!("Error processing chromosome: {}", e);
            }
        }
    }
    
    // Update genome with processed chromosomes and weighted average assembly index
    genome.chromosomes = updated_chromosomes;
    genome.assembly_index_avg = Some(total_assembly_index / total_weight);
    
    info!(
        "Genome processing complete. Average assembly index: {:.2}",
        genome.assembly_index_avg.unwrap()
    );
    
    // Store results in database
    if !args.skip_db {
        info!("Connecting to Neo4j database");
        let db_config = DatabaseConfig::new(
            &args.neo4j_uri,
            &args.neo4j_user,
            &args.neo4j_password,
            30, // 30 second timeout
        );
        
        let db = database::Database::new(&db_config).await?;
        
        // Test connection
        db.test_connection().await?;
        
        // Store genome
        db.store_genome(&genome).await?;
        
        // Store motif graphs for each chromosome
        for (i, graph) in chromosome_graphs.iter().enumerate() {
            let chromosome = &genome.chromosomes[i];
            db.store_motif_graph(&genome.id, &chromosome.id, graph).await?;
        }
        
        // Run analysis
        let analyzer = analysis::Analyzer::new(
            &args.neo4j_uri,
            &args.neo4j_user,
            &args.neo4j_password,
            &args.output,
        ).await?;
        
        // Generate comparative analysis
        analyzer.calculate_group_statistics().await?;
        analyzer.analyze_size_vs_index_correlation().await?;
        analyzer.export_chromosome_metrics_csv().await?;
        
        info!("Analysis complete and results saved to {}", args.output);
    } else {
        info!("Skipping database operations as requested");
    }
    
    info!("Orbweaver processing complete!");
    
    Ok(())
}
