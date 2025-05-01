use anyhow::{Context, Result};
use log::{debug, info, warn};
use neo4rs::{Graph, Node, query};
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::time::Duration;

use crate::genome::{Chromosome, Genome};
use crate::graph::{MotifGraph, MotifNode};

/// Database configuration
pub struct DatabaseConfig {
    /// Neo4j URI
    pub uri: String,
    
    /// Neo4j username
    pub username: String,
    
    /// Neo4j password
    pub password: String,
    
    /// Connection timeout in seconds
    pub timeout: u64,
}

impl DatabaseConfig {
    /// Create a new database configuration
    pub fn new(uri: &str, username: &str, password: &str, timeout: u64) -> Self {
        Self {
            uri: uri.to_string(),
            username: username.to_string(),
            password: password.to_string(),
            timeout,
        }
    }
}

/// Neo4j database connection and operations
pub struct Database {
    /// The Neo4j graph connection
    graph: Graph,
}

impl Database {
    /// Create a new database connection
    pub async fn new(config: &DatabaseConfig) -> Result<Self> {
        info!("Connecting to Neo4j at {}", config.uri);
        
        let graph = Graph::new(&config.uri, &config.username, &config.password)
            .await
            .context("Failed to connect to Neo4j database")?;
        
        info!("Successfully connected to Neo4j");
        
        Ok(Self { graph })
    }
    
    /// Test the database connection
    pub async fn test_connection(&self) -> Result<()> {
        let result = self.graph.run(query("RETURN 1")).await?;
        let exists = result.into_iter().next().is_some();
        
        if exists {
            debug!("Database connection test successful");
        } else {
            warn!("Database connection test failed: no results returned");
        }
        
        Ok(())
    }
    
    /// Store a genome in the database
    pub async fn store_genome(&self, genome: &Genome) -> Result<()> {
        info!("Storing genome {} in database", genome.id);
        
        // Create genome node
        let genome_query = query(
            "CREATE (g:Genome {
                genome_id: $genome_id,
                species_name: $species_name,
                related_group: $related_group,
                original_n_percent_avg: $original_n_percent_avg,
                total_processed_size_bp: $total_processed_size,
                assembly_index_avg: $assembly_index_avg
            }) RETURN g"
        )
        .param("genome_id", &genome.id)
        .param("species_name", &genome.species_name)
        .param("related_group", &genome.related_group)
        .param("original_n_percent_avg", genome.original_n_percent_avg)
        .param("total_processed_size", genome.total_processed_size as i64)
        .param("assembly_index_avg", genome.assembly_index_avg.unwrap_or(0.0));
        
        let result = self.graph.run(genome_query).await?;
        let genome_node = result.into_iter().next()
            .context("Failed to create genome node")?
            .get::<Node>("g")
            .context("Failed to retrieve created genome node")?;
        
        debug!("Created genome node with ID: {}", genome_node.id());
        
        // Create chromosome nodes and link to genome
        for chromosome in &genome.chromosomes {
            // Skip chromosomes without assembly index (not processed)
            if chromosome.assembly_index.is_none() {
                warn!("Skipping chromosome {} as it has no assembly index", chromosome.id);
                continue;
            }
            
            self.store_chromosome(genome, chromosome).await?;
        }
        
        info!("Successfully stored genome {} with {} chromosomes", genome.id, genome.chromosomes.len());
        
        Ok(())
    }
    
    /// Store a chromosome and its motifs in the database
    async fn store_chromosome(&self, genome: &Genome, chromosome: &Chromosome) -> Result<()> {
        info!("Storing chromosome {} in database", chromosome.id);
        
        // Create chromosome node
        let chromosome_query = query(
            "MATCH (g:Genome {genome_id: $genome_id})
            CREATE (c:Chromosome {
                chromosome_id: $chromosome_id,
                processed_length_bp: $processed_length,
                original_n_percent: $original_n_percent,
                assembly_index: $assembly_index,
                graph_depth: $graph_depth,
                avg_motif_length: $avg_motif_length
            })-[:PART_OF]->(g)
            RETURN c"
        )
        .param("genome_id", &genome.id)
        .param("chromosome_id", &chromosome.id)
        .param("processed_length", chromosome.processed_length as i64)
        .param("original_n_percent", chromosome.original_n_percent)
        .param("assembly_index", chromosome.assembly_index.unwrap_or(0) as i64)
        .param("graph_depth", chromosome.graph_depth.unwrap_or(0) as i64)
        .param("avg_motif_length", chromosome.avg_motif_length.unwrap_or(0.0));
        
        let result = self.graph.run(chromosome_query).await?;
        let chromosome_node = result.into_iter().next()
            .context("Failed to create chromosome node")?
            .get::<Node>("c")
            .context("Failed to retrieve created chromosome node")?;
        
        debug!("Created chromosome node with ID: {}", chromosome_node.id());
        
        Ok(())
    }
    
    /// Store a motif graph in the database
    pub async fn store_motif_graph(
        &self, 
        genome_id: &str, 
        chromosome_id: &str, 
        graph: &MotifGraph
    ) -> Result<()> {
        info!("Storing motif graph for chromosome {} of genome {}", chromosome_id, genome_id);
        
        // Create motif nodes
        let mut node_map = HashMap::new();
        
        for node_idx in graph.graph.node_indices() {
            let motif = &graph.graph[node_idx];
            
            // Create motif node
            let motif_query = query(
                "CREATE (m:Motif {
                    sequence: $sequence,
                    length: $length,
                    creation_iteration: $creation_iteration
                }) RETURN m"
            )
            .param("sequence", &motif.sequence)
            .param("length", motif.length as i64)
            .param("creation_iteration", motif.creation_iteration as i64);
            
            let result = self.graph.run(motif_query).await?;
            let motif_node = result.into_iter().next()
                .context("Failed to create motif node")?
                .get::<Node>("m")
                .context("Failed to retrieve created motif node")?;
            
            node_map.insert(node_idx.index(), motif_node.id());
            
            // Link motif to chromosome with frequency
            let link_query = query(
                "MATCH 
                    (m:Motif), 
                    (c:Chromosome {chromosome_id: $chromosome_id})-[:PART_OF]->(:Genome {genome_id: $genome_id})
                WHERE id(m) = $motif_id
                CREATE (m)-[:FOUND_IN {frequency: $frequency}]->(c)"
            )
            .param("chromosome_id", chromosome_id)
            .param("genome_id", genome_id)
            .param("motif_id", motif_node.id())
            .param("frequency", motif.frequency as i64);
            
            self.graph.run(link_query).await?;
        }
        
        // Create edges for motif derivations
        for edge_idx in graph.graph.edge_indices() {
            let (source, target) = graph.graph.edge_endpoints(edge_idx).unwrap();
            
            let source_id = node_map.get(&source.index()).unwrap();
            let target_id = node_map.get(&target.index()).unwrap();
            
            let edge_query = query(
                "MATCH (source:Motif), (target:Motif)
                WHERE id(source) = $source_id AND id(target) = $target_id
                CREATE (source)-[:DERIVED_FROM]->(target)"
            )
            .param("source_id", *source_id)
            .param("target_id", *target_id);
            
            self.graph.run(edge_query).await?;
        }
        
        info!(
            "Successfully stored motif graph with {} nodes and {} edges",
            graph.node_count, graph.edge_count
        );
        
        Ok(())
    }
} 