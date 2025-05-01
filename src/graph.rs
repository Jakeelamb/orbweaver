use anyhow::Result;
use log::{debug, info};
use petgraph::graph::{DiGraph, NodeIndex};
use petgraph::visit::Dfs;
use std::collections::HashMap;
use std::rc::Rc;

use crate::slp::{Motif, SLP};

/// Represents a node in the motif graph
#[derive(Debug, Clone)]
pub struct MotifNode {
    /// The sequence of the motif
    pub sequence: String,
    
    /// Frequency of the motif in the chromosome
    pub frequency: usize,
    
    /// Length of the motif
    pub length: usize,
    
    /// Iteration when this motif was created (0 for base nucleotides)
    pub creation_iteration: usize,
}

/// Edge types in the motif graph
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EdgeType {
    /// Indicates that one motif is derived from two other motifs
    DerivedFrom,
}

/// Represents a graph of motifs, based on an SLP
pub struct MotifGraph {
    /// The directed graph structure
    pub graph: DiGraph<MotifNode, EdgeType>,
    
    /// Map from motif sequences to node indices
    pub motif_to_node: HashMap<String, NodeIndex>,
    
    /// The depth of the graph (longest path from a leaf to root)
    pub depth: usize,
    
    /// Number of nodes in the graph
    pub node_count: usize,
    
    /// Number of edges in the graph
    pub edge_count: usize,
    
    /// Average motif length
    pub avg_motif_length: f64,
}

impl MotifGraph {
    /// Create a new motif graph from an SLP
    pub fn from_slp(slp: &SLP) -> Result<Self> {
        let mut graph = DiGraph::<MotifNode, EdgeType>::new();
        let mut motif_to_node = HashMap::new();
        
        info!("Converting SLP to motif graph with {} rules", slp.rules.len());
        
        // First, add all motifs as nodes
        for (sequence, motif) in &slp.motifs {
            let node = MotifNode {
                sequence: sequence.clone(),
                frequency: motif.frequency,
                length: motif.length,
                creation_iteration: motif.creation_iteration,
            };
            
            let node_idx = graph.add_node(node);
            motif_to_node.insert(sequence.clone(), node_idx);
        }
        
        // Then, add all rules as edges
        for rule in &slp.rules {
            let lhs_idx = motif_to_node.get(&rule.lhs.sequence).unwrap();
            let rhs1_idx = motif_to_node.get(&rule.rhs1.sequence).unwrap();
            let rhs2_idx = motif_to_node.get(&rule.rhs2.sequence).unwrap();
            
            // Add edges
            graph.add_edge(*lhs_idx, *rhs1_idx, EdgeType::DerivedFrom);
            graph.add_edge(*lhs_idx, *rhs2_idx, EdgeType::DerivedFrom);
        }
        
        // Calculate graph metrics
        let node_count = graph.node_count();
        let edge_count = graph.edge_count();
        
        // Calculate average motif length
        let total_length: usize = graph.node_weights().map(|n| n.length).sum();
        let avg_motif_length = if node_count > 0 {
            total_length as f64 / node_count as f64
        } else {
            0.0
        };
        
        // Calculate graph depth
        let depth = Self::calculate_graph_depth(&graph, &motif_to_node);
        
        debug!(
            "Created motif graph: {} nodes, {} edges, depth: {}, avg motif length: {:.2}",
            node_count, edge_count, depth, avg_motif_length
        );
        
        Ok(Self {
            graph,
            motif_to_node,
            depth,
            node_count,
            edge_count,
            avg_motif_length,
        })
    }
    
    /// Calculate the depth of the graph (longest path from a leaf to root)
    fn calculate_graph_depth(
        graph: &DiGraph<MotifNode, EdgeType>,
        motif_to_node: &HashMap<String, NodeIndex>,
    ) -> usize {
        // Find the nodes with no outgoing edges (leaf nodes)
        let leaf_nodes: Vec<NodeIndex> = graph
            .node_indices()
            .filter(|&n| graph.neighbors_directed(n, petgraph::Direction::Outgoing).count() == 0)
            .collect();
        
        // For each leaf node, find the longest path to any node with no incoming edges
        let mut max_depth = 0;
        
        for &start_node in &leaf_nodes {
            let mut visited = vec![false; graph.node_count()];
            let mut depths = vec![0; graph.node_count()];
            let mut stack = vec![(start_node, 0)];
            
            while let Some((node, depth)) = stack.pop() {
                let node_idx = node.index();
                
                if visited[node_idx] {
                    continue;
                }
                
                visited[node_idx] = true;
                depths[node_idx] = depth;
                
                for neighbor in graph.neighbors_directed(node, petgraph::Direction::Incoming) {
                    stack.push((neighbor, depth + 1));
                }
            }
            
            // Find the maximum depth reached
            if let Some(&node_depth) = depths.iter().max() {
                max_depth = max_depth.max(node_depth);
            }
        }
        
        max_depth
    }
    
    /// Calculate the in-degree distribution of the graph
    pub fn calculate_in_degree_distribution(&self) -> HashMap<usize, usize> {
        let mut distribution = HashMap::new();
        
        for node in self.graph.node_indices() {
            let in_degree = self.graph.neighbors_directed(node, petgraph::Direction::Incoming).count();
            *distribution.entry(in_degree).or_insert(0) += 1;
        }
        
        distribution
    }
    
    /// Calculate the out-degree distribution of the graph
    pub fn calculate_out_degree_distribution(&self) -> HashMap<usize, usize> {
        let mut distribution = HashMap::new();
        
        for node in self.graph.node_indices() {
            let out_degree = self.graph.neighbors_directed(node, petgraph::Direction::Outgoing).count();
            *distribution.entry(out_degree).or_insert(0) += 1;
        }
        
        distribution
    }
} 