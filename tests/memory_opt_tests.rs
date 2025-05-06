#[cfg(test)]
mod memory_optimization_tests {
    use orbweaver::encode::dna_2bit::EncodedBase;
    use orbweaver::grammar::builder::GrammarBuilder;
    
    use std::time::Instant;
    
    fn generate_test_sequence(size: usize) -> Vec<EncodedBase> {
        // Generate a random-ish DNA sequence
        let mut sequence = Vec::with_capacity(size);
        for i in 0..size {
            // Simple pattern: cycle through ACGT
            let base = match i % 4 {
                0 => EncodedBase(0), // A
                1 => EncodedBase(1), // C
                2 => EncodedBase(2), // G
                3 => EncodedBase(3), // T
                _ => unreachable!(),
            };
            sequence.push(base);
        }
        sequence
    }
    
    // Test that the rule eviction strategy properly evicts rules
    #[test]
    fn test_rule_eviction() {
        // Generate a sequence with many patterns that could create rules
        let sequence = vec![
            EncodedBase(0), EncodedBase(1), // AC
            EncodedBase(0), EncodedBase(1), // AC
            EncodedBase(2), EncodedBase(3), // GT
            EncodedBase(2), EncodedBase(3), // GT
            EncodedBase(0), EncodedBase(2), // AG
            EncodedBase(0), EncodedBase(2), // AG
            EncodedBase(1), EncodedBase(3), // CT
            EncodedBase(1), EncodedBase(3), // CT
            EncodedBase(0), EncodedBase(3), // AT
            EncodedBase(0), EncodedBase(3), // AT
        ];
        
        // Create a builder with a small max rule count (4) and min usage (2)
        let min_rule_usage = 2;
        let max_rule_count = 4;
        let mut builder = GrammarBuilder::new(min_rule_usage, false)
            .with_max_rules(max_rule_count);
        
        // Build the grammar
        builder.build_grammar(&sequence).expect("Failed to build grammar");
        
        // Get the resulting grammar
        let (_, rules) = builder.get_grammar();
        
        // Verify that we have at most max_rule_count rules
        assert!(rules.len() <= max_rule_count, 
                "Expected at most {} rules, got {}", max_rule_count, rules.len());
        
        // Each rule should have usage_count >= min_rule_usage
        for (rule_id, rule) in rules {
            assert!(rule.usage_count >= min_rule_usage,
                    "Rule {} has usage count {} which is less than minimum {}", 
                    rule_id, rule.usage_count, min_rule_usage);
        }
        
        println!("Rule eviction test passed with {} rules (max was {})", 
                 rules.len(), max_rule_count);
    }
    
    // Test the streaming processing capability
    #[test]
    fn test_streaming_mode() {
        // Generate a larger sequence
        let sequence_size = 10000;
        let sequence = generate_test_sequence(sequence_size);
        
        // Define chunk size for streaming
        let chunk_size = 1000;
        
        // Create a builder with streaming mode enabled
        let mut streaming_builder = GrammarBuilder::new(2, false)
            .enable_streaming_mode();
        
        // Process the sequence in chunks
        let start_time = Instant::now();
        for chunk_start in (0..sequence.len()).step_by(chunk_size) {
            let end = std::cmp::min(chunk_start + chunk_size, sequence.len());
            let chunk = &sequence[chunk_start..end];
            
            streaming_builder.process_sequence_chunk(chunk)
                .expect("Failed to process chunk");
        }
        
        // Finalize the grammar
        streaming_builder.finalize_grammar().expect("Failed to finalize grammar");
        let streaming_time = start_time.elapsed();
        
        // Get the resulting grammar
        let (streaming_sequence, streaming_rules) = streaming_builder.get_grammar();
        
        // Create a builder without streaming for comparison
        let mut standard_builder = GrammarBuilder::new(2, false);
        
        // Process the entire sequence at once
        let start_time = Instant::now();
        standard_builder.build_grammar(&sequence).expect("Failed to build grammar");
        let standard_time = start_time.elapsed();
        
        // Get the resulting grammar
        let (standard_sequence, standard_rules) = standard_builder.get_grammar();
        
        // Print results for comparison
        println!("Streaming mode processed {} bases in {:?}", sequence_size, streaming_time);
        println!("Standard mode processed {} bases in {:?}", sequence_size, standard_time);
        println!("Streaming generated {} rules", streaming_rules.len());
        println!("Standard generated {} rules", standard_rules.len());
        
        // The results may not be identical due to different processing order,
        // but both should produce valid grammars
        assert!(!streaming_rules.is_empty(), "Streaming mode produced no rules");
        assert!(!streaming_sequence.is_empty(), "Streaming mode produced empty sequence");
    }
    
    // Test the 2-bit encoding memory savings
    #[test]
    fn test_2bit_encoding_memory_savings() {
        // Generate a sequence
        let sequence_size = 1000;
        
        // Create raw ASCII sequence
        let raw_sequence: Vec<u8> = (0..sequence_size)
            .map(|i| match i % 4 {
                0 => b'A',
                1 => b'C',
                2 => b'G',
                3 => b'T',
                _ => unreachable!(),
            })
            .collect();
        
        // Create 2-bit encoded sequence
        let encoded_sequence: Vec<EncodedBase> = raw_sequence.iter()
            .filter_map(|&b| EncodedBase::from_base(b))
            .collect();
        
        // Calculate memory usage
        let raw_memory = std::mem::size_of_val(&raw_sequence[0]) * raw_sequence.len();
        let encoded_memory = std::mem::size_of_val(&encoded_sequence[0]) * encoded_sequence.len();
        let saving_percentage = 100.0 * (1.0 - (encoded_memory as f64 / raw_memory as f64));
        
        println!("Raw memory usage: {} bytes", raw_memory);
        println!("Encoded memory usage: {} bytes", encoded_memory);
        println!("Memory saving: {:.1}%", saving_percentage);
        
        // Verify we're saving memory
        assert!(encoded_memory < raw_memory, 
                "2-bit encoding should use less memory than raw ASCII");
    }
} 