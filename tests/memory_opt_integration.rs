#[cfg(test)]
mod memory_optimization_integration_tests {
    use orbweaver::grammar::builder::GrammarBuilder;
    use orbweaver::encode::dna_2bit::EncodedBase;
    use orbweaver::encode::bitvec;
    use std::time::Instant;
    
    // Helper function to generate a test sequence
    fn generate_repeating_sequence(size: usize) -> Vec<EncodedBase> {
        let pattern = vec![
            // ACGTACGTACGT pattern
            EncodedBase(0), EncodedBase(1), EncodedBase(2), EncodedBase(3),
            EncodedBase(0), EncodedBase(1), EncodedBase(2), EncodedBase(3),
            EncodedBase(0), EncodedBase(1), EncodedBase(2), EncodedBase(3),
        ];
        
        let mut sequence = Vec::with_capacity(size);
        while sequence.len() < size {
            let remaining = size - sequence.len();
            let to_add = if remaining < pattern.len() {
                &pattern[0..remaining]
            } else {
                &pattern
            };
            sequence.extend_from_slice(to_add);
        }
        
        sequence
    }
    
    #[test]
    fn test_memory_optimized_processing() {
        // Generate a test sequence
        let sequence_size = 50000;
        println!("Generating test sequence of {} bases...", sequence_size);
        let sequence = generate_repeating_sequence(sequence_size);
        
        // Test standard processing
        let standard_start = Instant::now();
        let mut standard_builder = GrammarBuilder::new(2, false);
        standard_builder.build_grammar(&sequence).expect("Standard grammar building failed");
        let standard_duration = standard_start.elapsed();
        let (standard_seq, standard_rules) = standard_builder.get_grammar();
        
        println!("Standard processing: {:?}, rules: {}", standard_duration, standard_rules.len());
        
        // Test streaming processing
        let streaming_start = Instant::now();
        let mut streaming_builder = GrammarBuilder::new(2, false)
            .enable_streaming_mode()
            .with_max_rules(1000); // Limit rule count
        
        // Process in chunks
        let chunk_size = 5000;
        for chunk_start in (0..sequence.len()).step_by(chunk_size) {
            let end = std::cmp::min(chunk_start + chunk_size, sequence.len());
            let chunk = &sequence[chunk_start..end];
            
            streaming_builder.process_sequence_chunk(chunk)
                .expect("Streaming processing failed");
        }
        
        streaming_builder.finalize_grammar().expect("Failed to finalize streaming grammar");
        let streaming_duration = streaming_start.elapsed();
        let (streaming_seq, streaming_rules) = streaming_builder.get_grammar();
        
        println!("Streaming processing: {:?}, rules: {}", streaming_duration, streaming_rules.len());
        
        // Verify that both approaches produced valid results
        assert!(!standard_seq.is_empty(), "Standard sequence should not be empty");
        assert!(!standard_rules.is_empty(), "Standard rules should not be empty");
        assert!(!streaming_seq.is_empty(), "Streaming sequence should not be empty");
        assert!(!streaming_rules.is_empty(), "Streaming rules should not be empty");
        
        // Check that rule eviction worked (streaming should have <= max_rules)
        assert!(streaming_rules.len() <= 1000,
                "Streaming should have at most 1000 rules, has {}", streaming_rules.len());
        
        // Calculate memory savings from 2-bit encoding
        let raw_size = sequence_size;
        let encoded_size = std::mem::size_of::<EncodedBase>() * sequence.len();
        let (saving_pct, _) = bitvec::estimate_memory_savings(raw_size, encoded_size);
        
        println!("Memory savings from 2-bit encoding: {:.1}%", saving_pct);
        
        // Memory savings should be positive
        assert!(saving_pct > 0.0, "Memory savings should be positive");
    }
} 