//! Integration tests for the minimizer module.
//!
//! These tests verify the minimizer functionality from an external user's perspective,
//! testing only the public API.

use chloro_index::minimizer::{get_canonical, MinimizerIterator};
use chloro_index::encode_dna;

#[test]
fn test_canonical_form_public_api() {
    // Test that AAAA and TTTT (reverse complements) have the same canonical form
    let aaaa = encode_dna("AAAA").unwrap();
    let tttt = encode_dna("TTTT").unwrap();
    
    assert_eq!(get_canonical(aaaa, 4), get_canonical(tttt, 4));
}

#[test]
fn test_minimizer_iterator_public_api() {
    let seq = "ACGTACGTACGT";
    let iter = MinimizerIterator::new(seq, 4, 3);
    
    let results: Vec<_> = iter.collect();
    assert!(!results.is_empty());
    
    // Verify all results are valid
    for (pos, hash) in &results {
        assert!(*pos < seq.len());
        assert!(*hash < u64::MAX);
    }
}

#[test]
fn test_minimizer_real_world_sequence() {
    // Test with a longer, more realistic sequence
    let seq = "ACGTACGTACGTACGTACGTACGTACGTACGT";
    let iter = MinimizerIterator::new(seq, 8, 5);
    
    let results: Vec<_> = iter.collect();
    
    // Should extract multiple minimizers
    assert!(!results.is_empty());
    
    // Verify no duplicate consecutive minimizers
    let mut prev_hash = None;
    for (_, hash) in &results {
        if let Some(prev) = prev_hash {
            assert_ne!(prev, *hash, "Consecutive minimizers should be different");
        }
        prev_hash = Some(*hash);
    }
}

#[test]
fn test_minimizer_edge_cases() {
    // Test with minimum valid parameters
    let seq = "ACGT";
    let iter = MinimizerIterator::new(seq, 4, 1);
    let results: Vec<_> = iter.collect();
    assert_eq!(results.len(), 1);
    
    // Test with sequence exactly matching window size
    let seq2 = "ACGTACGT";
    let iter2 = MinimizerIterator::new(seq2, 4, 2);
    let results2: Vec<_> = iter2.collect();
    assert!(!results2.is_empty());
}

