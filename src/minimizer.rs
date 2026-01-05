//! Minimizer extraction algorithm for genomic sequences.
//!
//! This module provides functionality for extracting minimizers from DNA sequences,
//! which are used in bioinformatics for efficient sequence comparison and indexing.

use crate::{decode_dna, encode_dna, reverse_complement};

/// Computes the canonical form of a k-mer.
///
/// The canonical form is the lexicographically smaller value between the k-mer
/// and its reverse complement. This ensures that a k-mer and its reverse
/// complement map to the same canonical representation.
///
/// # Arguments
/// * `kmer` - Encoded k-mer as u64
/// * `k` - Length of the k-mer (number of bases)
///
/// # Returns
/// The canonical form of the k-mer (the smaller of the original and reverse complement)
///
/// # Examples
/// ```
/// use chloro_index::minimizer::get_canonical;
/// use chloro_index::encode_dna;
///
/// let aaaa = encode_dna("AAAA").unwrap();
/// let tttt = encode_dna("TTTT").unwrap();
/// assert_eq!(get_canonical(aaaa, 4), get_canonical(tttt, 4));
/// ```
pub fn get_canonical(kmer: u64, k: u8) -> u64 {
    // Decode the k-mer to a string
    let kmer_str = decode_dna(kmer, k as usize);
    
    // Compute the reverse complement
    let revcomp_str = reverse_complement(&kmer_str);
    
    // Encode the reverse complement
    let revcomp = encode_dna(&revcomp_str).unwrap_or(u64::MAX);
    
    // Return the smaller value
    kmer.min(revcomp)
}

/// Iterator that extracts minimizers from a DNA sequence.
///
/// A minimizer is the lexicographically smallest k-mer in a window of `w`
/// consecutive k-mers. This iterator slides a window of size `w` over the
/// sequence and yields the minimizer for each window position.
///
/// # Examples
/// ```
/// use chloro_index::minimizer::MinimizerIterator;
///
/// let seq = "ACGTACGTACGT";
/// let mut iter = MinimizerIterator::new(seq, 4, 3);
/// // Iterate over minimizers...
/// ```
pub struct MinimizerIterator {
    sequence: String,
    k: u8,
    w: usize,
    position: usize,
    last_minimizer: Option<u64>,
}

impl MinimizerIterator {
    /// Creates a new minimizer iterator.
    ///
    /// # Arguments
    /// * `sequence` - DNA sequence string
    /// * `k` - K-mer size (number of bases per k-mer)
    /// * `w` - Window size (number of consecutive k-mers in a window)
    ///
    /// # Panics
    /// Panics if `k` is 0 or greater than 32, or if `w` is 0.
    pub fn new(sequence: impl Into<String>, k: u8, w: usize) -> Self {
        assert!(k > 0 && k <= 32, "k must be between 1 and 32");
        assert!(w > 0, "w must be greater than 0");
        
        MinimizerIterator {
            sequence: sequence.into(),
            k,
            w,
            position: 0,
            last_minimizer: None,
        }
    }
    
    /// Finds the minimizer in a window starting at the given position.
    ///
    /// Returns `None` if the window extends beyond the sequence length.
    fn find_minimizer_in_window(&self, start_pos: usize) -> Option<(usize, u64)> {
        let seq_len = self.sequence.len();
        let k = self.k as usize;
        
        // Check if we have enough sequence for a full window
        let window_end = start_pos + (self.w - 1) + k;
        if window_end > seq_len {
            return None;
        }
        
        let mut minimizer = None;
        let mut minimizer_pos = 0;
        let mut minimizer_canonical = u64::MAX;
        
        // Iterate through all k-mers in the window
        for i in 0..self.w {
            let kmer_start = start_pos + i;
            let kmer_end = kmer_start + k;
            
            if kmer_end > seq_len {
                break;
            }
            
            // Extract and encode the k-mer
            if let Ok(kmer) = encode_dna(&self.sequence[kmer_start..kmer_end]) {
                let canonical = get_canonical(kmer, self.k);
                
                // Update minimizer if this is smaller
                if canonical < minimizer_canonical {
                    minimizer_canonical = canonical;
                    minimizer = Some(kmer);
                    minimizer_pos = kmer_start;
                }
            }
        }
        
        minimizer.map(|m| (minimizer_pos, get_canonical(m, self.k)))
    }
}

impl Iterator for MinimizerIterator {
    type Item = (usize, u64);
    
    /// Returns the next minimizer as `(position, canonical_hash)`.
    ///
    /// The position is 0-based and indicates where the minimizer k-mer starts
    /// in the original sequence. Only yields minimizers that differ from the
    /// previous one to avoid duplicates.
    fn next(&mut self) -> Option<Self::Item> {
        loop {
            if let Some((pos, canonical)) = self.find_minimizer_in_window(self.position) {
                // Only yield if different from the last minimizer
                if self.last_minimizer != Some(canonical) {
                    self.last_minimizer = Some(canonical);
                    self.position += 1;
                    return Some((pos, canonical));
                }
                // Skip duplicate minimizers
                self.position += 1;
            } else {
                // No more windows
                return None;
            }
        }
    }
}
