//! Genomic sequence index using B-Tree data structure.
//!
//! This module provides a B-Tree implementation optimized for indexing
//! genomic k-mers and their locations. DNA sequences are encoded as u64
//! values (2 bits per base) for efficient storage and lookup.

use serde::{Deserialize, Serialize};
use std::fs::File;
use std::io::{self, BufReader, BufWriter};
use std::mem;

// ============================================================================
// Public Types
// ============================================================================

/// Represents a genomic location with chromosome ID and position.
///
/// Positions are 1-based (following biological convention).
#[derive(Debug, Serialize, Deserialize, Clone, PartialEq, Eq)]
pub struct GenomicLocation {
    /// Chromosome identifier (typically 1-24 for human genome)
    pub chromosome: u8,
    /// 1-based position on the chromosome
    pub position: u32,
}

/// B-Tree index for genomic k-mer lookup.
///
/// This B-Tree stores encoded DNA sequences (k-mers) as keys and their
/// genomic locations as values. The tree maintains balance through node
/// splitting when nodes exceed capacity.
#[derive(Serialize, Deserialize)]
pub struct BTree {
    root: BTreeNode,
    /// Minimum degree of the B-Tree (t). Each node has at most 2t-1 keys.
    t: usize,
}

// ============================================================================
// Internal Types
// ============================================================================

/// Internal node structure for the B-Tree.
#[derive(Debug, Serialize, Deserialize, Clone)]
struct BTreeNode {
    /// Sorted array of encoded DNA sequence keys
    keys: Vec<u64>,
    /// Corresponding genomic locations for each key
    values: Vec<Vec<GenomicLocation>>,
    /// Child nodes (only used for non-leaf nodes)
    children: Vec<BTreeNode>,
    /// Whether this node is a leaf (has no children)
    is_leaf: bool,
}

impl BTreeNode {
    /// Creates a new B-Tree node.
    fn new(is_leaf: bool) -> Self {
        BTreeNode {
            keys: Vec::new(),
            values: Vec::new(),
            children: Vec::new(),
            is_leaf,
        }
    }
}

impl BTree {
    // ========================================================================
    // Public API
    // ========================================================================

    /// Creates a new B-Tree with the specified minimum degree.
    ///
    /// # Arguments
    /// * `t` - Minimum degree. Each node can have at most 2t-1 keys.
    ///          Must be at least 2 for a valid B-Tree.
    pub fn new(t: usize) -> Self {
        BTree {
            root: BTreeNode::new(true),
            t,
        }
    }

    /// Searches for a key in the B-Tree.
    ///
    /// # Arguments
    /// * `key` - Encoded DNA sequence to search for
    ///
    /// # Returns
    /// * `Some(Vec<GenomicLocation>)` - All genomic locations where this k-mer appears
    /// * `None` - Key not found
    pub fn search(&self, key: u64) -> Option<Vec<GenomicLocation>> {
        self.search_node(&self.root, key)
    }

    /// Inserts a key-value pair into the B-Tree.
    ///
    /// If the key already exists, the new locations are appended to the
    /// existing list of locations.
    ///
    /// # Arguments
    /// * `key` - Encoded DNA sequence
    /// * `locations` - Genomic locations where this k-mer appears
    pub fn insert(&mut self, key: u64, locations: Vec<GenomicLocation>) {
        if self.root.keys.len() == 2 * self.t - 1 {
            let mut new_root = BTreeNode::new(false);
            let old_root = mem::replace(&mut self.root, BTreeNode::new(true));
            new_root.children.push(old_root);
            
            Self::split_child(&mut new_root, 0, self.t);
            self.root = new_root;
        }
        
        Self::insert_non_full(&mut self.root, key, locations, self.t);
    }

    /// Indexes a DNA sequence by extracting all k-mers using a sliding window.
    ///
    /// This method processes the entire sequence, extracts all k-mers of
    /// length `K_MER_SIZE`, encodes them, and inserts them into the index
    /// with their genomic positions.
    ///
    /// # Arguments
    /// * `chr` - Chromosome identifier
    /// * `seq` - DNA sequence string (must contain only A, C, G, T)
    ///
    /// # Returns
    /// Number of k-mers successfully indexed
    pub fn bulk_insert_sequence(&mut self, chr: u8, seq: &str) -> usize {
        const K_MER_SIZE: usize = 32;
        
        if seq.len() < K_MER_SIZE {
            return 0;
        }

        let mut count = 0;
        for (i, window_bytes) in seq.as_bytes().windows(K_MER_SIZE).enumerate() {
            if let Ok(fragment) = std::str::from_utf8(window_bytes) {
                if let Ok(key) = encode_dna(fragment) {
                    let loc = GenomicLocation {
                        chromosome: chr,
                        position: (i + 1) as u32, // 1-based position
                    };
                    self.insert(key, vec![loc]);
                    count += 1;
                }
            }
        }
        count
    }

    /// Performs fuzzy search for a DNA query, finding exact matches and all 1-mismatch variants.
    ///
    /// This method first searches for an exact match of the query. Then it generates
    /// all possible 1-mismatch variants by substituting each position with A, C, G, or T
    /// (excluding the original character) and searches for each variant.
    ///
    /// # Arguments
    /// * `query` - DNA sequence string to search for (must contain only A, C, G, T)
    ///
    /// # Returns
    /// Vector of tuples `(sequence, locations)` where:
    /// * `sequence` - The matched DNA sequence (exact match or 1-mismatch variant)
    /// * `locations` - All genomic locations where this sequence appears
    ///
    /// The results may include duplicates if the same variant appears multiple times
    /// in the search process. Invalid sequences are silently skipped.
    pub fn search_fuzzy(&self, query: &str) -> Vec<(String, Vec<GenomicLocation>)> {
        let mut results = Vec::new();
        
        // First, search for exact match
        if let Ok(key) = encode_dna(query) {
            if let Some(locations) = self.search(key) {
                results.push((query.to_string(), locations));
            }
        }
        
        // Generate all 1-mismatch variants
        let query_chars: Vec<char> = query.chars().collect();
        let bases = ['A', 'C', 'G', 'T'];
        
        for i in 0..query_chars.len() {
            let original_char = query_chars[i].to_ascii_uppercase();
            
            // Try each base except the original
            for &base in &bases {
                if base != original_char {
                    // Create mutated query
                    let mut mutated_chars = query_chars.clone();
                    mutated_chars[i] = base;
                    let mutated_query: String = mutated_chars.iter().collect();
                    
                    // Encode and search
                    if let Ok(key) = encode_dna(&mutated_query) {
                        if let Some(locations) = self.search(key) {
                            results.push((mutated_query, locations));
                        }
                    }
                }
            }
        }
        
        results
    }

    /// Saves the B-Tree to a file using binary serialization.
    ///
    /// # Arguments
    /// * `filename` - Path to the output file
    ///
    /// # Errors
    /// Returns `io::Error` if file creation or serialization fails.
    pub fn save_to_file(&self, filename: &str) -> io::Result<()> {
        let file = File::create(filename)?;
        let writer = BufWriter::new(file);
        
        bincode::serialize_into(writer, self)
            .map_err(|e| io::Error::new(io::ErrorKind::Other, e))?;
            
        Ok(())
    }

    /// Loads a B-Tree from a file.
    ///
    /// # Arguments
    /// * `filename` - Path to the input file
    ///
    /// # Errors
    /// Returns `io::Error` if file opening or deserialization fails.
    pub fn load_from_file(filename: &str) -> io::Result<Self> {
        let file = File::open(filename)?;
        let reader = BufReader::new(file);

        let tree: BTree = bincode::deserialize_from(reader)
            .map_err(|e| io::Error::new(io::ErrorKind::Other, e))?;

        Ok(tree)
    }

    /// Prints the tree structure for debugging purposes.
    pub fn print_tree(&self) {
        self.print_node(&self.root, 0);
    }

    // ========================================================================
    // Internal Implementation
    // ========================================================================

    /// Recursively searches for a key starting from the given node.
    fn search_node(&self, node: &BTreeNode, key: u64) -> Option<Vec<GenomicLocation>> {
        let mut i = 0;
        while i < node.keys.len() && key > node.keys[i] {
            i += 1;
        }

        if i < node.keys.len() && key == node.keys[i] {
            return Some(node.values[i].clone());
        }

        if node.is_leaf {
            return None;
        }

        self.search_node(&node.children[i], key)
    }

    /// Inserts a key into a non-full node.
    ///
    /// This is the core insertion logic that handles both leaf and internal nodes.
    /// For existing keys, locations are merged rather than replaced.
    fn insert_non_full(node: &mut BTreeNode, key: u64, locations: Vec<GenomicLocation>, t: usize) {
        let mut i = 0;
        while i < node.keys.len() && key > node.keys[i] {
            i += 1;
        }

        if node.is_leaf {
            if i < node.keys.len() && key == node.keys[i] {
                // Key exists: merge locations
                node.values[i].extend(locations);
            } else {
                // Key doesn't exist: insert new entry
                node.keys.insert(i, key);
                node.values.insert(i, locations);
            }
        } else {
            // Internal node
            if i < node.keys.len() && key == node.keys[i] {
                // Key exists at this level: merge locations
                node.values[i].extend(locations);
                return;
            }

            // Check if child needs splitting before recursion
            if node.children[i].keys.len() == 2 * t - 1 {
                Self::split_child(node, i, t);
                // After split, adjust index if necessary
                if key > node.keys[i] {
                    i += 1;
                } else if key == node.keys[i] {
                    node.values[i].extend(locations);
                    return;
                }
            }

            Self::insert_non_full(&mut node.children[i], key, locations, t);
        }
    }

    /// Splits a full child node, promoting the middle key to the parent.
    ///
    /// This maintains the B-Tree property that nodes have between t-1 and 2t-1 keys.
    fn split_child(parent: &mut BTreeNode, i: usize, t: usize) {
        let full_child = &mut parent.children[i];
        let mut new_child = BTreeNode::new(full_child.is_leaf);
        
        // Move upper half of keys and values to new child
        new_child.keys.extend(full_child.keys.drain(t..));
        new_child.values.extend(full_child.values.drain(t..));
        
        // Move upper half of children if not a leaf
        if !full_child.is_leaf {
            new_child.children.extend(full_child.children.drain(t..));
        }

        // Promote middle key to parent
        let mid_key = full_child.keys.pop().unwrap();
        let mid_value = full_child.values.pop().unwrap();

        parent.keys.insert(i, mid_key);
        parent.values.insert(i, mid_value);
        parent.children.insert(i + 1, new_child);
    }

    /// Recursively prints node contents for debugging.
    fn print_node(&self, node: &BTreeNode, level: usize) {
        let indent = "  ".repeat(level);
        let pairs: Vec<String> = node.keys.iter()
            .zip(node.values.iter())
            .map(|(k, locs)| {
                let loc_strs: Vec<String> = locs.iter()
                    .map(|loc| format!("chr{}:{}", loc.chromosome, loc.position))
                    .collect();
                format!("{}: [{}]", k, loc_strs.join(", "))
            })
            .collect();
            
        println!("{}Level {}: {:?}", indent, level, pairs);
        
        if !node.is_leaf {
            for child in &node.children {
                self.print_node(child, level + 1);
            }
        }
    }
}

// ============================================================================
// DNA Encoding/Decoding
// ============================================================================

/// Maximum length of DNA sequence that can be encoded into a u64.
///
/// With 2 bits per base, a u64 can store up to 32 bases.
pub const MAX_DNA_LENGTH: usize = 32;

const BITS_PER_BASE: u32 = 2;

/// Encodes a DNA sequence into a u64 integer.
///
/// Each base is encoded using 2 bits: A=00, C=01, G=10, T=11.
/// The sequence is packed from left to right, with the first base
/// occupying the most significant bits.
///
/// # Arguments
/// * `dna` - DNA sequence string (case-insensitive)
///
/// # Returns
/// * `Ok(u64)` - Encoded sequence
/// * `Err(String)` - Error message if sequence is invalid or too long
///
/// # Examples
/// ```
/// use chloro_index::encode_dna;
/// assert_eq!(encode_dna("ACGT").unwrap(), 0b00011011);
/// ```
pub fn encode_dna(dna: &str) -> Result<u64, String> {
    if dna.is_empty() {
        return Err("Empty DNA sequence".to_string());
    }
    
    if dna.len() > MAX_DNA_LENGTH {
        return Err(format!(
            "Sequence too long: {} bases (maximum: {})",
            dna.len(),
            MAX_DNA_LENGTH
        ));
    }
    
    let mut encoded: u64 = 0;
    for base in dna.chars() {
        encoded <<= BITS_PER_BASE;
        encoded |= encode_base(base)? as u64;
    }
    
    Ok(encoded)
}

/// Decodes a u64 back into a DNA sequence string.
///
/// The length parameter is required because leading zeros in the
/// encoded value are indistinguishable from 'A' bases.
///
/// # Arguments
/// * `value` - Encoded DNA sequence
/// * `length` - Original sequence length
///
/// # Returns
/// Decoded DNA sequence string
pub fn decode_dna(mut value: u64, length: usize) -> String {
    if length == 0 {
        return String::new();
    }
    
    let mut sequence = String::with_capacity(length);
    for _ in 0..length {
        let base_bits = (value & 0b11) as u8;
        sequence.push(decode_base(base_bits));
        value >>= BITS_PER_BASE;
    }
    sequence.chars().rev().collect()
}

/// Encodes a single DNA base character to its 2-bit representation.
///
/// # Arguments
/// * `base` - DNA base character (A, C, G, or T, case-insensitive)
///
/// # Returns
/// * `Ok(u8)` - 2-bit encoding (0-3)
/// * `Err(String)` - Error if character is not a valid DNA base
fn encode_base(base: char) -> Result<u8, String> {
    match base.to_ascii_uppercase() {
        'A' => Ok(0b00),
        'C' => Ok(0b01),
        'G' => Ok(0b10),
        'T' => Ok(0b11),
        _ => Err(format!("Invalid DNA base: '{}' (expected A, C, G, or T)", base)),
    }
}

/// Decodes a 2-bit value back to its DNA base character.
///
/// # Arguments
/// * `bits` - 2-bit encoding (only lower 2 bits are used)
///
/// # Returns
/// DNA base character (A, C, G, or T)
fn decode_base(bits: u8) -> char {
    match bits & 0b11 {
        0b00 => 'A',
        0b01 => 'C',
        0b10 => 'G',
        0b11 => 'T',
        _ => unreachable!("2-bit value cannot exceed 0b11"),
    }
}

/// Computes the reverse complement of a DNA sequence.
///
/// This function reverses the input string and maps each character to its
/// complement: A->T, T->A, C->G, G->C, N->N.
///
/// Uses SIMD optimizations when available for better performance on long sequences.
///
/// # Arguments
/// * `dna` - DNA sequence string
///
/// # Returns
/// The reverse complement of the input sequence
pub fn reverse_complement(dna: &str) -> String {
    #[cfg(target_arch = "x86_64")]
    {
        if is_x86_feature_detected!("ssse3") && dna.len() >= 16 {
            return unsafe { reverse_complement_simd(dna) };
        }
    }
    
    // Fallback to scalar implementation
    reverse_complement_scalar(dna)
}

/// Scalar implementation of reverse complement (fallback).
fn reverse_complement_scalar(dna: &str) -> String {
    dna.chars()
        .rev()
        .map(|c| match c.to_ascii_uppercase() {
            'A' => 'T',
            'T' => 'A',
            'C' => 'G',
            'G' => 'C',
            'N' => 'N',
            _ => c, // Keep other characters as-is
        })
        .collect()
}

/// SIMD-optimized reverse complement using SSSE3 shuffle instructions.
///
/// This implementation uses SIMD to:
/// 1. Process 16 bytes at a time using SSSE3
/// 2. Use pshufb (_mm_shuffle_epi8) to reverse byte order within each chunk
/// 3. Use lookup table for complement mapping
///
/// # Safety
/// This function is unsafe because it uses SIMD intrinsics.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "ssse3")]
unsafe fn reverse_complement_simd(dna: &str) -> String {
    use std::arch::x86_64::*;
    
    let bytes = dna.as_bytes();
    let len = bytes.len();
    let mut result = Vec::with_capacity(len);
    
    // Complement lookup table: maps ASCII to complement
    // A(65) -> T(84), T(84) -> A(65), C(67) -> G(71), G(71) -> C(67), N(78) -> N(78)
    // For lowercase: a(97) -> t(116), t(116) -> a(97), c(99) -> g(103), g(103) -> c(99), n(110) -> n(110)
    let mut complement_table = [0u8; 256];
    complement_table[b'A' as usize] = b'T';
    complement_table[b'T' as usize] = b'A';
    complement_table[b'C' as usize] = b'G';
    complement_table[b'G' as usize] = b'C';
    complement_table[b'N' as usize] = b'N';
    complement_table[b'a' as usize] = b't';
    complement_table[b't' as usize] = b'a';
    complement_table[b'c' as usize] = b'g';
    complement_table[b'g' as usize] = b'c';
    complement_table[b'n' as usize] = b'n';
    // For other characters, map to themselves (identity)
    for i in 0..256 {
        if complement_table[i] == 0 {
            complement_table[i] = i as u8;
        }
    }
    
    // Process 16 bytes at a time using SSSE3
    const CHUNK_SIZE: usize = 16;
    let chunks = len / CHUNK_SIZE;
    let remainder = len % CHUNK_SIZE;
    
    // Reverse shuffle mask for 16 bytes (reverses byte order: [15, 14, 13, ..., 0])
    let reverse_mask = _mm_set_epi8(
        15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0
    );
    
    // Process chunks in reverse order (from end to beginning)
    for i in 0..chunks {
        let chunk_idx = chunks - 1 - i;
        let chunk_ptr = unsafe { bytes.as_ptr().add(chunk_idx * CHUNK_SIZE) };
        
        // Load 16 bytes
        let chunk = unsafe { _mm_loadu_si128(chunk_ptr as *const __m128i) };
        
        // Reverse bytes within the chunk using pshufb
        let reversed = _mm_shuffle_epi8(chunk, reverse_mask);
        
        // Store reversed bytes to temporary array
        let mut temp = [0u8; 16];
        unsafe {
            _mm_storeu_si128(temp.as_mut_ptr() as *mut __m128i, reversed);
        }
        
        // Apply complement mapping using lookup table
        for j in 0..16 {
            temp[j] = complement_table[temp[j] as usize];
        }
        
        result.extend_from_slice(&temp);
    }
    
    // Handle remainder using scalar code (process from end of sequence)
    if remainder > 0 {
        for i in (0..remainder).rev() {
            let idx = chunks * CHUNK_SIZE + i;
            result.push(complement_table[bytes[idx] as usize]);
        }
    }
    
    unsafe { String::from_utf8_unchecked(result) }
}

