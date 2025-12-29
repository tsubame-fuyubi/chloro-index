use serde::{Deserialize, Serialize};
use std::fs::File;
use std::io::{self, BufReader};
use std::mem;

/// B-Tree node structure
#[derive(Debug, Serialize, Deserialize, Clone)]
struct BTreeNode {
    keys: Vec<u64>,
    values: Vec<String>,
    children: Vec<BTreeNode>,
    is_leaf: bool,
}

impl BTreeNode {
    fn new(is_leaf: bool) -> Self {
        BTreeNode {
            keys: Vec::new(),
            values: Vec::new(),
            children: Vec::new(),
            is_leaf,
        }
    }
}

/// B-Tree data structure
#[derive(Serialize, Deserialize)]
pub struct BTree {
    root: BTreeNode,
    /// Minimum degree (t)
    t: usize,
}

impl BTree {
    /// Creates a new B-Tree instance with the given minimum degree
    pub fn new(t: usize) -> Self {
        BTree {
            root: BTreeNode::new(true),
            t,
        }
    }

    /// Searches for a key and returns its associated value if found
    pub fn search(&self, key: u64) -> Option<String> {
        self.search_node(&self.root, key)
    }

    fn search_node(&self, node: &BTreeNode, key: u64) -> Option<String> {
        let mut i = 0;
        while i < node.keys.len() && key > node.keys[i] {
            i += 1;
        }

        // Key found, return the associated value
        if i < node.keys.len() && key == node.keys[i] {
            return Some(node.values[i].clone());
        }

        // If leaf node and key not found, it doesn't exist
        if node.is_leaf {
            return None;
        }

        // Recursively search in the appropriate child node
        self.search_node(&node.children[i], key)
    }

    /// Inserts a key-value pair into the B-Tree
    pub fn insert(&mut self, key: u64, value: String) {
        // If root is full, split it and create a new root
        if self.root.keys.len() == 2 * self.t - 1 {
            let mut new_root = BTreeNode::new(false);
            let old_root = mem::replace(&mut self.root, BTreeNode::new(true));
            new_root.children.push(old_root);
            
            Self::split_child(&mut new_root, 0, self.t);
            self.root = new_root;
        }
        
        Self::insert_non_full(&mut self.root, key, value, self.t);
    }

    fn insert_non_full(node: &mut BTreeNode, key: u64, value: String, t: usize) {
        let mut i = node.keys.len();

        if node.is_leaf {
            // Find insertion position in leaf node
            while i > 0 && key < node.keys[i - 1] {
                i -= 1;
            }
            node.keys.insert(i, key);
            node.values.insert(i, value);
        } else {
            // Find the appropriate child node
            while i > 0 && key < node.keys[i - 1] {
                i -= 1;
            }

            // If child is full, split it first
            if node.children[i].keys.len() == 2 * t - 1 {
                Self::split_child(node, i, t);
                if key > node.keys[i] {
                    i += 1;
                }
            }
            
            Self::insert_non_full(&mut node.children[i], key, value, t);
        }
    }

    fn split_child(parent: &mut BTreeNode, i: usize, t: usize) {
        let full_child = &mut parent.children[i];
        
        // Create new node to store the second half of full_child
        let mut new_child = BTreeNode::new(full_child.is_leaf);
        
        // Move keys and values from index t onwards
        new_child.keys.extend(full_child.keys.drain(t..));
        new_child.values.extend(full_child.values.drain(t..));
        
        // If not a leaf, also move child pointers
        if !full_child.is_leaf {
            new_child.children.extend(full_child.children.drain(t..));
        }

        // Promote the middle key to parent
        // Safe to unwrap: full_child is guaranteed to have at least t keys
        let mid_key = full_child.keys.pop().unwrap();
        let mid_value = full_child.values.pop().unwrap();

        parent.keys.insert(i, mid_key);
        parent.values.insert(i, mid_value);
        parent.children.insert(i + 1, new_child);
    }

    /// Prints the tree structure (for debugging)
    pub fn print_tree(&self) {
        self.print_node(&self.root, 0);
    }

    fn print_node(&self, node: &BTreeNode, level: usize) {
        let indent = "  ".repeat(level);
        let pairs: Vec<String> = node.keys.iter()
            .zip(node.values.iter())
            .map(|(k, v)| format!("{}:{}", k, v))
            .collect();
            
        println!("{}Level {}: {:?}", indent, level, pairs);
        
        if !node.is_leaf {
            for child in &node.children {
                self.print_node(child, level + 1);
            }
        }
    }
    
    /// Persistence functions

    /// Saves the B-Tree to a file
    pub fn save_to_file(&self, filename: &str) -> io::Result<()> {
        let file = File::create(filename)?;
        serde_json::to_writer_pretty(file, self)
            .map_err(|e| io::Error::new(io::ErrorKind::Other, e))?;
        Ok(())
    }

    /// Loads a B-Tree from a file
    pub fn load_from_file(filename: &str) -> io::Result<Self> {
        let file = File::open(filename)?;
        let reader = BufReader::new(file);
        let tree = serde_json::from_reader(reader)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
        Ok(tree)
    }    
}

// --- DNA Encoding/Decoding Module ---

/// Maximum DNA sequence length that can be encoded in a u64.
/// Each base uses 2 bits, so we can store up to 32 bases (64 bits / 2 bits per base).
pub const MAX_DNA_LENGTH: usize = 32;

/// Bits per base pair in the encoding scheme.
const BITS_PER_BASE: u32 = 2;

/// Encodes a DNA sequence into a u64 integer.
/// 
/// Each base is encoded using 2 bits: A=00, C=01, G=10, T=11.
/// Bases are packed from left to right (most significant bits first).
/// 
/// # Arguments
/// 
/// * `dna` - A string containing only A, C, G, T (case-insensitive)
/// 
/// # Returns
/// 
/// * `Ok(u64)` - The encoded sequence as a u64 integer
/// * `Err(String)` - An error message if the sequence is invalid or too long
/// 
/// # Examples
/// 
/// ```
/// use chloro_index::encode_dna;
/// 
/// let encoded = encode_dna("ACGT").unwrap();
/// assert_eq!(encoded, 0b00011011); // A(00) C(01) G(10) T(11)
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
        encoded <<= BITS_PER_BASE; // Shift left to make room for the next base
        encoded |= encode_base(base)? as u64;
    }
    
    Ok(encoded)
}

/// Decodes a u64 integer back into a DNA sequence.
/// 
/// The length parameter is required because leading zeros in the encoded value
/// cannot be distinguished from encoded 'A' bases (which are also 00).
/// 
/// # Arguments
/// 
/// * `value` - The encoded u64 integer
/// * `length` - The original length of the DNA sequence
/// 
/// # Returns
/// 
/// * `String` - The decoded DNA sequence
/// 
/// # Examples
/// 
/// ```
/// use chloro_index::{encode_dna, decode_dna};
/// 
/// let encoded = encode_dna("ACGT").unwrap();
/// let decoded = decode_dna(encoded, 4);
/// assert_eq!(decoded, "ACGT");
/// ```
pub fn decode_dna(mut value: u64, length: usize) -> String {
    if length == 0 {
        return String::new();
    }
    
    let mut sequence = String::with_capacity(length);
    
    // Extract bases from least significant bits (right to left)
    for _ in 0..length {
        let base_bits = (value & 0b11) as u8;
        sequence.push(decode_base(base_bits));
        value >>= BITS_PER_BASE;
    }
    
    // Reverse because we decoded from right to left
    sequence.chars().rev().collect()
}

/// Encodes a single DNA base to its 2-bit representation.
/// 
/// * A/a -> 0b00 (0)
/// * C/c -> 0b01 (1)
/// * G/g -> 0b10 (2)
/// * T/t -> 0b11 (3)
fn encode_base(base: char) -> Result<u8, String> {
    match base.to_ascii_uppercase() {
        'A' => Ok(0b00),
        'C' => Ok(0b01),
        'G' => Ok(0b10),
        'T' => Ok(0b11),
        _ => Err(format!("Invalid DNA base: '{}' (expected A, C, G, or T)", base)),
    }
}

/// Decodes a 2-bit value back to its corresponding DNA base character.
fn decode_base(bits: u8) -> char {
    match bits & 0b11 {
        0b00 => 'A',
        0b01 => 'C',
        0b10 => 'G',
        0b11 => 'T',
        _ => unreachable!("2-bit value cannot exceed 0b11"),
    }
}

