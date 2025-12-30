use serde::{Deserialize, Serialize};
use std::fs::File;
use std::io::{self, BufReader};
use std::mem;

/// Genomic location with chromosome ID and position
#[derive(Debug, Serialize, Deserialize, Clone, PartialEq)]
pub struct GenomicLocation {
    pub chromosome: u8,
    pub position: u32,
}

#[derive(Debug, Serialize, Deserialize, Clone)]
struct BTreeNode {
    keys: Vec<u64>,
    values: Vec<Vec<GenomicLocation>>,
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

#[derive(Serialize, Deserialize)]
pub struct BTree {
    root: BTreeNode,
    t: usize,
}

impl BTree {
    pub fn new(t: usize) -> Self {
        BTree {
            root: BTreeNode::new(true),
            t,
        }
    }

    pub fn search(&self, key: u64) -> Option<Vec<GenomicLocation>> {
        self.search_node(&self.root, key)
    }

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

    fn insert_non_full(node: &mut BTreeNode, key: u64, locations: Vec<GenomicLocation>, t: usize) {
        let mut i = 0;
        while i < node.keys.len() && key > node.keys[i] {
            i += 1;
        }

        if node.is_leaf {
            if i < node.keys.len() && key == node.keys[i] {
                node.values[i].extend(locations);
            } else {
                node.keys.insert(i, key);
                node.values.insert(i, locations);
            }
        } else {
            if i < node.keys.len() && key == node.keys[i] {
                node.values[i].extend(locations);
                return;
            }

            if node.children[i].keys.len() == 2 * t - 1 {
                Self::split_child(node, i, t);
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

    fn split_child(parent: &mut BTreeNode, i: usize, t: usize) {
        let full_child = &mut parent.children[i];
        let mut new_child = BTreeNode::new(full_child.is_leaf);
        
        new_child.keys.extend(full_child.keys.drain(t..));
        new_child.values.extend(full_child.values.drain(t..));
        
        if !full_child.is_leaf {
            new_child.children.extend(full_child.children.drain(t..));
        }

        let mid_key = full_child.keys.pop().unwrap();
        let mid_value = full_child.values.pop().unwrap();

        parent.keys.insert(i, mid_key);
        parent.values.insert(i, mid_value);
        parent.children.insert(i + 1, new_child);
    }

    pub fn print_tree(&self) {
        self.print_node(&self.root, 0);
    }

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
    pub fn save_to_file(&self, filename: &str) -> io::Result<()> {
        let file = File::create(filename)?;
        serde_json::to_writer_pretty(file, self)
            .map_err(|e| io::Error::new(io::ErrorKind::Other, e))?;
        Ok(())
    }

    pub fn load_from_file(filename: &str) -> io::Result<Self> {
        let file = File::open(filename)?;
        let reader = BufReader::new(file);
        let tree = serde_json::from_reader(reader)
            .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
        Ok(tree)
    }    
}

pub const MAX_DNA_LENGTH: usize = 32;
const BITS_PER_BASE: u32 = 2;

/// Encodes a DNA sequence into a u64 (2 bits per base: A=00, C=01, G=10, T=11)
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

/// Decodes a u64 back into a DNA sequence (length required to distinguish leading zeros from 'A' bases)
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

fn encode_base(base: char) -> Result<u8, String> {
    match base.to_ascii_uppercase() {
        'A' => Ok(0b00),
        'C' => Ok(0b01),
        'G' => Ok(0b10),
        'T' => Ok(0b11),
        _ => Err(format!("Invalid DNA base: '{}' (expected A, C, G, or T)", base)),
    }
}

fn decode_base(bits: u8) -> char {
    match bits & 0b11 {
        0b00 => 'A',
        0b01 => 'C',
        0b10 => 'G',
        0b11 => 'T',
        _ => unreachable!("2-bit value cannot exceed 0b11"),
    }
}

