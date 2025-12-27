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
