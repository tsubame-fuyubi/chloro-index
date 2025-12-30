//! Command-line interface for ChloroDB genomic sequence index.
//!
//! This program provides an interactive shell for importing FASTA files
//! and querying genomic k-mer locations.

use std::fs::File;
use std::io::{self, BufRead, BufReader, Write};
use std::time::Instant;

use chloro_index::{BTree, encode_dna};

const DEFAULT_DB_FILE: &str = "chloro.db";
const MAX_DISPLAY_RESULTS: usize = 5;

fn main() {
    let filename = DEFAULT_DB_FILE;
    
    // Initialize or load database
    let mut tree = match BTree::load_from_file(filename) {
        Ok(t) => {
            println!("Database loaded from {}.", filename);
            t
        }
        Err(_) => {
            println!("Initialized new B-Tree database.");
            BTree::new(4)
        }
    };

    println!("=== ChloroDB: Genomic K-mer Index ===");
    println!("Commands: import <file>, find <seq>, exit");

    // Main command loop
    loop {
        print!("> ");
        io::stdout().flush().unwrap();
        
        let mut input = String::new();
        if io::stdin().read_line(&mut input).is_err() {
            break;
        }
        
        let parts: Vec<&str> = input.trim().split_whitespace().collect();
        if parts.is_empty() {
            continue;
        }

        match parts[0] {
            "import" => handle_import(&mut tree, &parts),
            "find" => handle_find(&tree, &parts),
            "exit" => {
                if let Err(e) = tree.save_to_file(filename) {
                    eprintln!("Error saving database: {}", e);
                } else {
                    println!("Database saved to {}.", filename);
                }
                break;
            }
            _ => println!("Unknown command. Available: import, find, exit"),
        }
    }
}

/// Handles the import command for FASTA files.
fn handle_import(tree: &mut BTree, parts: &[&str]) {
    if parts.len() < 2 {
        eprintln!("Usage: import <fasta_file_path>");
        return;
    }
    
    let path = parts[1];
    let start = Instant::now();
    println!("Reading file: {}", path);

    let file = match File::open(path) {
        Ok(f) => f,
        Err(e) => {
            eprintln!("Error opening file: {}", e);
            return;
        }
    };

    match import_fasta(tree, BufReader::new(file)) {
        Ok(total_count) => {
            println!(
                "Import complete in {:.2?}. Total k-mers indexed: {}",
                start.elapsed(),
                total_count
            );
        }
        Err(e) => {
            eprintln!("Error importing FASTA file: {}", e);
        }
    }
}

/// Parses a FASTA file and indexes all sequences.
///
/// This function implements a simple FASTA parser that extracts chromosome
/// IDs from header lines (lines starting with '>') and indexes the
/// corresponding DNA sequences.
///
/// # Arguments
/// * `tree` - B-Tree index to populate
/// * `reader` - Buffered reader for the FASTA file
///
/// # Returns
/// Total number of k-mers indexed, or an error if parsing fails
fn import_fasta<R: BufRead>(tree: &mut BTree, reader: R) -> io::Result<usize> {
    let mut current_chr: Option<u8> = None;
    let mut sequence_buffer = String::new();
    let mut total_count = 0;

    for line_result in reader.lines() {
        let line = line_result?;
        let line = line.trim();

        if line.starts_with('>') {
            // Process previous sequence before starting new one
            if let Some(chr) = current_chr {
                if !sequence_buffer.is_empty() {
                    let count = tree.bulk_insert_sequence(chr, &sequence_buffer);
                    println!("  Indexed Chr{}: {} k-mers", chr, count);
                    total_count += count;
                    sequence_buffer.clear();
                }
            }
            
            // Extract chromosome ID from header (e.g., ">chr1" -> 1)
            current_chr = parse_chromosome_id(line);
            if current_chr.is_none() {
                println!(
                    "  Warning: skipping header '{}' (ID must be numeric)",
                    line
                );
            }
        } else if !line.is_empty() {
            // Accumulate sequence data
            sequence_buffer.push_str(line);
        }
    }

    // Process final sequence
    if let Some(chr) = current_chr {
        if !sequence_buffer.is_empty() {
            let count = tree.bulk_insert_sequence(chr, &sequence_buffer);
            println!("  Indexed Chr{}: {} k-mers", chr, count);
            total_count += count;
        }
    }

    Ok(total_count)
}

/// Extracts chromosome ID from a FASTA header line.
///
/// This is a simple parser that extracts numeric characters from the
/// header. For example, ">chr1" or ">chromosome1" would both yield 1.
///
/// # Arguments
/// * `header` - FASTA header line (must start with '>')
///
/// # Returns
/// Chromosome ID if numeric characters are found, None otherwise
fn parse_chromosome_id(header: &str) -> Option<u8> {
    let id_str: String = header.chars().filter(|c| c.is_ascii_digit()).collect();
    id_str.parse::<u8>().ok()
}

/// Handles the find command for sequence lookup.
fn handle_find(tree: &BTree, parts: &[&str]) {
    if parts.len() < 2 {
        eprintln!("Usage: find <dna_sequence>");
        return;
    }

    let sequence = parts[1];
    
    // Early return on encoding error
    let key = match encode_dna(sequence) {
        Ok(k) => k,
        Err(e) => {
            eprintln!("Invalid DNA sequence: {}", e);
            return;
        }
    };
    
    // Early return if not found
    let locations = match tree.search(key) {
        Some(locs) => locs,
        None => {
            println!("Sequence not found in index.");
            return;
        }
    };
    
    // Display results
    println!("Found {} occurrence(s):", locations.len());
    for (i, loc) in locations.iter().enumerate() {
        if i >= MAX_DISPLAY_RESULTS {
            println!(
                "  ... ({} more)",
                locations.len() - MAX_DISPLAY_RESULTS
            );
            break;
        }
        println!("  - chr{}:{}", loc.chromosome, loc.position);
    }
}
