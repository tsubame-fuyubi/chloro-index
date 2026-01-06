//! Chloro-Index: High-performance DNA Indexer CLI
//!
//! This is the main entry point for the Chloro-Index v2.0-alpha CLI toolchain.
//! The legacy CLI demo is available via: cargo run --bin demo_v1_simd

use anyhow::{Context, Result};
use chloro_index::{BTree, GenomicLocation, encode_dna, minimizer::MinimizerIterator};
use clap::{Parser, Subcommand};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::PathBuf;

/// High-performance DNA Indexer
#[derive(Parser)]
#[command(name = "chloro")]
#[command(version)]
#[command(about = "High-performance DNA Indexer")]
#[command(propagate_version = true)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

/// Available subcommands
#[derive(Subcommand)]
enum Commands {
    /// Build an index from a FASTA file
    Index {
        /// Input FASTA file
        #[arg(short = 'i', long = "input", required = true, value_name = "FASTA")]
        input: PathBuf,
        
        /// Output index file (.db)
        #[arg(short = 'o', long = "output", required = true, value_name = "DB")]
        output: PathBuf,
        
        /// Use minimizer-based indexing strategy
        #[arg(short = 'm', long = "minimizer")]
        minimizer: bool,
    },
    
    /// Search for a DNA sequence in the index
    Search {
        /// Index database file
        #[arg(short = 'd', long = "db", required = true, value_name = "DB")]
        db: PathBuf,
        
        /// DNA sequence to search
        #[arg(short = 'q', long = "query", required = true)]
        query: String,
        
        /// Use minimizer-based search strategy
        #[arg(short, long)]
        minimizer: bool,
    },
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    
    match cli.command {
        Commands::Index { input, output, minimizer } => {
            index_fasta(&input, &output, minimizer)?;
        }
        Commands::Search { db, query, minimizer } => {
            search_index(&db, &query, minimizer)?;
        }
    }
    
    Ok(())
}

/// Indexes a FASTA file into a B-Tree database.
fn index_fasta(input: &PathBuf, output: &PathBuf, use_minimizer: bool) -> Result<()> {
    println!("Indexing {:?} to {:?}...", input, output);
    if use_minimizer {
        println!("Using minimizer-based indexing strategy (k=15, w=10)");
    } else {
        println!("Using standard 32-mer indexing strategy");
    }
    
    // Open and read FASTA file
    let file = File::open(input)
        .with_context(|| format!("Failed to open input file: {:?}", input))?;
    let reader = BufReader::new(file);
    
    // Create B-Tree with minimum degree 2
    let mut tree = BTree::new(2);
    
    // Parse FASTA and index sequences
    let total_keys = parse_and_index_fasta(reader, &mut tree, use_minimizer)
        .with_context(|| "Failed to parse and index FASTA file")?;
    
    // Save the index
    let output_str = output.to_str()
        .ok_or_else(|| anyhow::anyhow!("Output path contains invalid UTF-8: {:?}", output))?;
    tree.save_to_file(output_str)
        .with_context(|| format!("Failed to save index to: {:?}", output))?;
    
    // Get file size
    let file_size = std::fs::metadata(output)
        .with_context(|| format!("Failed to get file size for: {:?}", output))?
        .len();
    
    // Print statistics
    println!("✓ Indexed {} keys", total_keys);
    println!("✓ Saved index to {:?} ({} bytes)", output, file_size);
    
    Ok(())
}

/// Parses a FASTA file and indexes all sequences into the B-Tree.
fn parse_and_index_fasta<R: BufRead>(
    reader: R,
    tree: &mut BTree,
    use_minimizer: bool,
) -> Result<usize> {
    let mut total_keys = 0;
    let mut current_chr: Option<u32> = None;
    let mut current_seq = String::new();
    let mut chr_counter: u32 = 1;
    
    for line in reader.lines() {
        let line = line?;
        let line = line.trim();
        
        if line.is_empty() {
            continue;
        }
        
        if line.starts_with('>') {
            // Process previous sequence if any
            if !current_seq.is_empty() {
                if let Some(chr) = current_chr {
                    let keys = index_sequence(&current_seq, chr, tree, use_minimizer)?;
                    total_keys += keys;
                }
                current_seq.clear();
            }
            
            // Start new sequence
            current_chr = Some(chr_counter);
            chr_counter = chr_counter.saturating_add(1);
            if chr_counter == 0 {
                chr_counter = 1; // Avoid 0, keep it in 1-based range
            }
        } else {
            // Accumulate sequence lines
            current_seq.push_str(line);
        }
    }
    
    // Process last sequence
    if !current_seq.is_empty() {
        if let Some(chr) = current_chr {
            let keys = index_sequence(&current_seq, chr, tree, use_minimizer)?;
            total_keys += keys;
        }
    }
    
    Ok(total_keys)
}

/// Indexes a single DNA sequence using either standard or minimizer mode.
fn index_sequence(
    seq: &str,
    chr: u32,
    tree: &mut BTree,
    use_minimizer: bool,
) -> Result<usize> {
    // Normalize sequence: convert to uppercase and filter invalid characters
    let normalized: String = seq
        .chars()
        .filter_map(|c| {
            let upper = c.to_ascii_uppercase();
            match upper {
                'A' | 'C' | 'G' | 'T' => Some(upper),
                'N' | _ => None, // Skip N and other invalid characters
            }
        })
        .collect();
    
    if normalized.is_empty() {
        return Ok(0);
    }
    
    let mut count = 0;
    
    if use_minimizer {
        // Minimizer mode: k=15, w=10
        const K: u8 = 15;
        const W: usize = 10;
        
        let iter = MinimizerIterator::new(&normalized, K, W);
        
        for (pos, canonical_hash) in iter {
            let loc = GenomicLocation {
                chromosome: chr,
                position: (pos + 1) as u32, // 1-based position
            };
            tree.insert(canonical_hash, vec![loc]);
            count += 1;
        }
    } else {
        // Standard mode: 32-mer sliding window
        const K_MER_SIZE: usize = 32;
        
        if normalized.len() < K_MER_SIZE {
            return Ok(0);
        }
        
        for (i, window_bytes) in normalized.as_bytes().windows(K_MER_SIZE).enumerate() {
            if let Ok(fragment) = std::str::from_utf8(window_bytes) {
                if let Ok(key) = encode_dna(fragment) {
                    let loc = GenomicLocation {
                        chromosome: chr,
                        position: (i + 1) as u32, // 1-based position
                    };
                    tree.insert(key, vec![loc]);
                    count += 1;
                }
                // Skip invalid k-mers (e.g., containing N) - encode_dna will fail
            }
        }
    }
    
    Ok(count)
}

/// Searches for a DNA sequence in the index database.
fn search_index(db: &PathBuf, query: &str, use_minimizer: bool) -> Result<()> {
    println!("Searching for query '{}' in database {:?}...", query, db);
    
    // Load the index
    let db_str = db.to_str()
        .ok_or_else(|| anyhow::anyhow!("Database path contains invalid UTF-8: {:?}", db))?;
    let tree = BTree::load_from_file(db_str)
        .with_context(|| format!("Failed to load index from: {:?}", db))?;
    
    // Normalize query: convert to uppercase and filter invalid characters
    let normalized: String = query
        .chars()
        .filter_map(|c| {
            let upper = c.to_ascii_uppercase();
            match upper {
                'A' | 'C' | 'G' | 'T' => Some(upper),
                'N' | _ => None, // Skip N and other invalid characters
            }
        })
        .collect();
    
    if normalized.is_empty() {
        println!("No matches found. (Query contains no valid DNA bases)");
        return Ok(());
    }
    
    let mut found_any = false;
    
    if use_minimizer {
        // Minimizer mode: k=15, w=10 (must match Index parameters)
        const K: u8 = 15;
        const W: usize = 10;
        
        println!("Using minimizer-based search (k={}, w={})", K, W);
        
        let iter = MinimizerIterator::new(&normalized, K, W);
        let minimizers: Vec<(usize, u64)> = iter.collect();
        
        if minimizers.is_empty() {
            println!("No matches found. (Query too short for minimizer extraction)");
            return Ok(());
        }
        
        for (pos, canonical_hash) in minimizers {
            if let Some(locations) = tree.search(canonical_hash) {
                found_any = true;
                println!("Hit found for minimizer [{}] at query position {}:", canonical_hash, pos + 1);
                for loc in locations {
                    println!("  Chromosome {}: Position {}", loc.chromosome, loc.position);
                }
            }
        }
    } else {
        // Standard mode: use first 32 bases
        const K_MER_SIZE: usize = 32;
        
        if normalized.len() < K_MER_SIZE {
            println!("No matches found. (Query must be at least {} bases for standard search)", K_MER_SIZE);
            return Ok(());
        }
        
        println!("Using standard 32-mer search");
        
        // Take the first 32 bases
        let query_32mer = &normalized[..K_MER_SIZE];
        
        match encode_dna(query_32mer) {
            Ok(encoded_key) => {
                if let Some(locations) = tree.search(encoded_key) {
                    found_any = true;
                    println!("Match found for query '{}':", query_32mer);
                    for loc in locations {
                        println!("  Chromosome {}: Position {}", loc.chromosome, loc.position);
                    }
                }
            }
            Err(e) => {
                return Err(anyhow::anyhow!("Failed to encode query: {}", e));
            }
        }
    }
    
    if !found_any {
        println!("No matches found.");
    }
    
    Ok(())
}
