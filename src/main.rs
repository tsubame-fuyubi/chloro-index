use std::io::{self, Write};
use chloro_index::{BTree, GenomicLocation};
use chloro_index::encode_dna; 

fn main() {
    let filename = "chloro.db";

    // Try to load database, create new one if it doesn't exist
    let mut tree = match BTree::load_from_file(filename) {
        Ok(t) => {
            println!("Database loaded successfully.");
            t
        },
        Err(_) => {
            println!("No existing database found. Initializing new B-Tree.");
            BTree::new(4)
        }
    };

    println!("=== ChloroDB Index CLI ===");
    println!("Type 'help' for instructions.");

    loop {
        print!("> ");
        if io::stdout().flush().is_err() {
            break;
        }

        let mut input = String::new();
        if io::stdin().read_line(&mut input).is_err() {
            break;
        }
        
        let input = input.trim();
        let parts: Vec<&str> = input.split_whitespace().collect();

        if parts.is_empty() {
            continue;
        }

        match parts[0] {
            "dna_add" => {
                if parts.len() < 4 {
                    eprintln!("Error: Missing arguments. Usage: dna_add <ACGT> <chr> <pos>");
                    continue;
                }

                let seq = parts[1];
                match encode_dna(seq) {
                    Ok(key) => {
                        if let (Ok(chr), Ok(pos)) = (parts[2].parse::<u8>(), parts[3].parse::<u32>()) {
                            let location = GenomicLocation { chromosome: chr, position: pos };
                            println!("Encoded sequence: '{}' -> Key: {}", seq, key);
                            tree.insert(key, vec![location.clone()]);
                            println!("Inserted: {} -> chr{}:{}", key, location.chromosome, location.position);
                        } else {
                            eprintln!("Error: Invalid chromosome or position. Usage: dna_add <ACGT> <chr> <pos>");
                        }
                    }
                    Err(e) => eprintln!("Error: Invalid DNA sequence: {}", e),
                }
            }

            "dna_find" => {
                if parts.len() < 2 {
                    eprintln!("Error: Missing sequence. Usage: dna_find <ACGT>");
                    continue;
                }

                let seq = parts[1];
                match encode_dna(seq) {
                    Ok(key) => {
                        println!("Searching for sequence: '{}' (Key: {})", seq, key);
                        match tree.search(key) {
                            Some(locations) => {
                                let loc_strs: Vec<String> = locations.iter()
                                    .map(|loc| format!("chr{}:{}", loc.chromosome, loc.position))
                                    .collect();
                                println!("Found: {} => [{}]", key, loc_strs.join(", "));
                            }
                            None => println!("Not found: Key {}", key),
                        }
                    }
                    Err(e) => eprintln!("Error: Invalid DNA sequence: {}", e),
                }
            }

            "insert" => {
                if parts.len() < 4 {
                    eprintln!("Error: Missing arguments. Usage: insert <key> <chr> <pos>");
                    continue;
                }

                if let (Ok(key), Ok(chr), Ok(pos)) = (parts[1].parse::<u64>(), parts[2].parse::<u8>(), parts[3].parse::<u32>()) {
                    let location = GenomicLocation { chromosome: chr, position: pos };
                    tree.insert(key, vec![location.clone()]);
                    println!("Inserted: {} -> chr{}:{}", key, location.chromosome, location.position);
                } else {
                    eprintln!("Error: Invalid arguments. Usage: insert <key> <chr> <pos>");
                }
            }

            "search" => {
                if parts.len() < 2 {
                    eprintln!("Error: Missing key. Usage: search <key>");
                    continue;
                }

                if let Ok(key) = parts[1].parse::<u64>() {
                    match tree.search(key) {
                        Some(locations) => {
                            let loc_strs: Vec<String> = locations.iter()
                                .map(|loc| format!("chr{}:{}", loc.chromosome, loc.position))
                                .collect();
                            println!("Found: {} => [{}]", key, loc_strs.join(", "));
                        }
                        None => println!("Not found: Key {}", key),
                    }
                } else {
                    eprintln!("Error: Key must be a valid number.");
                }
            }

            "show" => {
                tree.print_tree();
            }

            "help" => {
                println!("Available commands:");
                println!("  dna_add <seq> <chr> <pos>  - Insert DNA sequence with location");
                println!("  dna_find <seq>             - Search for DNA sequence");
                println!("  insert <key> <chr> <pos>   - Insert key with location");
                println!("  search <key>               - Search by key");
                println!("  show                       - Display B-Tree structure");
                println!("  exit                       - Save database and exit");
            }

            "exit" => {
                print!("Saving database...");
                match tree.save_to_file(filename) {
                    Ok(_) => println!(" OK."),
                    Err(e) => eprintln!(" Failed: {}", e),
                }
                break;
            }

            _ => {
                eprintln!("Unknown command: '{}'. Type 'help' for a list of commands.", parts[0]);
            }
        }
    }
}