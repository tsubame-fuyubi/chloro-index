use std::io::{self, Write};
use chloro_index::BTree; 

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
            "insert" => {
                if parts.len() < 2 {
                    eprintln!("Error: Missing arguments. Usage: insert <key> [value]");
                    continue;
                }

                if let Ok(key) = parts[1].parse::<u64>() {
                    // Get value, generate default if not provided
                    let value = if parts.len() > 2 {
                        parts[2].to_string()
                    } else {
                        format!("val_{}", key)
                    };

                    println!("Inserted: {} -> {}", key, value);
                    tree.insert(key, value);
                } else {
                    eprintln!("Error: Key must be a valid number.");
                }
            }
            "search" => {
                if parts.len() < 2 {
                    eprintln!("Error: Missing key. Usage: search <key>");
                    continue;
                }

                if let Ok(key) = parts[1].parse::<u64>() {
                    match tree.search(key) {
                        Some(val) => {
                            println!("Found: {} => {}", key, val);
                        },
                        None => {
                            println!("Not found: Key {}", key);
                        }
                    }
                } else {
                    eprintln!("Error: Key must be a valid number.");
                }
            }
            "show" => {
                tree.print_tree();
            }
            "help" => {
                println!("Commands:");
                println!("  insert <key> [value]  - Insert a key-value pair");
                println!("  search <key>          - Search for a key");
                println!("  show                  - Display the tree structure");
                println!("  exit                  - Save and exit");
            }
            "exit" => {
                print!("Saving database...");
                match tree.save_to_file(filename) {
                    Ok(_) => println!(" Done."),
                    Err(e) => eprintln!(" Error saving file: {}", e),
                }
                break;
            }
            _ => {
                eprintln!("Unknown command: '{}'. Type 'help' for list.", parts[0]);
            }
        }
    }
}