#!/usr/bin/env python3
"""
Generate a smoke test FASTA file for chloro-index CLI testing.

This script creates a small FASTA file with known sequences that can be
used to verify the indexing and search functionality.
"""

import random
import os
from pathlib import Path

# Golden key for testing (32-mer)
GOLDEN_KEY = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"  # 32 A's
SECOND_KEY = "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"  # 32 C's

# Configuration
OUTPUT_DIR = Path("data")
OUTPUT_FILE = OUTPUT_DIR / "smoke_test.fasta"
SEQUENCE_LENGTH = 1000  # Base pairs per sequence (excluding the key)


def generate_random_dna(length: int) -> str:
    """Generate a random DNA sequence of given length."""
    bases = ['A', 'C', 'G', 'T']
    return ''.join(random.choice(bases) for _ in range(length))


def generate_fasta_content() -> str:
    """Generate FASTA file content with golden keys embedded."""
    lines = []
    
    # Sequence 1: chr1 with GOLDEN_KEY in the middle
    seq1_prefix = generate_random_dna(SEQUENCE_LENGTH // 2)
    seq1_suffix = generate_random_dna(SEQUENCE_LENGTH // 2)
    seq1 = seq1_prefix + GOLDEN_KEY + seq1_suffix
    
    lines.append(">chr1")
    # Write sequence in chunks of 80 characters (FASTA standard)
    for i in range(0, len(seq1), 80):
        lines.append(seq1[i:i+80])
    
    # Sequence 2: chr2 with SECOND_KEY in the middle
    seq2_prefix = generate_random_dna(SEQUENCE_LENGTH // 2)
    seq2_suffix = generate_random_dna(SEQUENCE_LENGTH // 2)
    seq2 = seq2_prefix + SECOND_KEY + seq2_suffix
    
    lines.append(">chr2")
    for i in range(0, len(seq2), 80):
        lines.append(seq2[i:i+80])
    
    # Sequence 3: chr3 with random sequence (no special key)
    seq3 = generate_random_dna(SEQUENCE_LENGTH)
    lines.append(">chr3")
    for i in range(0, len(seq3), 80):
        lines.append(seq3[i:i+80])
    
    return '\n'.join(lines)


def main():
    """Generate the smoke test FASTA file and print test commands."""
    # Create output directory if it doesn't exist
    OUTPUT_DIR.mkdir(exist_ok=True)
    
    # Generate FASTA content
    print(f"Generating smoke test FASTA file...")
    fasta_content = generate_fasta_content()
    
    # Write to file
    with open(OUTPUT_FILE, 'w') as f:
        f.write(fasta_content)
    
    file_size = OUTPUT_FILE.stat().st_size
    print(f"[OK] Generated {OUTPUT_FILE} ({file_size} bytes)")
    print()
    
    # Print test instructions
    print("=" * 70)
    print("SMOKE TEST COMMANDS")
    print("=" * 70)
    print()
    print("1. Index the FASTA file (Standard Mode):")
    print(f'   cargo run --bin chloro-index -- index -i {OUTPUT_FILE} -o data/smoke_test.db')
    print()
    print("2. Search for the golden key (Standard Mode):")
    print(f'   cargo run --bin chloro-index -- search -d data/smoke_test.db -q "{GOLDEN_KEY}"')
    print()
    print("3. Index the FASTA file (Minimizer Mode):")
    print(f'   cargo run --bin chloro-index -- index -i {OUTPUT_FILE} -o data/smoke_test_minimizer.db -m')
    print()
    print("4. Search using minimizer mode:")
    print(f'   cargo run --bin chloro-index -- search -d data/smoke_test_minimizer.db -q "{GOLDEN_KEY}" -m')
    print()
    print("=" * 70)
    print()
    print(f"Golden Key (32-mer): {GOLDEN_KEY}")
    print(f"Second Key (32-mer): {SECOND_KEY}")
    print()
    print("Expected results:")
    print(f"  - Standard search for '{GOLDEN_KEY}' should find chr1")
    print(f"  - Standard search for '{SECOND_KEY}' should find chr2")
    print("  - Minimizer search should find multiple hits")


if __name__ == "__main__":
    # Set random seed for reproducibility
    random.seed(42)
    main()

