# chloro-index
A genomic sequence indexing tool based on B-Tree data structure for fast location of DNA k-mers in genomes.

### Features
- Import FASTA files and index all 32-mers
- Query DNA sequence locations in genomes
- Compress DNA sequences to u64 using 2-bit encoding (A=00, C=01, G=10, T=11)
- Support persistent storage to disk

### Usage
```bash
cargo run
> import data/test.fasta
> find ACGTACGTACGTACGTACGTACGTACGTAC
> exit
```
