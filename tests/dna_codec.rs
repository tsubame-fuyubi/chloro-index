use chloro_index::{encode_dna, decode_dna, MAX_DNA_LENGTH};

#[test]
fn test_dna_encoding_basic() {
    let seq = "ACGT";
    let encoded = encode_dna(seq).unwrap();
    // A(00) C(01) G(10) T(11) -> 00011011 (binary) = 27 (decimal)
    assert_eq!(encoded, 0b00011011);
    
    let decoded = decode_dna(encoded, 4);
    assert_eq!(decoded, "ACGT");
}

#[test]
fn test_dna_encoding_case_insensitive() {
    let seq1 = "ACGT";
    let seq2 = "acgt";
    let seq3 = "AcGt";
    
    let encoded1 = encode_dna(seq1).unwrap();
    let encoded2 = encode_dna(seq2).unwrap();
    let encoded3 = encode_dna(seq3).unwrap();
    
    assert_eq!(encoded1, encoded2);
    assert_eq!(encoded1, encoded3);
    
    assert_eq!(decode_dna(encoded1, 4), "ACGT");
}

#[test]
fn test_invalid_dna_base() {
    assert!(encode_dna("ACZ").is_err());
    assert!(encode_dna("ACGTX").is_err());
}

#[test]
fn test_empty_sequence() {
    assert!(encode_dna("").is_err());
    assert_eq!(decode_dna(0, 0), "");
}

#[test]
fn test_max_length_sequence() {
    let seq = "A".repeat(MAX_DNA_LENGTH);
    assert!(encode_dna(&seq).is_ok());
    
    let seq_too_long = "A".repeat(MAX_DNA_LENGTH + 1);
    assert!(encode_dna(&seq_too_long).is_err());
}

#[test]
fn test_single_base() {
    for (base, expected_bits) in [('A', 0b00), ('C', 0b01), ('G', 0b10), ('T', 0b11)] {
        let encoded = encode_dna(&base.to_string()).unwrap();
        assert_eq!(encoded, expected_bits as u64);
        assert_eq!(decode_dna(encoded, 1), base.to_string());
    }
}

#[test]
fn test_round_trip() {
    let sequences = vec![
        "A",
        "AC",
        "ACGT",
        "AAAA",
        "TTTT",
        "ATCGATCG",
        "ACGTACGTACGTACGT",
    ];

    for seq in sequences {
        let encoded = encode_dna(seq).unwrap();
        let decoded = decode_dna(encoded, seq.len());
        assert_eq!(decoded, seq.to_uppercase());
    }
}

