use ahash::HashMap;
use dame::sort::{fill_hap, get_pieces_info, rc, read_primers, read_tags, run, Hap, SortArgs};
use std::io::Write;
use tempfile::tempdir;

// ── rc ────────────────────────────────────────────────────────────────────────

#[test]
fn test_rc_palindrome() {
    // ACGT reversed = TGCA, complement = ACGT
    assert_eq!(rc("ACGT"), "ACGT");
}

#[test]
fn test_rc_all_a() {
    assert_eq!(rc("AAAA"), "TTTT");
}

#[test]
fn test_rc_mixed() {
    // ATCG reversed = GCTA, complement each → CGAT
    assert_eq!(rc("ATCG"), "CGAT");
}

#[test]
fn test_rc_ambiguous_n() {
    assert_eq!(rc("N"), "N");
}

#[test]
fn test_rc_bytes_palindrome() {
    use dame::sort::rc_bytes;
    assert_eq!(rc_bytes(b"ACGT"), b"ACGT");
}

#[test]
fn test_rc_bytes_all_a() {
    use dame::sort::rc_bytes;
    assert_eq!(rc_bytes(b"AAAA"), b"TTTT");
}

#[test]
fn test_rc_bytes_iupac_full_table() {
    use dame::sort::rc_bytes;
    // Every IUPAC code reversed and complemented
    // Input:  ACGTMRWSYKVHDB
    // Reverse: BDHVKYSWRMTGCA
    // Complement: VHDBMRSWYKACGT
    assert_eq!(rc_bytes(b"ACGTMRWSYKVHDB"), b"VHDBMRSWYKACGT");
    // N stays N
    assert_eq!(rc_bytes(b"NNN"), b"NNN");
}

// ── iupac_matches ─────────────────────────────────────────────────────────────

#[test]
fn test_iupac_matches_exact() {
    use dame::sort::iupac_matches;
    assert!(iupac_matches(b'A', b'A'));
    assert!(iupac_matches(b'C', b'C'));
    assert!(iupac_matches(b'G', b'G'));
    assert!(iupac_matches(b'T', b'T'));
    assert!(!iupac_matches(b'A', b'C'));
    assert!(!iupac_matches(b'G', b'T'));
}

#[test]
fn test_iupac_matches_ambiguous() {
    use dame::sort::iupac_matches;
    // R = A or G
    assert!(iupac_matches(b'R', b'A'));
    assert!(iupac_matches(b'R', b'G'));
    assert!(!iupac_matches(b'R', b'C'));
    assert!(!iupac_matches(b'R', b'T'));
    // Y = C or T
    assert!(iupac_matches(b'Y', b'C'));
    assert!(iupac_matches(b'Y', b'T'));
    assert!(!iupac_matches(b'Y', b'A'));
    // N = anything
    assert!(iupac_matches(b'N', b'A'));
    assert!(iupac_matches(b'N', b'C'));
    assert!(iupac_matches(b'N', b'G'));
    assert!(iupac_matches(b'N', b'T'));
    // S = C or G
    assert!(iupac_matches(b'S', b'C'));
    assert!(iupac_matches(b'S', b'G'));
    assert!(!iupac_matches(b'S', b'A'));
    // W = A or T
    assert!(iupac_matches(b'W', b'A'));
    assert!(iupac_matches(b'W', b'T'));
    assert!(!iupac_matches(b'W', b'G'));
    // K = G or T
    assert!(iupac_matches(b'K', b'G'));
    assert!(iupac_matches(b'K', b'T'));
    assert!(!iupac_matches(b'K', b'A'));
    // M = A or C
    assert!(iupac_matches(b'M', b'A'));
    assert!(iupac_matches(b'M', b'C'));
    assert!(!iupac_matches(b'M', b'T'));
    // B = C, G, T
    assert!(iupac_matches(b'B', b'C'));
    assert!(iupac_matches(b'B', b'G'));
    assert!(iupac_matches(b'B', b'T'));
    assert!(!iupac_matches(b'B', b'A'));
    // D = A, G, T
    assert!(iupac_matches(b'D', b'A'));
    assert!(iupac_matches(b'D', b'G'));
    assert!(iupac_matches(b'D', b'T'));
    assert!(!iupac_matches(b'D', b'C'));
    // H = A, C, T
    assert!(iupac_matches(b'H', b'A'));
    assert!(iupac_matches(b'H', b'C'));
    assert!(iupac_matches(b'H', b'T'));
    assert!(!iupac_matches(b'H', b'G'));
    // V = A, C, G
    assert!(iupac_matches(b'V', b'A'));
    assert!(iupac_matches(b'V', b'C'));
    assert!(iupac_matches(b'V', b'G'));
    assert!(!iupac_matches(b'V', b'T'));
}

#[test]
fn test_iupac_matches_non_acgt_read_byte() {
    use dame::sort::iupac_matches;
    // Lowercase or non-ACGT read bytes are not valid FASTQ bases; wildcard returns false
    assert!(!iupac_matches(b'A', b'a'));
    assert!(!iupac_matches(b'N', b'a'));
    assert!(!iupac_matches(b'N', b'n'));
}

// ── find_primer ───────────────────────────────────────────────────────────────

#[test]
fn test_find_primer_exact() {
    use dame::sort::find_primer;
    // primer ACGT in XXXXACGTXXXX — should find at position 4
    let seq = b"XXXXACGTXXXX";
    let primer = b"ACGT";
    let result = find_primer(primer, seq);
    assert_eq!(result, Some((4, 8)));
}

#[test]
fn test_find_primer_iupac() {
    use dame::sort::find_primer;
    // primer with R (= A or G): GCRTGC matches GCATGC
    let seq = b"TTTTGCATGCTTTT";
    let primer = b"GCRTGC";
    let result = find_primer(primer, seq);
    assert_eq!(result, Some((4, 10)));
}

#[test]
fn test_find_primer_iupac_second_option() {
    use dame::sort::find_primer;
    // primer GCRTGC also matches GCGTGC (R = G)
    let seq = b"TTTTGCGTGCTTTT";
    let primer = b"GCRTGC";
    let result = find_primer(primer, seq);
    assert_eq!(result, Some((4, 10)));
}

#[test]
fn test_find_primer_not_found() {
    use dame::sort::find_primer;
    let seq = b"AAAAAAAAAA";
    let primer = b"GCATGC";
    assert_eq!(find_primer(primer, seq), None);
}

#[test]
fn test_find_primer_leftmost() {
    use dame::sort::find_primer;
    // primer appears twice — must return leftmost
    let seq = b"ACGTXXXXACGT";
    let primer = b"ACGT";
    let result = find_primer(primer, seq);
    assert_eq!(result, Some((0, 4)));
}

#[test]
fn test_find_primer_primer_longer_than_seq() {
    use dame::sort::find_primer;
    let seq = b"AC";
    let primer = b"ACGT";
    assert_eq!(find_primer(primer, seq), None);
}

// ── read_tags ─────────────────────────────────────────────────────────────────

#[test]
fn test_read_tags() {
    let dir = tempdir().unwrap();
    let path = dir.path().join("tags.txt");
    let mut f = std::fs::File::create(&path).unwrap();
    writeln!(f, "ACGT\tTag1\nTTTT\tTag2").unwrap();

    let tags = read_tags(path.to_str().unwrap()).unwrap();
    assert!(tags.contains_key("Tag1"));
    assert_eq!(tags["Tag1"][0], "ACGT");
    // rc("ACGT") == "ACGT" (palindrome)
    assert_eq!(tags["Tag1"][1], rc("ACGT"));
    assert!(tags.contains_key("Tag2"));
    assert_eq!(tags["Tag2"][0], "TTTT");
    assert_eq!(tags["Tag2"][1], "AAAA");
}

// ── read_primers ──────────────────────────────────────────────────────────────

#[test]
fn test_read_primers() {
    let dir = tempdir().unwrap();
    let path = dir.path().join("primers.txt");
    let mut f = std::fs::File::create(&path).unwrap();
    writeln!(f, "CO1\tACGT\tTTTT").unwrap();

    let primers = read_primers(path.to_str().unwrap()).unwrap();
    assert!(primers.contains_key("CO1"));
    let co1 = &primers["CO1"];
    // start_primers[0] = F = b"ACGT", start_primers[1] = R = b"TTTT"
    assert_eq!(co1.start_primers.len(), 2);
    assert_eq!(co1.start_primers[0], b"ACGT");
    assert_eq!(co1.start_primers[1], b"TTTT");
    // end_primers[0] = RC(F) = RC("ACGT") = "ACGT" (palindrome)
    // end_primers[1] = RC(R) = RC("TTTT") = "AAAA"
    assert_eq!(co1.end_primers.len(), 2);
    assert_eq!(co1.end_primers[0], b"ACGT"); // rc("ACGT") == "ACGT"
    assert_eq!(co1.end_primers[1], b"AAAA"); // rc("TTTT") == "AAAA"
}

// ── fill_hap ──────────────────────────────────────────────────────────────────

#[test]
fn test_fill_hap_new_entry() {
    let mut hap: Hap = indexmap::IndexMap::new();
    fill_hap(&mut hap, "Tag1", "Tag2", "CO1", "ACGTACGT");
    assert!(hap.contains_key("Tag1_Tag2"));
    let entry = &hap["Tag1_Tag2"];
    assert_eq!(entry.tag1, "Tag1");
    assert_eq!(entry.tag2, "Tag2");
    assert!(entry.seqs.contains_key("ACGTACGT"));
    assert_eq!(entry.seqs["ACGTACGT"].count, 1);
    assert_eq!(entry.seqs["ACGTACGT"].primer_name, "CO1");
}

#[test]
fn test_fill_hap_increment_count() {
    let mut hap: Hap = indexmap::IndexMap::new();
    fill_hap(&mut hap, "Tag1", "Tag2", "CO1", "ACGTACGT");
    fill_hap(&mut hap, "Tag1", "Tag2", "CO1", "ACGTACGT");
    assert_eq!(hap["Tag1_Tag2"].seqs["ACGTACGT"].count, 2);
}

#[test]
fn test_fill_hap_multiple_seqs() {
    let mut hap: Hap = indexmap::IndexMap::new();
    fill_hap(&mut hap, "Tag1", "Tag2", "CO1", "AAAA");
    fill_hap(&mut hap, "Tag1", "Tag2", "CO1", "CCCC");
    assert_eq!(hap["Tag1_Tag2"].seqs.len(), 2);
}

// ── get_pieces_info ───────────────────────────────────────────────────────────

/// Build a minimal tags + primers map for testing.
/// Tags: AAAA=Tag1, CCCC=Tag2, GGGG=Tag3, TTTT=Tag4
/// Primer CO1: F=ACGT, R=TGCA
fn make_test_tags() -> HashMap<String, Vec<String>> {
    let dir = tempdir().unwrap();
    let path = dir.path().join("tags.txt");
    let mut f = std::fs::File::create(&path).unwrap();
    writeln!(f, "AAAA\tTag1\nCCCC\tTag2\nGGGG\tTag3\nTTTT\tTag4").unwrap();
    // dir must stay alive until after read_tags; we shadow it
    let tags = read_tags(path.to_str().unwrap()).unwrap();
    tags
}

fn make_test_primers() -> indexmap::IndexMap<String, dame::sort::PrimerEntry> {
    let dir = tempdir().unwrap();
    let path = dir.path().join("primers.txt");
    let mut f = std::fs::File::create(&path).unwrap();
    writeln!(f, "CO1\tACGT\tTGCA").unwrap();
    read_primers(path.to_str().unwrap()).unwrap()
}

#[test]
fn test_get_pieces_info_forward_read() {
    // Read: AAAA ACGT ATATAT TGCA GGGG
    // tag1=AAAA (Tag1 forward), primer_F=ACGT, barcode=ATATAT, RC(R)=TGCA, tag2=GGGG
    // GGGG = RC(CCCC) → Tag2's RC seq → tag2=Tag2 (Tags[Tag2][1]=="GGGG")
    let tags = make_test_tags();
    let primers = make_test_primers();
    let line = "AAAAACGTATATATTGCAGGGG";
    let info = get_pieces_info(line, &primers, &tags, false);
    assert!(info.is_some(), "Expected Some(PieceInfo) for a valid forward read");
    let info = info.unwrap();
    assert_eq!(info.tag1, "Tag1");
    assert_eq!(info.tag2, "Tag2");
    assert_eq!(info.between, "ATATAT");
    assert_eq!(info.primer_name, "CO1");
}

#[test]
fn test_get_pieces_info_error_read() {
    let tags = make_test_tags();
    let primers = make_test_primers();
    let line = "NNNNNNNNNNNNNNNNNNNN";
    let info = get_pieces_info(line, &primers, &tags, false);
    assert!(info.is_none(), "Expected None for an error/ambiguous read");
}

// ── run (integration) ─────────────────────────────────────────────────────────

#[test]
fn test_run_sort_produces_output_files() {
    // fixtures relative to rust/ (CARGO_MANIFEST_DIR)
    let fixtures = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .unwrap()
        .join("tests/fixtures");

    let fastq = fixtures.join("sample.fastq");
    let primers = fixtures.join("Primers.txt");
    let tags = fixtures.join("Tags.txt");

    let tmp = tempdir().unwrap();
    let original_dir = std::env::current_dir().unwrap();
    std::env::set_current_dir(tmp.path()).unwrap();

    let result = run(SortArgs {
        fastq: fastq.to_str().unwrap().to_string(),
        primers: primers.to_str().unwrap().to_string(),
        tags: tags.to_str().unwrap().to_string(),
        keep_primers_seq: false,
    });

    // Restore cwd regardless of outcome
    std::env::set_current_dir(&original_dir).unwrap();

    result.unwrap();

    // SummaryCounts.txt should exist
    let summary = tmp.path().join("SummaryCounts.txt");
    assert!(summary.exists(), "SummaryCounts.txt should exist");

    // Tag1_Tag2.txt should exist (AAAA=Tag1 forward, GGGG=Tag2 RC → Tag2)
    let tag_file = tmp.path().join("Tag1_Tag2.txt");
    assert!(tag_file.exists(), "Tag1_Tag2.txt should exist");

    let contents = std::fs::read_to_string(&tag_file).unwrap();

    // The barcode "ATATATATAT" should appear (read1 + read2, count=2)
    assert!(
        contents.contains("ATATATATAT"),
        "Tag1_Tag2.txt should contain ATATATATAT, got:\n{contents}"
    );

    // The barcode "GCGCGCGCGCGC" should appear (read3 + read4, count=2)
    assert!(
        contents.contains("GCGCGCGCGCGC"),
        "Tag1_Tag2.txt should contain GCGCGCGCGCGC, got:\n{contents}"
    );

    // Verify counts: each barcode should have count 2
    for line in contents.lines() {
        let parts: Vec<&str> = line.split('\t').collect();
        assert_eq!(parts.len(), 5, "Each line should have 5 tab-separated fields");
        let count: u32 = parts[3].parse().expect("count field should be numeric");
        assert_eq!(count, 2, "Each barcode should have count 2, got: {line}");
    }

    // Check SummaryCounts.txt content
    let summary_contents = std::fs::read_to_string(&summary).unwrap();
    let lines: Vec<&str> = summary_contents.lines().collect();
    assert!(lines[0].starts_with("#tagName1"), "First line should be header");
    assert_eq!(lines.len(), 2, "Should have header + 1 tag combo row");
    let row: Vec<&str> = lines[1].split('\t').collect();
    assert_eq!(row[0], "Tag1");
    assert_eq!(row[1], "Tag2");
    assert_eq!(row[2], "2", "Should have 2 unique sequences");
    assert_eq!(row[3], "4", "Sum of counts should be 4 (2+2)");
}

#[test]
fn test_get_pieces_info_no_panic_on_inverted_primers() {
    // Both primers present but in wrong order: TGCA then ACGT
    // Without the guard, prim_ini_prim > prim_fin_prim and &line[a..b] panics.
    // With the guard, this returns None safely.
    let dir = tempdir().unwrap();
    let tag_file = dir.path().join("tags.txt");
    std::fs::write(&tag_file, "AAAA\tTag1\nCCCC\tTag2\n").unwrap();
    let prim_file = dir.path().join("primers.txt");
    std::fs::write(&prim_file, "CO1\tACGT\tTGCA\n").unwrap();
    let tags = read_tags(tag_file.to_str().unwrap()).unwrap();
    let primers = read_primers(prim_file.to_str().unwrap()).unwrap();
    // RC(R)=TGCA comes first, F primer ACGT comes second — inverted orientation
    // In forward branch: finds ACGT at pos 8, then RC(R)=TGCA would match at pos 4 (before start)
    // → prim_ini_prim=12 > prim_fin_prim=4 → guard triggers → returns None, no panic
    let bad_read = "AAAATGCAACGTCCCC"; // TGCA at [4,8), ACGT at [8,12) — end primer before start
    let result = get_pieces_info(bad_read, &primers, &tags, false);
    // Guard should return None, not panic
    assert!(result.is_none());
}
