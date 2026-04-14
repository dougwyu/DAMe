use anyhow::{Context, Result};
use clap::Args;
use indexmap::IndexMap;
use ahash::HashMap;
use std::fs::{File, OpenOptions};
use std::io::{BufRead, BufReader, Write};

#[derive(Args)]
pub struct SortArgs {
    #[arg(long = "fq")]
    pub fastq: String,
    #[arg(long = "primers")]
    pub primers: String,
    #[arg(long = "tags")]
    pub tags: String,
    #[arg(long = "keep-primers-seq")]
    pub keep_primers_seq: bool,
}

/// Reverse complement of a DNA sequence.
/// Complement mapping: A↔T, C↔G, G↔C, T↔A, M↔K, R↔Y, W↔W, S↔S, Y↔R, K↔M, V↔B, H↔D, D↔H, B↔V
pub fn rc(seq: &str) -> String {
    seq.chars()
        .rev()
        .map(|c| match c {
            'A' => 'T',
            'T' => 'A',
            'C' => 'G',
            'G' => 'C',
            'M' => 'K',
            'K' => 'M',
            'R' => 'Y',
            'Y' => 'R',
            'W' => 'W',
            'S' => 'S',
            'V' => 'B',
            'B' => 'V',
            'H' => 'D',
            'D' => 'H',
            _ => c,
        })
        .collect()
}

/// Byte-level reverse complement for DNA sequences.
/// Equivalent to `rc()` but operates on `&[u8]` directly — no UTF-8 decode, single pass.
/// Input bytes should be uppercase ASCII (standard FASTQ).
pub fn rc_bytes(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| match b {
            b'A' => b'T',
            b'T' => b'A',
            b'C' => b'G',
            b'G' => b'C',
            b'M' => b'K',
            b'K' => b'M',
            b'R' => b'Y',
            b'Y' => b'R',
            b'W' => b'W',
            b'S' => b'S',
            b'V' => b'B',
            b'B' => b'V',
            b'H' => b'D',
            b'D' => b'H',
            _ => b,
        })
        .collect()
}

/// Returns true if `primer_byte` (an IUPAC code) is compatible with `read_byte` (A/C/G/T).
/// Both bytes are expected to be uppercase ASCII.
pub fn iupac_matches(primer_byte: u8, read_byte: u8) -> bool {
    match primer_byte {
        b'A' => read_byte == b'A',
        b'C' => read_byte == b'C',
        b'G' => read_byte == b'G',
        b'T' => read_byte == b'T',
        b'R' => matches!(read_byte, b'A' | b'G'),
        b'Y' => matches!(read_byte, b'C' | b'T'),
        b'S' => matches!(read_byte, b'C' | b'G'),
        b'W' => matches!(read_byte, b'A' | b'T'),
        b'K' => matches!(read_byte, b'G' | b'T'),
        b'M' => matches!(read_byte, b'A' | b'C'),
        b'B' => matches!(read_byte, b'C' | b'G' | b'T'),
        b'D' => matches!(read_byte, b'A' | b'G' | b'T'),
        b'H' => matches!(read_byte, b'A' | b'C' | b'T'),
        b'V' => matches!(read_byte, b'A' | b'C' | b'G'),
        b'N' => matches!(read_byte, b'A' | b'C' | b'G' | b'T'),
        _ => false,
    }
}

/// Find the leftmost occurrence of `primer` in `seq` using IUPAC matching.
/// Returns `Some((start, end))` where `end = start + primer.len()`, or `None` if not found.
pub fn find_primer(primer: &[u8], seq: &[u8]) -> Option<(usize, usize)> {
    let plen = primer.len();
    let slen = seq.len();
    if plen > slen {
        return None;
    }
    for i in 0..=(slen - plen) {
        if primer.iter().zip(&seq[i..i + plen]).all(|(&p, &s)| iupac_matches(p, s)) {
            return Some((i, i + plen));
        }
    }
    None
}

pub struct PrimerEntry {
    /// [F_bytes, R_bytes] — start-side primers as raw bytes for find_primer
    pub start_primers: Vec<Vec<u8>>,
    /// [RC(F)_bytes, RC(R)_bytes] — end-side primers as raw bytes for find_primer
    pub end_primers: Vec<Vec<u8>>,
}

pub struct SeqEntry {
    pub count: u32,
    pub primer_name: String,
}

pub struct HapEntry {
    pub tag1: String,
    pub tag2: String,
    pub seqs: IndexMap<String, SeqEntry>,
}

pub type Hap = IndexMap<String, HapEntry>;

pub struct PieceInfo {
    pub tag1: String,
    pub tag2: String,
    pub primer_name: String,
    pub between: String,
}

/// Pre-built O(1) reverse lookup for tag sequences.
/// `by_fwd` maps forward tag bytes → tag name; `by_rc` maps RC tag bytes → tag name.
pub struct TagLookup {
    pub by_fwd: HashMap<Vec<u8>, String>,
    pub by_rc: HashMap<Vec<u8>, String>,
}

/// Read a Tags file (TagSeq\tTagName per line) and build O(1) reverse lookup maps.
pub fn read_tags(path: &str) -> Result<TagLookup> {
    let file = File::open(path).with_context(|| format!("Cannot open tags file: {path}"))?;
    let reader = BufReader::new(file);
    let mut by_fwd: HashMap<Vec<u8>, String> = HashMap::default();
    let mut by_rc: HashMap<Vec<u8>, String> = HashMap::default();

    for line in reader.lines() {
        let line = line?;
        let line = line.trim();
        if line.is_empty() {
            continue;
        }
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() < 2 {
            continue;
        }
        let seq = parts[0];
        let name = parts[1];
        by_fwd.insert(seq.as_bytes().to_vec(), name.to_string());
        by_rc.insert(rc(seq).as_bytes().to_vec(), name.to_string());
    }

    Ok(TagLookup { by_fwd, by_rc })
}

/// Read a Primers file (Name\tFwdSeq\tRevSeq per line).
/// Returns IndexMap<Name, PrimerEntry>
pub fn read_primers(path: &str) -> Result<IndexMap<String, PrimerEntry>> {
    let file = File::open(path).with_context(|| format!("Cannot open primers file: {path}"))?;
    let reader = BufReader::new(file);
    let mut primers: IndexMap<String, PrimerEntry> = IndexMap::new();

    for line in reader.lines() {
        let line = line?;
        let line = line.trim();
        if line.is_empty() {
            continue;
        }
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() < 3 {
            continue;
        }
        let name = parts[0];
        let f_raw = parts[1];
        let r_raw = parts[2];

        let frc = rc(f_raw);
        let rrc = rc(r_raw);

        primers.insert(
            name.to_string(),
            PrimerEntry {
                start_primers: vec![f_raw.as_bytes().to_vec(), r_raw.as_bytes().to_vec()],
                end_primers: vec![frc.as_bytes().to_vec(), rrc.as_bytes().to_vec()],
            },
        );
    }

    Ok(primers)
}

/// Fill the HAP map with a tag combination and barcode sequence.
pub fn fill_hap(hap: &mut Hap, tag1: &str, tag2: &str, primer_name: &str, between: &str) {
    let key = format!("{}_{}", tag1, tag2);
    let entry = hap.entry(key).or_insert_with(|| HapEntry {
        tag1: tag1.to_string(),
        tag2: tag2.to_string(),
        seqs: IndexMap::new(),
    });
    let seq_entry = entry.seqs.entry(between.to_string()).or_insert_with(|| SeqEntry {
        count: 0,
        primer_name: primer_name.to_string(),
    });
    seq_entry.count += 1;
}

/// Try to extract tag info from a sequence line.
/// Returns Some(PieceInfo) on success, None on failure/error.
pub fn get_pieces_info(
    line: &str,
    primers: &IndexMap<String, PrimerEntry>,
    tags: &TagLookup,
    keep_primers_seq: bool,
) -> Option<PieceInfo> {
    let seq = line.as_bytes();

    for (key, primer) in primers {
        // Try forward orientation: start_primers[0] (F) at left, end_primers[1] (RC(R)) at right
        if let Some((fwd_start, fwd_end)) = find_primer(&primer.start_primers[0], seq) {
            let (prim_ini_prim, prim_ini_tags) = if keep_primers_seq {
                (fwd_start, fwd_start)
            } else {
                (fwd_end, fwd_start)
            };
            if let Some((rev_start, rev_end)) = find_primer(&primer.end_primers[1], seq) {
                let (prim_fin_prim, prim_fin_tags) = if keep_primers_seq {
                    (rev_end, rev_end)
                } else {
                    (rev_start, rev_end)
                };
                if prim_ini_prim >= prim_fin_prim {
                    return None;
                }
                let between = &line[prim_ini_prim..prim_fin_prim];
                if between.is_empty() {
                    return None;
                }
                let tag1_str = &line[..prim_ini_tags];
                let tag2_str = &line[prim_fin_tags..];
                let tag_name1 = tags.by_fwd.get(tag1_str.as_bytes()).cloned();
                let tag_name2 = tags.by_rc.get(tag2_str.as_bytes()).cloned();
                if let (Some(tn1), Some(tn2)) = (tag_name1, tag_name2) {
                    return Some(PieceInfo {
                        tag1: tn1,
                        tag2: tn2,
                        primer_name: key.clone(),
                        between: between.to_string(),
                    });
                }
                return None;
            }
            // Forward start primer found but end primer not found → error
            return None;
        } else {
            // Try reverse orientation: start_primers[1] (R) at left, end_primers[0] (RC(F)) at right
            if let Some((fwd_start, fwd_end)) = find_primer(&primer.start_primers[1], seq) {
                let (prim_ini_prim, prim_ini_tags) = if keep_primers_seq {
                    (fwd_start, fwd_start)
                } else {
                    (fwd_end, fwd_start)
                };
                if let Some((rev_start, rev_end)) = find_primer(&primer.end_primers[0], seq) {
                    let (prim_fin_prim, prim_fin_tags) = if keep_primers_seq {
                        (rev_end, rev_end)
                    } else {
                        (rev_start, rev_end)
                    };
                    if prim_ini_prim >= prim_fin_prim {
                        return None;
                    }
                    let between_raw = &line[prim_ini_prim..prim_fin_prim];
                    if between_raw.is_empty() {
                        return None;
                    }
                    // rc_bytes is byte-level and single-pass; between_raw is ASCII DNA
                    let between_bytes = rc_bytes(between_raw.as_bytes());
                    let between = String::from_utf8(between_bytes)
                        .expect("rc_bytes output is valid ASCII by construction");
                    let tag1_str = &line[..prim_ini_tags];
                    let tag2_str = &line[prim_fin_tags..];
                    // In reverse orientation: tag roles are swapped
                    let tag_name2 = tags.by_fwd.get(tag1_str.as_bytes()).cloned();
                    let tag_name1 = tags.by_rc.get(tag2_str.as_bytes()).cloned();
                    if let (Some(tn1), Some(tn2)) = (tag_name1, tag_name2) {
                        return Some(PieceInfo {
                            tag1: tn1,
                            tag2: tn2,
                            primer_name: key.clone(),
                            between,
                        });
                    }
                    return None;
                }
                return None;
            }
        }
    }
    None
}

fn print_sorted_collapsed_counted_seqs(hap: &Hap) -> Result<()> {
    for (tag_comb, entry) in hap {
        let filename = format!("{}.txt", tag_comb);
        let mut out = OpenOptions::new()
            .write(true)
            .create(true)
            .truncate(true)
            .open(&filename)
            .with_context(|| format!("Cannot create output file: {filename}"))?;
        for (seq, seq_entry) in &entry.seqs {
            writeln!(
                out,
                "{}\t{}\t{}\t{}\t{}",
                seq_entry.primer_name, entry.tag1, entry.tag2, seq_entry.count, seq
            )?;
        }
    }
    Ok(())
}

fn print_summary_file(hap: &Hap) -> Result<()> {
    let mut out = OpenOptions::new()
        .write(true)
        .create(true)
        .truncate(true)
        .open("SummaryCounts.txt")
        .context("Cannot create SummaryCounts.txt")?;
    writeln!(out, "#tagName1\ttagName2\tNumUniqSeqs\tSumTotalFreq")?;
    for (_tag_comb, entry) in hap {
        let num_uniq = entry.seqs.len();
        let sum_freq: u32 = entry.seqs.values().map(|s| s.count).sum();
        writeln!(
            out,
            "{}\t{}\t{}\t{}",
            entry.tag1, entry.tag2, num_uniq, sum_freq
        )?;
    }
    Ok(())
}

pub fn run(args: SortArgs) -> Result<()> {
    let tags = read_tags(&args.tags)?;
    let primers = read_primers(&args.primers)?;
    let mut hap: Hap = IndexMap::new();
    let mut count_errors: u32 = 0;

    let mut reader = needletail::parse_fastx_file(&args.fastq)
        .with_context(|| format!("Cannot open FASTQ: {}", args.fastq))?;

    while let Some(record) = reader.next() {
        let seqrec = record.context("Error reading FASTQ record")?;
        let seq_bytes = seqrec.seq();
        let seq = std::str::from_utf8(&seq_bytes)
            .context("Non-UTF-8 bytes in FASTQ sequence — file may be corrupt")?;
        if seq.is_empty() {
            continue;
        }
        match get_pieces_info(seq, &primers, &tags, args.keep_primers_seq) {
            Some(info) => {
                fill_hap(&mut hap, &info.tag1, &info.tag2, &info.primer_name, &info.between);
            }
            None => {
                count_errors += 1;
            }
        }
    }

    print_sorted_collapsed_counted_seqs(&hap)?;
    print_summary_file(&hap)?;

    println!(
        "Number of erroneous sequences (with errors in the sequence of primer or tags, or no barcode amplified): {}",
        count_errors
    );

    Ok(())
}
