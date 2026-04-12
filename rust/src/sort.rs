use anyhow::{Context, Result};
use clap::Args;
use indexmap::IndexMap;
use regex::Regex;
use std::collections::HashMap;
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

/// Expand an ambiguous nucleotide character to its regex pattern string.
fn ambig_expand(seq: &str) -> String {
    seq.chars()
        .map(|c| match c {
            'A' => "A",
            'B' => "[CGT]",
            'C' => "C",
            'D' => "[AGT]",
            'G' => "G",
            'H' => "[ACT]",
            'K' => "[GT]",
            'M' => "[AC]",
            'N' => "[ACGT]",
            'R' => "[AG]",
            'S' => "[CG]",
            'T' => "T",
            'V' => "[ACG]",
            'W' => "[AT]",
            'Y' => "[CT]",
            _ => ".",
        })
        .collect()
}

pub struct PrimerEntry {
    pub a_side: Vec<String>,    // [F_pattern, R_pattern] as regex strings
    pub b_side: Vec<String>,    // [RC(F)_pattern, RC(R)_pattern] as regex strings
    pub a_side_re: Vec<Regex>,  // compiled a_side regexes
    pub b_side_re: Vec<Regex>,  // compiled b_side regexes
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

/// Read a Tags file (TagSeq\tTagName per line).
/// Returns HashMap<TagName, [forward_seq, RC(forward_seq)]>
pub fn read_tags(path: &str) -> Result<HashMap<String, Vec<String>>> {
    let file = File::open(path).with_context(|| format!("Cannot open tags file: {path}"))?;
    let reader = BufReader::new(file);
    let mut tags: HashMap<String, Vec<String>> = HashMap::new();

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
        let entry = tags.entry(name.to_string()).or_default();
        entry.push(seq.to_string());
        entry.push(rc(seq));
    }

    Ok(tags)
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

        // a_side: [F_pattern, R_pattern]
        let f_pat = ambig_expand(f_raw);
        let r_pat = ambig_expand(r_raw);
        // b_side: [RC(F)_pattern, RC(R)_pattern]
        let frc = rc(f_raw);
        let rrc = rc(r_raw);
        let frc_pat = ambig_expand(&frc);
        let rrc_pat = ambig_expand(&rrc);

        let a_side = vec![f_pat.clone(), r_pat.clone()];
        let b_side = vec![frc_pat.clone(), rrc_pat.clone()];
        let a_side_re = vec![
            Regex::new(&f_pat).with_context(|| format!("Bad forward regex: {f_pat}"))?,
            Regex::new(&r_pat).with_context(|| format!("Bad reverse regex: {r_pat}"))?,
        ];
        let b_side_re = vec![
            Regex::new(&frc_pat).with_context(|| format!("Bad RC(F) regex: {frc_pat}"))?,
            Regex::new(&rrc_pat).with_context(|| format!("Bad RC(R) regex: {rrc_pat}"))?,
        ];

        primers.insert(
            name.to_string(),
            PrimerEntry {
                a_side,
                b_side,
                a_side_re,
                b_side_re,
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
    tags: &HashMap<String, Vec<String>>,
    keep_primers_seq: bool,
) -> Option<PieceInfo> {
    for (key, primer) in primers {
        // Try forward orientation: find a_side[0] (F pattern) at start region
        if let Some(m_start) = primer.a_side_re[0].find(line) {
            let (prim_ini_prim, prim_ini_tags) = if keep_primers_seq {
                (m_start.start(), m_start.start())
            } else {
                (m_start.end(), m_start.start())
            };
            // Now find b_side[1] (RC(R) pattern) — the end anchor
            if let Some(m_end) = primer.b_side_re[1].find(line) {
                let (prim_fin_prim, prim_fin_tags) = if keep_primers_seq {
                    (m_end.end(), m_end.end())
                } else {
                    (m_end.start(), m_end.end())
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
                // tag1: forward seq matches tag1_str (TAGS[t][0] == tag1_str)
                let tag_name1 = tags.iter().find(|(_, v)| v[0] == tag1_str).map(|(k, _)| k.clone());
                // tag2: RC seq matches tag2_str (TAGS[t][1] == tag2_str)
                let tag_name2 = tags.iter().find(|(_, v)| v[1] == tag2_str).map(|(k, _)| k.clone());
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
            // forward primer found but end primer not found → error, don't try reverse
            return None;
        } else {
            // Try reverse orientation: find a_side[1] (R pattern) at start
            if let Some(m_start) = primer.a_side_re[1].find(line) {
                let (prim_ini_prim, prim_ini_tags) = if keep_primers_seq {
                    (m_start.start(), m_start.start())
                } else {
                    (m_start.end(), m_start.start())
                };
                // Find b_side[0] (RC(F) pattern) at end
                if let Some(m_end) = primer.b_side_re[0].find(line) {
                    let (prim_fin_prim, prim_fin_tags) = if keep_primers_seq {
                        (m_end.end(), m_end.end())
                    } else {
                        (m_end.start(), m_end.end())
                    };
                    if prim_ini_prim >= prim_fin_prim {
                        return None;
                    }
                    let between_raw = &line[prim_ini_prim..prim_fin_prim];
                    if between_raw.is_empty() {
                        return None;
                    }
                    let between = rc(between_raw);
                    let tag1_str = &line[..prim_ini_tags];
                    let tag2_str = &line[prim_fin_tags..];
                    // In reverse: tagName2 matches forward, tagName1 matches RC
                    let tag_name2 = tags.iter().find(|(_, v)| v[0] == tag1_str).map(|(k, _)| k.clone());
                    let tag_name1 = tags.iter().find(|(_, v)| v[1] == tag2_str).map(|(k, _)| k.clone());
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

    let file = File::open(&args.fastq).with_context(|| format!("Cannot open FASTQ: {}", args.fastq))?;
    let mut reader = BufReader::new(file);
    let mut line_buf = String::new();

    loop {
        // Read header line
        line_buf.clear();
        let n = reader.read_line(&mut line_buf)?;
        if n == 0 {
            break; // EOF
        }
        // Read seq line
        line_buf.clear();
        let n = reader.read_line(&mut line_buf)?;
        if n == 0 {
            break;
        }
        let seq = line_buf.trim_end_matches('\n').trim_end_matches('\r').to_string();

        // Read "+" line
        line_buf.clear();
        reader.read_line(&mut line_buf)?;
        // Read qual line
        line_buf.clear();
        reader.read_line(&mut line_buf)?;

        if seq.is_empty() {
            continue;
        }

        match get_pieces_info(&seq, &primers, &tags, args.keep_primers_seq) {
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
