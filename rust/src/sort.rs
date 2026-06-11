use anyhow::{Context, Result};
use clap::Args;
use indexmap::IndexMap;
use ahash::HashMap;
use std::collections::HashSet;
use std::fs::{File, OpenOptions};
use std::io::{BufRead, BufReader, Write};

/// Demultiplex pooled FASTQ reads by primer and tag sequences.
///
/// Reads with structure [fwd_tag][fwd_primer][amplicon][rc(rev_primer)][rc(rev_tag)]
/// are classified into per-tag-pair collapsed output files.  Both forward and
/// reverse-complement orientations are detected automatically.  IUPAC ambiguity
/// codes are supported in primer sequences.
#[derive(Args)]
pub struct SortArgs {
    /// Input FASTQ file (plain or gzip)
    #[arg(long = "fq")]
    pub fastq: String,

    /// Primers file: tab-separated columns (Name, ForwardSeq, ReverseSeq), one primer
    /// set per line.  IUPAC ambiguity codes (R, Y, S, W, K, M, B, D, H, V, N) are
    /// supported and matched correctly in both exact and mismatch modes.
    #[arg(short = 'p', long = "primers")]
    pub primers: String,

    /// Tags file: tab-separated columns (TagSequence, TagName), one tag per line
    #[arg(short = 't', long = "tags")]
    pub tags: String,

    /// Retain primer sequences in the output amplicon.
    /// By default primer sequences are stripped from each read before writing.
    #[arg(long = "keep-primers-seq")]
    pub keep_primers_seq: bool,

    /// PCRsetsInfo file (same tab-separated format as `filter --ps-info`:
    /// SampleName, FwdTagName, RevTagName, PoolNumber).
    /// When supplied, writes SplitSummary.txt categorising reads into five classes:
    /// valid pair / same-tag pair / different-pool pair / one tag only / no tags found.
    #[arg(long = "ps-info")]
    pub psinfo: Option<String>,

    /// Maximum Hamming mismatches allowed per primer (applied independently to the
    /// forward and reverse primers).  When non-zero, the mismatch-tolerant classifier
    /// is used: tag candidates are pre-screened with `--mt` first, then primer
    /// mismatches are counted only for passing candidates.  Ambiguous reads (two
    /// equally-good matches) are written to the error output rather than classified.
    /// [default: 0]
    #[arg(short = 'm', long = "m", default_value_t = 0)]
    pub max_primer_mm: u32,

    /// Maximum Hamming mismatches allowed per tag (applied independently to tag1
    /// and tag2).  [default: 0]
    #[arg(long = "mt", default_value_t = 0)]
    pub max_tag_mm: u32,
}

// ── Reverse complement ─────────────────────────────────────────────────────────

/// Reverse complement of a DNA sequence (string version, used during file loading).
pub fn rc(seq: &str) -> String {
    seq.chars()
        .rev()
        .map(|c| match c {
            'A' => 'T', 'T' => 'A', 'C' => 'G', 'G' => 'C',
            'M' => 'K', 'K' => 'M', 'R' => 'Y', 'Y' => 'R',
            'W' => 'W', 'S' => 'S', 'V' => 'B', 'B' => 'V',
            'H' => 'D', 'D' => 'H', _ => c,
        })
        .collect()
}

/// Byte-level reverse complement — single pass over `&[u8]`, no UTF-8 decode.
pub fn rc_bytes(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| match b {
            b'A' => b'T', b'T' => b'A', b'C' => b'G', b'G' => b'C',
            b'M' => b'K', b'K' => b'M', b'R' => b'Y', b'Y' => b'R',
            b'W' => b'W', b'S' => b'S', b'V' => b'B', b'B' => b'V',
            b'H' => b'D', b'D' => b'H', _ => b,
        })
        .collect()
}

// ── IUPAC helpers ──────────────────────────────────────────────────────────────

/// Returns true if `primer_byte` (IUPAC code) is compatible with `read_byte` (A/C/G/T).
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
/// Returns `Some((start, end))` or `None`.
pub fn find_primer(primer: &[u8], seq: &[u8]) -> Option<(usize, usize)> {
    let plen = primer.len();
    let slen = seq.len();
    if plen > slen { return None; }
    for i in 0..=(slen - plen) {
        if primer.iter().zip(&seq[i..i + plen]).all(|(&p, &s)| iupac_matches(p, s)) {
            return Some((i, i + plen));
        }
    }
    None
}

/// IUPAC-aware Hamming distance between a primer and an equal-length read region.
pub fn hamming_iupac_primer(primer: &[u8], region: &[u8]) -> u32 {
    primer.iter().zip(region.iter())
        .map(|(&p, &r)| if iupac_matches(p, r) { 0u32 } else { 1u32 })
        .sum()
}

/// Exact Hamming distance between two equal-length byte slices (for plain-ACGT tags).
pub fn hamming_exact(a: &[u8], b: &[u8]) -> u32 {
    a.iter().zip(b.iter())
        .map(|(&x, &y)| if x == y { 0u32 } else { 1u32 })
        .sum()
}

// ── Data structures ────────────────────────────────────────────────────────────

pub struct PrimerEntry {
    /// [F_bytes, R_bytes] — start-side primers
    pub start_primers: Vec<Vec<u8>>,
    /// [RC(F)_bytes, RC(R)_bytes] — end-side primers
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

/// Outcome of attempting to classify one read.
pub enum SortOutcome {
    /// Both primers and both tags matched.
    Match(PieceInfo),
    /// Both primers matched but at least one tag did not (partial match → HAP_err).
    PartialMatch {
        tag1: Option<String>,
        tag2: Option<String>,
        primer_name: String,
        between: String,
    },
    /// No primers found (no match → no_tags_seqs).
    NoMatch,
}

/// One entry in the tag list, holding both orientations for Hamming candidate search.
pub struct TagEntry {
    pub name: String,
    pub fwd: Vec<u8>,
    pub rc: Vec<u8>,
}

/// Reverse-lookup maps for O(1) exact tag matching, plus an ordered list for
/// Hamming-distance candidate enumeration.
pub struct TagLookup {
    pub by_fwd: HashMap<Vec<u8>, String>,
    pub by_rc:  HashMap<Vec<u8>, String>,
    pub list:   Vec<TagEntry>,
}

// ── File readers ───────────────────────────────────────────────────────────────

/// Read a Tags file (TagSeq\tTagName per line).
pub fn read_tags(path: &str) -> Result<TagLookup> {
    let file = File::open(path).with_context(|| format!("Cannot open tags file: {path}"))?;
    let mut by_fwd: HashMap<Vec<u8>, String> = HashMap::default();
    let mut by_rc:  HashMap<Vec<u8>, String> = HashMap::default();
    let mut list:   Vec<TagEntry>            = Vec::new();

    for line in BufReader::new(file).lines() {
        let line = line?;
        let line = line.trim();
        if line.is_empty() { continue; }
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() < 2 { continue; }
        let seq  = parts[0];
        let name = parts[1];
        let rc_seq = rc(seq);
        by_fwd.insert(seq.as_bytes().to_vec(), name.to_string());
        by_rc .insert(rc_seq.as_bytes().to_vec(), name.to_string());
        list.push(TagEntry {
            name: name.to_string(),
            fwd:  seq.as_bytes().to_vec(),
            rc:   rc_seq.as_bytes().to_vec(),
        });
    }
    Ok(TagLookup { by_fwd, by_rc, list })
}

/// Read a Primers file (Name\tFwdSeq\tRevSeq per line).
pub fn read_primers(path: &str) -> Result<IndexMap<String, PrimerEntry>> {
    let file = File::open(path).with_context(|| format!("Cannot open primers file: {path}"))?;
    let mut primers: IndexMap<String, PrimerEntry> = IndexMap::new();

    for line in BufReader::new(file).lines() {
        let line = line?;
        let line = line.trim();
        if line.is_empty() { continue; }
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() < 3 { continue; }
        let name  = parts[0];
        let f_raw = parts[1];
        let r_raw = parts[2];
        let frc = rc(f_raw);
        let rrc = rc(r_raw);
        primers.insert(
            name.to_string(),
            PrimerEntry {
                start_primers: vec![f_raw.as_bytes().to_vec(), r_raw.as_bytes().to_vec()],
                end_primers:   vec![frc.as_bytes().to_vec(),   rrc.as_bytes().to_vec()],
            },
        );
    }
    Ok(primers)
}

/// Read a PCRsetsInfo file and return the set of valid (tag1_lc, tag2_lc) pairs
/// for the given pool number.
/// Format: SampleName\tTag1\tTag2\tPoolNumber
pub fn read_valid_pairs(path: &str, pool_num: u32) -> Result<HashSet<(String, String)>> {
    let file = File::open(path).with_context(|| format!("Cannot open PSinfo file: {path}"))?;
    let mut valid: HashSet<(String, String)> = HashSet::new();
    for line in BufReader::new(file).lines() {
        let line = line?;
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() < 4 { continue; }
        if let Ok(n) = parts[3].parse::<u32>() {
            if n == pool_num {
                valid.insert((parts[1].to_lowercase(), parts[2].to_lowercase()));
            }
        }
    }
    Ok(valid)
}

// ── HAP operations ─────────────────────────────────────────────────────────────

/// Insert or increment a sequence in the HAP map.
pub fn fill_hap(hap: &mut Hap, tag1: &str, tag2: &str, primer_name: &str, between: &str) {
    let key   = format!("{}_{}", tag1, tag2);
    let entry = hap.entry(key).or_insert_with(|| HapEntry {
        tag1: tag1.to_string(),
        tag2: tag2.to_string(),
        seqs: IndexMap::new(),
    });
    if let Some(se) = entry.seqs.get_mut(between) {
        se.count += 1;
    } else {
        entry.seqs.insert(
            between.to_string(),
            SeqEntry { count: 1, primer_name: primer_name.to_string() },
        );
    }
}

// ── Read classification functions ──────────────────────────────────────────────

/// Exact-match classification (m=0, mt=0 path).
pub fn get_pieces_info(
    line: &str,
    primers: &IndexMap<String, PrimerEntry>,
    tags: &TagLookup,
    keep_primers_seq: bool,
) -> SortOutcome {
    let seq = line.as_bytes();

    for (key, primer) in primers {
        // ── forward: F primer on left, RC(R) on right ──
        if let Some((fwd_start, fwd_end)) = find_primer(&primer.start_primers[0], seq) {
            let (amp_ini, tag1_end) = if keep_primers_seq {
                (fwd_start, fwd_start)
            } else {
                (fwd_end, fwd_start)
            };
            if let Some((rev_start, rev_end)) = find_primer(&primer.end_primers[1], seq) {
                let (amp_fin, tag2_start) = if keep_primers_seq {
                    (rev_end, rev_end)
                } else {
                    (rev_start, rev_end)
                };
                if amp_ini >= amp_fin { return SortOutcome::NoMatch; }
                let between = &line[amp_ini..amp_fin];
                if between.is_empty() { return SortOutcome::NoMatch; }
                let tn1 = tags.by_fwd.get(&seq[..tag1_end]).cloned();
                let tn2 = tags.by_rc .get(&seq[tag2_start..]).cloned();
                return match (tn1, tn2) {
                    (Some(t1), Some(t2)) => SortOutcome::Match(PieceInfo {
                        tag1: t1, tag2: t2,
                        primer_name: key.clone(),
                        between: between.to_string(),
                    }),
                    (t1, t2) => SortOutcome::PartialMatch {
                        tag1: t1, tag2: t2,
                        primer_name: key.clone(),
                        between: between.to_string(),
                    },
                };
            }
            return SortOutcome::NoMatch;
        }

        // ── reverse: R primer on left, RC(F) on right ──
        if let Some((fwd_start, fwd_end)) = find_primer(&primer.start_primers[1], seq) {
            let (amp_ini, tag1_end) = if keep_primers_seq {
                (fwd_start, fwd_start)
            } else {
                (fwd_end, fwd_start)
            };
            if let Some((rev_start, rev_end)) = find_primer(&primer.end_primers[0], seq) {
                let (amp_fin, tag2_start) = if keep_primers_seq {
                    (rev_end, rev_end)
                } else {
                    (rev_start, rev_end)
                };
                if amp_ini >= amp_fin { return SortOutcome::NoMatch; }
                let between_raw = &line[amp_ini..amp_fin];
                if between_raw.is_empty() { return SortOutcome::NoMatch; }
                let between = String::from_utf8(rc_bytes(between_raw.as_bytes()))
                    .expect("rc_bytes output is valid ASCII");
                // In reverse orientation tag roles are swapped
                let tn2 = tags.by_fwd.get(&seq[..tag1_end]).cloned();
                let tn1 = tags.by_rc .get(&seq[tag2_start..]).cloned();
                return match (tn1, tn2) {
                    (Some(t1), Some(t2)) => SortOutcome::Match(PieceInfo {
                        tag1: t1, tag2: t2,
                        primer_name: key.clone(),
                        between,
                    }),
                    (t1, t2) => SortOutcome::PartialMatch {
                        tag1: t1, tag2: t2,
                        primer_name: key.clone(),
                        between,
                    },
                };
            }
            return SortOutcome::NoMatch;
        }
    }
    SortOutcome::NoMatch
}

/// Tag-anchored mismatch classification (m>0 or mt>0 path).
///
/// Precomputes tag candidates once per read (O(n_tags) scan), then checks primer
/// mismatches only for those candidates — typically 1–2 per end.
pub fn get_pieces_info_mismatch(
    line: &str,
    primers: &IndexMap<String, PrimerEntry>,
    tags: &TagLookup,
    keep_primers_seq: bool,
    max_primer_mm: u32,
    max_tag_mm: u32,
) -> SortOutcome {
    let seq = line.as_bytes();
    let read_len = seq.len();

    // Scan all tags once per read.
    let mut t1_cands: Vec<(String, u32, usize)> = Vec::new();
    let mut t2_cands: Vec<(String, u32, usize)> = Vec::new();

    for entry in &tags.list {
        let t1_len = entry.fwd.len();
        if t1_len <= read_len {
            let mm = hamming_exact(&entry.fwd, &seq[..t1_len]);
            if mm <= max_tag_mm {
                t1_cands.push((entry.name.clone(), mm, t1_len));
            }
        }
        let t2_len = entry.rc.len();
        if t2_len <= read_len {
            let mm = hamming_exact(&entry.rc, &seq[read_len - t2_len..]);
            if mm <= max_tag_mm {
                t2_cands.push((entry.name.clone(), mm, t2_len));
            }
        }
    }

    if t1_cands.is_empty() || t2_cands.is_empty() {
        return SortOutcome::NoMatch;
    }

    let mut best_mm: Option<u32> = None;
    let mut best: Option<PieceInfo> = None;
    let mut ambiguous = false;

    for (key, primer) in primers {
        let f_raw   = &primer.start_primers[0];
        let r_raw   = &primer.start_primers[1];
        let frc_raw = &primer.end_primers[0];
        let rrc_raw = &primer.end_primers[1];

        // ── forward: [tag1_fwd][F][amplicon][RC(R)][tag2_rc] ──
        for (t1_name, t1_mm, t1_len) in &t1_cands {
            let f_end = t1_len + f_raw.len();
            if f_end > read_len { continue; }
            let f_mm = hamming_iupac_primer(f_raw, &seq[*t1_len..f_end]);
            if f_mm > max_primer_mm { continue; }

            for (t2_name, t2_mm, t2_len) in &t2_cands {
                let rrc_end = read_len - t2_len;
                if rrc_raw.len() > rrc_end { continue; }
                let rrc_start = rrc_end - rrc_raw.len();
                if rrc_start < f_end { continue; }
                let rrc_mm = hamming_iupac_primer(rrc_raw, &seq[rrc_start..rrc_end]);
                if rrc_mm > max_primer_mm { continue; }
                let (amp_start, amp_end) = if keep_primers_seq {
                    (*t1_len, rrc_end)
                } else {
                    (f_end, rrc_start)
                };
                if amp_start >= amp_end { continue; }
                let between = &line[amp_start..amp_end];
                if between.is_empty() { continue; }
                let total = t1_mm + f_mm + rrc_mm + t2_mm;
                update_best(total, t1_name, t2_name, key, between,
                            &mut best_mm, &mut best, &mut ambiguous);
            }
        }

        // ── reverse: [tag1_fwd][R][RC(amplicon)][RC(F)][tag2_rc] ──
        for (t1_name, t1_mm, t1_len) in &t1_cands {
            let r_end = t1_len + r_raw.len();
            if r_end > read_len { continue; }
            let r_mm = hamming_iupac_primer(r_raw, &seq[*t1_len..r_end]);
            if r_mm > max_primer_mm { continue; }

            for (t2_name, t2_mm, t2_len) in &t2_cands {
                let frc_end = read_len - t2_len;
                if frc_raw.len() > frc_end { continue; }
                let frc_start = frc_end - frc_raw.len();
                if frc_start < r_end { continue; }
                let frc_mm = hamming_iupac_primer(frc_raw, &seq[frc_start..frc_end]);
                if frc_mm > max_primer_mm { continue; }
                let (amp_start, amp_end) = if keep_primers_seq {
                    (*t1_len, frc_end)
                } else {
                    (r_end, frc_start)
                };
                if amp_start >= amp_end { continue; }
                let between_raw = &line[amp_start..amp_end];
                if between_raw.is_empty() { continue; }
                let between = String::from_utf8(rc_bytes(between_raw.as_bytes()))
                    .expect("rc_bytes output is valid ASCII");
                let total = t1_mm + r_mm + frc_mm + t2_mm;
                update_best(total, t1_name, t2_name, key, &between,
                            &mut best_mm, &mut best, &mut ambiguous);
            }
        }
    }

    if ambiguous { SortOutcome::NoMatch } else { best.map(SortOutcome::Match).unwrap_or(SortOutcome::NoMatch) }
}

/// Update the best-match tracker: replace on lower total, flag ambiguous on tie.
#[inline]
fn update_best(
    total: u32,
    t1_name: &str, t2_name: &str, primer_name: &str, between: &str,
    best_mm: &mut Option<u32>,
    best: &mut Option<PieceInfo>,
    ambiguous: &mut bool,
) {
    match *best_mm {
        None => {
            *best_mm   = Some(total);
            *best      = Some(PieceInfo {
                tag1: t1_name.to_string(), tag2: t2_name.to_string(),
                primer_name: primer_name.to_string(), between: between.to_string(),
            });
            *ambiguous = false;
        }
        Some(prev) if total < prev => {
            *best_mm   = Some(total);
            *best      = Some(PieceInfo {
                tag1: t1_name.to_string(), tag2: t2_name.to_string(),
                primer_name: primer_name.to_string(), between: between.to_string(),
            });
            *ambiguous = false;
        }
        Some(prev) if total == prev => { *ambiguous = true; }
        _ => {}
    }
}

// ── Output functions ───────────────────────────────────────────────────────────

fn print_sorted_collapsed_counted_seqs(hap: &Hap) -> Result<()> {
    for (tag_comb, entry) in hap {
        let filename = format!("{}.txt", tag_comb);
        let mut out = OpenOptions::new().write(true).create(true).truncate(true)
            .open(&filename)
            .with_context(|| format!("Cannot create output file: {filename}"))?;
        for (seq, se) in &entry.seqs {
            writeln!(out, "{}\t{}\t{}\t{}\t{}", se.primer_name, entry.tag1, entry.tag2, se.count, seq)?;
        }
    }
    Ok(())
}

fn print_summary_file(hap: &Hap) -> Result<()> {
    let mut out = OpenOptions::new().write(true).create(true).truncate(true)
        .open("SummaryCounts.txt")
        .context("Cannot create SummaryCounts.txt")?;
    writeln!(out, "#tagName1\ttagName2\tNumUniqSeqs\tSumTotalFreq")?;
    for entry in hap.values() {
        let num_uniq  = entry.seqs.len();
        let sum_freq: u32 = entry.seqs.values().map(|s| s.count).sum();
        writeln!(out, "{}\t{}\t{}\t{}", entry.tag1, entry.tag2, num_uniq, sum_freq)?;
    }
    Ok(())
}

fn pct(num: u64, denom: u64) -> String {
    if denom == 0 || num == 0 { "0.0".to_string() }
    else { format!("{:.2}", num as f64 / denom as f64 * 100.0) }
}

/// Total and unique sequence counts for a slice of HapEntry references.
fn entry_stats(entries: &[&HapEntry]) -> (u64, u64) {
    let total: u64 = entries.iter().flat_map(|e| e.seqs.values()).map(|s| s.count as u64).sum();
    let uniq:  u64 = entries.iter().map(|e| e.seqs.len() as u64).sum();
    (total, uniq)
}

fn write_detail_rows(entries: &[&HapEntry], out: &mut impl Write) -> Result<()> {
    for e in entries {
        let n_uniq  = e.seqs.len();
        let n_total: u64 = e.seqs.values().map(|s| s.count as u64).sum();
        writeln!(out, "{}\t{}\t{}\t{}", e.tag1, e.tag2, n_uniq, n_total)?;
    }
    Ok(())
}

/// Write `splitSummaryByPSInfo_{prefix}_{pool}.txt`, matching the Python output format.
pub fn print_split_summary_file(
    hap: &Hap,
    hap_err: &Hap,
    no_tags_seqs: &HashMap<String, u32>,
    prefix: &str,
    pool: &str,
    valid_pairs: Option<&HashSet<(String, String)>>,
) -> Result<()> {
    // Split HAP into "same" (valid pairs) and "diff" (everything else) using references.
    let is_valid = |e: &&HapEntry| -> bool {
        match valid_pairs {
            Some(vp) => vp.contains(&(e.tag1.to_lowercase(), e.tag2.to_lowercase())),
            None     => e.tag1.to_lowercase() == e.tag2.to_lowercase(),
        }
    };
    let same: Vec<&HapEntry> = hap.values().filter(|e| is_valid(e)).collect();
    let diff: Vec<&HapEntry> = hap.values().filter(|e| !is_valid(e)).collect();
    let err_entries: Vec<&HapEntry> = hap_err.values().collect();

    let (same_total,   same_uniq)   = entry_stats(&same);
    let (diff_total,   diff_uniq)   = entry_stats(&diff);
    let (onetag_total, onetag_uniq) = entry_stats(&err_entries);
    let notag_total: u64 = no_tags_seqs.values().map(|&v| v as u64).sum();
    let notag_uniq:  u64 = no_tags_seqs.len() as u64;

    let grand_total = same_total + diff_total + onetag_total + notag_total;
    let grand_uniq  = same_uniq  + diff_uniq  + onetag_uniq  + notag_uniq;

    let filename = format!("splitSummaryByPSInfo_{}_{}.txt", prefix, pool);
    let mut out = OpenOptions::new().write(true).create(true).truncate(true)
        .open(&filename)
        .with_context(|| format!("Cannot create {filename}"))?;

    // Header row
    writeln!(out, "{:<52}\tTotal seqs\tTotal unique seqs\t% Total seqs\t%Total unique seqs", "")?;

    // Summary rows
    writeln!(out, "{:<52}\t{}\t{}\t{}\t{}",
        "Tag combinations where the tag pair was used",
        same_total, same_uniq, pct(same_total, grand_total), pct(same_uniq, grand_uniq))?;
    writeln!(out, "Tag combinations where both tags used")?;
    writeln!(out, "{:<52}\t{}\t{}\t{}\t{}",
        "but not in this combination",
        diff_total, diff_uniq, pct(diff_total, grand_total), pct(diff_uniq, grand_uniq))?;
    writeln!(out, "{:<52}\t{}\t{}\t{}\t{}",
        "Tag combinations where only one of the tags was used",
        onetag_total, onetag_uniq, pct(onetag_total, grand_total), pct(onetag_uniq, grand_uniq))?;
    writeln!(out, "{:<52}\t{}\t{}\t{}\t{}",
        "Tag combinations where neither tag was used",
        notag_total, notag_uniq, pct(notag_total, grand_total), pct(notag_uniq, grand_uniq))?;
    writeln!(out, "{:<52}\t{}\t{}\t100.00\t100.00", "Total", grand_total, grand_uniq)?;
    writeln!(out)?;

    // Detail sections
    writeln!(out, "Tag combinations where the tag pair was used.")?;
    writeln!(out, "---------------------------------------------")?;
    writeln!(out, "Tag1Name\tTag2Name\tNumUniqSeqs\tSumTotalFreq")?;
    write_detail_rows(&same, &mut out)?;
    writeln!(out)?;

    writeln!(out, "Tag combinations where both tags used - but not in this combination.")?;
    writeln!(out, "--------------------------------------------------------------------")?;
    writeln!(out, "Tag1Name\tTag2Name\tNumUniqSeqs\tSumTotalFreq")?;
    write_detail_rows(&diff, &mut out)?;
    writeln!(out)?;

    writeln!(out, "Tag combinations where only one of the tags was used.")?;
    writeln!(out, "-----------------------------------------------------")?;
    writeln!(out, "Tag1Name\tTag2Name\tNumUniqSeqs\tSumTotalFreq")?;
    write_detail_rows(&err_entries, &mut out)?;
    writeln!(out)?;

    writeln!(out, "Tag combinations where neither tag was used.")?;
    writeln!(out, "--------------------------------------------")?;
    writeln!(out, "Tag1Name\tTag2Name\tNumUniqSeqs\tSumTotalFreq")?;

    Ok(())
}

// ── Entry point ────────────────────────────────────────────────────────────────

pub fn run(args: SortArgs) -> Result<()> {
    let tags    = read_tags(&args.tags)?;
    let primers = read_primers(&args.primers)?;

    let mut hap:          Hap = IndexMap::new();
    let mut hap_err:      Hap = IndexMap::new();
    let mut no_tags_seqs: HashMap<String, u32> = HashMap::default();
    let mut count_errors: u32 = 0;

    let mut reader = needletail::parse_fastx_file(&args.fastq)
        .with_context(|| format!("Cannot open FASTQ: {}", args.fastq))?;

    while let Some(record) = reader.next() {
        let seqrec = record.context("Error reading FASTQ record")?;
        let seq_bytes = seqrec.seq();
        let seq = std::str::from_utf8(&seq_bytes)
            .context("Non-UTF-8 bytes in FASTQ sequence")?;
        if seq.is_empty() { continue; }

        let outcome = if args.max_primer_mm == 0 && args.max_tag_mm == 0 {
            get_pieces_info(seq, &primers, &tags, args.keep_primers_seq)
        } else {
            get_pieces_info_mismatch(
                seq, &primers, &tags,
                args.keep_primers_seq,
                args.max_primer_mm,
                args.max_tag_mm,
            )
        };

        match outcome {
            SortOutcome::Match(info) => {
                fill_hap(&mut hap, &info.tag1, &info.tag2, &info.primer_name, &info.between);
            }
            SortOutcome::PartialMatch { tag1, tag2, primer_name, between } => {
                let t1 = tag1.as_deref().unwrap_or("");
                let t2 = tag2.as_deref().unwrap_or("");
                if !t1.is_empty() || !t2.is_empty() {
                    fill_hap(&mut hap_err,
                             if t1.is_empty() { "none" } else { t1 },
                             if t2.is_empty() { "none" } else { t2 },
                             &primer_name, &between);
                } else {
                    *no_tags_seqs.entry(between).or_insert(0) += 1;
                }
                count_errors += 1;
            }
            SortOutcome::NoMatch => {
                *no_tags_seqs.entry(seq.to_string()).or_insert(0) += 1;
                count_errors += 1;
            }
        }
    }

    // Derive prefix (first '_'-delimited token of FASTQ basename) and pool (cwd name).
    let fq_path = std::path::Path::new(&args.fastq);
    let fq_abs  = std::fs::canonicalize(fq_path).unwrap_or_else(|_| fq_path.to_path_buf());
    let mut basename = fq_abs.file_name()
        .and_then(|n| n.to_str()).unwrap_or("").to_string();
    for ext in &[".fastq.gz", ".fq.gz", ".fastq", ".fq"] {
        if basename.ends_with(ext) {
            basename.truncate(basename.len() - ext.len());
            break;
        }
    }
    let prefix = basename.split('_').next().unwrap_or("").to_string();
    let pool   = std::env::current_dir().ok()
        .and_then(|p| p.file_name().map(|n| n.to_string_lossy().to_string()))
        .unwrap_or_default();

    // Resolve valid pairs from PCRsetsInfo if provided.
    let valid_pairs_opt: Option<HashSet<(String, String)>> = if let Some(ref psinfo_path) = args.psinfo {
        let pool_digits: String = pool.chars().filter(|c| c.is_ascii_digit()).collect();
        if let Ok(pool_num) = pool_digits.parse::<u32>() {
            Some(read_valid_pairs(psinfo_path, pool_num)?)
        } else {
            None
        }
    } else {
        None
    };

    print_sorted_collapsed_counted_seqs(&hap)?;
    print_summary_file(&hap)?;

    if args.psinfo.is_some() {
        print_split_summary_file(
            &hap, &hap_err, &no_tags_seqs,
            &prefix, &pool,
            valid_pairs_opt.as_ref(),
        )?;
    }

    println!(
        "Number of erroneous sequences (with errors in the sequence of primer or tags, or no barcode amplified): {}",
        count_errors
    );

    Ok(())
}
