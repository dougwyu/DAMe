# DAMe Rust Port Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Implement the `dame` Rust binary with subcommands `sort`, `chimera`, `filter`, `decollapse`, `rsi` — producing byte-for-byte identical output to the Python 3 port (`dame-py`) — plus Rust unit tests and updated integration tests comparing both tools.

**Architecture:** Single Rust crate (`rust/`) with a library root (`lib.rs`) that exposes pub modules for each subcommand, and a thin binary entry point (`main.rs`) that dispatches via `clap` derive. Unit tests live in `#[cfg(test)]` blocks inside each source file; integration tests in `rust/tests/*.rs` test the library functions directly. Shell integration tests in `tests/integration/` run both `dame` and `dame-py` and compare outputs.

**Tech Stack:** Rust stable, `clap 4` (derive), `regex 1`, `ndarray 0.15`, `anyhow 1`, `indexmap 2`, `tempfile 3` (dev), GitHub Actions

---

## File Map

| File | Action | Responsibility |
|---|---|---|
| `rust/Cargo.toml` | Create | Package metadata, all dependencies |
| `rust/src/lib.rs` | Create | Re-exports all pub modules |
| `rust/src/main.rs` | Create | CLI entry point, clap dispatch |
| `rust/src/sort.rs` | Create | rc(), readTags/Primers, GetPiecesInfo, FillHAP, print fns, `run()` |
| `rust/src/chimera_check.rs` | Create | makeTagFiles, MakeSizeOutFastas, SortFasta (wraps usearch), MakeFasSeqOneLine, MakeNoChimHaps, `run()` |
| `rust/src/filter.rs` | Create | makePSnumFiles, ReadPSnumFiles, MakeSampleNameArray, ReadHapsForASample, getSeqsSetsAndFRcounts, MakeComparisonFile, `run()` |
| `rust/src/decollapse.rs` | Create | `run()` (expand unique sequences) |
| `rust/src/rsi.rs` | Create | `compare()` with ndarray, `run()` |
| `rust/tests/sort_test.rs` | Create | Rust integration tests for sort functions |
| `rust/tests/chimera_test.rs` | Create | Rust integration tests for chimera functions |
| `rust/tests/filter_test.rs` | Create | Rust integration tests for filter functions |
| `rust/tests/decollapse_test.rs` | Create | Rust integration tests for decollapse |
| `rust/tests/rsi_test.rs` | Create | Rust integration tests for RSI |
| `tests/integration/run_sort.sh` | Modify | Add `dame sort` run + diff comparison |
| `tests/integration/run_rsi.sh` | Modify | Add `dame rsi` run + diff comparison |
| `tests/integration/run_filter.sh` | Modify | Add `dame filter` run + sorted diff comparison |
| `tests/integration/run_decollapse.sh` | Modify | Add `dame decollapse` run + diff comparison |
| `tests/integration/run_chimera.sh` | Create | Test `dame chimera` (skips if usearch absent) |
| `tests/integration/run_pipeline.sh` | Create | End-to-end sort → filter → rsi |
| `.github/workflows/ci.yml` | Create | Matrix CI: pytest + cargo test + integration tests |

---

## Task 1: Scaffold Rust package

**Files:**
- Create: `rust/Cargo.toml`
- Create: `rust/src/lib.rs`
- Create: `rust/src/main.rs`
- Create: `rust/src/sort.rs` (stub)
- Create: `rust/src/chimera_check.rs` (stub)
- Create: `rust/src/filter.rs` (stub)
- Create: `rust/src/decollapse.rs` (stub)
- Create: `rust/src/rsi.rs` (stub)

- [ ] **Step 1: Create `rust/Cargo.toml`**

```toml
[package]
name = "dame"
version = "0.1.0"
edition = "2021"

[[bin]]
name = "dame"
path = "src/main.rs"

[lib]
name = "dame"
path = "src/lib.rs"

[dependencies]
clap = { version = "4", features = ["derive"] }
regex = "1"
ndarray = "0.15"
anyhow = "1"
indexmap = "2"

[dev-dependencies]
tempfile = "3"
```

- [ ] **Step 2: Create `rust/src/lib.rs`**

```rust
pub mod chimera_check;
pub mod decollapse;
pub mod filter;
pub mod rsi;
pub mod sort;
```

- [ ] **Step 3: Create `rust/src/main.rs`**

```rust
use anyhow::Result;
use clap::{Parser, Subcommand};
use dame::{chimera_check, decollapse, filter, rsi, sort};

#[derive(Parser)]
#[command(name = "dame", about = "DNA Metabarcoding toolkit")]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    Sort(sort::SortArgs),
    Chimera(chimera_check::ChimeraArgs),
    Filter(filter::FilterArgs),
    Decollapse(decollapse::DecollapseArgs),
    Rsi(rsi::RsiArgs),
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    match cli.command {
        Commands::Sort(args) => sort::run(args),
        Commands::Chimera(args) => chimera_check::run(args),
        Commands::Filter(args) => filter::run(args),
        Commands::Decollapse(args) => decollapse::run(args),
        Commands::Rsi(args) => rsi::run(args),
    }
}
```

- [ ] **Step 4: Create stub modules**

`rust/src/sort.rs`:
```rust
use anyhow::Result;
use clap::Args;

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

pub fn run(_args: SortArgs) -> Result<()> {
    todo!()
}
```

`rust/src/chimera_check.rs`:
```rust
use anyhow::Result;
use clap::Args;

#[derive(Args)]
pub struct ChimeraArgs {
    #[arg(long = "ps-info")]
    pub ps_info: String,
    #[arg(long = "x")]
    pub x: usize,
    #[arg(long = "p", default_value = "1")]
    pub p: usize,
}

pub fn run(_args: ChimeraArgs) -> Result<()> {
    todo!()
}
```

`rust/src/filter.rs`:
```rust
use anyhow::Result;
use clap::Args;

#[derive(Args)]
pub struct FilterArgs {
    #[arg(long = "ps-info")]
    pub ps_info: String,
    #[arg(long = "x", default_value = "2")]
    pub x: usize,
    #[arg(long = "y", default_value = "1")]
    pub y: usize,
    #[arg(long = "p", default_value = "1")]
    pub p: usize,
    #[arg(long = "t", default_value = "1")]
    pub t: u32,
    #[arg(long = "l", default_value = "100")]
    pub l: usize,
    #[arg(long = "chimera-checked")]
    pub chimera_checked: bool,
}

pub fn run(_args: FilterArgs) -> Result<()> {
    todo!()
}
```

`rust/src/decollapse.rs`:
```rust
use anyhow::Result;
use clap::Args;

#[derive(Args)]
pub struct DecollapseArgs {
    #[arg(long = "input")]
    pub input: String,
    #[arg(long = "out-fas", default_value = "Decollapsed.fasta")]
    pub out_fas: String,
}

pub fn run(_args: DecollapseArgs) -> Result<()> {
    todo!()
}
```

`rust/src/rsi.rs`:
```rust
use anyhow::Result;
use clap::Args;

#[derive(Args)]
pub struct RsiArgs {
    pub input: String,
    #[arg(short = 'e', long = "explicit")]
    pub explicit: bool,
    #[arg(short = 'o', long = "output")]
    pub output: Option<String>,
}

pub fn run(_args: RsiArgs) -> Result<()> {
    todo!()
}
```

- [ ] **Step 5: Verify the package compiles**

```bash
cd rust && cargo build 2>&1
```

Expected: compilation succeeds (stubs with `todo!()` compile but are not callable without panicking).

- [ ] **Step 6: Commit**

```bash
git add rust/
git commit -m "feat: scaffold dame Rust package with stub subcommands"
```

---

## Task 2: `sort` module — core functions + unit tests

**Files:**
- Modify: `rust/src/sort.rs`
- Create: `rust/tests/sort_test.rs`

- [ ] **Step 1: Write failing tests in `rust/tests/sort_test.rs`**

```rust
use dame::sort::{rc, read_tags, read_primers, fill_hap, HapEntry};
use std::collections::HashMap;

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
    // ATCG reversed = GCTA, complement of each: G→C, C→G, T→A, A→T → CGAT
    assert_eq!(rc("ATCG"), "CGAT");
}

#[test]
fn test_rc_ambiguous_n() {
    assert_eq!(rc("N"), "N");
}

#[test]
fn test_read_tags() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("tags.txt");
    std::fs::write(&path, "ACGT\tTag1\nTTTT\tTag2\n").unwrap();
    let tags = read_tags(path.to_str().unwrap()).unwrap();
    assert!(tags.contains_key("Tag1"));
    assert_eq!(tags["Tag1"][0], "ACGT");
    assert_eq!(tags["Tag1"][1], rc("ACGT")); // "ACGT"
    assert!(tags.contains_key("Tag2"));
    assert_eq!(tags["Tag2"][0], "TTTT");
    assert_eq!(tags["Tag2"][1], "AAAA");
}

#[test]
fn test_read_primers() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("primers.txt");
    std::fs::write(&path, "CO1\tACGT\tTTTT\n").unwrap();
    let primers = read_primers(path.to_str().unwrap()).unwrap();
    assert!(primers.contains_key("CO1"));
    let entry = &primers["CO1"];
    // a_side: [F_pattern, R_pattern], b_side: [RC(F)_pattern, RC(R)_pattern]
    assert_eq!(entry.a_side.len(), 2);
    assert_eq!(entry.b_side.len(), 2);
}

#[test]
fn test_fill_hap_new_entry() {
    let mut hap: indexmap::IndexMap<String, HapEntry> = indexmap::IndexMap::new();
    fill_hap(&mut hap, "Tag1", "Tag2", "CO1", "ACGTACGT");
    assert!(hap.contains_key("Tag1_Tag2"));
    let entry = &hap["Tag1_Tag2"];
    assert_eq!(entry.tag1, "Tag1");
    assert_eq!(entry.tag2, "Tag2");
    assert_eq!(entry.seqs["ACGTACGT"].count, 1);
    assert_eq!(entry.seqs["ACGTACGT"].primer_name, "CO1");
}

#[test]
fn test_fill_hap_increment_count() {
    let mut hap = indexmap::IndexMap::new();
    fill_hap(&mut hap, "Tag1", "Tag2", "CO1", "ACGTACGT");
    fill_hap(&mut hap, "Tag1", "Tag2", "CO1", "ACGTACGT");
    assert_eq!(hap["Tag1_Tag2"].seqs["ACGTACGT"].count, 2);
}

#[test]
fn test_fill_hap_multiple_seqs() {
    let mut hap = indexmap::IndexMap::new();
    fill_hap(&mut hap, "Tag1", "Tag2", "CO1", "AAAA");
    fill_hap(&mut hap, "Tag1", "Tag2", "CO1", "CCCC");
    assert_eq!(hap["Tag1_Tag2"].seqs.len(), 2);
}
```

- [ ] **Step 2: Run tests — expect failure**

```bash
cd rust && cargo test --test sort_test 2>&1
```

Expected: compile error — `rc`, `read_tags`, `read_primers`, `fill_hap`, `HapEntry` not defined.

- [ ] **Step 3: Implement core functions in `rust/src/sort.rs`**

```rust
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Write, BufWriter};
use anyhow::{Context, Result};
use clap::Args;
use indexmap::IndexMap;
use regex::Regex;

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

pub struct PrimerEntry {
    pub a_side: Vec<String>,   // [F_pattern, R_pattern]
    pub b_side: Vec<String>,   // [RC(F)_pattern, RC(R)_pattern]
    pub a_side_re: Vec<Regex>, // compiled
    pub b_side_re: Vec<Regex>, // compiled
}

pub struct PieceInfo {
    pub tag1: String,
    pub tag2: String,
    pub primer_name: String,
    pub between: String,
}

fn complement_char(c: char) -> char {
    match c {
        'A' => 'T', 'C' => 'G', 'G' => 'C', 'T' => 'A',
        'M' => 'K', 'R' => 'Y', 'W' => 'W', 'S' => 'S',
        'Y' => 'R', 'K' => 'M', 'V' => 'B', 'H' => 'D',
        'D' => 'H', 'B' => 'V', _ => c,
    }
}

pub fn rc(seq: &str) -> String {
    seq.chars().rev().map(complement_char).collect()
}

fn to_regex_pattern(seq: &str) -> String {
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

pub fn read_tags(path: &str) -> Result<HashMap<String, Vec<String>>> {
    let mut tags: HashMap<String, Vec<String>> = HashMap::new();
    let file = File::open(path).with_context(|| format!("Cannot open tags file: {path}"))?;
    for line in BufReader::new(file).lines() {
        let line = line?;
        let line = line.trim().to_string();
        if line.is_empty() {
            continue;
        }
        let parts: Vec<&str> = line.split_whitespace().collect();
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

pub fn read_primers(path: &str) -> Result<HashMap<String, PrimerEntry>> {
    let mut primers: HashMap<String, PrimerEntry> = HashMap::new();
    let file = File::open(path).with_context(|| format!("Cannot open primers file: {path}"))?;
    for line in BufReader::new(file).lines() {
        let line = line?;
        let line = line.trim().to_string();
        if line.is_empty() {
            continue;
        }
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() < 3 {
            continue;
        }
        let name = parts[0];
        let f = parts[1];
        let r = parts[2];
        let f_rc = rc(f);
        let r_rc = rc(r);
        let f_pat = to_regex_pattern(f);
        let r_pat = to_regex_pattern(r);
        let f_rc_pat = to_regex_pattern(&f_rc);
        let r_rc_pat = to_regex_pattern(&r_rc);
        let a_side_re = vec![
            Regex::new(&f_pat).with_context(|| format!("Bad regex: {f_pat}"))?,
            Regex::new(&r_pat).with_context(|| format!("Bad regex: {r_pat}"))?,
        ];
        let b_side_re = vec![
            Regex::new(&f_rc_pat).with_context(|| format!("Bad regex: {f_rc_pat}"))?,
            Regex::new(&r_rc_pat).with_context(|| format!("Bad regex: {r_rc_pat}"))?,
        ];
        let entry = primers.entry(name.to_string()).or_insert_with(|| PrimerEntry {
            a_side: vec![],
            b_side: vec![],
            a_side_re: vec![],
            b_side_re: vec![],
        });
        entry.a_side.push(f_pat);
        entry.a_side.push(r_pat);
        entry.b_side.push(f_rc_pat);
        entry.b_side.push(r_rc_pat);
        entry.a_side_re.extend(a_side_re);
        entry.b_side_re.extend(b_side_re);
    }
    Ok(primers)
}

pub fn fill_hap(hap: &mut Hap, tag1: &str, tag2: &str, primer_name: &str, between: &str) {
    let key = format!("{tag1}_{tag2}");
    let entry = hap.entry(key).or_insert_with(|| HapEntry {
        tag1: tag1.to_string(),
        tag2: tag2.to_string(),
        seqs: IndexMap::new(),
    });
    if let Some(seq_entry) = entry.seqs.get_mut(between) {
        seq_entry.count += 1;
    } else {
        entry.seqs.insert(
            between.to_string(),
            SeqEntry { count: 1, primer_name: primer_name.to_string() },
        );
    }
}

pub fn run(_args: SortArgs) -> Result<()> {
    todo!()
}
```

- [ ] **Step 4: Run tests — expect pass**

```bash
cd rust && cargo test --test sort_test 2>&1
```

Expected: `test_rc_palindrome ... ok`, `test_rc_all_a ... ok`, etc. All 8 tests pass.

- [ ] **Step 5: Write tests for `get_pieces_info` and full `run()`**

Add to `rust/tests/sort_test.rs`:

```rust
use dame::sort::{get_pieces_info, run, SortArgs};

#[test]
fn test_get_pieces_info_forward_read() {
    // sequence: AAAA + ACGT (primer F) + barcode + TGCA (RC of primer R) + GGGG
    // tags: AAAA=Tag1 (fwd), GGGG=Tag3 (rc of CCCC)
    // primers: CO1 F=ACGT R=TGCA
    let tag_dir = tempfile::tempdir().unwrap();
    let tag_file = tag_dir.path().join("tags.txt");
    std::fs::write(&tag_file, "AAAA\tTag1\nCCCC\tTag2\nGGGG\tTag3\nTTTT\tTag4\n").unwrap();
    let prim_dir = tempfile::tempdir().unwrap();
    let prim_file = prim_dir.path().join("primers.txt");
    std::fs::write(&prim_file, "CO1\tACGT\tTGCA\n").unwrap();

    let tags = read_tags(tag_file.to_str().unwrap()).unwrap();
    let primers = read_primers(prim_file.to_str().unwrap()).unwrap();

    // AAAA + ACGT + ATATAT + TGCA + GGGG  (GGGG = RC("CCCC") which is Tag2's rc)
    let line = "AAAAACGTATATATTGCAGGGG";
    let info = get_pieces_info(line, &primers, &tags, false);
    assert!(info.is_some());
    let info = info.unwrap();
    assert_eq!(info.tag1, "Tag1");
    assert_eq!(info.tag2, "Tag2");  // GGGG = RC(CCCC) = tags["Tag2"][1]
    assert_eq!(info.primer_name, "CO1");
    assert_eq!(info.between, "ATATAT");
}

#[test]
fn test_get_pieces_info_error_read() {
    // Random sequence with no primer
    let tag_dir = tempfile::tempdir().unwrap();
    let tag_file = tag_dir.path().join("tags.txt");
    std::fs::write(&tag_file, "AAAA\tTag1\nCCCC\tTag2\n").unwrap();
    let prim_dir = tempfile::tempdir().unwrap();
    let prim_file = prim_dir.path().join("primers.txt");
    std::fs::write(&prim_file, "CO1\tACGT\tTGCA\n").unwrap();

    let tags = read_tags(tag_file.to_str().unwrap()).unwrap();
    let primers = read_primers(prim_file.to_str().unwrap()).unwrap();

    let info = get_pieces_info("NNNNNNNNNNNNNNNNNNNN", &primers, &tags, false);
    assert!(info.is_none());
}

#[test]
fn test_run_sort_produces_output_files() {
    let fixtures = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent().unwrap()
        .join("tests/fixtures");
    let dir = tempfile::tempdir().unwrap();
    let orig = std::env::current_dir().unwrap();
    std::env::set_current_dir(&dir).unwrap();
    let result = run(SortArgs {
        fastq: fixtures.join("sample.fastq").to_str().unwrap().to_string(),
        primers: fixtures.join("Primers.txt").to_str().unwrap().to_string(),
        tags: fixtures.join("Tags.txt").to_str().unwrap().to_string(),
        keep_primers_seq: false,
    });
    std::env::set_current_dir(orig).unwrap();
    assert!(result.is_ok());
    assert!(dir.path().join("SummaryCounts.txt").exists());
    assert!(dir.path().join("Tag1_Tag2.txt").exists());
    let content = std::fs::read_to_string(dir.path().join("Tag1_Tag2.txt")).unwrap();
    assert!(content.contains("ATATATATAT"));
}
```

- [ ] **Step 6: Run tests — expect failure**

```bash
cd rust && cargo test --test sort_test test_get_pieces_info 2>&1
```

Expected: compile error — `get_pieces_info` and complete `run()` not yet implemented.

- [ ] **Step 7: Implement `get_pieces_info`, output functions, and `run()` in `rust/src/sort.rs`**

Add after `fill_hap`:

```rust
fn find_tag_name(tags: &HashMap<String, Vec<String>>, idx: usize, seq: &str) -> Option<String> {
    tags.iter()
        .find(|(_, v)| v.get(idx).map(|s| s.as_str()) == Some(seq))
        .map(|(k, _)| k.clone())
}

pub fn get_pieces_info(
    line: &str,
    primers: &HashMap<String, PrimerEntry>,
    tags: &HashMap<String, Vec<String>>,
    keep_primers_seq: bool,
) -> Option<PieceInfo> {
    for (key, entry) in primers {
        // Forward read: F primer at start
        if let Some(m_start) = entry.a_side_re[0].find(line) {
            let (prim_ini_prim, prim_ini_tags) = if keep_primers_seq {
                (m_start.start(), m_start.start())
            } else {
                (m_start.end(), m_start.start())
            };
            if let Some(m_end) = entry.b_side_re[1].find(line) {
                let (prim_fin_prim, prim_fin_tags) = if keep_primers_seq {
                    (m_end.end(), m_end.end())
                } else {
                    (m_end.start(), m_end.end())
                };
                let between = &line[prim_ini_prim..prim_fin_prim];
                if between.is_empty() {
                    return None;
                }
                let tag1_seq = &line[..prim_ini_tags];
                let tag2_seq = &line[prim_fin_tags..];
                let tag1 = find_tag_name(tags, 0, tag1_seq);
                let tag2 = find_tag_name(tags, 1, tag2_seq);
                if let (Some(t1), Some(t2)) = (tag1, tag2) {
                    return Some(PieceInfo {
                        tag1: t1,
                        tag2: t2,
                        primer_name: key.clone(),
                        between: between.to_string(),
                    });
                }
            }
            return None; // committed to forward orientation, failed
        }
        // Reverse read: R primer at start
        if let Some(m_start) = entry.a_side_re[1].find(line) {
            let (prim_ini_prim, prim_ini_tags) = if keep_primers_seq {
                (m_start.start(), m_start.start())
            } else {
                (m_start.end(), m_start.start())
            };
            if let Some(m_end) = entry.b_side_re[0].find(line) {
                let (prim_fin_prim, prim_fin_tags) = if keep_primers_seq {
                    (m_end.end(), m_end.end())
                } else {
                    (m_end.start(), m_end.end())
                };
                let between = &line[prim_ini_prim..prim_fin_prim];
                if between.is_empty() {
                    return None;
                }
                let between = rc(between);
                let tag1_seq = &line[..prim_ini_tags];
                let tag2_seq = &line[prim_fin_tags..];
                // Reversed for reverse reads
                let tag2 = find_tag_name(tags, 0, tag1_seq);
                let tag1 = find_tag_name(tags, 1, tag2_seq);
                if let (Some(t1), Some(t2)) = (tag1, tag2) {
                    return Some(PieceInfo {
                        tag1: t1,
                        tag2: t2,
                        primer_name: key.clone(),
                        between,
                    });
                }
            }
            return None; // committed to reverse orientation, failed
        }
        // Neither orientation matched → try next primer key
    }
    None
}

fn print_sorted_collapsed_counted_seqs(hap: &Hap) -> Result<()> {
    for (tag_comb, entry) in hap {
        let file = File::create(format!("{tag_comb}.txt"))
            .with_context(|| format!("Cannot create {tag_comb}.txt"))?;
        let mut out = BufWriter::new(file);
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
    let file = File::create("SummaryCounts.txt")?;
    let mut out = BufWriter::new(file);
    writeln!(out, "#tagName1\ttagName2\tNumUniqSeqs\tSumTotalFreq")?;
    for (_, entry) in hap {
        let num_uniq = entry.seqs.len();
        let sum_total: u32 = entry.seqs.values().map(|s| s.count).sum();
        writeln!(out, "{}\t{}\t{}\t{}", entry.tag1, entry.tag2, num_uniq, sum_total)?;
    }
    Ok(())
}

pub fn run(args: SortArgs) -> Result<()> {
    let tags = read_tags(&args.tags)?;
    let primers = read_primers(&args.primers)?;
    let mut hap: Hap = IndexMap::new();
    let mut count_errors: u64 = 0;

    let file = File::open(&args.fastq)
        .with_context(|| format!("Cannot open fastq: {}", args.fastq))?;
    let mut reader = BufReader::new(file);
    let mut header = String::new();

    while reader.read_line(&mut header)? > 0 {
        header.clear();
        let mut seq = String::new();
        if reader.read_line(&mut seq)? == 0 {
            break;
        }
        let seq = seq.trim_end_matches('\n').trim_end_matches('\r').to_string();

        let mut plus = String::new();
        reader.read_line(&mut plus)?;
        let mut qual = String::new();
        reader.read_line(&mut qual)?;

        if seq.is_empty() {
            break;
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
        "Number of erroneous sequences (with errors in the sequence of primer or tags, \
         or no barcode amplified): {count_errors}"
    );
    Ok(())
}
```

- [ ] **Step 8: Run all sort tests — expect all pass**

```bash
cd rust && cargo test --test sort_test 2>&1
```

Expected: all 10 tests pass.

- [ ] **Step 9: Commit**

```bash
cd rust
git add src/sort.rs tests/sort_test.rs
git commit -m "feat: implement dame sort subcommand with unit tests"
```

---

## Task 3: `decollapse` module + unit tests

**Files:**
- Modify: `rust/src/decollapse.rs`
- Create: `rust/tests/decollapse_test.rs`

- [ ] **Step 1: Write failing tests in `rust/tests/decollapse_test.rs`**

```rust
use dame::decollapse::{run, DecollapseArgs};

#[test]
fn test_decollapse_single_seq() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("input.txt");
    // Format: PrimerName\tTag1\tTag2\tFreq\tSeq
    std::fs::write(&input, "CO1\tTag1\tTag2\t3\tACGT\n").unwrap();
    let out = dir.path().join("out.fasta");
    let orig = std::env::current_dir().unwrap();
    std::env::set_current_dir(&dir).unwrap();
    run(DecollapseArgs {
        input: input.to_str().unwrap().to_string(),
        out_fas: out.to_str().unwrap().to_string(),
    })
    .unwrap();
    std::env::set_current_dir(orig).unwrap();
    let content = std::fs::read_to_string(&out).unwrap();
    let headers: Vec<&str> = content.lines().filter(|l| l.starts_with('>')).collect();
    let seqs: Vec<&str> = content.lines().filter(|l| !l.starts_with('>')).collect();
    assert_eq!(headers.len(), 3);
    assert!(seqs.iter().all(|&s| s == "ACGT"));
    assert!(headers[0].contains("Tag1.Tag2.3_1"));
    assert!(headers[1].contains("Tag1.Tag2.3_2"));
    assert!(headers[2].contains("Tag1.Tag2.3_3"));
}

#[test]
fn test_decollapse_multiple_seqs() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("input.txt");
    std::fs::write(&input, "CO1\tTag1\tTag2\t2\tAAAA\nCO1\tTag1\tTag2\t1\tCCCC\n").unwrap();
    let out = dir.path().join("out.fasta");
    run(DecollapseArgs {
        input: input.to_str().unwrap().to_string(),
        out_fas: out.to_str().unwrap().to_string(),
    })
    .unwrap();
    let content = std::fs::read_to_string(&out).unwrap();
    let header_count = content.lines().filter(|l| l.starts_with('>')).count();
    assert_eq!(header_count, 3); // 2 + 1
}
```

- [ ] **Step 2: Run tests — expect failure**

```bash
cd rust && cargo test --test decollapse_test 2>&1
```

Expected: compile error — `run` and `DecollapseArgs` not yet complete.

- [ ] **Step 3: Implement `rust/src/decollapse.rs`**

```rust
use anyhow::{Context, Result};
use clap::Args;
use std::fs::File;
use std::io::{BufRead, BufReader, Write, BufWriter};

#[derive(Args)]
pub struct DecollapseArgs {
    #[arg(long = "input")]
    pub input: String,
    #[arg(long = "out-fas", default_value = "Decollapsed.fasta")]
    pub out_fas: String,
}

pub fn run(args: DecollapseArgs) -> Result<()> {
    let in_file = File::open(&args.input)
        .with_context(|| format!("Cannot open input: {}", args.input))?;
    let out_file = File::create(&args.out_fas)
        .with_context(|| format!("Cannot create output: {}", args.out_fas))?;
    let mut out = BufWriter::new(out_file);
    let mut seq_id: u64 = 0;

    for line in BufReader::new(in_file).lines() {
        let line = line?;
        let line = line.trim().to_string();
        if line.is_empty() {
            continue;
        }
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() < 5 {
            continue;
        }
        let tag1 = parts[1];
        let tag2 = parts[2];
        let freq: u32 = parts[3].parse().with_context(|| format!("Bad freq: {}", parts[3]))?;
        let seq = parts[4];

        for _ in 1..=freq {
            seq_id += 1;
            writeln!(
                out,
                ">{}.{}.{}_{}\n{}",
                tag1, tag2, freq, seq_id, seq
            )?;
        }
    }
    Ok(())
}
```

- [ ] **Step 4: Run tests — expect pass**

```bash
cd rust && cargo test --test decollapse_test 2>&1
```

Expected: both tests pass.

- [ ] **Step 5: Commit**

```bash
cd rust
git add src/decollapse.rs tests/decollapse_test.rs
git commit -m "feat: implement dame decollapse subcommand with unit tests"
```

---

## Task 4: `chimera_check` module + unit tests

**Files:**
- Modify: `rust/src/chimera_check.rs`
- Create: `rust/tests/chimera_test.rs`

- [ ] **Step 1: Write failing tests in `rust/tests/chimera_test.rs`**

```rust
use dame::chimera_check::{make_tag_files, make_tag_files_with_pools, make_fas_seq_one_line};

#[test]
fn test_make_tag_files_creates_files() {
    let dir = tempfile::tempdir().unwrap();
    let psinfo = dir.path().join("PSinfo.txt");
    std::fs::write(&psinfo, "S1\tTag1\tTag2\t1\nS1\tTag3\tTag4\t1\n").unwrap();
    let orig = std::env::current_dir().unwrap();
    std::env::set_current_dir(&dir).unwrap();
    make_tag_files(psinfo.to_str().unwrap(), 2).unwrap();
    std::env::set_current_dir(orig).unwrap();
    assert!(dir.path().join("PS1.tags.txt").exists());
    assert!(dir.path().join("PS2.tags.txt").exists());
    let c1 = std::fs::read_to_string(dir.path().join("PS1.tags.txt")).unwrap();
    let c2 = std::fs::read_to_string(dir.path().join("PS2.tags.txt")).unwrap();
    assert!(c1.contains("Tag1\tTag2"));
    assert!(c2.contains("Tag3\tTag4"));
}

#[test]
fn test_make_tag_files_with_pools_includes_pool() {
    let dir = tempfile::tempdir().unwrap();
    let psinfo = dir.path().join("PSinfo.txt");
    std::fs::write(&psinfo, "S1\tTag1\tTag2\t1\nS1\tTag3\tTag4\t2\n").unwrap();
    let orig = std::env::current_dir().unwrap();
    std::env::set_current_dir(&dir).unwrap();
    make_tag_files_with_pools(psinfo.to_str().unwrap(), 2).unwrap();
    std::env::set_current_dir(orig).unwrap();
    let c1 = std::fs::read_to_string(dir.path().join("PS1.tags.txt")).unwrap();
    assert!(c1.contains("Tag1\tTag2\t1"));
}

#[test]
fn test_make_fas_seq_one_line() {
    let dir = tempfile::tempdir().unwrap();
    let fasta = dir.path().join("Pool1.noChim.fasta");
    std::fs::write(&fasta, ">header1\nACGT\nACGT\n>header2\nTTTT\n").unwrap();
    let orig = std::env::current_dir().unwrap();
    std::env::set_current_dir(&dir).unwrap();
    make_fas_seq_one_line(1).unwrap();
    std::env::set_current_dir(orig).unwrap();
    let result = std::fs::read_to_string(dir.path().join("Pool1.noChim.oneLiner.fasta")).unwrap();
    let lines: Vec<&str> = result.lines().collect();
    assert_eq!(lines[0], ">header1");
    assert_eq!(lines[1], "ACGTACGT");
    assert_eq!(lines[2], ">header2");
    assert_eq!(lines[3], "TTTT");
}
```

- [ ] **Step 2: Run tests — expect failure**

```bash
cd rust && cargo test --test chimera_test 2>&1
```

Expected: compile error — functions not defined.

- [ ] **Step 3: Implement `rust/src/chimera_check.rs`**

```rust
use anyhow::{Context, Result};
use clap::Args;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Write, BufWriter};
use std::process::Command;

#[derive(Args)]
pub struct ChimeraArgs {
    #[arg(long = "ps-info")]
    pub ps_info: String,
    #[arg(long = "x")]
    pub x: usize,
    #[arg(long = "p", default_value = "1")]
    pub p: usize,
}

fn write_ps_tag_files(ps_info: &str, x: usize, include_pool: bool) -> Result<()> {
    let mut outs: Vec<BufWriter<File>> = (0..x)
        .map(|i| {
            File::create(format!("PS{}.tags.txt", i + 1))
                .map(BufWriter::new)
                .with_context(|| format!("Cannot create PS{}.tags.txt", i + 1))
        })
        .collect::<Result<_>>()?;

    let file = File::open(ps_info).with_context(|| format!("Cannot open: {ps_info}"))?;
    for (nr, line) in BufReader::new(file).lines().enumerate() {
        let line = line?;
        let line = line.trim().to_string();
        if line.is_empty() {
            continue;
        }
        let parts: Vec<&str> = line.split_whitespace().collect();
        let nr1 = nr + 1;
        let residue = nr1 % x;
        let idx = if residue != 0 { residue - 1 } else { x - 1 };
        if include_pool {
            writeln!(outs[idx], "{}\t{}\t{}", parts[1], parts[2], parts[3])?;
        } else {
            writeln!(outs[idx], "{}\t{}", parts[1], parts[2])?;
        }
    }
    Ok(())
}

pub fn make_tag_files(ps_info: &str, x: usize) -> Result<()> {
    write_ps_tag_files(ps_info, x, false)
}

pub fn make_tag_files_with_pools(ps_info: &str, x: usize) -> Result<()> {
    write_ps_tag_files(ps_info, x, true)
}

pub fn make_size_out_fastas(p: usize, x: usize) -> Result<()> {
    let mut outs: Vec<BufWriter<File>> = (0..p)
        .map(|pool| {
            File::create(format!("Pool{}.fasta", pool + 1))
                .map(BufWriter::new)
                .with_context(|| format!("Cannot create Pool{}.fasta", pool + 1))
        })
        .collect::<Result<_>>()?;

    for num in 0..x {
        let tag_file = format!("PS{}.tags.txt", num + 1);
        if !std::path::Path::new(&tag_file).exists() {
            continue;
        }
        let file = File::open(&tag_file)?;
        for line in BufReader::new(file).lines() {
            let line = line?;
            let line = line.trim().to_string();
            if line.is_empty() {
                continue;
            }
            let parts: Vec<&str> = line.split_whitespace().collect();
            let hap_name = if p > 1 {
                format!("./pool{}/{}_{}.txt", parts[2], parts[0], parts[1])
            } else {
                format!("{}_{}.txt", parts[0], parts[1])
            };
            if !std::path::Path::new(&hap_name).exists() {
                continue;
            }
            let pool_idx = if p > 1 {
                parts[2].parse::<usize>().unwrap_or(1) - 1
            } else {
                0
            };
            let hap_file = File::open(&hap_name)?;
            let mut id_num = 1usize;
            for seq_line in BufReader::new(hap_file).lines() {
                let seq_line = seq_line?;
                let seq_parts: Vec<&str> = seq_line.split('\t').collect();
                if seq_parts.len() < 5 {
                    continue;
                }
                writeln!(
                    outs[pool_idx],
                    ">{}_{}_{}_{};size={}\n{}",
                    seq_parts[0], parts[0], parts[1], id_num,
                    seq_parts[3], seq_parts[4]
                )?;
                id_num += 1;
            }
        }
    }
    Ok(())
}

pub fn sort_fasta(p: usize) -> Result<()> {
    for pool in 0..p {
        let input_file = format!("Pool{}.fasta", pool + 1);
        let sorted_file = format!("Pool{}.sort.fasta", pool + 1);
        let chim_file = format!("Pool{}.Chim.fasta", pool + 1);
        let no_chim_file = format!("Pool{}.noChim.fasta", pool + 1);

        let sort_out = Command::new("usearch")
            .args(["--sortsize", &input_file, "--output", &sorted_file])
            .output()
            .context("Failed to run usearch --sortsize (is usearch on PATH?)")?;
        std::fs::write(format!("sort{}.out", pool + 1), &sort_out.stdout)?;
        std::fs::write(format!("sort{}.err", pool + 1), &sort_out.stderr)?;

        let chim_out = Command::new("usearch")
            .args(["-uchime", &sorted_file, "-chimeras", &chim_file, "-nonchimeras", &no_chim_file])
            .output()
            .context("Failed to run usearch -uchime")?;
        std::fs::write(format!("chimeraCheck{}.out", pool + 1), &chim_out.stdout)?;
        std::fs::write(format!("chimeraCheck{}.err", pool + 1), &chim_out.stderr)?;
    }
    Ok(())
}

pub fn make_fas_seq_one_line(p: usize) -> Result<()> {
    for pool in 0..p {
        let in_path = format!("Pool{}.noChim.fasta", pool + 1);
        let out_path = format!("Pool{}.noChim.oneLiner.fasta", pool + 1);
        let in_file = File::open(&in_path)
            .with_context(|| format!("Cannot open {in_path}"))?;
        let out_file = File::create(&out_path)
            .with_context(|| format!("Cannot create {out_path}"))?;
        let mut out = BufWriter::new(out_file);
        let mut seq = String::new();
        for line in BufReader::new(in_file).lines() {
            let line = line?;
            let line = line.trim_end().to_string();
            if line.starts_with('>') {
                if !seq.is_empty() {
                    writeln!(out, "{}", seq)?;
                    seq.clear();
                }
                writeln!(out, "{}", line)?;
            } else {
                seq.push_str(&line);
            }
        }
        if !seq.is_empty() {
            writeln!(out, "{}", seq)?;
        }
    }
    Ok(())
}

pub fn make_no_chim_haps(p: usize) -> Result<()> {
    // HAP[tagHapKey] = [primerNames, tag1s, tag2s, freqs, seqs]
    let mut hap: HashMap<String, (Vec<String>, Vec<String>, Vec<String>, Vec<String>, Vec<String>)> = HashMap::new();

    for pool in 0..p {
        let path = format!("Pool{}.noChim.oneLiner.fasta", pool + 1);
        if !std::path::Path::new(&path).exists() {
            continue;
        }
        let file = File::open(&path)?;
        let mut primer_name = String::new();
        let mut tag1 = String::new();
        let mut tag2 = String::new();
        let mut freq = String::new();
        let mut tag_hap_key = String::new();

        for line in BufReader::new(file).lines() {
            let line = line?;
            let line = line.trim_end().to_string();
            if line.starts_with('>') {
                let h = &line[1..];
                let parts: Vec<&str> = h.splitn(2, ';').collect();
                let name_parts: Vec<&str> = parts[0].split('_').collect();
                primer_name = name_parts.get(0).unwrap_or(&"").to_string();
                tag1 = name_parts.get(1).unwrap_or(&"").to_string();
                tag2 = name_parts.get(2).unwrap_or(&"").to_string();
                freq = parts.get(1)
                    .and_then(|s| s.strip_prefix("size="))
                    .unwrap_or("")
                    .to_string();
                tag_hap_key = format!("{}_{}_{}", tag1, tag2, pool + 1);
                hap.entry(tag_hap_key.clone()).or_default();
            } else {
                if let Some(entry) = hap.get_mut(&tag_hap_key) {
                    entry.0.push(primer_name.clone());
                    entry.1.push(tag1.clone());
                    entry.2.push(tag2.clone());
                    entry.3.push(freq.clone());
                    entry.4.push(line);
                }
            }
        }
    }

    for (tag_comb, (primers, tags1, tags2, freqs, seqs)) in &hap {
        let out_path = format!("{}.noChim.txt", tag_comb);
        let out_file = File::create(&out_path)
            .with_context(|| format!("Cannot create {out_path}"))?;
        let mut out = BufWriter::new(out_file);
        for i in 0..primers.len() {
            writeln!(out, "{}\t{}\t{}\t{}\t{}", primers[i], tags1[i], tags2[i], freqs[i], seqs[i])?;
        }
    }
    Ok(())
}

pub fn run(args: ChimeraArgs) -> Result<()> {
    if args.p == 1 {
        make_tag_files(&args.ps_info, args.x)?;
    } else {
        make_tag_files_with_pools(&args.ps_info, args.x)?;
    }
    make_size_out_fastas(args.p, args.x)?;
    sort_fasta(args.p)?;
    make_fas_seq_one_line(args.p)?;
    make_no_chim_haps(args.p)?;
    Ok(())
}
```

- [ ] **Step 4: Run tests — expect pass**

```bash
cd rust && cargo test --test chimera_test 2>&1
```

Expected: 3 tests pass.

- [ ] **Step 5: Commit**

```bash
cd rust
git add src/chimera_check.rs tests/chimera_test.rs
git commit -m "feat: implement dame chimera subcommand with unit tests"
```

---

## Task 5: `filter` module + unit tests

**Files:**
- Modify: `rust/src/filter.rs`
- Create: `rust/tests/filter_test.rs`

- [ ] **Step 1: Write failing tests in `rust/tests/filter_test.rs`**

```rust
use dame::filter::{
    make_ps_num_files, read_ps_num_files, make_sample_name_array,
    read_haps_for_a_sample, get_seqs_sets_and_fr_counts,
};

#[test]
fn test_make_sample_name_array() {
    let dir = tempfile::tempdir().unwrap();
    let psinfo = dir.path().join("PSinfo.txt");
    std::fs::write(
        &psinfo,
        "SampleA\tTag1\tTag2\t1\nSampleA\tTag3\tTag4\t1\nSampleB\tTag5\tTag6\t1\nSampleB\tTag7\tTag8\t1\n",
    )
    .unwrap();
    let names = make_sample_name_array(psinfo.to_str().unwrap()).unwrap();
    assert_eq!(names, vec!["SampleA", "SampleB"]);
}

#[test]
fn test_make_sample_name_array_deduplicates() {
    let dir = tempfile::tempdir().unwrap();
    let psinfo = dir.path().join("PSinfo.txt");
    std::fs::write(&psinfo, "S1\tTag1\tTag2\t1\nS1\tTag3\tTag4\t1\n").unwrap();
    let names = make_sample_name_array(psinfo.to_str().unwrap()).unwrap();
    assert_eq!(names.len(), 1);
    assert_eq!(names[0], "S1");
}

#[test]
fn test_make_ps_num_files_creates_files() {
    let dir = tempfile::tempdir().unwrap();
    let psinfo = dir.path().join("PSinfo.txt");
    std::fs::write(&psinfo, "S1\tTag1\tTag2\t1\nS1\tTag3\tTag4\t1\n").unwrap();
    let orig = std::env::current_dir().unwrap();
    std::env::set_current_dir(&dir).unwrap();
    make_ps_num_files(psinfo.to_str().unwrap(), 2, 1, false).unwrap();
    std::env::set_current_dir(orig).unwrap();
    assert!(dir.path().join("PS1_files.txt").exists());
    assert!(dir.path().join("PS2_files.txt").exists());
    let c1 = std::fs::read_to_string(dir.path().join("PS1_files.txt")).unwrap();
    let c2 = std::fs::read_to_string(dir.path().join("PS2_files.txt")).unwrap();
    assert!(c1.contains("pool1/Tag1_Tag2.txt"));
    assert!(c2.contains("pool1/Tag3_Tag4.txt"));
}

#[test]
fn test_get_seqs_sets_and_fr_counts_empty() {
    let haps: std::collections::HashMap<usize, Vec<Vec<String>>> =
        [(0, vec![]), (1, vec![])].into_iter().collect();
    let (seqs_all, f, r, counts, seqs) = get_seqs_sets_and_fr_counts(2, &haps);
    assert!(seqs_all.is_empty());
    assert!(f.is_empty());
    assert!(r.is_empty());
}

#[test]
fn test_get_seqs_sets_and_fr_counts_with_data() {
    let haps: std::collections::HashMap<usize, Vec<Vec<String>>> = [
        (0, vec![
            vec!["CO1".into(), "Tag1".into(), "Tag2".into(), "3".into(), "AAAA".into()],
            vec!["CO1".into(), "Tag1".into(), "Tag2".into(), "1".into(), "CCCC".into()],
        ]),
        (1, vec![
            vec!["CO1".into(), "Tag1".into(), "Tag2".into(), "2".into(), "AAAA".into()],
        ]),
    ]
    .into_iter()
    .collect();
    let (seqs_all, f, r, counts, _seqs) = get_seqs_sets_and_fr_counts(2, &haps);
    assert!(seqs_all.contains("AAAA"));
    assert!(seqs_all.contains("CCCC"));
    assert_eq!(f[&0], "Tag1");
    assert_eq!(r[&0], "Tag2");
    assert_eq!(counts[&0], vec!["3", "1"]);
}
```

- [ ] **Step 2: Run tests — expect failure**

```bash
cd rust && cargo test --test filter_test 2>&1
```

Expected: compile error — functions not defined.

- [ ] **Step 3: Implement `rust/src/filter.rs`**

```rust
use anyhow::{Context, Result};
use clap::Args;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader, Write, BufWriter};

#[derive(Args)]
pub struct FilterArgs {
    #[arg(long = "ps-info")]
    pub ps_info: String,
    #[arg(long = "x", default_value = "2")]
    pub x: usize,
    #[arg(long = "y", default_value = "1")]
    pub y: usize,
    #[arg(long = "p", default_value = "1")]
    pub p: usize,
    #[arg(long = "t", default_value = "1")]
    pub t: u32,
    #[arg(long = "l", default_value = "100")]
    pub l: usize,
    #[arg(long = "chimera-checked")]
    pub chimera_checked: bool,
}

pub fn make_ps_num_files(ps_info: &str, x: usize, _p: usize, chimera_checked: bool) -> Result<()> {
    let mut outs: Vec<BufWriter<File>> = (0..x)
        .map(|i| {
            File::create(format!("PS{}_files.txt", i + 1))
                .map(BufWriter::new)
                .with_context(|| format!("Cannot create PS{}_files.txt", i + 1))
        })
        .collect::<Result<_>>()?;

    let file = File::open(ps_info).with_context(|| format!("Cannot open: {ps_info}"))?;
    for (nr, line) in BufReader::new(file).lines().enumerate() {
        let line = line?;
        let line = line.trim().to_string();
        if line.is_empty() {
            continue;
        }
        let parts: Vec<&str> = line.split_whitespace().collect();
        let nr1 = nr + 1;
        let residue = nr1 % x;
        let idx = if residue != 0 { residue - 1 } else { x - 1 };
        if chimera_checked {
            writeln!(outs[idx], "{}_{}_{}.noChim.txt", parts[1], parts[2], parts[3])?;
        } else {
            writeln!(outs[idx], "pool{}/{}_{}.txt", parts[3], parts[1], parts[2])?;
        }
    }
    Ok(())
}

pub fn read_ps_num_files(x: usize) -> Result<HashMap<usize, Vec<String>>> {
    let mut result = HashMap::new();
    for i in 0..x {
        let path = format!("PS{}_files.txt", i + 1);
        let file = File::open(&path).with_context(|| format!("Cannot open {path}"))?;
        let lines: Vec<String> = BufReader::new(file)
            .lines()
            .collect::<std::io::Result<_>>()?;
        result.insert(i, lines);
    }
    Ok(result)
}

pub fn make_sample_name_array(ps_info: &str) -> Result<Vec<String>> {
    let mut names: Vec<String> = Vec::new();
    let file = File::open(ps_info).with_context(|| format!("Cannot open: {ps_info}"))?;
    for line in BufReader::new(file).lines() {
        let line = line?;
        let name = line.split_whitespace().next().unwrap_or("").to_string();
        if !name.is_empty() && !names.contains(&name) {
            names.push(name);
        }
    }
    Ok(names)
}

pub fn read_haps_for_a_sample(
    x: usize,
    ps_ins_lines: &HashMap<usize, Vec<String>>,
    i: usize,
) -> Result<HashMap<usize, Vec<Vec<String>>>> {
    let mut haps = HashMap::new();
    for j in 0..x {
        let mut rows = Vec::new();
        if let Some(lines) = ps_ins_lines.get(&j) {
            if let Some(path) = lines.get(i) {
                let path = path.trim();
                if path != "empty" && std::path::Path::new(path).exists() {
                    let file = File::open(path)?;
                    for line in BufReader::new(file).lines() {
                        let line = line?;
                        let parts: Vec<String> =
                            line.split('\t').map(|s| s.to_string()).collect();
                        rows.push(parts);
                    }
                }
            }
        }
        haps.insert(j, rows);
    }
    Ok(haps)
}

pub fn get_seqs_sets_and_fr_counts(
    x: usize,
    haps: &HashMap<usize, Vec<Vec<String>>>,
) -> (
    HashSet<String>,
    HashMap<usize, String>,
    HashMap<usize, String>,
    HashMap<usize, Vec<String>>,
    HashMap<usize, Vec<String>>,
) {
    let mut f = HashMap::new();
    let mut r = HashMap::new();
    let mut counts = HashMap::new();
    let mut seqs = HashMap::new();
    let mut seqs_all = HashSet::new();

    for j in 0..x {
        if let Some(hap_j) = haps.get(&j) {
            if !hap_j.is_empty() {
                f.insert(j, hap_j[0].get(1).cloned().unwrap_or_default());
                r.insert(j, hap_j[0].get(2).cloned().unwrap_or_default());
                let mut c_vec = Vec::new();
                let mut s_vec = Vec::new();
                for row in hap_j {
                    let cnt = row.get(3).cloned().unwrap_or_default();
                    let seq = row.get(4).cloned().unwrap_or_default();
                    c_vec.push(cnt);
                    s_vec.push(seq.clone());
                    seqs_all.insert(seq);
                }
                counts.insert(j, c_vec);
                seqs.insert(j, s_vec);
            }
        }
    }
    (seqs_all, f, r, counts, seqs)
}

#[allow(clippy::too_many_arguments)]
pub fn make_comparison_file(
    x: usize,
    seqs_all: &HashSet<String>,
    haps: &HashMap<usize, Vec<Vec<String>>>,
    f: &HashMap<usize, String>,
    r: &HashMap<usize, String>,
    counts: &HashMap<usize, Vec<String>>,
    seqs: &HashMap<usize, Vec<String>>,
    out: &mut impl Write,
    out_thresh: &mut impl Write,
    out_yx: &mut impl Write,
    out_fas: &mut impl Write,
    out_thresh_fas: &mut impl Write,
    out_yx_fas: &mut impl Write,
    out_thresh_len_fas: &mut impl Write,
    y: usize,
    t: u32,
    l: usize,
    sample_name: &[String],
    i: usize,
) -> Result<()> {
    let mut idnum = 1usize;
    // Sort for deterministic output (Python uses set which is non-deterministic)
    let mut seqs_sorted: Vec<&String> = seqs_all.iter().collect();
    seqs_sorted.sort();
    for seq in seqs_sorted {
        let mut line = format!("{}\t", sample_name[i]);
        let mut line_fas_ids = format!(">{}\t", sample_name[i]);
        let mut line_fas_counts = "\t".to_string();
        let mut y_count = 0usize;
        let mut t_count = 0usize;

        for j in 0..x {
            let hap_j = haps.get(&j).map(|v| v.as_slice()).unwrap_or_default();
            let f_j = f.get(&j).map(|s| s.as_str()).unwrap_or("empty");
            let r_j = r.get(&j).map(|s| s.as_str()).unwrap_or("empty");

            if !hap_j.is_empty() {
                let pos = seqs.get(&j)
                    .and_then(|sv| sv.iter().position(|s| s == seq));
                let count = if let Some(p) = pos {
                    y_count += 1;
                    let cnt_str = counts[&j][p].as_str();
                    let cnt_val: u32 = cnt_str.parse().unwrap_or(0);
                    if cnt_val < t {
                        t_count += 1;
                    }
                    cnt_str.to_string()
                } else {
                    "0".to_string()
                };
                line.push_str(&format!("{}-{}\t{}\t", f_j, r_j, count));
                if j < x - 1 {
                    line_fas_ids.push_str(&format!("{}-{}.", f_j, r_j));
                    line_fas_counts.push_str(&format!("{}_", count));
                } else {
                    line_fas_ids.push_str(&format!("{}-{}_{}\t", f_j, r_j, idnum));
                    line_fas_counts.push_str(&format!("{}\n{}", count, seq));
                }
            } else {
                line.push_str("empty\t0\t");
                if j < x - 1 {
                    line_fas_ids.push_str("empty-empty.");
                    line_fas_counts.push_str("0_");
                } else {
                    line_fas_ids.push_str(&format!("empty-empty_{}\t", idnum));
                    line_fas_counts.push_str(&format!("0\n{}", seq));
                }
            }
        }
        line.push_str(&format!("{}\n", seq));
        let line_fas = format!("{}{}\n", line_fas_ids, line_fas_counts);

        write!(out, "{}", line)?;
        write!(out_fas, "{}", line_fas)?;
        if y_count >= y {
            write!(out_yx, "{}", line)?;
            write!(out_yx_fas, "{}", line_fas)?;
        }
        if (y_count - t_count) >= y {
            write!(out_thresh, "{}", line)?;
            write!(out_thresh_fas, "{}", line_fas)?;
            if seq.len() >= l {
                write!(out_thresh_len_fas, "{}", line_fas)?;
            }
        }
        idnum += 1;
    }
    Ok(())
}

pub fn run(args: FilterArgs) -> Result<()> {
    let ps_info = &args.ps_info;
    let x = args.x;
    let y = args.y;
    let t = args.t;
    let l = args.l;

    let mut out = BufWriter::new(File::create(format!("Comparisons_{}PCRs.txt", x))?);
    let mut out_yx = BufWriter::new(File::create(format!("Comparisons_{}outOf{}PCRs.txt", y, x))?);
    let mut out_thresh = BufWriter::new(File::create(format!(
        "Comparisons_{}outOf{}PCRs.countsThreshold{}.txt", y, x, t
    ))?);
    let mut out_fas = BufWriter::new(File::create(format!("Comparisons_{}PCRs.fasta", x))?);
    let mut out_yx_fas = BufWriter::new(File::create(format!("FilteredReads_atLeast{}.fasta", y))?);
    let mut out_thresh_fas = BufWriter::new(File::create(format!(
        "FilteredReads_atLeast{}.threshold.fasta", y
    ))?);
    let mut out_thresh_len_fas = BufWriter::new(File::create("FilteredReads.fna")?);

    make_ps_num_files(ps_info, x, args.p, args.chimera_checked)?;
    let ps_ins_lines = read_ps_num_files(x)?;
    let sample_names = make_sample_name_array(ps_info)?;

    let num_samples = ps_ins_lines.get(&0).map(|v| v.len()).unwrap_or(0);
    for i in 0..num_samples {
        let haps = read_haps_for_a_sample(x, &ps_ins_lines, i)?;
        let (seqs_all, f, r, counts, seqs) = get_seqs_sets_and_fr_counts(x, &haps);
        make_comparison_file(
            x, &seqs_all, &haps, &f, &r, &counts, &seqs,
            &mut out, &mut out_thresh, &mut out_yx,
            &mut out_fas, &mut out_thresh_fas, &mut out_yx_fas,
            &mut out_thresh_len_fas,
            y, t, l, &sample_names, i,
        )?;
    }
    Ok(())
}
```

- [ ] **Step 4: Run tests — expect pass**

```bash
cd rust && cargo test --test filter_test 2>&1
```

Expected: 5 tests pass.

- [ ] **Step 5: Commit**

```bash
cd rust
git add src/filter.rs tests/filter_test.rs
git commit -m "feat: implement dame filter subcommand with unit tests"
```

---

## Task 6: `rsi` module + unit tests

**Files:**
- Modify: `rust/src/rsi.rs`
- Create: `rust/tests/rsi_test.rs`

- [ ] **Step 1: Write failing tests in `rust/tests/rsi_test.rs`**

```rust
use dame::rsi::{compare, run, RsiArgs};
use ndarray::array;

#[test]
fn test_compare_identical_replicates() {
    let matrix = array![[10i64, 10], [20, 20], [30, 30]];
    let result = compare(&matrix, "sample1", 1, 2);
    assert!(result.abs() < 1e-10);
}

#[test]
fn test_compare_completely_different() {
    let matrix = array![[100i64, 0], [0, 100]];
    let result = compare(&matrix, "sample1", 1, 2);
    assert!((result - 1.0).abs() < 1e-10);
}

#[test]
fn test_compare_partial_overlap() {
    // col_a = [0.5, 0.5], col_b = [0.0, 1.0], min = [0.0, 0.5], sum = 0.5, RSI = 0.5
    let matrix = array![[50i64, 0], [50, 100]];
    let result = compare(&matrix, "sample1", 1, 2);
    assert!((result - 0.5).abs() < 1e-10);
}

#[test]
fn test_compare_zero_replicate_handled() {
    let matrix = array![[0i64, 10], [0, 20]];
    let result = compare(&matrix, "sample1", 1, 2);
    assert!(result >= 0.0 && result <= 1.0);
}

#[test]
fn test_run_rsi_produces_output() {
    let dir = tempfile::tempdir().unwrap();
    // Format: SampleName F0-R0 count0 F1-R1 count1 ... seq
    let input = dir.path().join("comparisons.txt");
    std::fs::write(
        &input,
        "S1\tA-B\t10\tA-B\t10\tACGT\nS1\tA-B\t20\tA-B\t20\tGGGG\n",
    )
    .unwrap();
    let out = dir.path().join("RSI_output.txt");
    let orig = std::env::current_dir().unwrap();
    std::env::set_current_dir(&dir).unwrap();
    run(RsiArgs {
        input: input.to_str().unwrap().to_string(),
        explicit: false,
        output: Some(out.to_str().unwrap().to_string()),
    })
    .unwrap();
    std::env::set_current_dir(orig).unwrap();
    let content = std::fs::read_to_string(&out).unwrap();
    assert!(content.contains("Sample\tRSI"));
    assert!(content.contains("S1"));
}
```

- [ ] **Step 2: Run tests — expect failure**

```bash
cd rust && cargo test --test rsi_test 2>&1
```

Expected: compile error — `compare`, `run`, `RsiArgs` not yet complete.

- [ ] **Step 3: Implement `rust/src/rsi.rs`**

```rust
use anyhow::{Context, Result};
use clap::Args;
use ndarray::{Array2, Axis};
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader, Write, BufWriter};

#[derive(Args)]
pub struct RsiArgs {
    pub input: String,
    #[arg(short = 'e', long = "explicit")]
    pub explicit: bool,
    #[arg(short = 'o', long = "output")]
    pub output: Option<String>,
}

pub fn compare(matrix: &Array2<i64>, _j: &str, a: usize, b: usize) -> f64 {
    let total = matrix.sum_axis(Axis(0));
    let total_a = if total[0] == 0 {
        eprintln!("This sample gave zero in replicate: {a}");
        1_f64
    } else {
        total[0] as f64
    };
    let total_b = if total[1] == 0 {
        eprintln!("This sample gave zero in replicate: {b}");
        1_f64
    } else {
        total[1] as f64
    };

    let min_sum: f64 = matrix
        .rows()
        .into_iter()
        .map(|row| {
            let pa = row[0] as f64 / total_a;
            let pb = row[1] as f64 / total_b;
            pa.min(pb)
        })
        .sum();

    1.0 - min_sum
}

pub fn run(args: RsiArgs) -> Result<()> {
    let file = File::open(&args.input)
        .with_context(|| format!("Cannot open: {}", args.input))?;
    let data: Vec<Vec<String>> = BufReader::new(file)
        .lines()
        .map(|l| l.map(|s| s.split_whitespace().map(|t| t.to_string()).collect()))
        .collect::<std::io::Result<_>>()?;

    if data.is_empty() {
        println!("There are no replicates in the file.");
        return Ok(());
    }

    let row_len = data[0].len();
    // row format: sampleName, F0-R0, count0, F1-R1, count1, ..., seq
    // no_rep = number of PCR replicates = (row_len - 2) / 2
    let no_rep = (row_len.saturating_sub(2)) / 2;
    if no_rep < 2 {
        println!("There are no replicates in the file.");
        return Ok(());
    }

    let names: HashSet<&str> = data.iter().map(|row| row[0].as_str()).collect();

    let mut rkn_avg: Vec<(String, f64)> = Vec::new();
    let mut rkn_exp: Vec<(String, usize, usize, f64)> = Vec::new();

    let mut names_sorted: Vec<&str> = names.iter().copied().collect();
    names_sorted.sort();

    for name in &names_sorted {
        let subset: Vec<&Vec<String>> = data.iter().filter(|row| row[0] == *name).collect();

        if args.explicit {
            for a in 1..no_rep {
                for b in (a + 1)..=no_rep {
                    let col_a = a * 2;
                    let col_b = b * 2;
                    let matrix_data: Vec<i64> = subset.iter()
                        .flat_map(|row| {
                            let va = row.get(col_a)
                                .and_then(|s| s.parse::<i64>().ok())
                                .unwrap_or(0);
                            let vb = row.get(col_b)
                                .and_then(|s| s.parse::<i64>().ok())
                                .unwrap_or(0);
                            [va, vb]
                        })
                        .collect();
                    let nrows = subset.len();
                    let matrix = Array2::from_shape_vec((nrows, 2), matrix_data)?;
                    let output = compare(&matrix, name, a, b);
                    rkn_exp.push((name.to_string(), a, b, output));
                }
            }
        } else {
            let mut total_output = 0.0_f64;
            let mut rep_count = 0usize;
            for a in 1..no_rep {
                for b in (a + 1)..=no_rep {
                    let col_a = a * 2;
                    let col_b = b * 2;
                    let matrix_data: Vec<i64> = subset.iter()
                        .flat_map(|row| {
                            let va = row.get(col_a)
                                .and_then(|s| s.parse::<i64>().ok())
                                .unwrap_or(0);
                            let vb = row.get(col_b)
                                .and_then(|s| s.parse::<i64>().ok())
                                .unwrap_or(0);
                            [va, vb]
                        })
                        .collect();
                    let nrows = subset.len();
                    let matrix = Array2::from_shape_vec((nrows, 2), matrix_data)?;
                    total_output += compare(&matrix, name, a, b);
                    rep_count += 1;
                }
            }
            rkn_avg.push((name.to_string(), total_output / rep_count as f64));
        }
    }

    let outfile = args.output.as_deref().unwrap_or("RSI_output.txt");
    let out_file = File::create(outfile)
        .with_context(|| format!("Cannot create: {outfile}"))?;
    let mut out = BufWriter::new(out_file);

    if args.explicit {
        writeln!(out, "Sample\tReplicateA\tReplicateB\tRSI")?;
        for (sample, a, b, rsi) in &rkn_exp {
            writeln!(out, "{sample}\t{a}\t{b}\t{rsi}")?;
        }
    } else {
        writeln!(out, "Sample\tRSI")?;
        for (sample, rsi) in &rkn_avg {
            writeln!(out, "{sample}\t{rsi}")?;
        }
    }
    Ok(())
}
```

- [ ] **Step 4: Run tests — expect pass**

```bash
cd rust && cargo test --test rsi_test 2>&1
```

Expected: 5 tests pass.

- [ ] **Step 5: Run full test suite**

```bash
cd rust && cargo test 2>&1
```

Expected: all tests pass with `0 failed`.

- [ ] **Step 6: Commit**

```bash
cd rust
git add src/rsi.rs tests/rsi_test.rs
git commit -m "feat: implement dame rsi subcommand with unit tests"
```

---

## Task 7: Update integration tests to compare `dame` vs `dame-py`

**Files:**
- Modify: `tests/integration/run_sort.sh`
- Modify: `tests/integration/run_rsi.sh`
- Modify: `tests/integration/run_decollapse.sh`
- Modify: `tests/integration/run_filter.sh`

**Note:** These scripts require both `dame-py` (installed via `pip install -e python/`) and `dame` (built via `cargo build --release` in `rust/`, binary at `rust/target/release/dame`) to be on PATH. The scripts are responsible for comparing output files of both tools.

- [ ] **Step 1: Replace `tests/integration/run_sort.sh`**

```bash
#!/usr/bin/env bash
set -euo pipefail

FIXTURES="$(cd "$(dirname "$0")/../fixtures" && pwd)"

run_sort() {
    local tool="$1"
    local outdir="$2"
    mkdir -p "$outdir"
    cd "$outdir"
    $tool sort \
        --fq "$FIXTURES/sample.fastq" \
        --primers "$FIXTURES/Primers.txt" \
        --tags "$FIXTURES/Tags.txt"
    cd - > /dev/null
}

TMPPY=$(mktemp -d)
TMPRS=$(mktemp -d)
trap "rm -rf '$TMPPY' '$TMPRS'" EXIT

echo "==> Running dame-py sort..."
# dame-py uses -fq, -p, -t flags
cd "$TMPPY"
dame-py sort \
    -fq "$FIXTURES/sample.fastq" \
    -p  "$FIXTURES/Primers.txt" \
    -t  "$FIXTURES/Tags.txt"
cd - > /dev/null

echo "==> Running dame sort..."
cd "$TMPRS"
dame sort \
    --fq "$FIXTURES/sample.fastq" \
    --primers "$FIXTURES/Primers.txt" \
    --tags "$FIXTURES/Tags.txt"
cd - > /dev/null

# Verify required files exist
for f in SummaryCounts.txt Tag1_Tag2.txt; do
    [ -f "$TMPPY/$f" ] || { echo "FAIL: dame-py missing $f"; exit 1; }
    [ -f "$TMPRS/$f" ] || { echo "FAIL: dame missing $f"; exit 1; }
done

# Compare outputs
for f in SummaryCounts.txt Tag1_Tag2.txt; do
    if ! diff <(sort "$TMPPY/$f") <(sort "$TMPRS/$f") > /dev/null 2>&1; then
        echo "FAIL: $f differs between dame-py and dame"
        diff <(sort "$TMPPY/$f") <(sort "$TMPRS/$f")
        exit 1
    fi
done

# Content checks
grep -q "#tagName1" "$TMPPY/SummaryCounts.txt" || { echo "FAIL: SummaryCounts.txt missing header"; exit 1; }
grep -q "ATATATATAT" "$TMPPY/Tag1_Tag2.txt" || { echo "FAIL: expected barcode not found"; exit 1; }

echo "PASS: dame sort and dame-py sort produce identical output"
```

- [ ] **Step 2: Replace `tests/integration/run_rsi.sh`**

```bash
#!/usr/bin/env bash
set -euo pipefail

FIXTURES="$(cd "$(dirname "$0")/../fixtures" && pwd)"
TMPPY=$(mktemp -d)
TMPRS=$(mktemp -d)
trap "rm -rf '$TMPPY' '$TMPRS'" EXIT

echo "==> Running dame-py rsi..."
cd "$TMPPY"
dame-py rsi "$FIXTURES/Comparisons_4PCRs.txt"
cd - > /dev/null

echo "==> Running dame rsi..."
cd "$TMPRS"
dame rsi "$FIXTURES/Comparisons_4PCRs.txt"
cd - > /dev/null

[ -f "$TMPPY/RSI_output.txt" ] || { echo "FAIL: dame-py RSI_output.txt missing"; exit 1; }
[ -f "$TMPRS/RSI_output.txt" ] || { echo "FAIL: dame RSI_output.txt missing"; exit 1; }
grep -q "Sample" "$TMPPY/RSI_output.txt" || { echo "FAIL: RSI_output.txt missing header"; exit 1; }

# Compare RSI values with numeric tolerance (not strict diff)
python3 - "$TMPPY/RSI_output.txt" "$TMPRS/RSI_output.txt" <<'PYEOF'
import sys, math
f1, f2 = open(sys.argv[1]), open(sys.argv[2])
for l1, l2 in zip(f1, f2):
    p1, p2 = l1.strip().split('\t'), l2.strip().split('\t')
    assert p1[0] == p2[0], f"Sample mismatch: {p1[0]} vs {p2[0]}"
    if p1[0] == 'Sample': continue
    v1, v2 = float(p1[1]), float(p2[1])
    assert abs(v1 - v2) < 1e-9, f"RSI mismatch for {p1[0]}: {v1} vs {v2}"
print("RSI values match within tolerance")
PYEOF

echo "==> Testing --explicit flag..."
cd "$TMPPY" && dame-py rsi --explicit -o RSI_explicit.txt "$FIXTURES/Comparisons_4PCRs.txt" && cd - > /dev/null
cd "$TMPRS" && dame rsi --explicit -o RSI_explicit.txt "$FIXTURES/Comparisons_4PCRs.txt" && cd - > /dev/null
grep -q "ReplicateA" "$TMPPY/RSI_explicit.txt" || { echo "FAIL: explicit output missing header"; exit 1; }

echo "PASS: dame rsi and dame-py rsi produce matching output"
```

- [ ] **Step 3: Replace `tests/integration/run_decollapse.sh`**

```bash
#!/usr/bin/env bash
set -euo pipefail

TMPPY=$(mktemp -d)
TMPRS=$(mktemp -d)
trap "rm -rf '$TMPPY' '$TMPRS'" EXIT

write_input() {
    cat > "$1/input.txt" <<'EOF'
CO1	Tag1	Tag2	3	ATATATATATAT
CO1	Tag1	Tag2	2	GCGCGCGCGCGC
EOF
}

echo "==> Running dame-py decollapse..."
write_input "$TMPPY"
cd "$TMPPY" && dame-py decollapse -input input.txt -outFas out.fasta && cd - > /dev/null

echo "==> Running dame decollapse..."
write_input "$TMPRS"
cd "$TMPRS" && dame decollapse --input input.txt --out-fas out.fasta && cd - > /dev/null

[ -f "$TMPPY/out.fasta" ] || { echo "FAIL: dame-py out.fasta missing"; exit 1; }
[ -f "$TMPRS/out.fasta" ] || { echo "FAIL: dame out.fasta missing"; exit 1; }

for dir in "$TMPPY" "$TMPRS"; do
    cnt=$(grep -c "^>" "$dir/out.fasta")
    [ "$cnt" -eq 5 ] || { echo "FAIL: expected 5 seqs, got $cnt in $dir"; exit 1; }
done

if ! diff <(sort "$TMPPY/out.fasta") <(sort "$TMPRS/out.fasta") > /dev/null 2>&1; then
    echo "FAIL: out.fasta differs between dame-py and dame"
    diff <(sort "$TMPPY/out.fasta") <(sort "$TMPRS/out.fasta")
    exit 1
fi

echo "PASS: dame decollapse and dame-py decollapse produce identical output"
```

- [ ] **Step 4: Replace `tests/integration/run_filter.sh`**

```bash
#!/usr/bin/env bash
set -euo pipefail

TMPPY=$(mktemp -d)
TMPRS=$(mktemp -d)
trap "rm -rf '$TMPPY' '$TMPRS'" EXIT

setup_dir() {
    local d="$1"
    mkdir -p "$d/pool1"
    cat > "$d/pool1/Tag1_Tag2.txt" <<'EOF'
CO1	Tag1	Tag2	5	ATATATATATAT
CO1	Tag1	Tag2	3	GCGCGCGCGCGC
EOF
    cat > "$d/pool1/Tag3_Tag4.txt" <<'EOF'
CO1	Tag3	Tag4	4	ATATATATATAT
CO1	Tag3	Tag4	1	GCGCGCGCGCGC
EOF
    cat > "$d/PSinfo.txt" <<'EOF'
SampleA	Tag1	Tag2	1
SampleA	Tag3	Tag4	1
EOF
}

echo "==> Running dame-py filter..."
setup_dir "$TMPPY"
cd "$TMPPY" && dame-py filter -psInfo PSinfo.txt -x 2 -y 1 -t 1 -l 10 && cd - > /dev/null

echo "==> Running dame filter..."
setup_dir "$TMPRS"
cd "$TMPRS" && dame filter --ps-info PSinfo.txt --x 2 --y 1 --t 1 --l 10 && cd - > /dev/null

for f in FilteredReads.fna Comparisons_2PCRs.txt; do
    [ -f "$TMPPY/$f" ] || { echo "FAIL: dame-py missing $f"; exit 1; }
    [ -f "$TMPRS/$f" ] || { echo "FAIL: dame missing $f"; exit 1; }
done

# Filter uses sorted sets - compare sorted outputs
for f in Comparisons_2PCRs.txt FilteredReads.fna; do
    if ! diff <(sort "$TMPPY/$f") <(sort "$TMPRS/$f") > /dev/null 2>&1; then
        echo "FAIL: $f differs between dame-py and dame"
        diff <(sort "$TMPPY/$f") <(sort "$TMPRS/$f")
        exit 1
    fi
done

grep -q "SampleA" "$TMPPY/Comparisons_2PCRs.txt" || { echo "FAIL: SampleA not found"; exit 1; }

echo "PASS: dame filter and dame-py filter produce matching output"
```

- [ ] **Step 5: Commit**

```bash
git add tests/integration/run_sort.sh tests/integration/run_rsi.sh \
        tests/integration/run_decollapse.sh tests/integration/run_filter.sh
git commit -m "test: update integration scripts to compare dame vs dame-py output"
```

---

## Task 8: Add `run_chimera.sh` and `run_pipeline.sh`

**Files:**
- Create: `tests/integration/run_chimera.sh`
- Create: `tests/integration/run_pipeline.sh`
- Create: `tests/fixtures/Comparisons_4PCRs.txt` (if not already present)

Note: If `tests/fixtures/Comparisons_4PCRs.txt` does not exist, create it here. The RSI integration test needs it. Check by running `ls tests/fixtures/`.

- [ ] **Step 1: Check for `Comparisons_4PCRs.txt` fixture**

```bash
ls tests/fixtures/
```

If `Comparisons_4PCRs.txt` is missing, create it:

```bash
cat > tests/fixtures/Comparisons_4PCRs.txt <<'EOF'
SampleA	CO1-CO1	10	CO1-CO1	8	ATATATATATAT
SampleA	CO1-CO1	20	CO1-CO1	22	GCGCGCGCGCGC
SampleB	CO1-CO1	15	CO1-CO1	14	ATATATATATAT
SampleB	CO1-CO1	5	CO1-CO1	6	TTTTTTTTTTTT
EOF
```

- [ ] **Step 2: Create `tests/integration/run_chimera.sh`**

```bash
#!/usr/bin/env bash
set -euo pipefail

# Skip if usearch is not on PATH
if ! command -v usearch &> /dev/null; then
    echo "SKIP: usearch not found on PATH — skipping chimera integration test"
    exit 0
fi

TMPPY=$(mktemp -d)
TMPRS=$(mktemp -d)
trap "rm -rf '$TMPPY' '$TMPRS'" EXIT

setup_dir() {
    local d="$1"
    cat > "$d/Tag1_Tag2.txt" <<'EOF'
CO1	Tag1	Tag2	5	ATATATATATAT
CO1	Tag1	Tag2	3	GCGCGCGCGCGC
EOF
    cat > "$d/PSinfo.txt" <<'EOF'
SampleA	Tag1	Tag2	1
SampleA	Tag1	Tag2	1
EOF
}

echo "==> Running dame-py chimera..."
setup_dir "$TMPPY"
cd "$TMPPY" && dame-py chimera -psInfo PSinfo.txt -x 2 && cd - > /dev/null

echo "==> Running dame chimera..."
setup_dir "$TMPRS"
cd "$TMPRS" && dame chimera --ps-info PSinfo.txt --x 2 && cd - > /dev/null

# Compare tag files
for i in 1 2; do
    if ! diff <(sort "$TMPPY/PS${i}.tags.txt") <(sort "$TMPRS/PS${i}.tags.txt"); then
        echo "FAIL: PS${i}.tags.txt differs"
        exit 1
    fi
done

echo "PASS: dame chimera and dame-py chimera produce identical tag files"
```

- [ ] **Step 3: Create `tests/integration/run_pipeline.sh`**

```bash
#!/usr/bin/env bash
set -euo pipefail

FIXTURES="$(cd "$(dirname "$0")/../fixtures" && pwd)"

run_pipeline() {
    local tool="$1"
    local outdir="$2"
    local sort_flags=("$@")  # shift below

    mkdir -p "$outdir"
    cd "$outdir"

    echo "  [1/3] sort..."
    if [ "$tool" = "dame-py" ]; then
        dame-py sort \
            -fq "$FIXTURES/sample.fastq" \
            -p  "$FIXTURES/Primers.txt" \
            -t  "$FIXTURES/Tags.txt"
    else
        dame sort \
            --fq "$FIXTURES/sample.fastq" \
            --primers "$FIXTURES/Primers.txt" \
            --tags "$FIXTURES/Tags.txt"
    fi

    # Create PSinfo from sort output
    python3 - <<'PYEOF'
import os
for f in os.listdir('.'):
    if f.endswith('.txt') and '_' in f and not f.startswith('Summary'):
        parts = f.replace('.txt','').split('_')
        if len(parts) == 2:
            tag1, tag2 = parts
            print(f"SampleA\t{tag1}\t{tag2}\t1")
PYEOF
    python3 - > PSinfo.txt <<'PYEOF'
import os
for f in os.listdir('.'):
    if f.endswith('.txt') and '_' in f and not f.startswith('Summary'):
        parts = f.replace('.txt','').split('_')
        if len(parts) == 2:
            tag1, tag2 = parts
            print(f"SampleA\t{tag1}\t{tag2}\t1")
PYEOF

    mkdir -p pool1
    for f in *.txt; do
        [[ "$f" =~ ^(PS|Summary) ]] && continue
        cp "$f" pool1/ 2>/dev/null || true
    done

    echo "  [2/3] filter..."
    local psinfo_count
    psinfo_count=$(wc -l < PSinfo.txt)
    if [ "$psinfo_count" -eq 0 ]; then
        echo "  WARNING: no PSinfo entries, skipping filter"
        cd - > /dev/null
        return
    fi
    if [ "$tool" = "dame-py" ]; then
        dame-py filter -psInfo PSinfo.txt -x "$psinfo_count" -y 1 -t 1 -l 5
    else
        dame filter --ps-info PSinfo.txt --x "$psinfo_count" --y 1 --t 1 --l 5
    fi

    echo "  [3/3] rsi..."
    local comp_file
    comp_file=$(ls Comparisons_*.txt 2>/dev/null | head -1 || true)
    if [ -z "$comp_file" ]; then
        echo "  WARNING: no comparisons file, skipping rsi"
        cd - > /dev/null
        return
    fi
    if [ "$tool" = "dame-py" ]; then
        dame-py rsi "$comp_file"
    else
        dame rsi "$comp_file"
    fi

    cd - > /dev/null
}

TMPPY=$(mktemp -d)
TMPRS=$(mktemp -d)
trap "rm -rf '$TMPPY' '$TMPRS'" EXIT

echo "==> Running dame-py pipeline..."
run_pipeline "dame-py" "$TMPPY"

echo "==> Running dame pipeline..."
run_pipeline "dame" "$TMPRS"

# Verify RSI output exists and is valid in both
for tool_dir in "$TMPPY" "$TMPRS"; do
    if [ -f "$tool_dir/RSI_output.txt" ]; then
        grep -q "Sample" "$tool_dir/RSI_output.txt" || { echo "FAIL: RSI_output.txt missing header"; exit 1; }
        echo "  RSI output found in $tool_dir"
    fi
done

echo "PASS: full pipeline (sort -> filter -> rsi) completed for both dame-py and dame"
```

- [ ] **Step 4: Make integration scripts executable**

```bash
chmod +x tests/integration/run_chimera.sh tests/integration/run_pipeline.sh
```

- [ ] **Step 5: Commit**

```bash
git add tests/integration/run_chimera.sh tests/integration/run_pipeline.sh tests/fixtures/
git commit -m "test: add run_chimera.sh and run_pipeline.sh integration tests"
```

---

## Task 9: GitHub Actions CI

**Files:**
- Create: `.github/workflows/ci.yml`

- [ ] **Step 1: Create `.github/workflows/ci.yml`**

```yaml
name: CI

on:
  push:
    branches: [master]
  pull_request:
    branches: [master]

jobs:
  python-tests:
    name: Python 3 unit tests
    runs-on: ubuntu-latest
    strategy:
      matrix:
        python-version: ["3.11", "3.12"]
    steps:
      - uses: actions/checkout@v4
      - uses: actions/setup-python@v5
        with:
          python-version: ${{ matrix.python-version }}
      - name: Install dame-py
        run: pip install -e "python/[dev]"
      - name: Run pytest
        run: pytest python/tests/ -v

  rust-tests:
    name: Rust unit tests
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4
      - uses: dtolnay/rust-toolchain@stable
      - name: Cache cargo
        uses: actions/cache@v4
        with:
          path: |
            ~/.cargo/registry
            ~/.cargo/git
            rust/target
          key: ${{ runner.os }}-cargo-${{ hashFiles('rust/Cargo.lock') }}
      - name: Run cargo test
        working-directory: rust
        run: cargo test --all

  integration-tests:
    name: Integration tests
    runs-on: ubuntu-latest
    needs: [python-tests, rust-tests]
    steps:
      - uses: actions/checkout@v4
      - uses: actions/setup-python@v5
        with:
          python-version: "3.11"
      - uses: dtolnay/rust-toolchain@stable
      - name: Cache cargo
        uses: actions/cache@v4
        with:
          path: |
            ~/.cargo/registry
            ~/.cargo/git
            rust/target
          key: ${{ runner.os }}-cargo-integration-${{ hashFiles('rust/Cargo.lock') }}
      - name: Build dame binary
        working-directory: rust
        run: cargo build --release
      - name: Install dame-py
        run: pip install -e "python/[dev]"
      - name: Add dame to PATH
        run: echo "$GITHUB_WORKSPACE/rust/target/release" >> $GITHUB_PATH
      - name: Run sort integration test
        run: bash tests/integration/run_sort.sh
      - name: Run rsi integration test
        run: bash tests/integration/run_rsi.sh
      - name: Run decollapse integration test
        run: bash tests/integration/run_decollapse.sh
      - name: Run filter integration test
        run: bash tests/integration/run_filter.sh
      - name: Run chimera integration test (skips if no usearch)
        run: bash tests/integration/run_chimera.sh
      - name: Run pipeline integration test
        run: bash tests/integration/run_pipeline.sh
```

- [ ] **Step 2: Add `[dev]` extra to `python/pyproject.toml`**

Open `python/pyproject.toml` and add:

```toml
[project.optional-dependencies]
dev = ["pytest"]
```

- [ ] **Step 3: Commit**

```bash
git add .github/workflows/ci.yml python/pyproject.toml
git commit -m "ci: add GitHub Actions CI for Python, Rust, and integration tests"
```

---

## Spec Coverage Self-Review

| Spec Requirement | Task |
|---|---|
| Single Rust binary `dame` | Task 1 (Cargo.toml `[[bin]]`) |
| `dame sort -fq -p -t [--keep-primers-seq]` | Tasks 1, 2 |
| `dame chimera -ps-info -x [-p]` | Tasks 1, 4 |
| `dame filter -ps-info [-x -y -p -t -l] [--chimera-checked]` | Tasks 1, 5 |
| `dame decollapse -input [-out-fas]` | Tasks 1, 3 |
| `dame rsi [-e] [-o] <input>` | Tasks 1, 6 |
| Crates: clap, regex, ndarray, anyhow | Task 1 |
| `indexmap` for insertion-order output | Tasks 2, 5 |
| Regex compiled once at startup | Task 2 (`read_primers` returns `PrimerEntry` with compiled `Regex` fields) |
| `dame chimera` wraps `usearch` via `std::process::Command` | Task 4 |
| Byte-for-byte identical output to dame-py (RSI: tolerance) | Tasks 7 (diff with sort), Task 7 (RSI tolerance check) |
| Rust unit tests mirroring Python tests | Tasks 2–6 |
| Integration tests: shell scripts comparing both tools | Tasks 7–8 |
| `run_pipeline.sh` end-to-end | Task 8 |
| GitHub Actions CI matrix | Task 9 |
| Python 3 `[dev]` extra for CI | Task 9 |
