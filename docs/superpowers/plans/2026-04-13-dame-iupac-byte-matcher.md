# DAMe IUPAC Byte Matcher Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the regex-based IUPAC primer matcher in `rust/src/sort.rs` with a byte sliding-window approach, removing the `regex` crate dependency and bringing DAMe's Rust sort step into the same performance tier as Begum (~30–50× over Python).

**Architecture:** All changes are confined to `rust/src/sort.rs`. The public interface of every other module is untouched. `PrimerEntry` drops its four `Regex`-based fields and gains two `Vec<Vec<u8>>` byte-sequence fields. A new private `iupac_matches(primer_byte: u8, read_byte: u8) -> bool` lookup function and a `find_primer(primer: &[u8], seq: &[u8]) -> Option<(usize, usize)>` sliding-window function replace all four `Regex::find()` call sites in `get_pieces_info`. The `regex` crate is removed from `Cargo.toml`. The Python implementation (`dame-py`) is unchanged and continues to use regex — the integration test suite diffs their outputs to verify byte-for-byte agreement on the tutorial dataset.

**Tech Stack:** Rust stable, `cargo test`, existing `tempfile` dev-dependency, existing integration shell scripts in `tests/integration/`.

---

## Background: what currently happens

In `rust/src/sort.rs`, `read_primers()` compiles four `Regex` objects per primer (forward, reverse, RC-forward, RC-reverse) and stores them in `PrimerEntry`. In `get_pieces_info()`, called on every FASTQ read, four `Regex::find()` calls scan the full read sequence. This costs ~1,000–2,000 CPU cycles per `.find()` call, making the sort step ~10× slower than it needs to be for a Rust implementation.

The replacement stores each primer as raw bytes (`Vec<u8>`) and uses a left-to-right sliding window with an IUPAC lookup table (`iupac_matches`). This is semantically equivalent for exact IUPAC matching (no mismatches) and costs ~10–20 CPU cycles per position — the same algorithmic complexity but with a ~100× smaller constant.

## File map

| File | Action | What changes |
|------|--------|-------------|
| `rust/src/sort.rs` | Modify | Remove `ambig_expand`, `Regex` imports; replace `PrimerEntry` fields; add `iupac_matches` + `find_primer`; rewrite `read_primers` + `get_pieces_info` |
| `rust/Cargo.toml` | Modify | Remove `regex = "1"` dependency |
| `rust/tests/sort_test.rs` | Modify | Update `test_read_primers` assertions to match new `PrimerEntry` fields; add `test_iupac_matches` and `test_find_primer` unit tests |

No other files change.

---

## Task 1: Add `iupac_matches` and `find_primer` with unit tests (TDD first)

**Files:**
- Modify: `rust/src/sort.rs` (add two new `pub` functions near top, before `PrimerEntry`)
- Modify: `rust/tests/sort_test.rs` (add tests)

The IUPAC truth table used by `iupac_matches`:

| Primer byte | Matches read bytes |
|---|---|
| A | A |
| C | C |
| G | G |
| T | T |
| R | A, G |
| Y | C, T |
| S | C, G |
| W | A, T |
| K | G, T |
| M | A, C |
| B | C, G, T |
| D | A, G, T |
| H | A, C, T |
| V | A, C, G |
| N | A, C, G, T |

Input bytes are always uppercase ASCII (FASTQ sequences). Both `primer_byte` and `read_byte` are already uppercase — no case folding needed.

- [ ] **Step 1: Write the failing tests**

Add to `rust/tests/sort_test.rs` (after the existing `test_rc_ambiguous_n` test, before `test_read_tags`):

```rust
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
```

- [ ] **Step 2: Run tests — confirm they fail**

```bash
cd /path/to/DAMe
cargo test --manifest-path rust/Cargo.toml test_iupac_matches 2>&1 | grep -E "error|FAILED|not found"
```

Expected: compile error — `iupac_matches` and `find_primer` not yet defined.

- [ ] **Step 3: Add `iupac_matches` and `find_primer` to `sort.rs`**

In `rust/src/sort.rs`, add these two functions immediately after the `rc()` function (before the existing `ambig_expand` function). Mark them `pub` so the test crate can access them:

```rust
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
        if primer.iter().zip(&seq[i..]).all(|(&p, &s)| iupac_matches(p, s)) {
            return Some((i, i + plen));
        }
    }
    None
}
```

- [ ] **Step 4: Run new tests — confirm they pass**

```bash
cargo test --manifest-path rust/Cargo.toml test_iupac 2>&1 | tail -5
cargo test --manifest-path rust/Cargo.toml test_find_primer 2>&1 | tail -5
```

Expected: all 8 new tests pass.

- [ ] **Step 5: Commit**

```bash
git add rust/src/sort.rs rust/tests/sort_test.rs
git commit -m "feat: add iupac_matches and find_primer byte-level primer matcher"
```

---

## Task 2: Replace `PrimerEntry` and `read_primers` — store bytes, not Regex

**Files:**
- Modify: `rust/src/sort.rs` — change `PrimerEntry` struct and `read_primers` function
- Modify: `rust/tests/sort_test.rs` — update `test_read_primers` assertions

`PrimerEntry` currently holds four `Vec<String>` (pattern strings) and four compiled `Regex` objects. Replace all of these with two `Vec<Vec<u8>>` fields: one for the start-side primers and one for the end-side primers, as raw byte sequences ready for `find_primer`.

The four byte sequences per primer are:
- `start_primers[0]`: forward primer F bytes (used as start anchor in forward orientation)
- `start_primers[1]`: reverse primer R bytes (used as start anchor in reverse orientation)
- `end_primers[0]`: RC(F) bytes (used as end anchor in reverse orientation)
- `end_primers[1]`: RC(R) bytes (used as end anchor in forward orientation)

This preserves the exact same indexing logic as the current `a_side_re`/`b_side_re` fields.

- [ ] **Step 1: Update `test_read_primers`**

In `rust/tests/sort_test.rs`, replace:

```rust
fn test_read_primers() {
    let dir = tempdir().unwrap();
    let path = dir.path().join("primers.txt");
    let mut f = std::fs::File::create(&path).unwrap();
    writeln!(f, "CO1\tACGT\tTTTT").unwrap();

    let primers = read_primers(path.to_str().unwrap()).unwrap();
    assert!(primers.contains_key("CO1"));
    let co1 = &primers["CO1"];
    assert_eq!(co1.a_side.len(), 2);
    assert_eq!(co1.b_side.len(), 2);
    assert_eq!(co1.a_side_re.len(), 2);
    assert_eq!(co1.b_side_re.len(), 2);
}
```

With:

```rust
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
```

- [ ] **Step 2: Run — confirm test fails (field names don't exist yet)**

```bash
cargo test --manifest-path rust/Cargo.toml test_read_primers 2>&1 | grep -E "error|no field"
```

Expected: compile error — `start_primers` field not found.

- [ ] **Step 3: Replace `PrimerEntry` struct in `sort.rs`**

Find and replace the `PrimerEntry` struct (currently lines 70–75):

```rust
// REMOVE this entire struct:
pub struct PrimerEntry {
    pub a_side: Vec<String>,    // [F_pattern, R_pattern] as regex strings
    pub b_side: Vec<String>,    // [RC(F)_pattern, RC(R)_pattern] as regex strings
    pub a_side_re: Vec<Regex>,  // compiled a_side regexes
    pub b_side_re: Vec<Regex>,  // compiled b_side regexes
}

// REPLACE with:
pub struct PrimerEntry {
    /// [F_bytes, R_bytes] — start-side primers as raw bytes for find_primer
    pub start_primers: Vec<Vec<u8>>,
    /// [RC(F)_bytes, RC(R)_bytes] — end-side primers as raw bytes for find_primer
    pub end_primers: Vec<Vec<u8>>,
}
```

- [ ] **Step 4: Replace `read_primers` function in `sort.rs`**

The current `read_primers` function (lines 126–177) calls `ambig_expand()` to build regex strings, then compiles them with `Regex::new()`. Replace the body of `read_primers` entirely (keep the signature):

```rust
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
```

Also remove the now-unused `ambig_expand` function entirely (lines 46–68 in the original).

- [ ] **Step 5: Remove `use regex::Regex` import from `sort.rs`**

At the top of `sort.rs`, remove the line:
```rust
use regex::Regex;
```

- [ ] **Step 6: Remove `regex` from `Cargo.toml`**

In `rust/Cargo.toml`, remove:
```toml
regex = "1"
```

- [ ] **Step 7: Run — confirm test_read_primers passes (other tests will fail until Task 3)**

```bash
cargo test --manifest-path rust/Cargo.toml test_read_primers 2>&1 | tail -5
```

Expected: `test_read_primers` passes. Other tests using `get_pieces_info` will fail to compile because `a_side_re`/`b_side_re` no longer exist — that is expected and fixed in Task 3.

- [ ] **Step 8: Commit**

```bash
git add rust/src/sort.rs rust/Cargo.toml rust/tests/sort_test.rs
git commit -m "refactor: replace Regex fields in PrimerEntry with raw byte vecs, remove regex crate"
```

---

## Task 3: Rewrite `get_pieces_info` to use `find_primer`

**Files:**
- Modify: `rust/src/sort.rs` — rewrite `get_pieces_info` body only (signature unchanged)

The current function calls `primer.a_side_re[0].find(line)` (etc.) returning `Option<Match>` where `.start()` and `.end()` give byte positions. `find_primer` returns `Option<(usize, usize)>` with the same semantics. The replacement is a direct structural substitution.

The function signature does not change:
```rust
pub fn get_pieces_info(
    line: &str,
    primers: &IndexMap<String, PrimerEntry>,
    tags: &HashMap<String, Vec<String>>,
    keep_primers_seq: bool,
) -> Option<PieceInfo>
```

- [ ] **Step 1: Write the updated `get_pieces_info` body**

Replace the entire body of `get_pieces_info` (from the `for (key, primer) in primers {` line to the closing `None`) with:

```rust
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
                let tag_name1 = tags.iter().find(|(_, v)| v[0] == tag1_str).map(|(k, _)| k.clone());
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
                    let between = rc(between_raw);
                    let tag1_str = &line[..prim_ini_tags];
                    let tag2_str = &line[prim_fin_tags..];
                    // In reverse orientation: tag roles are swapped
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
```

- [ ] **Step 2: Run all existing unit tests**

```bash
cargo test --manifest-path rust/Cargo.toml 2>&1 | tail -20
```

Expected: all tests pass (previously 13, now 13 + the 8 new ones from Task 1 = 21 total).

- [ ] **Step 3: Run clippy**

```bash
cargo clippy --manifest-path rust/Cargo.toml -- -D warnings 2>&1
```

Expected: no warnings. Fix any that appear before proceeding.

- [ ] **Step 4: Commit**

```bash
git add rust/src/sort.rs
git commit -m "perf: replace regex find with byte sliding-window in get_pieces_info"
```

---

## Task 4: Run integration tests to verify output parity with dame-py

**Files:** None changed — this is a verification task only.

The integration scripts in `tests/integration/` run both `dame-py` and the new `dame` binary on the tutorial dataset and diff their outputs. This is the definitive check that removing the regex engine did not change any demultiplexing decisions.

- [ ] **Step 1: Build the new binary**

```bash
cargo build --release --manifest-path rust/Cargo.toml 2>&1 | tail -3
```

Expected: `Finished release profile`.

- [ ] **Step 2: Verify `dame` is on PATH (or use absolute path)**

```bash
# Option A: add to PATH
export PATH="$(pwd)/rust/target/release:$PATH"
dame --help | head -3

# Option B: verify absolute path works
./rust/target/release/dame --help | head -3
```

- [ ] **Step 3: Run sort integration test**

```bash
bash tests/integration/run_sort.sh
```

Expected: `PASS` — no diff between dame and dame-py on `SummaryCounts.txt` and `tag*.txt` files for the tutorial dataset.

- [ ] **Step 4: Run filter integration test**

```bash
bash tests/integration/run_filter.sh
```

Expected: `PASS`.

- [ ] **Step 5: Run RSI integration test**

```bash
bash tests/integration/run_rsi.sh
```

Expected: `PASS` — RSI values within 1×10⁻⁹ tolerance.

- [ ] **Step 6: Run decollapse integration test**

```bash
bash tests/integration/run_decollapse.sh
```

Expected: `PASS`.

- [ ] **Step 7: Run full pipeline integration test**

```bash
bash tests/integration/run_pipeline.sh
```

Expected: `PASS`.

- [ ] **Step 8: If any test fails, diagnose before proceeding**

A failure here means the byte matcher produced a different match position than the regex engine for some read. To diagnose:
1. Note which output file differs (e.g., `SummaryCounts.txt` count mismatch)
2. Diff the tag files: `diff <(sort dame_out/tag*.txt) <(sort damepy_out/tag*.txt)`
3. Find a read that the two tools classified differently by running both on a single-read FASTQ
4. Check if the discrepancy is in primer position detection or tag lookup

Do not proceed to Task 5 until all integration tests pass.

- [ ] **Step 9: Commit (no file changes — just note the verification passed)**

This task produces no code changes. The commit from Task 3 is sufficient. Document the verification in the PR description.

---

## Task 5: Benchmark and update README

**Files:**
- Modify: `README.md` — update Performance section with v2.2 results

Run the sort benchmark on the 196k-read synthetic dataset at `/tmp/dame_bench/` to measure the improvement. The benchmark script from the v2.1 work (`/tmp/run_bench2.sh`) can be reused.

- [ ] **Step 1: Run sort benchmark (3 runs each)**

```bash
DAME_RS="$(pwd)/rust/target/release/dame"
BENCH="/tmp/dame_bench"
PRIMERS="$BENCH/Primers.txt"
TAGS="$BENCH/Tags.txt"
POOL1="$BENCH/Pool1.fastq"
POOL2="$BENCH/Pool2.fastq"
ts() { python3 -c "import time; print(int(time.time()*1000))"; }

echo "=== dame 2.2 sort ==="
for i in 1 2 3; do
  T0=$(ts); "$DAME_RS" sort --fq "$POOL1" --primers "$PRIMERS" --tags "$TAGS" > /dev/null 2>&1; T1=$(ts); P1=$((T1-T0))
  rm -f "$BENCH"/tag*.txt "$BENCH"/SummaryCounts.txt
  T0=$(ts); "$DAME_RS" sort --fq "$POOL2" --primers "$PRIMERS" --tags "$TAGS" > /dev/null 2>&1; T2=$(ts); P2=$((T2-T0))
  rm -f "$BENCH"/tag*.txt "$BENCH"/SummaryCounts.txt
  echo "  Run $i: Pool1=${P1}ms Pool2=${P2}ms"
done

echo "=== dame-py sort ==="
for i in 1 2 3; do
  T0=$(ts); dame-py sort -fq "$POOL1" -p "$PRIMERS" -t "$TAGS" > /dev/null 2>&1; T1=$(ts); P1=$((T1-T0))
  rm -f "$BENCH"/tag*.txt "$BENCH"/SummaryCounts.txt
  T0=$(ts); dame-py sort -fq "$POOL2" -p "$PRIMERS" -t "$TAGS" > /dev/null 2>&1; T2=$(ts); P2=$((T2-T0))
  rm -f "$BENCH"/tag*.txt "$BENCH"/SummaryCounts.txt
  echo "  Run $i: Pool1=${P1}ms Pool2=${P2}ms"
done
```

- [ ] **Step 2: Update the Performance section in `README.md`**

Add a `Rust 2.2` column to the large-dataset table (alongside the existing 2.0 and 2.1 columns). Use the median of runs 2 and 3 (run 1 is the cold-cache outlier). Also update the explanatory text below the table. Also bump the title to `DAMe v2.2`.

- [ ] **Step 3: Update dev history item 7 in `README.md`**

Add a note to development history item 7 mentioning the byte matcher:

```
... (a) ... needletail ... (b) ... ahash ... (c) the regex-based IUPAC
primer matcher in `sort` was replaced with a byte sliding-window and IUPAC
lookup table, removing the `regex` crate dependency and achieving ~X× sort
speedup over Python on large pools.
```

Fill in the actual measured speedup.

- [ ] **Step 4: Bump `rust/Cargo.toml` version to `0.3.0`**

```toml
version = "0.3.0"
```

- [ ] **Step 5: Commit**

```bash
git add README.md rust/Cargo.toml
git commit -m "docs: update benchmark tables for v2.2, bump crate to 0.3.0"
```

---

## Self-review

**Spec coverage:**
- ✅ Remove regex crate: Task 2 removes `regex` from Cargo.toml and `use regex::Regex`
- ✅ `iupac_matches` with full IUPAC table: Task 1, Step 3
- ✅ `find_primer` sliding window: Task 1, Step 3
- ✅ `PrimerEntry` refactored: Task 2
- ✅ `read_primers` rewritten: Task 2
- ✅ `get_pieces_info` rewritten: Task 3
- ✅ `ambig_expand` removed: Task 2, Step 4
- ✅ Unit tests for new functions: Task 1
- ✅ `test_read_primers` updated: Task 2
- ✅ Integration tests verify parity with dame-py: Task 4
- ✅ Benchmark and README update: Task 5

**Placeholder scan:** No TBDs or "similar to above" references found. All code blocks are complete.

**Type consistency:**
- `find_primer(&[u8], &[u8]) -> Option<(usize, usize)>` defined in Task 1, used identically in Task 3
- `PrimerEntry.start_primers: Vec<Vec<u8>>` and `.end_primers: Vec<Vec<u8>>` defined in Task 2, referenced correctly in Task 2 test update and Task 3 implementation
- `iupac_matches(u8, u8) -> bool` defined in Task 1, used inside `find_primer` in Task 1
