# Tag Mismatches + Anchored Matching (sort) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a per-tag substitution tolerance (`-mt` / `--tag-mismatches`, default 0) to the `sort` step in both implementations, via a tag-anchored matcher used whenever primer or tag mismatches are allowed.

**Architecture:** `run()` keeps today's exact matcher for `m==0 && mt==0` (byte-identical) and dispatches to a new **anchored** matcher when `m>0 || mt>0`. The anchored matcher finds tag candidates at the read ends by IUPAC Hamming distance (≤ mt), verifies the primers at the expected offsets (≤ m), scores each valid assembly by total mismatches, keeps the unique minimum, and discards ties as ambiguous. A startup warning flags an unsafe `mt` relative to the tag set's minimum pairwise distance.

**Tech Stack:** Rust (`clap`, `needletail`, `ahash`, `indexmap`), Python 3 (stdlib), pytest, `cargo test`, bash integration scripts, GitHub Actions.

**Reference spec:** `docs/superpowers/specs/2026-06-11-tag-mismatches-design.md`

---

## File Structure

- `rust/src/sort.rs` — add `hamming_iupac`, extend `TagLookup` with ordered `fwd_list`/`rc_list`, add `get_pieces_info_anchored`, add `min_equal_length_tag_distance`, add `tag_mismatches` to `SortArgs`, dispatch + warning in `run`.
- `rust/tests/sort_test.rs` — unit tests for the above.
- `python/dame/modules_sort.py` — add `hamming_iupac`, `GetPiecesInfoMismatch`, `min_equal_length_tag_distance`.
- `python/dame/sort.py` — `-mt` flag, dispatch + warning in `run`.
- `python/tests/test_sort.py` — unit tests.
- `tests/fixtures/sample_tag_err.fastq` — new fixture (single tag-base error).
- `tests/integration/run_sort_tag_mismatch.sh` — new cross-impl test.
- `.github/workflows/ci.yml` — run the new script.
- `rust/Cargo.toml`, `README.md`, `tutorial/README.md` — docs + version bump.

**Parity note (applies throughout):** Rust and Python must produce byte-identical output. Both iterate primers and tags in file order; both score `tag1_mm + start_primer_mm + end_primer_mm + tag2_mm`, keep the unique minimum, and discard ties. Keep the two anchored matchers structurally identical.

---

## Task 1: Rust `hamming_iupac` + `TagLookup` ordered lists

**Files:**
- Modify: `rust/src/sort.rs`
- Test: `rust/tests/sort_test.rs`

- [ ] **Step 1: Write failing tests**

In `rust/tests/sort_test.rs`, add after the `iupac_matches` tests:

```rust
#[test]
fn test_hamming_iupac_exact() {
    use dame::sort::hamming_iupac;
    assert_eq!(hamming_iupac(b"ACGT", b"ACGT"), 0);
    assert_eq!(hamming_iupac(b"ACGT", b"ACGA"), 1);
    assert_eq!(hamming_iupac(b"ACGT", b"AAGA"), 2);
}

#[test]
fn test_hamming_iupac_respects_iupac() {
    use dame::sort::hamming_iupac;
    // R matches A or G at zero cost
    assert_eq!(hamming_iupac(b"GCRTGC", b"GCATGC"), 0);
    assert_eq!(hamming_iupac(b"GCRTGC", b"GCGTGC"), 0);
    assert_eq!(hamming_iupac(b"GCRTGC", b"GCCTGC"), 1);
}

#[test]
fn test_hamming_iupac_length_mismatch_is_max() {
    use dame::sort::hamming_iupac;
    assert_eq!(hamming_iupac(b"ACGT", b"ACG"), usize::MAX);
}

#[test]
fn test_read_tags_ordered_lists() {
    use dame::sort::read_tags;
    let dir = tempdir().unwrap();
    let path = dir.path().join("tags.txt");
    std::fs::write(&path, "AAAA\tTag1\nCCCC\tTag2\n").unwrap();
    let tags = read_tags(path.to_str().unwrap()).unwrap();
    // fwd_list preserves file order: (seq, name)
    assert_eq!(tags.fwd_list[0], (b"AAAA".to_vec(), "Tag1".to_string()));
    assert_eq!(tags.fwd_list[1], (b"CCCC".to_vec(), "Tag2".to_string()));
    // rc_list holds RC seqs in the same order: rc(AAAA)=TTTT, rc(CCCC)=GGGG
    assert_eq!(tags.rc_list[0], (b"TTTT".to_vec(), "Tag1".to_string()));
    assert_eq!(tags.rc_list[1], (b"GGGG".to_vec(), "Tag2".to_string()));
}
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `cargo test --manifest-path rust/Cargo.toml --test sort_test hamming_iupac`
Expected: FAIL — `hamming_iupac` not found; `fwd_list`/`rc_list` fields not found.

- [ ] **Step 3: Add `hamming_iupac`**

In `rust/src/sort.rs`, add immediately after the `iupac_matches` function:

```rust
/// Count positions where `region` fails the IUPAC constraint of `pattern`.
/// Returns `usize::MAX` when lengths differ, so an out-of-range region never
/// satisfies any mismatch budget.
pub fn hamming_iupac(pattern: &[u8], region: &[u8]) -> usize {
    if pattern.len() != region.len() {
        return usize::MAX;
    }
    pattern
        .iter()
        .zip(region)
        .filter(|(&p, &r)| !iupac_matches(p, r))
        .count()
}
```

- [ ] **Step 4: Extend `TagLookup` and populate the lists in `read_tags`**

In `rust/src/sort.rs`, replace the `TagLookup` struct:

```rust
/// Reverse lookups for tag sequences.
/// `by_fwd`/`by_rc` are O(1) maps used by the exact path; `fwd_list`/`rc_list`
/// are file-ordered, name-deduplicated lists iterated by the anchored matcher.
pub struct TagLookup {
    pub by_fwd: HashMap<Vec<u8>, String>,
    pub by_rc: HashMap<Vec<u8>, String>,
    pub fwd_list: Vec<(Vec<u8>, String)>,
    pub rc_list: Vec<(Vec<u8>, String)>,
}
```

Then in `read_tags`, replace the body's accumulator setup and loop tail. The function becomes:

```rust
pub fn read_tags(path: &str) -> Result<TagLookup> {
    let file = File::open(path).with_context(|| format!("Cannot open tags file: {path}"))?;
    let reader = BufReader::new(file);
    let mut by_fwd: HashMap<Vec<u8>, String> = HashMap::default();
    let mut by_rc: HashMap<Vec<u8>, String> = HashMap::default();
    let mut fwd_list: Vec<(Vec<u8>, String)> = Vec::new();
    let mut rc_list: Vec<(Vec<u8>, String)> = Vec::new();
    let mut seen: std::collections::HashSet<String> = std::collections::HashSet::new();

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
        // First occurrence per name only, mirroring the Python readTags dedup.
        if seen.insert(name.to_string()) {
            fwd_list.push((seq.as_bytes().to_vec(), name.to_string()));
            rc_list.push((rc(seq).as_bytes().to_vec(), name.to_string()));
        }
    }

    Ok(TagLookup { by_fwd, by_rc, fwd_list, rc_list })
}
```

- [ ] **Step 5: Run tests to verify they pass**

Run: `cargo test --manifest-path rust/Cargo.toml --test sort_test`
Expected: PASS (all existing + 4 new).

- [ ] **Step 6: Commit**

```bash
git add rust/src/sort.rs rust/tests/sort_test.rs
git commit -m "feat(rust): hamming_iupac + ordered tag lists for anchored matching

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

## Task 2: Rust anchored matcher `get_pieces_info_anchored`

**Files:**
- Modify: `rust/src/sort.rs`
- Test: `rust/tests/sort_test.rs`

Context: `PrimerEntry.start_primers = [F, R]`, `end_primers = [RC(F), RC(R)]` (byte vecs). `make_test_tags()` gives AAAA=Tag1, CCCC=Tag2, GGGG=Tag3, TTTT=Tag4; `make_test_primers()` gives CO1 F=ACGT R=TGCA. `rc_bytes` and `PieceInfo` already exist.

- [ ] **Step 1: Write failing tests**

In `rust/tests/sort_test.rs`, add after `test_get_pieces_info_lowercase_read`:

```rust
#[test]
fn test_anchored_forward_tag_mismatch() {
    use dame::sort::get_pieces_info_anchored;
    // Forward read AAAA|ACGT|ATATAT|TGCA|GGGG with tag1 AAAA miscalled AAAG.
    // tag1 Tag1 within 1 mismatch; everything else exact. mt=1 recovers; mt=0 rejects.
    let tags = make_test_tags();
    let primers = make_test_primers();
    let line = "AAAGACGTATATATTGCAGGGG";
    assert!(get_pieces_info_anchored(line, &primers, &tags, false, 0, 0).is_none());
    let info = get_pieces_info_anchored(line, &primers, &tags, false, 0, 1).unwrap();
    assert_eq!(info.tag1, "Tag1");
    assert_eq!(info.tag2, "Tag2");
    assert_eq!(info.between, "ATATAT");
    assert_eq!(info.primer_name, "CO1");
}

#[test]
fn test_anchored_primer_only_parity() {
    use dame::sort::get_pieces_info_anchored;
    // m=1, mt=0: forward primer ACGT miscalled ACGA — same result the leftmost
    // path produces for this fixture.
    let tags = make_test_tags();
    let primers = make_test_primers();
    let line = "AAAAACGAATATATTGCAGGGG";
    assert!(get_pieces_info_anchored(line, &primers, &tags, false, 0, 0).is_none());
    let info = get_pieces_info_anchored(line, &primers, &tags, false, 1, 0).unwrap();
    assert_eq!(info.tag1, "Tag1");
    assert_eq!(info.tag2, "Tag2");
    assert_eq!(info.between, "ATATAT");
}

#[test]
fn test_anchored_combined_primer_and_tag_mismatch() {
    use dame::sort::get_pieces_info_anchored;
    // tag1 AAAA->AAAG (1) and forward primer ACGT->ACGA (1): needs m>=1 AND mt>=1.
    let tags = make_test_tags();
    let primers = make_test_primers();
    let line = "AAAGACGAATATATTGCAGGGG";
    assert!(get_pieces_info_anchored(line, &primers, &tags, false, 0, 1).is_none());
    assert!(get_pieces_info_anchored(line, &primers, &tags, false, 1, 0).is_none());
    let info = get_pieces_info_anchored(line, &primers, &tags, false, 1, 1).unwrap();
    assert_eq!(info.tag1, "Tag1");
    assert_eq!(info.tag2, "Tag2");
    assert_eq!(info.between, "ATATAT");
}

#[test]
fn test_anchored_reverse_orientation_tag_mismatch() {
    use dame::sort::{get_pieces_info_anchored, read_tags, read_primers};
    // Non-palindromic primers F=AAAG R=CCGT so orientation is unambiguous.
    // Reverse read: CCCC|CCGT|GGGGGG|CTTT|TTTT, with end tag TTTT miscalled TTTA.
    // tag2-role tag (from the read end, RC side) within 1 mismatch.
    let dir = tempdir().unwrap();
    let tf = dir.path().join("t.txt");
    std::fs::write(&tf, "AAAA\tTag1\nCCCC\tTag2\nGGGG\tTag3\nTTTT\tTag4\n").unwrap();
    let pf = dir.path().join("p.txt");
    std::fs::write(&pf, "CO1\tAAAG\tCCGT\n").unwrap();
    let tags = read_tags(tf.to_str().unwrap()).unwrap();
    let primers = read_primers(pf.to_str().unwrap()).unwrap();
    let line = "CCCCCCGAGGGGGGCTTTTTTA"; // last base T->A in the trailing TTTT tag
    assert!(get_pieces_info_anchored(line, &primers, &tags, false, 0, 0).is_none());
    let info = get_pieces_info_anchored(line, &primers, &tags, false, 1, 1).unwrap();
    assert_eq!(info.tag1, "Tag1");
    assert_eq!(info.tag2, "Tag2");
    assert_eq!(info.between, "CCCCCC");
    assert_eq!(info.primer_name, "CO1");
}

#[test]
fn test_anchored_ambiguous_tag_is_discarded() {
    use dame::sort::{get_pieces_info_anchored, read_tags, read_primers};
    // Two tags equidistant from the read's start tag at mt=1 -> ambiguous -> None.
    // Tags AAAA and AATT are both Hamming-1 from AAAT.
    let dir = tempdir().unwrap();
    let tf = dir.path().join("t.txt");
    std::fs::write(&tf, "AAAA\tTagA\nAATT\tTagB\nGGGG\tTagR\n").unwrap();
    let pf = dir.path().join("p.txt");
    std::fs::write(&pf, "CO1\tACGT\tTGCA\n").unwrap();
    let tags = read_tags(tf.to_str().unwrap()).unwrap();
    let primers = read_primers(pf.to_str().unwrap()).unwrap();
    // start tag region AAAT (1 from AAAA, 1 from AATT); rest exact; end tag GGGG=rc(CCCC)? no.
    // Use end tag rc(TagR=GGGG)=CCCC so tag2 resolves to TagR unambiguously.
    let line = "AAATACGTATATATTGCACCCC";
    let res = get_pieces_info_anchored(line, &primers, &tags, false, 0, 1);
    assert!(res.is_none(), "equidistant tag1 candidates must be discarded as ambiguous");
}
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `cargo test --manifest-path rust/Cargo.toml --test sort_test anchored`
Expected: FAIL — `get_pieces_info_anchored` not found.

- [ ] **Step 3: Implement `get_pieces_info_anchored`**

In `rust/src/sort.rs`, add after `get_pieces_info`:

```rust
/// Anchored matcher used when primer or tag mismatches are allowed.
/// Finds tag candidates at the read ends (IUPAC Hamming <= max_tag_mm), checks
/// primers at the anchored offsets (<= max_primer_mm), scores each valid
/// assembly by total mismatches, returns the unique minimum, and returns None
/// when zero assemblies are valid or the minimum is tied (ambiguous).
pub fn get_pieces_info_anchored(
    line: &str,
    primers: &IndexMap<String, PrimerEntry>,
    tags: &TagLookup,
    keep_primers_seq: bool,
    max_primer_mm: usize,
    max_tag_mm: usize,
) -> Option<PieceInfo> {
    let line_upper = line.to_ascii_uppercase();
    let seq = line_upper.as_bytes();
    let slen = seq.len();

    let mut best: Option<PieceInfo> = None;
    let mut best_score: usize = usize::MAX;
    let mut tied = false;

    for (key, primer) in primers {
        for orientation in 0..2usize {
            // orientation 0 (forward): start=F, end=RC(R)
            // orientation 1 (reverse): start=R, end=RC(F)
            let (start_primer, end_primer) = if orientation == 0 {
                (&primer.start_primers[0], &primer.end_primers[1])
            } else {
                (&primer.start_primers[1], &primer.end_primers[0])
            };

            for (tag1_seq, tag1_name) in &tags.fwd_list {
                let t1l = tag1_seq.len();
                if t1l + start_primer.len() > slen {
                    continue;
                }
                let tag1_mm = hamming_iupac(tag1_seq, &seq[0..t1l]);
                if tag1_mm > max_tag_mm {
                    continue;
                }
                let p_start = t1l;
                let p_end = t1l + start_primer.len();
                let start_mm = hamming_iupac(start_primer, &seq[p_start..p_end]);
                if start_mm > max_primer_mm {
                    continue;
                }

                for (tag2_seq, tag2_name) in &tags.rc_list {
                    let t2l = tag2_seq.len();
                    if t2l + end_primer.len() > slen {
                        continue;
                    }
                    let tag2_start = slen - t2l;
                    let tag2_mm = hamming_iupac(tag2_seq, &seq[tag2_start..]);
                    if tag2_mm > max_tag_mm {
                        continue;
                    }
                    let ep_start = tag2_start - end_primer.len();
                    let ep_end = tag2_start;
                    if ep_start < p_end {
                        continue; // primers overlap / no amplicon room
                    }
                    let end_mm = hamming_iupac(end_primer, &seq[ep_start..ep_end]);
                    if end_mm > max_primer_mm {
                        continue;
                    }

                    let (b_start, b_end) = if keep_primers_seq {
                        (p_start, ep_end)
                    } else {
                        (p_end, ep_start)
                    };
                    if b_start >= b_end {
                        continue;
                    }
                    let between_raw = &line_upper[b_start..b_end];
                    let between = if orientation == 0 {
                        between_raw.to_string()
                    } else {
                        String::from_utf8(rc_bytes(between_raw.as_bytes()))
                            .expect("rc_bytes output is valid ASCII by construction")
                    };

                    // Role assignment matches the exact path: forward keeps
                    // (start->tag1, end->tag2); reverse swaps them.
                    let (tn1, tn2) = if orientation == 0 {
                        (tag1_name.clone(), tag2_name.clone())
                    } else {
                        (tag2_name.clone(), tag1_name.clone())
                    };

                    let score = tag1_mm + start_mm + end_mm + tag2_mm;
                    if score < best_score {
                        best_score = score;
                        best = Some(PieceInfo {
                            tag1: tn1,
                            tag2: tn2,
                            primer_name: key.clone(),
                            between,
                        });
                        tied = false;
                    } else if score == best_score {
                        tied = true;
                    }
                }
            }
        }
    }

    if tied {
        return None;
    }
    best
}
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `cargo test --manifest-path rust/Cargo.toml --test sort_test`
Expected: PASS (all existing + 5 new anchored tests).

- [ ] **Step 5: Commit**

```bash
git add rust/src/sort.rs rust/tests/sort_test.rs
git commit -m "feat(rust): anchored sort matcher with tag + primer mismatch tolerance

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

## Task 3: Rust tag-distance warning + CLI flag + dispatch

**Files:**
- Modify: `rust/src/sort.rs`
- Test: `rust/tests/sort_test.rs`

- [ ] **Step 1: Write failing test for `min_equal_length_tag_distance`**

In `rust/tests/sort_test.rs` add:

```rust
#[test]
fn test_min_equal_length_tag_distance() {
    use dame::sort::{min_equal_length_tag_distance, read_tags};
    let dir = tempdir().unwrap();
    let tf = dir.path().join("t.txt");
    // AAAA vs AATT = 2; AAAA vs AAAT = 1; AATT vs AAAT = 1 -> min 1
    std::fs::write(&tf, "AAAA\tT1\nAATT\tT2\nAAAT\tT3\n").unwrap();
    let tags = read_tags(tf.to_str().unwrap()).unwrap();
    assert_eq!(min_equal_length_tag_distance(&tags), Some(1));

    // single tag -> no pair -> None
    let tf2 = dir.path().join("t2.txt");
    std::fs::write(&tf2, "AAAA\tT1\n").unwrap();
    let tags2 = read_tags(tf2.to_str().unwrap()).unwrap();
    assert_eq!(min_equal_length_tag_distance(&tags2), None);
}
```

- [ ] **Step 2: Run to verify it fails**

Run: `cargo test --manifest-path rust/Cargo.toml --test sort_test min_equal_length`
Expected: FAIL — function not found.

- [ ] **Step 3: Implement `min_equal_length_tag_distance`**

In `rust/src/sort.rs`, add after `read_tags`:

```rust
/// Minimum pairwise Hamming distance among equal-length forward tag sequences.
/// Returns None if no two tags share a length.
pub fn min_equal_length_tag_distance(tags: &TagLookup) -> Option<usize> {
    let seqs: Vec<&Vec<u8>> = tags.fwd_list.iter().map(|(s, _)| s).collect();
    let mut best: Option<usize> = None;
    for i in 0..seqs.len() {
        for j in (i + 1)..seqs.len() {
            if seqs[i].len() != seqs[j].len() {
                continue;
            }
            let d = seqs[i]
                .iter()
                .zip(seqs[j].iter())
                .filter(|(a, b)| a != b)
                .count();
            best = Some(best.map_or(d, |b| b.min(d)));
        }
    }
    best
}
```

- [ ] **Step 4: Run to verify it passes**

Run: `cargo test --manifest-path rust/Cargo.toml --test sort_test min_equal_length`
Expected: PASS.

- [ ] **Step 5: Add `tag_mismatches` to `SortArgs`**

In `rust/src/sort.rs`, in `SortArgs`, add after `primer_mismatches`:

```rust
    #[arg(long = "tag-mismatches", default_value = "0")]
    pub tag_mismatches: usize,
```

- [ ] **Step 6: Update `SortArgs` literals in tests**

In `rust/tests/sort_test.rs`, the `test_run_sort_produces_output_files` builds `SortArgs { ... }`. Add `tag_mismatches: 0,` after `primer_mismatches: 0,`.

- [ ] **Step 7: Dispatch + warning in `run`**

In `rust/src/sort.rs` `run`, replace the matcher dispatch block. Currently:

```rust
        match get_pieces_info(seq, &primers, &tags, args.keep_primers_seq, args.primer_mismatches) {
```

Replace the lines from `let mut count_errors` setup through the `match` so that, after `let primers = read_primers(...)?;` and before the read loop, the warning is emitted; and the per-read `match` dispatches. Concretely, insert right after `let primers = read_primers(&args.primers)?;`:

```rust
    if args.tag_mismatches > 0 {
        if let Some(min_d) = min_equal_length_tag_distance(&tags) {
            if 2 * args.tag_mismatches + 1 > min_d {
                let safe = (min_d.saturating_sub(1)) / 2;
                eprintln!(
                    "WARNING: --tag-mismatches {} may misassign reads: minimum tag \
                     Hamming distance is {}, so the safe maximum is {}.",
                    args.tag_mismatches, min_d, safe
                );
            }
        }
    }
    let use_anchored = args.primer_mismatches > 0 || args.tag_mismatches > 0;
```

Then change the per-read `match` expression to:

```rust
        let info = if use_anchored {
            get_pieces_info_anchored(
                seq,
                &primers,
                &tags,
                args.keep_primers_seq,
                args.primer_mismatches,
                args.tag_mismatches,
            )
        } else {
            get_pieces_info(seq, &primers, &tags, args.keep_primers_seq, args.primer_mismatches)
        };
        match info {
```

(The `Some(info) => { ... } None => { count_errors += 1; }` arms stay unchanged.)

- [ ] **Step 8: Run the full Rust suite**

Run: `cargo test --manifest-path rust/Cargo.toml`
Expected: PASS (all targets).

- [ ] **Step 9: Commit**

```bash
git add rust/src/sort.rs rust/tests/sort_test.rs
git commit -m "feat(rust): --tag-mismatches flag, dispatch, and tag-distance warning

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

## Task 4: Python `hamming_iupac` + `min_equal_length_tag_distance`

**Files:**
- Modify: `python/dame/modules_sort.py`
- Test: `python/tests/test_sort.py`

- [ ] **Step 1: Write failing tests**

In `python/tests/test_sort.py`, update the import line to add the new names:

```python
from dame.modules_sort import (
    RC, readTags, readPrimers, FillHAP, GetPiecesInfo, iupac_matches, find_primer,
    hamming_iupac, GetPiecesInfoMismatch, min_equal_length_tag_distance,
)
```

Add tests at the end:

```python
def test_hamming_iupac():
    assert hamming_iupac("ACGT", "ACGT") == 0
    assert hamming_iupac("ACGT", "ACGA") == 1
    assert hamming_iupac("GCRTGC", "GCATGC") == 0   # R matches A
    assert hamming_iupac("GCRTGC", "GCCTGC") == 1
    # length mismatch -> sentinel greater than any real budget
    assert hamming_iupac("ACGT", "ACG") > 3


def test_min_equal_length_tag_distance():
    assert min_equal_length_tag_distance({"T1": ["AAAA", "TTTT"],
                                          "T2": ["AATT", "AATT"],
                                          "T3": ["AAAT", "ATTT"]}) == 1
    assert min_equal_length_tag_distance({"T1": ["AAAA", "TTTT"]}) is None
```

- [ ] **Step 2: Run to verify they fail**

Run: `pytest python/tests/test_sort.py -k "hamming_iupac or min_equal_length" -v`
Expected: FAIL — ImportError.

- [ ] **Step 3: Implement both functions**

In `python/dame/modules_sort.py`, add after `find_primer`:

```python
def hamming_iupac(pattern, region):
    """Count positions where region fails pattern's IUPAC constraint.
    Returns a sentinel greater than len(pattern) when lengths differ."""
    if len(pattern) != len(region):
        return len(pattern) + 1
    return sum(0 if iupac_matches(p, r) else 1 for p, r in zip(pattern, region))


def min_equal_length_tag_distance(TAGS):
    """Minimum pairwise Hamming distance among equal-length forward tag
    sequences (TAGS[name][0]). Returns None if no two tags share a length."""
    seqs = [TAGS[t][0] for t in TAGS]
    best = None
    for i in range(len(seqs)):
        for j in range(i + 1, len(seqs)):
            a, b = seqs[i], seqs[j]
            if len(a) != len(b):
                continue
            d = sum(1 for x, y in zip(a, b) if x != y)
            best = d if best is None else min(best, d)
    return best
```

- [ ] **Step 4: Run to verify they pass**

Run: `pytest python/tests/test_sort.py -k "hamming_iupac or min_equal_length" -v`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add python/dame/modules_sort.py python/tests/test_sort.py
git commit -m "feat(py): hamming_iupac + min tag-distance helper

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

## Task 5: Python anchored matcher `GetPiecesInfoMismatch`

**Files:**
- Modify: `python/dame/modules_sort.py`
- Test: `python/tests/test_sort.py`

Context: `PRIMERS[name] = [[F, R], [Frc, Rrc]]`. `TAGS[name] = [fwd_seq, rc_seq]`. `_mismatch_fixtures(tmp_path)` (already in the test file) returns `(PRIMERS, TAGS)` for CO1 F=ACGT R=TGCA and tags AAAA=Tag1/CCCC=Tag2/GGGG=Tag3/TTTT=Tag4.

- [ ] **Step 1: Write failing tests**

In `python/tests/test_sort.py`, add at the end:

```python
def test_anchored_forward_tag_mismatch(tmp_path):
    PRIMERS, TAGS = _mismatch_fixtures(tmp_path)
    line = "AAAGACGTATATATTGCAGGGG"  # tag1 AAAA->AAAG
    assert GetPiecesInfoMismatch(line, PRIMERS, TAGS, False, 0, 0) == [1]
    info = GetPiecesInfoMismatch(line, PRIMERS, TAGS, False, 0, 1)
    assert info == ["Tag1", "Tag2", "CO1", "ATATAT"]


def test_anchored_primer_only_parity(tmp_path):
    PRIMERS, TAGS = _mismatch_fixtures(tmp_path)
    line = "AAAAACGAATATATTGCAGGGG"  # primer ACGT->ACGA
    assert GetPiecesInfoMismatch(line, PRIMERS, TAGS, False, 0, 0) == [1]
    info = GetPiecesInfoMismatch(line, PRIMERS, TAGS, False, 1, 0)
    assert info == ["Tag1", "Tag2", "CO1", "ATATAT"]


def test_anchored_combined(tmp_path):
    PRIMERS, TAGS = _mismatch_fixtures(tmp_path)
    line = "AAAGACGAATATATTGCAGGGG"  # tag1 and primer each 1 off
    assert GetPiecesInfoMismatch(line, PRIMERS, TAGS, False, 0, 1) == [1]
    assert GetPiecesInfoMismatch(line, PRIMERS, TAGS, False, 1, 0) == [1]
    assert GetPiecesInfoMismatch(line, PRIMERS, TAGS, False, 1, 1) == ["Tag1", "Tag2", "CO1", "ATATAT"]


def test_anchored_reverse_orientation_tag_mismatch(tmp_path):
    tags_file = tmp_path / "t.txt"
    tags_file.write_text("AAAA\tTag1\nCCCC\tTag2\nGGGG\tTag3\nTTTT\tTag4\n")
    primers_file = tmp_path / "p.txt"
    primers_file.write_text("CO1\tAAAG\tCCGT\n")
    TAGS = readTags(str(tags_file), {})
    PRIMERS = readPrimers(str(primers_file), {})
    line = "CCCCCCGAGGGGGGCTTTTTTA"  # trailing tag TTTT->TTTA
    assert GetPiecesInfoMismatch(line, PRIMERS, TAGS, False, 0, 0) == [1]
    info = GetPiecesInfoMismatch(line, PRIMERS, TAGS, False, 1, 1)
    assert info == ["Tag1", "Tag2", "CO1", "CCCCCC"]


def test_anchored_ambiguous_tag_is_discarded(tmp_path):
    tags_file = tmp_path / "t.txt"
    tags_file.write_text("AAAA\tTagA\nAATT\tTagB\nGGGG\tTagR\n")
    primers_file = tmp_path / "p.txt"
    primers_file.write_text("CO1\tACGT\tTGCA\n")
    TAGS = readTags(str(tags_file), {})
    PRIMERS = readPrimers(str(primers_file), {})
    line = "AAATACGTATATATTGCACCCC"  # AAAT is Hamming-1 from both AAAA and AATT
    assert GetPiecesInfoMismatch(line, PRIMERS, TAGS, False, 0, 1) == [1]
```

- [ ] **Step 2: Run to verify they fail**

Run: `pytest python/tests/test_sort.py -k "anchored" -v`
Expected: FAIL — `GetPiecesInfoMismatch` not found.

- [ ] **Step 3: Implement `GetPiecesInfoMismatch`**

In `python/dame/modules_sort.py`, add after `GetPiecesInfo`:

```python
def GetPiecesInfoMismatch(line, PRIMERS, TAGS, keepPrimersSeq, max_primer_mm, max_tag_mm):
    """Anchored matcher: tolerate <= max_primer_mm per primer and <= max_tag_mm
    per tag. Returns [tag1, tag2, primer, between] for the unique best-scoring
    assembly, or [1] on no match / ambiguous tie. Mirrors the Rust
    get_pieces_info_anchored."""
    line = line.upper()
    slen = len(line)
    best = None
    best_score = None
    tied = False

    for key in PRIMERS:
        F, R = PRIMERS[key][0][0], PRIMERS[key][0][1]
        Frc, Rrc = PRIMERS[key][1][0], PRIMERS[key][1][1]
        for orientation in (0, 1):
            if orientation == 0:
                start_primer, end_primer = F, Rrc
            else:
                start_primer, end_primer = R, Frc

            for tag1_name in TAGS:
                tag1_seq = TAGS[tag1_name][0]
                t1l = len(tag1_seq)
                if t1l + len(start_primer) > slen:
                    continue
                tag1_mm = hamming_iupac(tag1_seq, line[0:t1l])
                if tag1_mm > max_tag_mm:
                    continue
                p_start = t1l
                p_end = t1l + len(start_primer)
                start_mm = hamming_iupac(start_primer, line[p_start:p_end])
                if start_mm > max_primer_mm:
                    continue

                for tag2_name in TAGS:
                    tag2_seq = TAGS[tag2_name][1]
                    t2l = len(tag2_seq)
                    if t2l + len(end_primer) > slen:
                        continue
                    tag2_start = slen - t2l
                    tag2_mm = hamming_iupac(tag2_seq, line[tag2_start:])
                    if tag2_mm > max_tag_mm:
                        continue
                    ep_start = tag2_start - len(end_primer)
                    ep_end = tag2_start
                    if ep_start < p_end:
                        continue
                    end_mm = hamming_iupac(end_primer, line[ep_start:ep_end])
                    if end_mm > max_primer_mm:
                        continue

                    if keepPrimersSeq:
                        b_start, b_end = p_start, ep_end
                    else:
                        b_start, b_end = p_end, ep_start
                    if b_start >= b_end:
                        continue
                    between = line[b_start:b_end]
                    if orientation == 1:
                        between = RC(between)

                    if orientation == 0:
                        tn1, tn2 = tag1_name, tag2_name
                    else:
                        tn1, tn2 = tag2_name, tag1_name

                    score = tag1_mm + start_mm + end_mm + tag2_mm
                    if best_score is None or score < best_score:
                        best_score = score
                        best = [tn1, tn2, key, between]
                        tied = False
                    elif score == best_score:
                        tied = True

    if tied or best is None:
        return [1]
    return best
```

- [ ] **Step 4: Run to verify they pass**

Run: `pytest python/tests/test_sort.py -k "anchored" -v`
Expected: PASS (5 new tests).

- [ ] **Step 5: Commit**

```bash
git add python/dame/modules_sort.py python/tests/test_sort.py
git commit -m "feat(py): anchored sort matcher with tag + primer mismatch tolerance

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

## Task 6: Python CLI flag + dispatch + warning

**Files:**
- Modify: `python/dame/sort.py`
- Test: `python/tests/test_sort.py`

- [ ] **Step 1: Write a failing dispatch/parse smoke test**

In `python/tests/test_sort.py`, add:

```python
def test_sort_argparser_accepts_tag_mismatches():
    import argparse
    import dame.sort as sortmod
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers()
    sortmod.register_subcommand(sub)
    args = parser.parse_args(["sort", "-fq", "x", "-p", "y", "-t", "z", "-mt", "1"])
    assert args.tag_mismatches == 1
    args2 = parser.parse_args(["sort", "--fq", "x", "--primers", "y", "--tags", "z",
                               "--tag-mismatches", "2"])
    assert args2.tag_mismatches == 2
```

- [ ] **Step 2: Run to verify it fails**

Run: `pytest python/tests/test_sort.py -k tag_mismatches -v`
Expected: FAIL — unrecognized arguments `-mt`.

- [ ] **Step 3: Add the flag, imports, and dispatch in `sort.py`**

In `python/dame/sort.py`, update the import block:

```python
from dame.modules_sort import (
    readTags, readPrimers, GetPiecesInfo, GetPiecesInfoMismatch, FillHAP,
    PrintSortedCollapsedCountedSeqs, PrintSummaryFile, min_equal_length_tag_distance,
)
from dame.utils import smart_open
```

Add `import sys` at the very top of the file (needed for the stderr warning):

```python
import sys
```

In `register_subcommand`, add after the `-m`/`--primer-mismatches` argument:

```python
    p.add_argument("-mt", "--tag-mismatches", dest="tag_mismatches", type=int, default=0,
                   help="Max substitutions allowed per tag (tag1 and tag2 independently) [default 0]")
```

In `run`, after `PRIMERS = readPrimers(args.p, PRIMERS)`, add the warning and dispatch flag:

```python
    if args.tag_mismatches > 0:
        min_d = min_equal_length_tag_distance(TAGS)
        if min_d is not None and 2 * args.tag_mismatches + 1 > min_d:
            safe = (min_d - 1) // 2
            print(f"WARNING: --tag-mismatches {args.tag_mismatches} may misassign reads: "
                  f"minimum tag Hamming distance is {min_d}, so the safe maximum is {safe}.",
                  file=sys.stderr)
    use_anchored = args.primer_mismatches > 0 or args.tag_mismatches > 0
```

Replace the `Info = GetPiecesInfo(...)` call in the read loop with:

```python
            if use_anchored:
                Info = GetPiecesInfoMismatch(line, PRIMERS, TAGS, args.keepPrimersSeq,
                                             args.primer_mismatches, args.tag_mismatches)
            else:
                Info = GetPiecesInfo(line, PRIMERS, TAGS, args.keepPrimersSeq, args.primer_mismatches)
```

- [ ] **Step 4: Run the full Python suite**

Run: `pytest python/tests/ -v`
Expected: PASS (all, including the new parse test).

- [ ] **Step 5: Commit**

```bash
git add python/dame/sort.py python/tests/test_sort.py
git commit -m "feat(py): --tag-mismatches flag, dispatch, and tag-distance warning

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

## Task 7: Cross-implementation integration test + CI

**Files:**
- Create: `tests/fixtures/sample_tag_err.fastq`
- Create: `tests/integration/run_sort_tag_mismatch.sh`
- Modify: `.github/workflows/ci.yml`

Read structure: tag1 `AAAA`(Tag1) → miscalled `AAAG`; primer `ACGT`; barcode `ATATATATATAT`; `TGCA`; tag2 `GGGG`(=RC(CCCC)=Tag2). Fixtures `Primers.txt` (`CO1<TAB>ACGT<TAB>TGCA`) and `Tags.txt` (AAAA=Tag1…) already exist.

- [ ] **Step 1: Create the fixture**

Create `tests/fixtures/sample_tag_err.fastq` with exactly:

```
@errtag1
AAAGACGTATATATATATATTGCAGGGG
+
IIIIIIIIIIIIIIIIIIIIIIIIIIII
@errtag2
AAAGACGTATATATATATATTGCAGGGG
+
IIIIIIIIIIIIIIIIIIIIIIIIIIII
```

- [ ] **Step 2: Create the integration script**

Create `tests/integration/run_sort_tag_mismatch.sh`:

```bash
#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
FIXTURES="$REPO_ROOT/tests/fixtures"
DAME_BIN="$REPO_ROOT/rust/target/release/dame"

if [ ! -f "$DAME_BIN" ]; then
    echo "SKIP: dame binary not found at $DAME_BIN (run: cd rust && cargo build --release)"
    exit 0
fi
if ! command -v dame-py >/dev/null 2>&1; then
    echo "SKIP: dame-py not found on PATH (run: pip install -e python/)"
    exit 0
fi

TMPPY=$(mktemp -d)
TMPRS=$(mktemp -d)
TMPN0=$(mktemp -d)
trap "rm -rf '$TMPPY' '$TMPRS' '$TMPN0'" EXIT

echo "==> dame-py sort -mt 1..."
cd "$TMPPY"
dame-py sort -fq "$FIXTURES/sample_tag_err.fastq" -p "$FIXTURES/Primers.txt" -t "$FIXTURES/Tags.txt" -mt 1

echo "==> dame sort --tag-mismatches 1..."
cd "$TMPRS"
"$DAME_BIN" sort --fq "$FIXTURES/sample_tag_err.fastq" --primers "$FIXTURES/Primers.txt" --tags "$FIXTURES/Tags.txt" --tag-mismatches 1

for D in "$TMPPY" "$TMPRS"; do
    if [ ! -f "$D/Tag1_Tag2.txt" ]; then
        echo "FAIL: $D did not produce Tag1_Tag2.txt at tag-mismatches=1"; exit 1
    fi
done

echo "==> Comparing Tag1_Tag2.txt (py vs rust)..."
if ! diff "$TMPPY/Tag1_Tag2.txt" "$TMPRS/Tag1_Tag2.txt"; then
    echo "FAIL: Tag1_Tag2.txt differs between dame-py and dame at tag-mismatches=1"; exit 1
fi
if ! grep -q "ATATATATATAT" "$TMPRS/Tag1_Tag2.txt"; then
    echo "FAIL: expected recovered barcode ATATATATATAT not found"; exit 1
fi

echo "==> Confirming rejection at tag-mismatches=0 (both tools)..."
cd "$TMPN0"
"$DAME_BIN" sort --fq "$FIXTURES/sample_tag_err.fastq" --primers "$FIXTURES/Primers.txt" --tags "$FIXTURES/Tags.txt"
if [ -f "$TMPN0/Tag1_Tag2.txt" ]; then
    echo "FAIL: rust read should NOT be recovered at tag-mismatches=0"; exit 1
fi
TMPN0PY=$(mktemp -d)
cd "$TMPN0PY"
dame-py sort -fq "$FIXTURES/sample_tag_err.fastq" -p "$FIXTURES/Primers.txt" -t "$FIXTURES/Tags.txt"
if [ -f "$TMPN0PY/Tag1_Tag2.txt" ]; then
    echo "FAIL: dame-py read should NOT be recovered at tag-mismatches=0"; rm -rf "$TMPN0PY"; exit 1
fi
rm -rf "$TMPN0PY"

echo "PASS: dame and dame-py sort agree at tag-mismatches=1 and reject at 0"
```

- [ ] **Step 3: Make it executable and build the binary**

```bash
chmod +x tests/integration/run_sort_tag_mismatch.sh
cargo build --release --manifest-path rust/Cargo.toml
```

- [ ] **Step 4: Run it**

Run: `bash tests/integration/run_sort_tag_mismatch.sh` (with `dame-py` installed)
Expected final line: `PASS: dame and dame-py sort agree at tag-mismatches=1 and reject at 0`
If outputs diverge, report BLOCKED with the diff — do NOT edit source to force a pass; a divergence means the two anchored matchers differ and must be reconciled.

- [ ] **Step 5: Wire into CI**

In `.github/workflows/ci.yml`, in `integration-tests`, after the `Run sort mismatch integration test` step, add:

```yaml
      - name: Run sort tag-mismatch integration test
        run: bash tests/integration/run_sort_tag_mismatch.sh
```

- [ ] **Step 6: Commit**

```bash
git add tests/fixtures/sample_tag_err.fastq tests/integration/run_sort_tag_mismatch.sh .github/workflows/ci.yml
git commit -m "test: cross-impl tag-mismatch integration test + CI

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

## Task 8: Docs + version bump

**Files:**
- Modify: `rust/Cargo.toml`, `README.md`, `tutorial/README.md`

- [ ] **Step 1: Bump crate version**

In `rust/Cargo.toml`, change `version = "0.5.0"` to `version = "0.6.0"`.

- [ ] **Step 2: README**

In `README.md`:

(a) Change the title heading `# DAMe v2.4: DNA Metabarcoding toolkit` to `# DAMe v2.5: DNA Metabarcoding toolkit`.

(b) In the "Pipeline overview" section, update the `dame sort` line to:

```
dame sort     -fq POOL.fastq --primers P.txt --tags T.txt [--primer-mismatches N] [--tag-mismatches N]
```

(c) Add a development-history entry after the existing item 10:

```markdown
11. **DAMe v2.5 — Tag mismatches + anchored matching.**  `sort` gained
    `--tag-mismatches N` (`-mt` in Python; default 0), a per-tag substitution
    tolerance.  When `--primer-mismatches` or `--tag-mismatches` is non-zero,
    sort uses a tag-anchored matcher: it finds tag candidates at the read ends
    by IUPAC Hamming distance, checks the primers at the expected offsets,
    scores each valid assembly by total mismatches, keeps the unique minimum,
    and discards ambiguous ties.  At the defaults (both 0) the original exact
    matcher runs unchanged (byte-identical).  A startup warning flags an unsafe
    `--tag-mismatches` relative to the tag set's minimum Hamming distance.
    Design adapted from a community PR (dougwyu/DAMe#1, @jiyinqiu).
```

- [ ] **Step 3: tutorial/README.md**

In `tutorial/README.md`, in the "Sort Options" table (where `--primer-mismatches` is documented), add a row:

```
| `--tag-mismatches N` | `-mt N` | Allow up to N substitutions per tag (tag1 and tag2 independently, IUPAC-aware). Triggers the anchored matcher; ambiguous reads are discarded. | 0 |
```

- [ ] **Step 4: Verify build + tests**

Run: `cargo build --release --manifest-path rust/Cargo.toml && cargo test --manifest-path rust/Cargo.toml && pytest python/tests/ -v`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add rust/Cargo.toml rust/Cargo.lock README.md tutorial/README.md
git commit -m "docs: document --tag-mismatches; bump to v2.5 (crate 0.6.0)

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

## Final verification

```bash
cargo test --manifest-path rust/Cargo.toml
pip install -e python/ && pytest python/tests/ -v
cargo build --release --manifest-path rust/Cargo.toml
export PATH="$PWD/rust/target/release:$PATH"
bash tests/integration/run_sort.sh            # m=0,mt=0 regression (byte-identical)
bash tests/integration/run_sort_mismatch.sh   # m=1 via anchored path
bash tests/integration/run_sort_tag_mismatch.sh  # mt=1
bash tests/integration/run_pipeline.sh
```

Expected: all PASS.
