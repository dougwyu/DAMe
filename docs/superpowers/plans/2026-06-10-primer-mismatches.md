# Primer Mismatches (sort) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a configurable, per-primer substitution tolerance (`-m`/`--primer-mismatches`, default 0) to the `sort` step in both the Rust and Python implementations, so reads with sequencing errors inside the primer region are no longer discarded.

**Architecture:** Both implementations gain a `find_primer(primer, seq, max_mismatches)` that scans left→right and returns the first window whose IUPAC-mismatch count is ≤ N (leftmost within budget). The budget is threaded through `get_pieces_info`/`GetPiecesInfo` to each of the four primer matches independently. N=0 reduces exactly to today's exact-match behavior, guaranteed byte-identical by the existing integration tests. Tags stay exact; only substitutions are tolerated.

**Tech Stack:** Rust (`clap`, `needletail`, `ahash`, `indexmap`), Python 3 (stdlib only — the `re` dependency is removed from `modules_sort.py`), pytest, `cargo test`, bash integration scripts, GitHub Actions.

**Reference spec:** `docs/superpowers/specs/2026-06-10-primer-mismatches-design.md`

---

## File Structure

- `rust/src/sort.rs` — `find_primer` gains a `max_mismatches` param; `get_pieces_info` and `SortArgs` thread it through.
- `rust/tests/sort_test.rs` — existing `find_primer`/`get_pieces_info`/`run` tests updated to the new signatures; new mismatch tests added.
- `python/dame/modules_sort.py` — add `iupac_matches` + `find_primer`; rewrite `readPrimers` (drop regex/`AMBIG`) and `GetPiecesInfo` (use `find_primer`, add `maxMismatches`); remove `import re`.
- `python/dame/sort.py` — add the `-m`/`--primer-mismatches` CLI arg; drop the `AMBIG` dict; update the `readPrimers`/`GetPiecesInfo` call sites.
- `python/tests/test_sort.py` — update `test_readPrimers` to the new signature; add `find_primer`/`iupac_matches`/`GetPiecesInfo` mismatch tests.
- `tests/fixtures/sample_primer_err.fastq` — new fixture: reads that only match with 1 primer mismatch.
- `tests/integration/run_sort_mismatch.sh` — new cross-implementation diff at `--primer-mismatches 1`.
- `.github/workflows/ci.yml` — run the new integration script.
- `rust/Cargo.toml`, `README.md`, `tutorial/README.md` — version bump + docs.

---

## Task 1: Rust `find_primer` gains a mismatch budget

**Files:**
- Modify: `rust/src/sort.rs` (the `find_primer` fn and its 4 call sites inside `get_pieces_info`)
- Test: `rust/tests/sort_test.rs`

- [ ] **Step 1: Update existing `find_primer` tests to the new 3-arg signature and add mismatch tests**

In `rust/tests/sort_test.rs`, every existing `find_primer` test calls `find_primer(primer, seq)`. Add `, 0` to each. The six existing tests become:

```rust
#[test]
fn test_find_primer_exact() {
    use dame::sort::find_primer;
    let seq = b"XXXXACGTXXXX";
    let primer = b"ACGT";
    assert_eq!(find_primer(primer, seq, 0), Some((4, 8)));
}

#[test]
fn test_find_primer_iupac() {
    use dame::sort::find_primer;
    let seq = b"TTTTGCATGCTTTT";
    let primer = b"GCRTGC";
    assert_eq!(find_primer(primer, seq, 0), Some((4, 10)));
}

#[test]
fn test_find_primer_iupac_second_option() {
    use dame::sort::find_primer;
    let seq = b"TTTTGCGTGCTTTT";
    let primer = b"GCRTGC";
    assert_eq!(find_primer(primer, seq, 0), Some((4, 10)));
}

#[test]
fn test_find_primer_not_found() {
    use dame::sort::find_primer;
    let seq = b"AAAAAAAAAA";
    let primer = b"GCATGC";
    assert_eq!(find_primer(primer, seq, 0), None);
}

#[test]
fn test_find_primer_leftmost() {
    use dame::sort::find_primer;
    let seq = b"ACGTXXXXACGT";
    let primer = b"ACGT";
    assert_eq!(find_primer(primer, seq, 0), Some((0, 4)));
}

#[test]
fn test_find_primer_primer_longer_than_seq() {
    use dame::sort::find_primer;
    let seq = b"AC";
    let primer = b"ACGT";
    assert_eq!(find_primer(primer, seq, 0), None);
}
```

Then add these new tests directly after them:

```rust
#[test]
fn test_find_primer_one_mismatch_found() {
    use dame::sort::find_primer;
    // ACGA differs from primer ACGT at the last base; N=1 accepts it.
    let seq = b"XXXXACGAXXXX";
    let primer = b"ACGT";
    assert_eq!(find_primer(primer, seq, 0), None);
    assert_eq!(find_primer(primer, seq, 1), Some((4, 8)));
}

#[test]
fn test_find_primer_two_mismatches_rejected_at_one() {
    use dame::sort::find_primer;
    // AAGT differs from ACGT at two positions; N=1 rejects, N=2 accepts.
    let seq = b"XXXXAAGAXXXX";
    let primer = b"ACGT";
    assert_eq!(find_primer(primer, seq, 1), None);
    assert_eq!(find_primer(primer, seq, 2), Some((4, 8)));
}

#[test]
fn test_find_primer_iupac_plus_mismatch() {
    use dame::sort::find_primer;
    // GCRTGC: R matches A or G (0 cost). Read GCATGA differs only at last base.
    let seq = b"TTTTGCATGATTTT";
    let primer = b"GCRTGC";
    assert_eq!(find_primer(primer, seq, 0), None);
    assert_eq!(find_primer(primer, seq, 1), Some((4, 10)));
}

#[test]
fn test_find_primer_leftmost_within_budget() {
    use dame::sort::find_primer;
    // GCTTGC (1 mismatch) sits before an exact GCATGC. Leftmost-within-budget
    // returns the earlier near-match, NOT the later exact match.
    // index: T0 T1 G2 C3 T4 T5 G6 C7 A8 T9 G10 C11
    let seq = b"TTGCTTGCATGC";
    let primer = b"GCATGC";
    assert_eq!(find_primer(primer, seq, 1), Some((2, 8)));
}
```

- [ ] **Step 2: Run the new tests to verify they fail to compile**

Run: `cargo test --manifest-path rust/Cargo.toml find_primer`
Expected: FAIL — compile error, `find_primer` takes 2 arguments but 3 were supplied.

- [ ] **Step 3: Add the `max_mismatches` parameter to `find_primer`**

In `rust/src/sort.rs`, replace the whole `find_primer` function (currently lines ~94-108):

```rust
/// Find the leftmost occurrence of `primer` in `seq` using IUPAC matching,
/// tolerating up to `max_mismatches` substitutions (positions where the read
/// base does not satisfy the primer's IUPAC code).
/// Returns the first window within budget as `Some((start, end))`, or `None`.
/// `max_mismatches == 0` is exact IUPAC matching.
pub fn find_primer(primer: &[u8], seq: &[u8], max_mismatches: usize) -> Option<(usize, usize)> {
    let plen = primer.len();
    let slen = seq.len();
    if plen > slen {
        return None;
    }
    for i in 0..=(slen - plen) {
        let mut mismatches = 0usize;
        let mut ok = true;
        for (&p, &s) in primer.iter().zip(&seq[i..i + plen]) {
            if !iupac_matches(p, s) {
                mismatches += 1;
                if mismatches > max_mismatches {
                    ok = false;
                    break;
                }
            }
        }
        if ok {
            return Some((i, i + plen));
        }
    }
    None
}
```

- [ ] **Step 4: Update the 4 `find_primer` call sites in `get_pieces_info` to pass a literal `0` (temporary)**

In `rust/src/sort.rs`, inside `get_pieces_info`, there are exactly four calls. Add `, 0` to each so the crate compiles (the budget is wired to a real parameter in Task 2):

- `find_primer(&primer.start_primers[0], seq)` → `find_primer(&primer.start_primers[0], seq, 0)`
- `find_primer(&primer.end_primers[1], seq)` → `find_primer(&primer.end_primers[1], seq, 0)`
- `find_primer(&primer.start_primers[1], seq)` → `find_primer(&primer.start_primers[1], seq, 0)`
- `find_primer(&primer.end_primers[0], seq)` → `find_primer(&primer.end_primers[0], seq, 0)`

- [ ] **Step 5: Run the tests to verify they pass**

Run: `cargo test --manifest-path rust/Cargo.toml find_primer`
Expected: PASS — all `find_primer` tests (existing + new) green.

- [ ] **Step 6: Commit**

```bash
git add rust/src/sort.rs rust/tests/sort_test.rs
git commit -m "feat(rust): find_primer tolerates up to N substitutions

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

## Task 2: Rust thread the budget through `get_pieces_info`, `SortArgs`, `run`

**Files:**
- Modify: `rust/src/sort.rs` (`get_pieces_info` signature + 4 call sites, `SortArgs`, `run`)
- Test: `rust/tests/sort_test.rs`

- [ ] **Step 1: Update existing `get_pieces_info` / `run` tests and add a mismatch test**

In `rust/tests/sort_test.rs`:

Update the three existing `get_pieces_info` call sites to the new 5-arg signature (add `, 0`):
- `get_pieces_info(line, &primers, &tags, false)` in `test_get_pieces_info_forward_read` → `get_pieces_info(line, &primers, &tags, false, 0)`
- `get_pieces_info(line, &primers, &tags, false)` in `test_get_pieces_info_error_read` → `get_pieces_info(line, &primers, &tags, false, 0)`
- `get_pieces_info(bad_read, &primers, &tags, false)` in `test_get_pieces_info_no_panic_on_inverted_primers` → `get_pieces_info(bad_read, &primers, &tags, false, 0)`

Update the `SortArgs` literal in `test_run_sort_produces_output_files` to add the new field:

```rust
    let result = run(SortArgs {
        fastq: fastq.to_str().unwrap().to_string(),
        primers: primers.to_str().unwrap().to_string(),
        tags: tags.to_str().unwrap().to_string(),
        keep_primers_seq: false,
        primer_mismatches: 0,
    });
```

Add this new test after `test_get_pieces_info_forward_read`:

```rust
#[test]
fn test_get_pieces_info_forward_read_with_one_mismatch() {
    // Same as the forward-read test, but the forward primer ACGT is miscalled
    // as ACGA. At N=0 the read is rejected; at N=1 it sorts to Tag1_Tag2.
    let tags = make_test_tags();
    let primers = make_test_primers();
    let line = "AAAAACGAATATATTGCAGGGG"; // primer region ACGA (was ACGT)
    assert!(get_pieces_info(line, &primers, &tags, false, 0).is_none());
    let info = get_pieces_info(line, &primers, &tags, false, 1).unwrap();
    assert_eq!(info.tag1, "Tag1");
    assert_eq!(info.tag2, "Tag2");
    assert_eq!(info.between, "ATATAT");
    assert_eq!(info.primer_name, "CO1");
}
```

- [ ] **Step 2: Run the tests to verify they fail to compile**

Run: `cargo test --manifest-path rust/Cargo.toml get_pieces_info`
Expected: FAIL — compile error, `get_pieces_info` takes 4 arguments but 5 were supplied; missing field `primer_mismatches` in `SortArgs`.

- [ ] **Step 3: Add `primer_mismatches` to `SortArgs`**

In `rust/src/sort.rs`, in the `SortArgs` struct, add the field after `keep_primers_seq`:

```rust
    #[arg(short = 'm', long = "primer-mismatches", default_value = "0")]
    pub primer_mismatches: usize,
```

- [ ] **Step 4: Add `max_mismatches` to `get_pieces_info` and pass it to the 4 calls**

In `rust/src/sort.rs`, change the `get_pieces_info` signature:

```rust
pub fn get_pieces_info(
    line: &str,
    primers: &IndexMap<String, PrimerEntry>,
    tags: &TagLookup,
    keep_primers_seq: bool,
    max_mismatches: usize,
) -> Option<PieceInfo> {
```

Then change the four `find_primer(..., 0)` calls added in Task 1 to use `max_mismatches`:
- `find_primer(&primer.start_primers[0], seq, 0)` → `find_primer(&primer.start_primers[0], seq, max_mismatches)`
- `find_primer(&primer.end_primers[1], seq, 0)` → `find_primer(&primer.end_primers[1], seq, max_mismatches)`
- `find_primer(&primer.start_primers[1], seq, 0)` → `find_primer(&primer.start_primers[1], seq, max_mismatches)`
- `find_primer(&primer.end_primers[0], seq, 0)` → `find_primer(&primer.end_primers[0], seq, max_mismatches)`

- [ ] **Step 5: Pass the arg in `run`**

In `rust/src/sort.rs`, in `run`, the call is currently:

```rust
        match get_pieces_info(seq, &primers, &tags, args.keep_primers_seq) {
```

Change it to:

```rust
        match get_pieces_info(seq, &primers, &tags, args.keep_primers_seq, args.primer_mismatches) {
```

- [ ] **Step 6: Run the full Rust test suite**

Run: `cargo test --manifest-path rust/Cargo.toml`
Expected: PASS — all unit + integration tests green (including `test_run_sort_produces_output_files`, the N=0 regression).

- [ ] **Step 7: Commit**

```bash
git add rust/src/sort.rs rust/tests/sort_test.rs
git commit -m "feat(rust): thread --primer-mismatches through sort

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

## Task 3: Python IUPAC matcher + `find_primer`

**Files:**
- Modify: `python/dame/modules_sort.py` (add two functions near the top, after `RC`)
- Test: `python/tests/test_sort.py`

- [ ] **Step 1: Write the failing tests**

In `python/tests/test_sort.py`, update the import line to include the new functions:

```python
from dame.modules_sort import RC, readTags, readPrimers, FillHAP, GetPiecesInfo, iupac_matches, find_primer
```

Add these tests at the end of the file:

```python
def test_iupac_matches_exact():
    assert iupac_matches("A", "A")
    assert not iupac_matches("A", "C")


def test_iupac_matches_ambiguous():
    # R = A or G
    assert iupac_matches("R", "A")
    assert iupac_matches("R", "G")
    assert not iupac_matches("R", "C")
    # N = anything
    assert iupac_matches("N", "A")
    assert iupac_matches("N", "T")


def test_iupac_matches_non_acgt_read():
    # A read base of N never satisfies any primer code (matches Rust).
    assert not iupac_matches("N", "N")
    assert not iupac_matches("A", "N")


def test_find_primer_exact():
    assert find_primer("ACGT", "XXXXACGTXXXX", 0) == (4, 8)
    assert find_primer("ACGT", "AAAAAAAA", 0) is None


def test_find_primer_leftmost():
    assert find_primer("ACGT", "ACGTXXXXACGT", 0) == (0, 4)


def test_find_primer_one_mismatch():
    # ACGA differs from ACGT at the last base.
    assert find_primer("ACGT", "XXXXACGAXXXX", 0) is None
    assert find_primer("ACGT", "XXXXACGAXXXX", 1) == (4, 8)


def test_find_primer_two_mismatches_rejected_at_one():
    assert find_primer("ACGT", "XXXXAAGAXXXX", 1) is None
    assert find_primer("ACGT", "XXXXAAGAXXXX", 2) == (4, 8)


def test_find_primer_leftmost_within_budget():
    # GCTTGC (1 mismatch) before exact GCATGC — leftmost-within-budget wins.
    assert find_primer("GCATGC", "TTGCTTGCATGC", 1) == (2, 8)


def test_find_primer_longer_than_seq():
    assert find_primer("ACGT", "AC", 0) is None
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `pytest python/tests/test_sort.py -k "iupac_matches or find_primer" -v`
Expected: FAIL — `ImportError: cannot import name 'iupac_matches'`.

- [ ] **Step 3: Add `iupac_matches` and `find_primer` to `modules_sort.py`**

In `python/dame/modules_sort.py`, add immediately after the `RC` function (after line 9):

```python
_IUPAC = {
    'A': frozenset('A'), 'C': frozenset('C'), 'G': frozenset('G'),
    'T': frozenset('T'), 'R': frozenset('AG'), 'Y': frozenset('CT'),
    'S': frozenset('CG'), 'W': frozenset('AT'), 'K': frozenset('GT'),
    'M': frozenset('AC'), 'B': frozenset('CGT'), 'D': frozenset('AGT'),
    'H': frozenset('ACT'), 'V': frozenset('ACG'), 'N': frozenset('ACGT'),
}


def iupac_matches(primer_base, read_base):
    """True if read_base (A/C/G/T) satisfies the primer's IUPAC code."""
    allowed = _IUPAC.get(primer_base)
    return allowed is not None and read_base in allowed


def find_primer(primer, seq, max_mismatches=0):
    """Leftmost window in seq matching primer with <= max_mismatches
    substitutions (IUPAC-aware). Returns (start, end) or None."""
    plen = len(primer)
    slen = len(seq)
    if plen > slen:
        return None
    for i in range(slen - plen + 1):
        mismatches = 0
        ok = True
        for p, s in zip(primer, seq[i:i + plen]):
            if not iupac_matches(p, s):
                mismatches += 1
                if mismatches > max_mismatches:
                    ok = False
                    break
        if ok:
            return (i, i + plen)
    return None
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `pytest python/tests/test_sort.py -k "iupac_matches or find_primer" -v`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add python/dame/modules_sort.py python/tests/test_sort.py
git commit -m "feat(py): add IUPAC find_primer with mismatch budget

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

## Task 4: Python rewrite `readPrimers` + `GetPiecesInfo`, add CLI flag

**Files:**
- Modify: `python/dame/modules_sort.py` (`readPrimers`, `GetPiecesInfo`, remove `import re`)
- Modify: `python/dame/sort.py` (CLI arg, drop `AMBIG`, update call sites)
- Test: `python/tests/test_sort.py` (`test_readPrimers`, new `GetPiecesInfo` mismatch test)

- [ ] **Step 1: Update `test_readPrimers` and add a `GetPiecesInfo` mismatch test**

In `python/tests/test_sort.py`, replace the existing `test_readPrimers` with the new 2-arg signature (no `AMBIG`, raw primer strings stored):

```python
def test_readPrimers(tmp_path):
    primers_file = tmp_path / "primers.txt"
    primers_file.write_text("CO1\tACGT\tTTTT\n")
    PRIMERS = {}
    result = readPrimers(str(primers_file), PRIMERS)
    assert "CO1" in result
    assert len(result["CO1"]) == 2          # [forward_list, rc_list]
    assert result["CO1"][0] == ["ACGT", "TTTT"]      # F, R (raw)
    assert result["CO1"][1] == ["ACGT", "AAAA"]      # RC(F), RC(R)
```

Add this test at the end of the file:

```python
def _mismatch_fixtures(tmp_path):
    tags_file = tmp_path / "tags.txt"
    tags_file.write_text("AAAA\tTag1\nCCCC\tTag2\nGGGG\tTag3\nTTTT\tTag4\n")
    primers_file = tmp_path / "primers.txt"
    primers_file.write_text("CO1\tACGT\tTGCA\n")
    TAGS = readTags(str(tags_file), {})
    PRIMERS = readPrimers(str(primers_file), {})
    return PRIMERS, TAGS


def test_GetPiecesInfo_one_mismatch(tmp_path):
    PRIMERS, TAGS = _mismatch_fixtures(tmp_path)
    # Forward primer ACGT miscalled as ACGA.
    line = "AAAAACGAATATATTGCAGGGG"
    # N=0 -> rejected (error sentinel [1])
    assert GetPiecesInfo(line, PRIMERS, TAGS, False, 0) == [1]
    # N=1 -> sorts to Tag1_Tag2 with barcode ATATAT
    info = GetPiecesInfo(line, PRIMERS, TAGS, False, 1)
    assert info[0] == "Tag1"
    assert info[1] == "Tag2"
    assert info[2] == "CO1"
    assert info[3] == "ATATAT"
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `pytest python/tests/test_sort.py -k "readPrimers or GetPiecesInfo" -v`
Expected: FAIL — `test_readPrimers` fails (`readPrimers` still requires the `AMBIG` arg / stores regex strings) and `GetPiecesInfo` rejects the 5th positional arg.

- [ ] **Step 3: Rewrite `readPrimers` (drop `AMBIG`) in `modules_sort.py`**

In `python/dame/modules_sort.py`, replace the whole `readPrimers` function (currently lines 25-46):

```python
def readPrimers(primers, PRIMERS):
    with open(primers) as f:
        for line in f:
            line = line.rstrip().split()
            if not line:
                continue
            if line[0] not in PRIMERS:
                PRIMERS[line[0]] = [[], []]
            F = line[1]
            R = line[2]
            Frc = RC(F)
            Rrc = RC(R)
            PRIMERS[line[0]][0].append(F)
            PRIMERS[line[0]][0].append(R)
            PRIMERS[line[0]][1].append(Frc)
            PRIMERS[line[0]][1].append(Rrc)
    return PRIMERS
```

- [ ] **Step 4: Rewrite `GetPiecesInfo` to use `find_primer` with a `maxMismatches` param**

In `python/dame/modules_sort.py`, replace the whole `GetPiecesInfo` function (currently lines 49-109):

```python
def GetPiecesInfo(line, PRIMERS, TAGS, keepPrimersSeq, maxMismatches=0):
    for key in PRIMERS:
        # Forward orientation: F at start, RC(R) at end
        primIni = find_primer(PRIMERS[key][0][0], line, maxMismatches)
        if primIni is not None:
            if keepPrimersSeq:
                primIniPosPrim = primIni[0]
                primIniPosTags = primIni[0]
            else:
                primIniPosPrim = primIni[1]
                primIniPosTags = primIni[0]
            primFin = find_primer(PRIMERS[key][1][1], line, maxMismatches)
            if primFin is not None:
                if keepPrimersSeq:
                    primFinPosPrim = primFin[1]
                    primFinPosTags = primFin[1]
                else:
                    primFinPosPrim = primFin[0]
                    primFinPosTags = primFin[1]
                PrimerName = key
                between = line[primIniPosPrim:primFinPosPrim]
                if len(between) == 0:
                    return [1]
                tag1 = line[:primIniPosTags]
                tag2 = line[primFinPosTags:]
                tagName1 = [t for t in TAGS if TAGS[t][0] == tag1]
                tagName2 = [t for t in TAGS if TAGS[t][1] == tag2]
                if len(tagName1) > 0 and len(tagName2) > 0:
                    return [tagName1[0], tagName2[0], PrimerName, between]
                return [1]
            return [1]
        else:
            # Reverse orientation: R at start, RC(F) at end
            primIni = find_primer(PRIMERS[key][0][1], line, maxMismatches)
            if primIni is not None:
                if keepPrimersSeq:
                    primIniPosPrim = primIni[0]
                    primIniPosTags = primIni[0]
                else:
                    primIniPosPrim = primIni[1]
                    primIniPosTags = primIni[0]
                primFin = find_primer(PRIMERS[key][1][0], line, maxMismatches)
                if primFin is not None:
                    if keepPrimersSeq:
                        primFinPosPrim = primFin[1]
                        primFinPosTags = primFin[1]
                    else:
                        primFinPosPrim = primFin[0]
                        primFinPosTags = primFin[1]
                    PrimerName = key
                    between = line[primIniPosPrim:primFinPosPrim]
                    if len(between) == 0:
                        return [1]
                    between = RC(between)
                    tag1 = line[:primIniPosTags]
                    tag2 = line[primFinPosTags:]
                    tagName2 = [t for t in TAGS if TAGS[t][0] == tag1]
                    tagName1 = [t for t in TAGS if TAGS[t][1] == tag2]
                    if len(tagName1) > 0 and len(tagName2) > 0:
                        return [tagName1[0], tagName2[0], PrimerName, between]
                    return [1]
                return [1]
    return [1]
```

- [ ] **Step 5: Remove the now-unused `import re` in `modules_sort.py`**

In `python/dame/modules_sort.py`, the first line is `import re`. Delete it (the file no longer uses regex; `RC` uses `str.maketrans`).

- [ ] **Step 6: Update `sort.py` — drop `AMBIG`, add the CLI flag, update call sites**

In `python/dame/sort.py`:

Delete the `AMBIG = { ... }` dict (lines 9-13).

In `register_subcommand`, add this argument before `p.set_defaults(func=run)`:

```python
    p.add_argument("-m", "--primer-mismatches", dest="primer_mismatches", type=int, default=0,
                   help="Max substitutions allowed per primer match [default 0]")
```

Change the `readPrimers` call in `run` from:

```python
    PRIMERS = readPrimers(args.p, PRIMERS, AMBIG)
```

to:

```python
    PRIMERS = readPrimers(args.p, PRIMERS)
```

Change the `GetPiecesInfo` call in `run` from:

```python
            Info = GetPiecesInfo(line, PRIMERS, TAGS, args.keepPrimersSeq)
```

to:

```python
            Info = GetPiecesInfo(line, PRIMERS, TAGS, args.keepPrimersSeq, args.primer_mismatches)
```

In the standalone `main()`, add the matching argument before `run(...)`:

```python
    parser.add_argument("-m", "--primer-mismatches", dest="primer_mismatches", type=int, default=0)
```

- [ ] **Step 7: Run the full Python test suite**

Run: `pytest python/tests/ -v`
Expected: PASS — all tests green, including the updated `test_readPrimers` and the new mismatch tests.

- [ ] **Step 8: Commit**

```bash
git add python/dame/modules_sort.py python/dame/sort.py python/tests/test_sort.py
git commit -m "feat(py): sort --primer-mismatches via manual IUPAC matcher

Replaces the regex primer matcher with a byte-level IUPAC sliding window
mirroring the Rust implementation; removes the re dependency.

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

## Task 5: Cross-implementation integration test + CI

**Files:**
- Create: `tests/fixtures/sample_primer_err.fastq`
- Create: `tests/integration/run_sort_mismatch.sh`
- Modify: `.github/workflows/ci.yml`

- [ ] **Step 1: Create the mismatch fixture**

Create `tests/fixtures/sample_primer_err.fastq` with exactly this content (the forward primer `ACGT` is miscalled as `ACGA`; barcode `ATATATATATAT`; tags `AAAA`=Tag1, `GGGG`=RC(`CCCC`)=Tag2). Two identical reads → count 2.

```
@errprimer1
AAAAACGAATATATATATATTGCAGGGG
+
IIIIIIIIIIIIIIIIIIIIIIIIIIII
@errprimer2
AAAAACGAATATATATATATTGCAGGGG
+
IIIIIIIIIIIIIIIIIIIIIIIIIIII
```

- [ ] **Step 2: Create the integration script**

Create `tests/integration/run_sort_mismatch.sh` (mirrors `run_sort.sh`, but runs both tools with one allowed mismatch and asserts the read is recovered):

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

TMPPY=$(mktemp -d)
TMPRS=$(mktemp -d)
TMPN0=$(mktemp -d)
trap "rm -rf '$TMPPY' '$TMPRS' '$TMPN0'" EXIT

echo "==> Running dame-py sort -m 1..."
cd "$TMPPY"
dame-py sort \
    -fq "$FIXTURES/sample_primer_err.fastq" \
    -p  "$FIXTURES/Primers.txt" \
    -t  "$FIXTURES/Tags.txt" \
    -m 1

echo "==> Running dame sort --primer-mismatches 1..."
cd "$TMPRS"
"$DAME_BIN" sort \
    --fq "$FIXTURES/sample_primer_err.fastq" \
    --primers "$FIXTURES/Primers.txt" \
    --tags "$FIXTURES/Tags.txt" \
    --primer-mismatches 1

# With 1 mismatch allowed, the read must be recovered into Tag1_Tag2.txt
for D in "$TMPPY" "$TMPRS"; do
    if [ ! -f "$D/Tag1_Tag2.txt" ]; then
        echo "FAIL: $D did not produce Tag1_Tag2.txt at -m 1"; exit 1
    fi
done

echo "==> Comparing Tag1_Tag2.txt (py vs rust, -m 1)..."
if ! diff <(sort "$TMPPY/Tag1_Tag2.txt") <(sort "$TMPRS/Tag1_Tag2.txt"); then
    echo "FAIL: Tag1_Tag2.txt differs between dame-py and dame at -m 1"
    exit 1
fi

# Sanity: the recovered barcode and count must be present
if ! grep -q "ATATATATATAT" "$TMPRS/Tag1_Tag2.txt"; then
    echo "FAIL: expected recovered barcode ATATATATATAT not found"; exit 1
fi

echo "==> Confirming the read is NOT recovered at -m 0 (both tools)..."
cd "$TMPN0"
"$DAME_BIN" sort \
    --fq "$FIXTURES/sample_primer_err.fastq" \
    --primers "$FIXTURES/Primers.txt" \
    --tags "$FIXTURES/Tags.txt"
if [ -f "$TMPN0/Tag1_Tag2.txt" ]; then
    echo "FAIL: read should NOT be recovered at -m 0"; exit 1
fi

echo "PASS: dame and dame-py sort agree at -m 1 and reject at -m 0"
```

- [ ] **Step 3: Make the script executable and build the Rust binary**

```bash
chmod +x tests/integration/run_sort_mismatch.sh
cargo build --release --manifest-path rust/Cargo.toml
```

- [ ] **Step 4: Run the integration script locally**

Run: `bash tests/integration/run_sort_mismatch.sh`
Expected: `PASS: dame and dame-py sort agree at -m 1 and reject at -m 0`
(Requires `dame-py` installed: `pip install -e python/`.)

- [ ] **Step 5: Wire the script into CI**

In `.github/workflows/ci.yml`, in the `integration-tests` job, add a step after the existing `Run sort integration test` step:

```yaml
      - name: Run sort mismatch integration test
        run: bash tests/integration/run_sort_mismatch.sh
```

- [ ] **Step 6: Commit**

```bash
git add tests/fixtures/sample_primer_err.fastq tests/integration/run_sort_mismatch.sh .github/workflows/ci.yml
git commit -m "test: cross-impl primer-mismatch integration test + CI

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

## Task 6: Docs + version bump

**Files:**
- Modify: `rust/Cargo.toml` (version → 0.5.0)
- Modify: `README.md`
- Modify: `tutorial/README.md`

- [ ] **Step 1: Bump the crate version**

In `rust/Cargo.toml`, change `version = "0.4.0"` to `version = "0.5.0"`.

- [ ] **Step 2: Document the flag in README.md**

In `README.md`:

In the Rust `dame sort` quick-start block, add the flag (and add the same `--primer-mismatches 1` line to the Python `dame-py sort` example — note the Python example flags are tracked separately as a known doc issue, so only append the new flag in the form already shown there):

```
dame sort \
  --fq Pool1.fastq \
  --primers Primers.txt \
  --tags Tags.txt \
  --primer-mismatches 1
```

In the "Pipeline overview" section, update the `dame sort` line to mention the option:

```
dame sort     -fq POOL.fastq --primers P.txt --tags T.txt [--primer-mismatches N]
```

Add a development-history entry after item 9:

```markdown
10. **DAMe v2.4 — Configurable primer mismatches.**  `sort` gained a
    `-m`/`--primer-mismatches N` option (default 0) that tolerates up to N
    substitutions per primer match, IUPAC-aware, using leftmost-within-budget
    selection.  The budget applies independently to each of the four primer
    sites (forward/reverse orientation × start/end).  Tags are still matched
    exactly.  At N=0 the output is byte-identical to v2.3, verified by the
    existing integration tests; a new `run_sort_mismatch.sh` checks both
    implementations agree at N=1.  The Python primer matcher was rewritten from
    `re` to a manual IUPAC sliding window mirroring Rust, removing the `re`
    dependency.
```

Bump the title heading from `# DAMe v2.3: DNA Metabarcoding toolkit` to `# DAMe v2.4: DNA Metabarcoding toolkit`.

- [ ] **Step 3: Document the flag in tutorial/README.md**

In `tutorial/README.md`, find the section describing `sort` flags and add a bullet/line documenting:

```
--primer-mismatches N   Allow up to N substitutions per primer match
                        (IUPAC-aware, default 0). Tags are always matched exactly.
```

(Place it alongside the existing `--keep-primers-seq` documentation; match the surrounding formatting.)

- [ ] **Step 4: Verify nothing broke**

Run: `cargo build --release --manifest-path rust/Cargo.toml && cargo test --manifest-path rust/Cargo.toml && pytest python/tests/ -v`
Expected: PASS — build + all tests green.

- [ ] **Step 5: Commit**

```bash
git add rust/Cargo.toml rust/Cargo.lock README.md tutorial/README.md
git commit -m "docs: document --primer-mismatches; bump to v2.4 (crate 0.5.0)

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>"
```

---

## Final verification

- [ ] **Run the whole suite end to end**

```bash
cargo test --manifest-path rust/Cargo.toml
pip install -e python/ && pytest python/tests/ -v
cargo build --release --manifest-path rust/Cargo.toml
pip install -e python/
export PATH="$PWD/rust/target/release:$PATH"
bash tests/integration/run_sort.sh
bash tests/integration/run_sort_mismatch.sh
bash tests/integration/run_pipeline.sh
```

Expected: all PASS. `run_sort.sh` (N=0) confirms byte-identical regression; `run_sort_mismatch.sh` confirms N=1 recovery and cross-implementation agreement.
