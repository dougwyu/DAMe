# DAMe Tag Mismatches + Anchored Matching (sort) Design

**Date:** 2026-06-11
**Status:** Approved

## Goal

Add a configurable, per-tag substitution tolerance (`-mt` / `--tag-mismatches`,
default 0) to the `sort` step in both implementations, so reads with sequencing
errors in the tag region are recovered rather than discarded. Delivering this
introduces a **tag-anchored** matcher for the tolerant path, which also makes
primer matching more robust (the primer is checked at the expected offset after
the tag, eliminating spurious early-match risk).

Design adapted from PR #1 (`dougwyu/DAMe#1`, by @jiyinqiu).

## Scope

- `sort` subcommand only.
- Substitutions only (Hamming distance). Indels are out of scope.
- `mt` applies per tag: tag1 and tag2 each tolerate up to N mismatches,
  independently — mirroring how `m`/`--primer-mismatches` is per primer.

### Out of scope (YAGNI)

Indels / edit distance; per-tag-different `mt`; making the anchored matcher the
path for `m=0 & mt=0`; the PR's `filter` auto-output-directory feature.

## CLI

| Implementation | New flag | Default |
|---|---|---|
| Rust `dame sort` | `--tag-mismatches <N>` (long only; `-m` and `-t` are taken) | 0 |
| Python `dame-py sort` | `-mt` / `--tag-mismatches <N>` | 0 |

Existing `-m` / `--primer-mismatches` is unchanged.

## Dispatch: two matching paths

`run()` chooses the matcher per invocation:

- **`m == 0 && mt == 0`** → the existing exact matcher (`get_pieces_info` /
  `GetPiecesInfo`, backed by leftmost `find_primer`). Unchanged. Output is
  byte-identical to the current release, guaranteed by the existing unit and
  integration regression tests (`run_sort.sh`).
- **`m > 0 || mt > 0`** → the new **anchored** matcher.

Consequence: `find_primer` (leftmost scan with budget) and its unit tests remain
as the exact-path matcher. The production *mismatch* path is now the anchored
matcher, so the `m > 0` behavior shifts from leftmost-scan to anchored (a
robustness improvement; for the existing `run_sort_mismatch.sh` fixture the
result is the same, so that test continues to pass and becomes an
anchored-path regression guard).

## Anchored matcher

A new function in each language:

- Rust: `get_pieces_info_anchored(line, primers, tags, keep_primers_seq, max_primer_mm, max_tag_mm) -> Option<PieceInfo>`
- Python: `GetPiecesInfoMismatch(line, PRIMERS, TAGS, keepPrimersSeq, max_primer_mm, max_tag_mm)`

Reads are uppercased first (consistent with the existing matcher). For each
primer set, both orientations are attempted.

### Forward orientation

Read structure: `[tag1_fwd][F][amplicon][RC(R)][tag2_rc]`.

1. **Tag1 candidates** — for every tag `t`, if
   `hamming_iupac(t.fwd, read[0 : len(t.fwd)]) <= mt`, record `(t, tag1_mm)`.
2. **F primer at anchored offset** — for each tag1 candidate, let
   `off = len(t.fwd)`; if `hamming_iupac(F, read[off : off+len(F)]) <= m`,
   record `F_mm`. (Fixed offset — not a scan.)
3. **Tag2 candidates** — for every tag `t`, if
   `hamming_iupac(t.rc, read[-len(t.rc):]) <= mt`, record `(t, tag2_mm)`.
4. **RC(R) primer at anchored offset** — for each tag2 candidate, check
   `RC(R)` in `read[-len(tag2)-len(RC_R) : -len(tag2)]` with `<= m` mismatches;
   record `RC_R_mm`.
5. **Assemble & resolve** — form every valid `(tag1, F, RC(R), tag2)`
   combination; score by `tag1_mm + F_mm + RC_R_mm + tag2_mm`. Keep the unique
   minimum. If two or more combinations tie for the minimum, **discard the read
   as ambiguous** (counts toward the error total). Zero valid combinations →
   discard.

Amplicon boundaries follow the existing `keepPrimersSeq` logic (trim primers
when False, include them when True), using the anchored offsets.

### Reverse orientation

Read structure: `[tag1_fwd][R][RC(amplicon)][RC(F)][tag2_rc]`. Same five steps
with `R` as the start primer and `RC(F)` as the end primer, the tag-role swap
(`tag1` resolved against the RC tag sequence, `tag2` against the forward tag
sequence — matching the exact path), and `between = RC(between)`.

### Return value

Identical shape to the exact path (`PieceInfo` in Rust;
`[tag1, tag2, primer, between]` / `[1]` on failure in Python), so `FillHAP`,
the per-tag-combination output files, and the error counter are unchanged.

## Shared helper

`hamming_iupac(pattern, region) -> usize`: count of positions where
`iupac_matches(pattern_base, region_base)` is false. Returns a large sentinel
(or the caller guards) when `len(region) != len(pattern)` so out-of-range slices
never match. Built on the existing `iupac_matches`; reused for both tag and
primer offset checks in both languages. (`find_primer` is unchanged.)

## Tag-distance safety warning

At startup, only when `mt > 0`: compute the minimum pairwise Hamming distance
among **equal-length** forward tag sequences. If `2 * mt + 1 > min_distance`,
print a non-fatal warning to stderr, e.g.:

```
WARNING: --tag-mismatches 2 may misassign reads: minimum tag Hamming distance
is 3, so the safe maximum is 1. Reads within 2 mismatches of the wrong tag will
be miscorrected.
```

If no two tags share a length, skip the warning (length differences themselves
disambiguate). Computed identically in both implementations. It is stderr only,
not output, so the byte-identical guarantee is unaffected.

## Rust changes

- `rust/src/sort.rs`:
  - Add `hamming_iupac(pattern: &[u8], region: &[u8]) -> usize`.
  - Extend `TagLookup` with ordered lists for iteration:
    `fwd_list: Vec<(Vec<u8>, String)>` and `rc_list: Vec<(Vec<u8>, String)>`
    (the existing `by_fwd` / `by_rc` maps stay for the exact path).
  - Add `get_pieces_info_anchored(...)`.
  - Add `tag_mismatches: usize` to `SortArgs`
    (`#[arg(long = "tag-mismatches", default_value = "0")]`).
  - Add `min_equal_length_tag_distance(tags) -> Option<usize>` and emit the
    warning in `run()` when `tag_mismatches > 0`.
  - `run()` dispatches: exact path when both budgets 0, else anchored.

## Python changes

- `python/dame/modules_sort.py`:
  - Add `hamming_iupac(pattern, region)`.
  - Add `GetPiecesInfoMismatch(...)`.
  - Add `min_equal_length_tag_distance(TAGS)`.
  - `readPrimers` already stores raw primer strings; `TAGS` (name → [fwd, rc])
    is already iterable — no structural change needed.
- `python/dame/sort.py`:
  - Add `-mt` / `--tag-mismatches` (dest `tag_mismatches`, type int, default 0)
    to `register_subcommand`.
  - In `run()`: emit the warning when `tag_mismatches > 0`; dispatch to
    `GetPiecesInfoMismatch` when `m > 0 or mt > 0`, else `GetPiecesInfo`.

## Testing

### Unit (Rust `rust/tests/sort_test.rs`, Python `python/tests/test_sort.py`)

- `hamming_iupac`: exact, IUPAC-satisfied (0 cost), substitution counts,
  length-mismatch guard.
- Anchored forward recovery with one tag mismatch (`m=0, mt=1`): recovered;
  rejected at `mt=0`.
- Anchored reverse-orientation recovery with one tag mismatch.
- Combined `m=1, mt=1` recovery.
- Ambiguity: a read equidistant (tie) from two tags within `mt` → discarded.
- `m>0, mt=0` through the anchored matcher reproduces the primer-only result
  (parity with the leftmost-path fixtures).
- `min_equal_length_tag_distance`: correct value; warning predicate
  (`2*mt+1 > min`) fires at the threshold and not below it.

### Integration (`tests/integration/`)

- New `run_sort_tag_mismatch.sh`: fixture with a single tag-base error; run both
  `dame-py sort -mt 1` and `dame sort --tag-mismatches 1`; the read is recovered
  and the two outputs diff-identical; confirm rejection at `-mt 0`.
- Existing `run_sort.sh` (m=0,mt=0) and `run_sort_mismatch.sh` (m=1) remain
  as regression guards; the latter now exercises the anchored path.
- Wire `run_sort_tag_mismatch.sh` into `.github/workflows/ci.yml`.

## Docs / versioning

- `README.md`: document `--tag-mismatches`, the per-tag semantics, the
  best-match/ambiguity-discard resolution, and the safety guidance; add a
  development-history entry crediting PR #1 (@jiyinqiu); bump title to **v2.5**.
- `tutorial/README.md`: document the flag alongside `--primer-mismatches`.
- `rust/Cargo.toml`: version → **0.6.0**.

## Implementation order

1. Shared `hamming_iupac` + `TagLookup` lists (Rust) and `hamming_iupac`
   (Python), with unit tests.
2. Anchored matcher in each language, with unit tests (forward, reverse,
   combined, ambiguity, primer-only parity).
3. `min_equal_length_tag_distance` + startup warning, with unit tests.
4. `SortArgs` / CLI flag + `run()` dispatch, both languages.
5. Integration fixture + `run_sort_tag_mismatch.sh` + CI wiring.
6. README / tutorial / version bump.
