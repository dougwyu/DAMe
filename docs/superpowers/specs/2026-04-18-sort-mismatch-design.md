# Sort Mismatch Tolerance Design

**Date:** 2026-04-18  
**Scope:** `python/dame/sort.py` and `python/dame/modules_sort.py` only

---

## Background

The current `dame-py sort` step finds primers using exact regex matching (IUPAC-aware but zero mismatches). Tag matching is also exact. In practice, sequencing errors can occur in both primer and tag regions, causing valid reads to be discarded.

Tag analysis (tutorial/Tags.txt, 48 tags, 8 bp each):
- Minimum pairwise Hamming distance: **3**
- Safe mismatch correction range: **mt ≤ 1** (mt=2 risks wrong-tag assignment for 84 tag pairs at distance 3)

---

## New Parameters

Both parameters are optional and default to 0 (preserving current behaviour exactly).

| Parameter | Type | Default | Meaning |
|-----------|------|---------|---------|
| `-m N` | int ≥ 0 | 0 | Max mismatches allowed in each primer (F and R independently, each ≤ N) |
| `-mt N` | int ≥ 0 | 0 | Max mismatches allowed in each tag (tag1 and tag2 independently, each ≤ N) |

Example:
```
dame-py sort -fq reads.fq -p Primers.txt -t Tags.txt -m 2 -mt 1
```

---

## Algorithm: Tag-Anchored Mismatch Search

When `m = 0` and `mt = 0`: existing code path (`GetPiecesInfo` with `re.finditer`) runs unchanged.

When `m > 0` or `mt > 0`: `GetPiecesInfoMismatch` is called instead.

### Read structure

Forward orientation:
```
[tag1_fwd][F_primer][amplicon][RC(R_primer)][tag2_rc]
```

Reverse orientation:
```
[tag1_fwd][R_primer][RC(amplicon)][RC(F_primer)][tag2_rc]
```

### Steps (forward orientation)

1. **Find tag1 candidates** — for each tag in TAGS, compute Hamming distance between `TAGS[t][0]` and `read[0:len(tag)]`. Collect all tags with distance ≤ mt.
2. **Check F primer at anchored position** — for each tag1 candidate, verify that `read[len(tag1) : len(tag1)+len(F)]` matches F primer with ≤ m mismatches (IUPAC-aware).
3. **Find tag2 candidates** — for each tag in TAGS, compute Hamming distance between `TAGS[t][1]` (RC sequence) and `read[-len(tag):]`. Collect all tags with distance ≤ mt.
4. **Check RC(R) primer at anchored position** — for each tag2 candidate, verify that `read[-len(tag2)-len(RC_R) : -len(tag2)]` matches RC(R) with ≤ m mismatches (IUPAC-aware).
5. **Extract amplicon** — once a valid (tag1, F_match, RC_R_match, tag2) is confirmed, the amplicon boundaries follow the same `keepPrimersSeq` logic as `GetPiecesInfo`: if `keepPrimersSeq` is False, trim primers from both ends; if True, include them.
6. **Combine and resolve** — collect all valid (tag1, tag2, primer combo) tuples:
   - Exactly one valid combination → extract amplicon, call `FillHAP`
   - Multiple combinations → pick the one with the lowest total mismatch count (tag1_mm + F_mm + RC_R_mm + tag2_mm). If still tied → discard as ambiguous (CountErrors)
   - Zero combinations → discard (CountErrors)

Reverse orientation uses the same steps with R primer and RC(F) primer.

### IUPAC mismatch counting

A new helper `iupac_match(primer_base, read_base) -> bool` checks whether `read_base` satisfies the IUPAC constraint of `primer_base` (e.g., `N` matches any base, `R` matches A or G). A position counts as a mismatch only when `iupac_match` returns False.

```python
IUPAC_SETS = {
    'A': {'A'}, 'C': {'C'}, 'G': {'G'}, 'T': {'T'},
    'M': {'A','C'}, 'R': {'A','G'}, 'W': {'A','T'},
    'S': {'C','G'}, 'Y': {'C','T'}, 'K': {'G','T'},
    'V': {'A','C','G'}, 'H': {'A','C','T'},
    'D': {'A','G','T'}, 'B': {'C','G','T'},
    'N': {'A','C','G','T'},
}

def iupac_match(p, r):
    return r.upper() in IUPAC_SETS.get(p.upper(), set())

def hamming_iupac(primer, region):
    # len(region) must equal len(primer)
    return sum(0 if iupac_match(p, r) else 1 for p, r in zip(primer, region))
```

---

## Code Changes

### `modules_sort.py`

1. Add `IUPAC_SETS`, `iupac_match()`, `hamming_iupac()` near top of file.
2. Modify `readPrimers()` to also store raw IUPAC sequences at indices 2 and 3:
   - `PRIMERS[name][2] = [F_raw, R_raw]` (original sequences before AMBIG substitution)
   - `PRIMERS[name][3] = [Frc_raw, Rrc_raw]`
   - Indices 0 and 1 (regex patterns) are unchanged — existing `GetPiecesInfo` continues to work.
3. Add `GetPiecesInfoMismatch(line, PRIMERS, TAGS, keepPrimersSeq, max_primer_mm, max_tag_mm)`.

### `sort.py`

1. Add `-m` and `-mt` arguments to both `register_subcommand()` and `main()`.
2. In `run()`, dispatch based on `m`/`mt`:
   ```python
   if args.m == 0 and args.mt == 0:
       Info = GetPiecesInfo(line, PRIMERS, TAGS, args.keepPrimersSeq)
   else:
       Info = GetPiecesInfoMismatch(line, PRIMERS, TAGS, args.keepPrimersSeq, args.m, args.mt)
   ```
   The return value format of `GetPiecesInfoMismatch` is identical to `GetPiecesInfo` so the downstream `FillHAP` / error-counting logic needs no changes.

---

## What Does NOT Change

- `FillHAP`, `PrintSortedCollapsedCountedSeqs`, `PrintSummaryFile`, `PrintSplitSummaryFile`
- Output file format (`.txt` per tag combination)
- `filter.py` and `modules_filter.py`
- v1 code (`DAMe_1.0/bin/`)
- Rust implementation
