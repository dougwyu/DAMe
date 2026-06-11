# DAMe Primer Mismatches (sort) Design

**Date:** 2026-06-10
**Status:** Approved

## Goal

Add a configurable, per-primer substitution tolerance to the `sort` step in both
implementations, so reads with sequencing errors inside the primer region are no
longer discarded. Today both implementations require an *exact* IUPAC primer
match, so a single miscalled base in the primer drops the whole read.

The default is `0`, which must remain **byte-for-byte identical** to current
`sort` output. This is enforced by the existing cross-implementation integration
tests (`tests/integration/run_sort.sh`), which become the N=0 regression guard.

## Scope

- Applies to the `sort` subcommand only.
- Primers only. **Tags are still matched exactly** — this feature does not
  introduce tag mismatch tolerance.
- Substitutions only (Hamming distance). Indels are explicitly out of scope.

### Out of scope (YAGNI)

Indels / edit distance, tag mismatch tolerance, per-primer-different N,
"best-match" (fewest-mismatch) window selection, and any change to the other
subcommands.

## Matching Semantics

- A **mismatch** is a position within the primer-width window where the read
  base does **not** satisfy the primer's IUPAC code. IUPAC ambiguity is still
  honored — an ambiguous primer base matching one of its allowed bases costs 0.
- **Leftmost within budget.** Scan the read left→right; return the first window
  whose mismatch count is ≤ N. The per-window comparison early-exits once the
  count exceeds N. At N=0 this reduces exactly to today's "break on first
  mismatch" behavior, so N=0 output is unchanged.
- **Per-primer budget.** N applies independently to each of the four primer
  matches: forward-orientation start `F` and end `RC(R)`; reverse-orientation
  start `R` and end `RC(F)`. A single read can absorb up to N mismatches on each
  primer site.
- All downstream slicing is unchanged. Tag extraction, the `keepPrimersSeq`
  flag, and the barcode boundary all consume the start/end positions the matcher
  returns, exactly as today.

### Selection example (leftmost within budget)

Forward primer `GCATGC`, N=1, read region:

```
index:   0 1 2 3 4 5 6 7 8 9 10 11
base:    T T G C T T G C A T G  C
window @ i=2:  G C T T G C  → 1 mismatch (pos 2: T≠A)  → accepted (first within budget)
window @ i=6:  G C A T G C  → 0 mismatches             → not reached
```

The matcher returns `i=2`. (The alternative "best within budget" policy, which
would return `i=6`, is deliberately not implemented — see Out of scope.)

### Edge cases

- **N ≥ primer length:** every window is trivially within budget, so the first
  window (`i=0`) matches. This is degenerate but not special-cased; documented as
  user error.
- **Primer longer than read:** returns no match (unchanged).
- **Empty / unmatched reads:** counted as errors (unchanged).

## CLI

| Implementation | Flag | Default |
|---|---|---|
| Rust `dame sort` | `-m` / `--primer-mismatches <N>` | `0` |
| Python `dame-py sort` | `-m` / `--primer-mismatches <N>` | `0` |

Realigning the rest of the Python CLI to `--long` style (a separate review
finding) is **out of scope** for this feature.

## Rust Changes (`rust/src/sort.rs`)

- `find_primer(primer: &[u8], seq: &[u8], max_mismatches: usize) -> Option<(usize, usize)>`
  — count mismatches per window using `iupac_matches`, with an inner early-exit
  once the count exceeds `max_mismatches`; return the first window within budget.
  N=0 collapses to the current break-on-first-mismatch loop.
- `get_pieces_info(..., max_mismatches: usize)` — thread the budget through to
  all four `find_primer` call sites.
- `SortArgs` — add `primer_mismatches: usize` (`#[arg(short = 'm', long = "primer-mismatches", default_value = "0")]`).
- `run` — pass `args.primer_mismatches` into `get_pieces_info`.

## Python Changes (`python/dame/modules_sort.py`, `python/dame/sort.py`)

Replace the stdlib-`re` matcher with a manual IUPAC sliding-window that mirrors
the Rust algorithm exactly, so the two implementations stay algorithmically
identical (the premise of the integration suite). Stdlib `re` cannot express
"≤ N substitutions."

- Add `iupac_matches(primer_base, read_base)` and
  `find_primer(primer, seq, max_mismatches)` helpers mirroring Rust.
- `readPrimers` stores raw primer strings (F, R, RC(F), RC(R)) without regex
  expansion. The `AMBIG` dict (`sort.py`) and `import re` (`modules_sort.py`) are
  removed.
- `GetPiecesInfo(line, PRIMERS, TAGS, keepPrimersSeq, maxMismatches)` uses
  `find_primer` instead of `re.finditer`. All slicing/tag logic is unchanged.
- `sort.py` (`register_subcommand` and standalone `main`) add the
  `-m` / `--primer-mismatches` argument (default 0, `type=int`).

## Testing

### Rust unit (`rust/tests/sort_test.rs`)

- `find_primer` with `max_mismatches` 0/1/2: exact still works; 1 mismatch found
  at N=1; 2 mismatches rejected at N=1; IUPAC code + 1 substitution combined;
  leftmost-within-budget tiebreak (the `GCTTGC` vs later exact `GCATGC` case);
  primer longer than read.
- `get_pieces_info` on a read with a single primer-base error: `Some` at N=1,
  `None` at N=0.

### Python unit (`python/tests/test_sort.py`)

- Mirror the Rust tests 1:1 (same inputs, same expected outputs).

### Integration (`tests/integration/`)

- Add a fixture read carrying a single primer-base error. Run both `dame-py sort`
  and `dame sort` with `--primer-mismatches 1` (`-m 1` for Python) and `diff`
  outputs — they must match each other.
- Existing N=0 `run_sort.sh` stays as the byte-identical regression guard.
- Wire the new mismatch check into `.github/workflows/ci.yml`.

## Docs / Versioning

- `README.md`: document `-m` / `--primer-mismatches` for both implementations,
  add to the pipeline overview, bump to **v2.4** (crate **0.5.0**), add a
  development-history entry.
- `tutorial/README.md`: note the new flag.

## Implementation Order

1. Rust `find_primer` + threading + `SortArgs` flag, with Rust unit tests.
2. Python manual matcher + flag, with Python unit tests.
3. Integration fixture + mismatch script + CI wiring.
4. README / tutorial / version bump.
