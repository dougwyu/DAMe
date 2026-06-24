# Spec: `convert` subcommand (v2.6)

**Date:** 2026-06-24  
**Branch:** `convert-subcommand`  
**Status:** Approved — proceed to implementation plan

---

## Background

`DAMe_1.0/bin/convertToUSearch.py` is a script from the original DAMe v1.0 that converts
`FilteredReads.fna` output from `dame filter` into the input format expected by two
downstream clustering tools: USEARCH and sumaclust. It was not committed to the repository
originally but is a required step in the DAMe amplicon-clustering workflow.

This spec ports and modernises it as `dame convert` / `dame-py convert`, adding it to the
maintained v2 toolkit (both Python and Rust implementations).

### Bugs in the original that this spec fixes

| Bug | Fix |
|-----|-----|
| `--maxLength` default `1e6` is a float; arithmetic on it raises `TypeError` | Default `None` (no cap); activates only when flag is passed |
| `SampleFastas/` directory not created before writing → crash | `os.makedirs("SampleFastas", exist_ok=True)` before first write |
| `args.lmax += 1` then pad to `args.lmax` → sequence padded one nt too long | Remove pre-increment; pad to exactly `--max-length` |
| `import numpy as np` unused | Dropped |
| Script not version-controlled | Committed to `DAMe_1.0/bin/` (frozen reference) + ported into maintained toolkit |

---

## Input format

`FilteredReads.fna` produced by `dame filter`. Each sequence is two lines:

```
>SampleA Tag1-Tag2.Tag3-Tag4_1 5_4
ATCGATCG...
```

Header has exactly 3 whitespace-delimited tokens:
- `toks[0]` — `>SampleName`
- `toks[1]` — tag-pair descriptor (not used by `convert`)
- `toks[2]` — underscore-delimited replicate counts, e.g. `5_4` means pool1=5, pool2=4

Derived values:
- `sample = toks[0][1:]` (strip leading `>`)
- `size = sum(int(x) for x in toks[2].split("_"))`

---

## Output formats

### Default: sumaclust

File: `FilteredReads.forsumaclust.fna`

```
>SampleA:1 count=9
ATCGATCG...
```

- `:N` is a per-record counter (starts at 1, increments with each passing sequence)
- `count=<size>` is the total read count across replicates
- No N-padding

### With `-u` / `--usearch`: USEARCH

File: `FilteredReads.forusearch.fna`

```
>SampleA;size=9
ATCGATCG...
```

- `;size=<size>` is the total read count
- When `--max-length N` is set, sequences are right-padded with `N` characters to
  exactly `N` bases. Sequences longer than `N` after filtering still pass unchanged (the
  filter already enforced the length bound).
- N-padding is only applied in USEARCH mode and only when `--max-length` is explicitly
  provided (opt-in, not default).

### With `-s` / `--sample-fastas`: per-sample files

Creates directory `SampleFastas/` (if it does not exist) and writes one file per sample:

```
SampleFastas/<SampleName>.fixed.fasta
```

Each file contains the same header and sequence lines as the main output file (sumaclust
or USEARCH format, whichever is active). This flag is independent of `-u` and can be
combined with it.

---

## CLI flags

### New canonical names (both Python and Rust)

| Flag | Short | Type | Default | Description |
|------|-------|------|---------|-------------|
| `--in-fasta` | `-i` | str | required | Input `FilteredReads.fna` |
| `--min-length` | | int | 0 | Drop sequences shorter than N |
| `--max-length` | | int | None | Drop sequences longer than N; pad to N in USEARCH mode |
| `--usearch` | `-u` | flag | false | Write USEARCH format (default: sumaclust) |
| `--sample-fastas` | `-s` | flag | false | Write per-sample fastas to `SampleFastas/` |

### Python backwards-compatible aliases

The Python CLI also accepts the original v1.0 spellings so that existing scripts
and muscle-memory are not broken:

| Original flag | Canonical alias |
|---------------|----------------|
| `--inFasta` | `--in-fasta` |
| `-lmin` / `--minLength` | `--min-length` |
| `-lmax` / `--maxLength` | `--max-length` |
| `--sampleFastas` | `--sample-fastas` |

Rust uses only the canonical names (it has no v1.0 compatibility obligation).

---

## Length filtering

A sequence passes if:

```
len(sequence) >= min_length  AND  (max_length is None OR len(sequence) <= max_length)
```

Both bounds are inclusive. No pre-increment/post-increment adjustments.

---

## Processing logic (pseudocode)

```
open output file (sumaclust or usearch name)
if --sample-fastas: makedirs("SampleFastas", exist_ok=True)

counter = 1
for each header+sequence pair in input:
    parse sample, size from header
    if len(sequence) < min_length: skip
    if max_length is not None and len(sequence) > max_length: skip

    if usearch:
        out_header = f">{sample};size={size}"
        out_seq = sequence.ljust(max_length, 'N') if max_length else sequence
    else:
        out_header = f">{sample}:{counter} count={size}"
        out_seq = sequence
        counter += 1

    write out_header + "\n" + out_seq + "\n" to main output
    if sample-fastas:
        open SampleFastas/<sample>.fixed.fasta (create if new, append if seen)
        write out_header + "\n" + out_seq + "\n"

close all files
```

---

## Files to create / modify

### New files

| File | Purpose |
|------|---------|
| `python/dame/convert.py` | Python `convert` subcommand |
| `python/tests/test_convert.py` | Python unit tests |
| `rust/src/convert.rs` | Rust `convert` subcommand |
| `rust/tests/convert_test.rs` | Rust unit tests |
| `tests/fixtures/FilteredReads_small.fna` | Shared fixture (3–4 records, 2 samples) |
| `tests/integration/run_convert.sh` | Cross-implementation parity test |

### Modified files

| File | Change |
|------|--------|
| `python/dame/__main__.py` | Register `convert` subcommand |
| `rust/src/main.rs` | Add `convert` subcommand dispatch |
| `rust/Cargo.toml` | Bump crate version `0.6.x` → `0.7.0` |
| `README.md` | Add `convert` to command table; bump version to v2.6 |
| `tutorial/README.md` | Add `convert` step after `filter` |
| `.github/workflows/ci.yml` | Add `run_convert.sh` step |
| `DAMe_1.0/bin/convertToUSearch.py` | Commit frozen original (currently untracked) |

---

## Integration test design (`run_convert.sh`)

1. Run `dame-py convert -i FilteredReads_small.fna` → `FilteredReads.forsumaclust.fna`
2. Run `dame convert -i FilteredReads_small.fna` → `FilteredReads.forsumaclust.fna`
3. `diff` the two outputs — must be byte-identical
4. Repeat for USEARCH mode (`-u`)
5. Repeat for USEARCH + `--max-length` (fixed-length padding)
6. Confirm `SampleFastas/` is created and non-empty with `-s`

---

## Test fixture (`FilteredReads_small.fna`)

```
>Sample1 Tag1-Tag2.Tag3-Tag4_1 5_4
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
>Sample1 Tag1-Tag2.Tag3-Tag4_2 3_2
GCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCG
>Sample2 Tag5-Tag6.Tag7-Tag8_1 10_8
TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
>Sample2 Tag5-Tag6.Tag7-Tag8_2 1_0
AAAA
```

Record 4 (`AAAA`, 4 nt) is used to test `--min-length 10` filtering.

---

## Parity invariant

Python and Rust must produce byte-identical output for all flag combinations. The
integration test enforces this.

---

## Version bump

- Repository: `v2.5` → `v2.6`
- Rust crate: `0.6.x` → `0.7.0`
- Python package: `2.5.x` → `2.6.0`
