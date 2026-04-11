# DAMe Python 3 + Rust Port Design

**Date:** 2026-04-11
**Status:** Approved

## Goals

Port DAMe (a DNA metabarcoding pipeline toolkit) from Python 2 to Python 3 and Rust.

- **Python 3 port:** Reference implementation and correctness oracle
- **Rust port:** Primary deliverable — single compiled binary, no Python dependency, native performance
- **Goals for Rust port:** Performance (large FASTQ files), distribution (single binary), modernization

## Repository Structure

```
DAMe/
  python/
    dame/
      __init__.py
      sort.py
      chimera_check.py
      filter.py
      decollapse.py
      rsi.py
    tests/
      test_sort.py
      test_chimera_check.py
      test_filter.py
      test_decollapse.py
      test_rsi.py
    pyproject.toml        ← `dame-py` package, entry point: `dame-py`
  rust/
    Cargo.toml            ← workspace root
    src/
      main.rs             ← CLI entry point, dispatches subcommands
      sort.rs
      chimera_check.rs
      filter.rs
      decollapse.rs
      rsi.rs
    tests/
      sort_test.rs
      chimera_check_test.rs
      filter_test.rs
      decollapse_test.rs
      rsi_test.rs
  tests/
    integration/
      run_sort.sh
      run_chimera.sh
      run_filter.sh
      run_decollapse.sh
      run_rsi.sh
      run_pipeline.sh     ← end-to-end: sort → chimera → filter → rsi
    fixtures/             ← symlink or copy of example/ data
  example/                ← existing example data (unchanged)
  bin/                    ← existing Python 2 originals (kept for reference)
  DAMe_Manual.pdf
  README.txt
```

The existing `bin/` and `example/` directories are left untouched to preserve git history and serve as Python 2 reference.

## Python 3 Port

### Approach

Mechanical Python 2→3 upgrade — no logic changes. The module structure is preserved (`sort.py` imports from `modules_sort.py` equivalents), reorganized into a package under `python/dame/`.

### Python 2→3 Changes

| Python 2 | Python 3 |
|---|---|
| `print x` | `print(x)` |
| `dict.has_key(k)` | `k in dict` |
| `string.maketrans(...)` | `str.maketrans(...)` |
| `optparse` (RSI.py) | `argparse` |
| Integer division behavior | explicit `//` or `float()` where needed |
| `file.readline()` loop | kept as-is (valid in Python 3) |

### Entry Point

`pyproject.toml` defines a `dame-py` console script that dispatches subcommands:
```
dame-py sort ...
dame-py chimera ...
dame-py filter ...
dame-py decollapse ...
dame-py rsi ...
```

### Dependencies

- `numpy` (already used in `rsi.py`, retain for Python 3 port)
- `pytest` (dev dependency)

## Rust Port

### CLI Interface

Single binary using `clap` (derive feature):

```
dame sort        -fq <file> -p <primers> -t <tags> [--keep-primers-seq]
dame chimera     -ps-info <file> -x <int> [-p <int>]
dame filter      -ps-info <file> [-x <int>] [-y <int>] [-p <int>] [-t <int>] [-l <int>] [--chimera-checked]
dame decollapse  -input <file> [-out-fas <file>]
dame rsi         [-e] [-o <file>] <input>
```

### Module Structure

Each subcommand is implemented in its own file, called from `main.rs`:

- `src/sort.rs` — FASTQ demultiplexing, tag/primer matching, collapse to unique sequences
- `src/chimera_check.rs` — file preparation for UCHIME chimera removal
- `src/filter.rs` — PCR replicate presence, abundance, and length filtering
- `src/decollapse.rs` — expand collapsed sequences to individual reads
- `src/rsi.rs` — Renkonen Similarity Index calculation

### Key Crates

| Need | Crate |
|---|---|
| CLI parsing | `clap` (derive feature) |
| Regex | `regex` |
| Matrix math (RSI) | `ndarray` (replaces numpy) |
| Error handling | `anyhow` |

### Performance

- `sort` is the hot path: reads FASTQ line-by-line via `BufReader`, avoids unnecessary heap allocation in the inner loop
- Regex patterns compiled once at startup, reused across all reads
- Same algorithm as Python — no algorithmic changes, just native execution speed

### External Dependencies

`dame chimera` wraps the external `usearch` binary (same as the Python 2 original), invoked via `std::process::Command`. `usearch` must be on `PATH` at runtime. The Rust port does not bundle or replace `usearch`.

### Output Format

Byte-for-byte identical to the Python 3 port for all text/FASTA outputs: same tab-separated fields, same filenames, same FASTA header formats. Exception: RSI floating point values are compared with a numeric tolerance (not strict `diff`) in both unit and integration tests.

## Testing Strategy

### Python 3 Unit Tests (`python/tests/`, pytest)

- One test file per module
- Test each function with small synthetic inputs
- Edge cases: empty tag match, reverse complement reads, zero-count sequences, single PCR replicate

### Rust Unit Tests (`rust/tests/`, `#[test]`)

- Mirror Python unit tests 1:1 — same inputs, same expected outputs
- Test each function in `sort.rs`, `filter.rs`, etc. independently

### Integration Tests (`tests/integration/`)

- Shell scripts run both `dame-py <subcommand>` and `dame <subcommand>` on the same fixture files
- `diff` every output file — any discrepancy fails the test
- One script per subcommand, plus `run_pipeline.sh` for end-to-end: `sort → chimera → filter → rsi`

### CI (GitHub Actions)

- Matrix: Python 3.11+, Rust stable
- Jobs: `pytest`, `cargo test`, integration tests
- Integration tests run last, after both unit suites pass

## Implementation Order

1. Python 3 port (establishes correctness baseline)
2. Integration test fixtures and diff harness (validates Python 3 port against example data)
3. Rust port, subcommand by subcommand, in pipeline order: `sort` → `chimera` → `filter` → `decollapse` → `rsi`
4. Rust unit tests added alongside each subcommand
5. CI configuration
