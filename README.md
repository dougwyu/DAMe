# DAMe v2.3: DNA Metabarcoding toolkit

DAMe demultiplexes pooled metabarcoding / eDNA FASTQ reads by primer and
tag sequences (**sort**), optionally removes chimeric sequences (**chimera**),
filters amplicons across PCR replicates (**filter**), computes the Renkonen
Similarity Index between replicates (**rsi**), and expands collapsed sequences
back to individual reads (**decollapse**).  It is available in two
implementations:

| | Python 3 (`python/`) | Rust (`rust/`) |
|---|---|---|
| Requirements | Python ≥ 3.11, numpy | Rust stable (cargo) |
| Install | `pip install -e python/` | `cd rust && cargo build --release` |
| Entry point | `dame-py` | `dame` |
| Input FASTQ | plain or gzip | plain or gzip |
| Chimera check | via `usearch` on PATH | via `usearch` on PATH |

Reads are expected to have the structure
`[fwd_tag][fwd_primer][amplicon][rc(rev_primer)][rc(rev_tag)]`.
Both forward and reverse-complement orientations are detected automatically
during sort.  IUPAC ambiguity codes are supported in primer sequences.

## Performance

Small dataset — tutorial (392 reads, 1 pool):

| Step | Python 3 | Rust 2.0 | Rust 2.1 | 2.1 vs Python |
|------|----------|----------|----------|---------------|
| sort | ~280 ms | ~38 ms | ~36 ms | ~8× |

Large dataset — synthetic benchmark (196,000 reads, 2 pools of ~100k reads each):

| Step | Python 3 | Rust 2.0 | Rust 2.1 | Rust 2.2 | Rust 2.3 | 2.3 vs Python |
|------|----------|----------|----------|----------|----------|---------------|
| sort (plain FASTQ) | ~500 ms/pool | ~97 ms/pool | ~89 ms/pool | ~111 ms/pool (†) | ~105 ms/pool | ~4.8× |
| sort (gzip FASTQ)  | — | not supported | ~102 ms/pool | ~111 ms/pool | ~105 ms/pool | new in 2.1 |
| filter | ~149 ms | ~39 ms | ~38 ms | ~38 ms | ~38 ms | ~4× |
| rsi | ~150 ms | ~33 ms | ~34 ms | ~34 ms | ~34 ms | ~4.4× |

(†) The v2.2 sort figure was remeasured under the current toolchain
(rustc 1.94) and machine-load conditions for an apples-to-apples comparison
with v2.3.  The previous v2.2 figure of ~59 ms/pool was recorded on a
different toolchain/machine state; the absolute numbers are not directly
comparable across measurement sessions, but the relative v2.3-vs-v2.2 gap
reported here is a true interleaved benchmark.

The v2.1 sort improvement (~9% over 2.0) comes from `needletail`'s faster
FASTQ parsing and `ahash` replacing `SipHash` for DNA-string key lookups.
Filter and RSI are I/O-bound on small collapsed-sequence files at this scale,
so the hasher change has negligible effect there.  The gzip overhead (~15%
in 2.1) is effectively eliminated in v2.2: the new binary handles both plain
and `.fastq.gz` input at the same speed.

The v2.2 byte-matcher change replaced four `Regex::find()` calls per read
with a hand-written IUPAC sliding-window (`iupac_matches` + `find_primer`) and
removed the `regex` crate dependency entirely.  On the CO1 tutorial primers
(no ambiguity codes), the `regex` crate's SIMD-compiled DFA remains
competitive with naive byte-by-byte scanning, so the sort time is similar to
v2.1.  The benefit of v2.2 is correctness assurance on highly ambiguous
primers (e.g. many N or degenerate positions) and a simpler dependency tree.

The v2.3 sort changes target three per-read constant-factor costs in the hot
loop: (a) a byte-level `rc_bytes(&[u8]) -> Vec<u8>` replaces the char-based
`rc(&str) -> String` in the reverse-orientation branch of `get_pieces_info`,
eliminating UTF-8 decode on every reverse-orientation read; (b) `read_tags`
now pre-builds two `HashMap<Vec<u8>, String>` reverse-lookup maps (forward
and RC) at startup, replacing the per-read `O(N_tags)` linear scan in
`get_pieces_info` with `O(1)` hash lookups; (c) `fill_hap` no longer calls
`between.to_string()` on reads that land on an already-seen barcode — the
common case at typical amplicon duplication rates.  On the CO1 tutorial
primers with 8 tags, the measured sort speedup is ~5% (the small N_tags
limits the benefit of the O(1) lookup, and the synthetic dataset is
parsing-bound).  Larger tag panels should benefit more.

## Quick start

### Python version

```bash
pip install -e python/

dame-py sort \
  --fq Pool1.fastq \
  --primers Primers.txt \
  --tags Tags.txt

dame-py filter \
  --ps-info PSinfo.txt \
  --x 2 --y 2 --t 2 --l 50

dame-py rsi Comparisons_2PCRs.txt
```

### Rust version

```bash
cd rust && cargo build --release
# binary: rust/target/release/dame

dame sort \
  --fq Pool1.fastq \
  --primers Primers.txt \
  --tags Tags.txt

dame filter \
  --ps-info PSinfo.txt \
  --x 2 --y 2 --t 2 --l 50

dame rsi Comparisons_2PCRs.txt
```

## Pipeline overview

```
dame sort     -fq POOL.fastq --primers P.txt --tags T.txt
              → TagA_TagB.txt (collapsed unique seqs + counts) per tag pair
              → SummaryCounts.txt

dame chimera  --ps-info PSinfo.txt --x 2          # requires usearch on PATH
              → TagA_TagB_Pool.noChim.txt

dame filter   --ps-info PSinfo.txt --x 2 --y 2 --t 2 --l 50
              → Comparisons_2PCRs.txt (all seqs, all replicates)
              → FilteredReads.fna   (passed all thresholds)

dame rsi      Comparisons_2PCRs.txt
              → RSI_output.txt

dame decollapse --input TagA_TagB.txt --out-fas decollapsed.fasta
```

Run `dame sort` once per sequencing pool (in separate directories), then
point `dame filter` at a PSinfo file that maps each tag-pair file to its
sample and pool.

## Development history

DAMe was originally written in Python 2 by Zepeda-Mendoza et al. (2016).
[Claude Code](https://claude.ai/code) was used to modernise and extend the
codebase:

1. **Python 2 → Python 3 port.**  Five compatibility issues were fixed:
   `str.maketrans` replacing the removed `string.maketrans`; `dict.has_key()`
   replaced throughout with `k in dict`; `optparse` replaced with `argparse`;
   integer division made explicit with `//`; and `subprocess.PIPE` output
   handled as `bytes` rather than `str`.  A bug in `filter` was also fixed: the
   0-indexed key `"0"` (not `"1"`) is the correct reference length for
   `PSinsLines`.

2. **Formal test suite.**  A pytest suite (`python/tests/`) was written
   covering all five subcommands with synthetic inputs generated via `tmp_path`
   fixtures, and shell integration tests (`tests/integration/`) that smoke-test
   the installed `dame-py` entry point end-to-end.

3. **Rust port.**  Claude Code wrote a design specification
   (`docs/superpowers/specs/`) and an implementation plan
   (`docs/superpowers/plans/`) then executed it using the
   [Superpowers plugin](https://github.com/superpowers-sh/superpowers) for
   Claude Code, specifically the `subagent-driven-development` skill — a fresh
   subagent per task with spec-compliance and code-quality review after each.
   The result is `rust/`: a full single-binary reimplementation using
   [clap](https://github.com/clap-rs/clap) for the CLI,
   [regex](https://github.com/rust-lang/regex) for IUPAC primer matching,
   [ndarray](https://github.com/rust-ndarray/ndarray) for RSI matrix arithmetic,
   and [indexmap](https://github.com/indexmap-rs/indexmap) to preserve
   insertion-order output identical to the Python 3 port.

4. **Integration tests comparing both implementations.**  Shell scripts in
   `tests/integration/` run both `dame-py` and `dame` on the same inputs and
   diff their outputs (sorted before comparison where set-based ordering could
   differ; RSI values compared with 1 × 10⁻⁹ floating-point tolerance).

5. **Tutorial and dataset.**  A synthetic dataset
   (`tutorial/generate_tutorial_data.py`) was generated to demonstrate IUPAC
   primer matching, both read orientations, and every filter outcome
   (pass, fail propPCRs, fail minCount, fail minLength).  The tutorial lives in
   `tutorial/README.md`.

6. **GitHub Actions CI.**  A matrix workflow (`.github/workflows/ci.yml`) runs
   pytest on Python 3.11 and 3.12, `cargo test` on Rust stable, and all six
   integration scripts on every push and pull request to `master`.

7. **DAMe v2.1 — Rust performance improvements.**  Two targeted optimisations
   were added to the Rust binary: (a) the manual four-line FASTQ reading loop
   in `sort` was replaced with
   [needletail 0.5](https://github.com/onecodex/needletail), which also adds
   transparent gzip input support; (b) all `HashMap` and `HashSet` instances
   across every module now use the
   [ahash 0.8](https://github.com/tkaitchuck/aHash) non-cryptographic hasher,
   which is substantially faster for DNA string keys.  A pre-existing fragile
   `HashMap` + manual order-tracking `Vec` pattern in `chimera_check` was also
   replaced with `IndexMap`, consistent with the convention in `sort`.

8. **DAMe v2.2 — IUPAC byte matcher.**  The `regex` crate was removed from
   the Rust binary.  The four `Regex::find()` calls per read in `sort` were
   replaced with a hand-written byte sliding-window: `iupac_matches(u8, u8)`
   looks up the full 15-code IUPAC truth table, and `find_primer(&[u8], &[u8])`
   scans the read bytes leftmost-first.  `PrimerEntry` now stores raw byte
   slices instead of compiled `Regex` objects; `read_primers` no longer calls
   `Regex::new()`.  All integration tests continue to produce byte-identical
   output to `dame-py`.  The sort throughput on CO1 tutorial primers is similar
   to v2.1 (the `regex` DFA is competitive with naive byte scanning on
   exact-match patterns); the main benefit is a simpler dependency tree and
   guaranteed correctness on highly degenerate primer sequences.

9. **DAMe v2.3 — Sort hot-loop constant-factor reductions.**  Three targeted
   changes to per-read work in `sort`: (a) a new `rc_bytes(&[u8]) -> Vec<u8>`
   replaces the char-based `rc(&str) -> String` on the reverse-orientation
   branch of `get_pieces_info`, eliminating UTF-8 decoding of every
   reverse-orientation read; (b) `read_tags` now returns a `TagLookup` struct
   holding two pre-built `HashMap<Vec<u8>, String>` reverse-lookup maps
   (forward and RC), replacing the `O(N_tags)` `tags.iter().find(...)` linear
   scan in `get_pieces_info` with `O(1)` hash lookup on the raw byte slice;
   (c) `fill_hap` no longer calls `between.to_string()` on reads that hit an
   already-seen barcode — the common case at typical amplicon duplication
   rates.  All integration tests continue to produce byte-identical output to
   `dame-py`.  On CO1 tutorial primers with 8 tags, the measured sort speedup
   is ~5%; larger tag panels should benefit more from the O(1) lookup.

## Documentation

See `tutorial/README.md` for a full walkthrough covering all input file
formats, every command-line flag, output file formats, and a worked example
dataset that demonstrates all filter outcomes.

The original DAMe v1.0 manual, scripts, and example data are preserved in
`DAMe_1.0/` (`bin/`, `example/`, `README.txt`, `DAMe_Manual.pdf`).

## Input file formats

**Primers.txt** — one primer set per line, tab-separated:
```
CO1	GCRTGC	CTGACT
```
(`Name`, `ForwardSeq`, `ReverseSeq`; IUPAC ambiguity codes supported)

**Tags.txt** — one tag per line, tab-separated:
```
AACCGGT	tag1
TTGGCCA	tag2
```
(`TagSequence`, `TagName`)

**PSinfo.txt** — one PCR replicate per line, tab-separated:
```
Sample1	tag1	tag2	1
Sample1	tag3	tag4	2
```
(`SampleName`, `FwdTagName`, `RevTagName`, `PoolNumber`)

## Testing

```bash
# Python unit tests
pytest python/tests/ -v

# Rust unit + integration tests
cargo test --manifest-path rust/Cargo.toml

# Shell integration tests (requires both dame-py and dame on PATH)
bash tests/integration/run_sort.sh
bash tests/integration/run_rsi.sh
bash tests/integration/run_filter.sh
bash tests/integration/run_decollapse.sh
bash tests/integration/run_chimera.sh   # skips if usearch not found
bash tests/integration/run_pipeline.sh
```

## Repository layout

```
python/                          Python 3 implementation (dame-py entry point)
  dame/
    __main__.py                  CLI entry point and subcommand dispatch
    sort.py / modules_sort.py    Demultiplex reads by tag+primer (sort subcommand)
    filter.py / modules_filter.py  Filter amplicons across PCR replicates
    chimera_check.py / modules_chimera_check.py  Chimera detection via usearch
    rsi.py                       Renkonen Similarity Index
    decollapse.py                Expand collapsed sequences back to reads
  tests/                         pytest unit test suite

rust/                            Rust implementation (dame binary)
  src/
    main.rs                      CLI entry point (clap dispatch)
    lib.rs                       Module declarations
    sort.rs                      Demultiplex reads by tag+primer
    filter.rs                    Filter amplicons across PCR replicates
    chimera_check.rs             Chimera detection via usearch
    rsi.rs                       Renkonen Similarity Index
    decollapse.rs                Expand collapsed sequences back to reads
  tests/                         Rust integration tests (one file per subcommand)
  Cargo.toml                     Crate manifest and dependencies

tests/
  fixtures/                      Shared test input files (FASTQ, primers, tags)
  integration/                   Shell scripts comparing dame vs dame-py output

tutorial/                        Synthetic dataset and step-by-step walkthrough
docs/                            Design specs and implementation plans
DAMe_1.0/                        Original DAMe v1.0 (Python 2 scripts, example data, manual)
```

## Citation

Zepeda-Mendoza, M.L., Bohmann, K., Carmona Baez, A. and Gilbert, M.T.P. (2016).
DAMe: A toolkit for the initial processing of datasets with PCR replicates of
double-tagged amplicons for DNA metabarcoding analyses. *BMC Research Notes*,
9(1), 255. https://doi.org/10.1186/s13104-016-2064-9
