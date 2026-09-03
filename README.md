# DAMe v3.0.0: DNA Metabarcoding toolkit

DAMe demultiplexes pooled metabarcoding / eDNA FASTQ reads by primer and
tag sequences (**sort**), optionally removes chimeric sequences (**chimera**),
filters amplicons across PCR replicates (**filter**), converts filtered reads
to USEARCH or sumaclust input (**convert**), computes the Renkonen Similarity
Index between replicates (**rsi**), and expands collapsed sequences back to
individual reads (**decollapse**).  It is available in two
implementations:

| | Python 3 (`python/`) | Rust (`rust/`) |
|---|---|---|
| Requirements | Python ≥ 3.11, numpy | Rust stable (cargo) |
| Install | `pip install -e python/` (or `uv tool install --editable ./python`) | `cd rust && cargo build --release` |
| Entry point | `dame-py` | `dame` |
| Input FASTQ | plain or gzip | plain or gzip |
| Chimera check | via `usearch` on PATH | via `usearch` on PATH |

Reads are expected to have the structure
`[fwd_tag][fwd_primer][amplicon][rc(rev_primer)][rc(rev_tag)]`.
Both forward and reverse-complement orientations are detected automatically
during sort. IUPAC ambiguity codes are supported in primer sequences. The
`--primer-mismatches` and `--tag-mismatches` options tolerate substitutions
only: insertions or deletions in a primer or tag are not tolerated, and such
reads normally remain unassigned. Amplicon length may vary, including through
insertions or deletions within the amplicon. `filter` does not perform mismatch
correction; it compares the amplicon sequences produced by `sort` exactly.

## Performance

DAMe v3 provides Python and Rust implementations developed in parallel.
The figures below are median wall-clock times for complete CLI invocations in
one Linux/aarch64 Docker environment, using two warmups and ten interleaved
measured runs. Outputs were checked for equivalence before accepting results.

### Sort

The default sort workload contained 98,000 synthetic reads:

| Implementation | Median | Throughput | Comparison |
|---|---:|---:|---:|
| Upstream DAMe v1 (Python 2.7) | 335.3 ms | 292.3k reads/s | baseline |
| DAMe Python v3 | 242.2 ms | 404.6k reads/s | **1.4× v1** |
| DAMe Rust v3 | 50.55 ms | 1.939M reads/s | **6.6× v1; 4.8× Python v3** |

Rust gains speed because the work repeated for every read runs as compiled
machine code rather than being stepped through by the Python interpreter.
`needletail` reads FASTQ records efficiently, and Rust compares DNA letters
as simple bytes instead of repeatedly creating and examining Python strings.
It also builds lookup tables for the tag sequences once at startup. That lets
it find a tag directly for each read instead of scanning the full tag list,
while avoiding many short-lived strings and other temporary objects. Python v3
now delays NumPy loading until RSI needs it, reducing startup work for `sort`.

The difference becomes much larger when primer or tag mismatches are allowed.
Python's exact path can use a regular-expression engine implemented in compiled
code, but its mismatch-tolerant path must loop through candidate primers, tags
and DNA positions in Python while counting differences. Python v3 now performs
each IUPAC Hamming comparison in one direct loop, avoiding a generator and a
function call for every DNA letter. Rust still performs the same checks with
compiled byte-level loops at fixed positions in each read. On this dataset it
was about 28× faster on the mismatch-enabled paths; the exact ratio can change
with the tag panel, read layout and permitted mismatch count. Upstream DAMe v1
has no mismatch capability and is therefore not included in this table.

| v3 sort mode | Python v3 | Rust v3 | Rust speedup |
|---|---:|---:|---:|
| Default / exact | 404.6k reads/s | 1.939M reads/s | 4.8× |
| One primer mismatch | 79.75k reads/s | 2.223M reads/s | 27.9× |
| One tag mismatch | 80.14k reads/s | 2.219M reads/s | 27.7× |

### Filter

The scaled filter workload contained 40,000 collapsed input records:

| Implementation | Median | Throughput | Comparison |
|---|---:|---:|---:|
| Upstream DAMe v1 (Python 2.7) | 296.6 ms | 134.9k records/s | baseline |
| DAMe Python v3 | 86.35 ms | 463.2k records/s | **3.4× v1; 1.2× Rust v3** |
| DAMe Rust v3 | 106.0 ms | 377.4k records/s | **2.8× v1** |

Filtering compares collapsed sequences across PCR replicates, applies count and
length thresholds, and writes several result files. Both v3 implementations now
index counts by sequence instead of repeatedly scanning each replicate. Python
also avoids importing NumPy for `filter`, so on this workload it finishes 1.2×
faster than Rust despite deterministic output sorting. Rust still benefits from
compiled loops and low per-record overhead. Tiny inputs are omitted because
process startup dominates them.

See the full [sort report](docs/performance-sort-v1-v3.md) and
[filter report](docs/performance-filter-v1-v3.md) for methodology and scope.

## Quick start

### Python version

`dame-py` accepts both the `--long` flags shown below (matching the Rust
binary) and the original single-dash forms (`-fq`, `-p`, `-t`, `-psInfo`,
`-x`, …) for backward compatibility.

Install with pip, or with [uv](https://docs.astral.sh/uv/) (which puts `dame-py`
in its own isolated, uv-managed environment available on your `PATH`):

```bash
pip install -e python/                 # standard pip (editable)
# or:
uv tool install --editable ./python    # via uv — run `dame-py` from anywhere
#   update deps later: uv tool upgrade dame-py
#   uninstall:         uv tool uninstall dame-py

dame-py sort \
  --fq Pool1.fastq \
  --primers Primers.txt \
  --tags Tags.txt \
  --primer-mismatches 1

dame-py filter \
  --ps-info PSinfo.txt \
  --x 2 --y 2 --t 2 --l 50

dame-py rsi Comparisons_2PCRs.txt
```

### Rust version

```bash
cd rust && cargo build --release
# binary is here: DAMe/rust/target/release/dame

dame sort \
  --fq Pool1.fastq \
  --primers Primers.txt \
  --tags Tags.txt \
  --primer-mismatches 1

dame filter \
  --ps-info PSinfo.txt \
  --x 2 --y 2 --t 2 --l 50

dame rsi Comparisons_2PCRs.txt
```

## Pipeline overview

```
dame sort     -fq POOL.fastq --primers P.txt --tags T.txt [--primer-mismatches N] [--tag-mismatches N]
              → TagA_TagB.txt (collapsed unique seqs + counts) per tag pair
              → SummaryCounts.txt

dame chimera  --ps-info PSinfo.txt --x 2          # requires usearch on PATH
              → TagA_TagB_Pool.noChim.txt

dame filter   --ps-info PSinfo.txt --x 2 --y 2 --t 2 --l 50
              → Comparisons_2PCRs.txt (all seqs, all replicates)
              → FilteredReads.fna   (passed all thresholds)

dame convert  -i FilteredReads.fna [-u] [--min-length N] [--max-length N] [-s]
              → FilteredReads.forsumaclust.fna  (sumaclust input, default)
              → FilteredReads.forusearch.fna    (USEARCH input, with -u)
              → SampleFastas/<Sample>.fixed.fasta  (per-sample, with -s)

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
   covering all six subcommands with synthetic inputs generated via `tmp_path`
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
   pytest on Python 3.11 and 3.12, `cargo test` on Rust stable, and all seven
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

10. **DAMe v2.4 — Configurable primer mismatches.**  `sort` gained a
    `-m`/`--primer-mismatches N` option (default 0) that tolerates up to N
    substitutions per primer match, IUPAC-aware, using leftmost-within-budget
    selection.  The budget applies independently to each of the four primer
    sites (forward/reverse orientation × start/end).  Tags are still matched
    exactly.  At N=0 the output is byte-identical to v2.3, verified by the
    existing integration tests; a new `run_sort_mismatch.sh` checks both
    implementations agree at N=1.  The Python *mismatch* matcher is a manual
    IUPAC sliding window mirroring Rust; the exact (no-mismatch) path keeps the
    faster compiled-regex matcher (see Performance).

11. **Polish and community-PR adoption.**  Several smaller improvements landed
    alongside v2.5, some adapted from a community PR (dougwyu/DAMe#1, @jiyinqiu):
    `dame-py` now reads gzip FASTQ transparently (`utils.smart_open`), matching
    the Rust binary; the Python CLI accepts `--long` flag aliases (`--fq`,
    `--primers`, `--ps-info`, …) alongside the original single-dash forms;
    `filter` output is deterministic (sequences are sorted, making the two
    implementations byte-identical and letting the integration test compare them
    in full); reads are uppercased before matching so soft-masked bases are not
    silently dropped; malformed `decollapse` rows are skipped instead of
    crashing; and the dead standalone `main()` entry points were removed.  A
    reproducible `benchmark/` harness was added; benchmarking it caught and
    fixed a ~10× regression in the Python `sort` default path (the regex matcher
    had been replaced wholesale for mismatch support — it is now used again for
    the no-mismatch path, see Performance).

12. **DAMe v2.5 — Tag mismatches + anchored matching.**  `sort` gained
    `--tag-mismatches N` (`-mt` in Python; default 0), a per-tag substitution
    tolerance.  When `--primer-mismatches` or `--tag-mismatches` is non-zero,
    sort uses a tag-anchored matcher: it finds tag candidates at the read ends
    by IUPAC Hamming distance, checks the primers at the expected offsets,
    scores each valid assembly by total mismatches, keeps the unique minimum,
    and discards ambiguous ties.  At the defaults (both 0) the original exact
    matcher runs unchanged (byte-identical).  A startup warning flags an unsafe
    `--tag-mismatches` relative to the tag set's minimum Hamming distance.
    Design adapted from a community PR (dougwyu/DAMe#1, @jiyinqiu).

13. **DAMe v2.6 — Convert subcommand.**  A new `convert` subcommand converts
    `FilteredReads.fna` (DAMe filter output) to USEARCH (`>Sample;size=N`) or
    sumaclust (`>Sample:N count=N`) input format, with optional length filtering
    (`--min-length`, `--max-length`), N-padding to a fixed width in USEARCH
    mode, and per-sample FASTA output (`--sample-fastas` → `SampleFastas/`).
    Fixes three bugs present in the original `convertToUSearch.py` v1.0 script:
    a float default for `lmax` (`1e6`) that broke integer comparisons, a missing
    `os.makedirs` call that crashed when `SampleFastas/` did not exist, and an
    off-by-one in the N-padding width.  The Python port accepts the original v1.0
    flag spellings (`--inFasta`, `-lmin`, `-lmax`, `--sampleFastas`) for backward
    compatibility.

14. **DAMe v3.0.0 — Input robustness and cross-implementation parity.**  A
    review of both implementations found six defects where they silently
    disagreed, each reproduced before and after the fix.  `dame-py` stopped
    reading a FASTQ at the first blank sequence line, discarding every
    remaining read while reporting success; on a 20-record file that lost 12 of
    14 reads.  Eleven Rust parsers split input files on tab alone where
    `dame-py` and DAMe v1.0 split on any whitespace, so a single stray space
    beside a tab in `Tags.txt` cut a 292-read sort to 61, silently and with
    exit code 0.  `dame rsi` looped forever on a file with fewer than two
    columns, a `usize` underflow defeating its own guard.  `dame-py` raised
    `IndexError` on a short `convert` header and on a blank line in `PSinfo`.
    Python's `rsi` row order came from a set, so it varied between runs of the
    same input and never matched the Rust output.

    Two **breaking changes**: `sort` now refuses a `Tags.txt` or `Primers.txt`
    that is malformed rather than resolving it arbitrarily.  Previously a
    repeated tag sequence made the assigned tag name, and so the output
    filename, depend on lookup order (`dame-py` took the first, `dame` the
    last); a repeated tag name was matched by `dame` on its exact path but not
    its anchored one, so a single binary gave different answers depending on
    whether `--primer-mismatches` was set; and a repeated primer set name meant
    `dame-py` sorted 4 reads where `dame` matched nothing at all.  Tag names,
    tag sequences and primer set names must now each be unique, and every line
    must carry all its fields.  Files that previously sorted with one of
    several arbitrary answers now fail with a clear message.

    `dame rsi` also formats values as Python does (`0.0`, not `0`; `1e-05`, not
    `1e-5`), verified over 21,424 values.  The Rust CLI gained help text for
    every flag in all six subcommands, and `filter --p` is documented as
    accepted but ignored — pooling is driven by column 4 of `PSinfo`.

    The root cause of the whole class was a test gap: every fixture in the repo
    was well formed, which is precisely where the two implementations agree.
    `tests/integration/run_malformed.sh` now feeds both binaries damaged input
    (a stray space, a blank line, a short header, a truncated record, a
    single-column file) and each case was verified to fail when its fix is
    reverted.  The suite grew from 59 pytest / 43 cargo tests to 74 / 84.

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

Every line must give all three fields and each primer set name must be unique;
`sort` refuses a Primers file that breaks either rule.  Two sets may share a
forward or reverse sequence, which is a legitimate multiplex design.

**Tags.txt** — one tag per line, tab-separated:
```
AACCGGT	tag1
TTGGCCA	tag2
```
(`TagSequence`, `TagName`)

Every line must give both fields, and names and sequences must be one to one:
`sort` refuses a Tags file with an incomplete entry, or one that repeats a name
or a sequence, naming the line and both conflicting entries.  A repeated sequence
would leave the tag a read is assigned to dependent on lookup order, and that
name becomes an output filename.  A repeated name is worse: only the first
sequence for it is ever matched, so reads carrying the others are silently
counted as erroneous.

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
bash tests/integration/run_sort_mismatch.sh       # --primer-mismatches parity
bash tests/integration/run_sort_tag_mismatch.sh   # --tag-mismatches parity
bash tests/integration/run_rsi.sh
bash tests/integration/run_filter.sh
bash tests/integration/run_convert.sh
bash tests/integration/run_malformed.sh        # damaged-input parity
bash tests/integration/run_decollapse.sh
bash tests/integration/run_chimera.sh   # skips if usearch not found
bash tests/integration/run_pipeline.sh

# Sort throughput benchmark (Python vs Rust)
bash benchmark/run_sort_benchmark.sh
```

## Repository layout

```
python/                          Python 3 implementation (dame-py entry point)
  dame/
    __main__.py                  CLI entry point and subcommand dispatch
    sort.py / modules_sort.py    Demultiplex reads by tag+primer (sort subcommand)
    filter.py / modules_filter.py  Filter amplicons across PCR replicates
    chimera_check.py / modules_chimera_check.py  Chimera detection via usearch
    convert.py                   Convert FilteredReads.fna to USEARCH/sumaclust format
    rsi.py                       Renkonen Similarity Index
    decollapse.py                Expand collapsed sequences back to reads
    utils.py                     Shared helpers (transparent gzip input)
  tests/                         pytest unit test suite

rust/                            Rust implementation (dame binary)
  src/
    main.rs                      CLI entry point (clap dispatch)
    lib.rs                       Module declarations
    sort.rs                      Demultiplex reads by tag+primer
    filter.rs                    Filter amplicons across PCR replicates
    chimera_check.rs             Chimera detection via usearch
    convert.rs                   Convert FilteredReads.fna to USEARCH/sumaclust format
    rsi.rs                       Renkonen Similarity Index
    decollapse.rs                Expand collapsed sequences back to reads
  tests/                         Rust integration tests (one file per subcommand)
  Cargo.toml                     Crate manifest and dependencies

tests/
  fixtures/                      Shared test input files (FASTQ, primers, tags)
  integration/                   Shell scripts comparing dame vs dame-py output

tutorial/                        Synthetic dataset and step-by-step walkthrough
benchmark/                       Reproducible sort throughput benchmark
docs/                            Design specs and implementation plans
DAMe_1.0/                        Original DAMe v1.0 (Python 2 scripts, example data, manual)
```

## Citation

Zepeda-Mendoza, M.L., Bohmann, K., Carmona Baez, A. and Gilbert, M.T.P. (2016).
DAMe: A toolkit for the initial processing of datasets with PCR replicates of
double-tagged amplicons for DNA metabarcoding analyses. *BMC Research Notes*,
9(1), 255. https://doi.org/10.1186/s13104-016-2064-9
