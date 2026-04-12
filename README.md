# DAMe v2.0: DNA Metabarcoding toolkit

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

Benchmarked on the tutorial dataset (392 single-end reads, 1 pool):

| Step | Python 3 | Rust | Speedup |
|------|----------|------|---------|
| sort | ~280 ms | ~38 ms | ~7× |

Rust startup overhead dominates at small scale; the speedup grows substantially
on datasets with tens of thousands of reads per pool.

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

DAMe was originally written in Python 2 by Zepeda Mendoza et al. (2016).
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
python/      Python 3 source (dame-py entry point)
rust/        Rust source (dame binary)
tests/       Developer and CI tests
tutorial/    Tutorial dataset and walkthrough
docs/        Design specs and implementation plans
DAMe_1.0/    Original DAMe v1.0 (Python 2 scripts, example data, manual)
```

## Citation

Zepeda Mendoza, M.L., Sicheritz-Pontén, T. and Gilbert, M.T.P. (2016).
Environmental genes and genomes: understanding the differences and challenges in
the approaches and software for their analyses. *Briefings in Bioinformatics*,
17(4), 745–754.
