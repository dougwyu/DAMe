# convert subcommand Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add `dame convert` / `dame-py convert` that converts `FilteredReads.fna` to USEARCH or sumaclust format, porting and fixing `DAMe_1.0/bin/convertToUSearch.py`.

**Architecture:** Single-file module per implementation (`convert.py` / `convert.rs`) with a testable pure function (`convert` / `process`) called by a thin `run(args)` dispatcher. Python registers backwards-compat aliases; Rust uses canonical names only.

**Tech Stack:** Python 3.11+ stdlib only; Rust (anyhow, clap derive, std::io). No new dependencies.

---

## File map

| File | Action | Responsibility |
|------|--------|---------------|
| `tests/fixtures/FilteredReads_small.fna` | Create | Shared test fixture (4 records, 2 samples) |
| `python/dame/convert.py` | Create | Python convert logic + CLI registration |
| `python/tests/test_convert.py` | Create | Python unit tests |
| `python/dame/__main__.py` | Modify | Register convert subcommand |
| `rust/src/convert.rs` | Create | Rust convert logic + inline unit tests |
| `rust/src/lib.rs` | Modify | Declare convert module |
| `rust/src/main.rs` | Modify | Add Convert subcommand dispatch |
| `rust/Cargo.toml` | Modify | Bump version 0.6.0 → 0.7.0 |
| `python/pyproject.toml` | Modify | Bump version 3.0.0 → 3.1.0 |
| `tests/integration/run_convert.sh` | Create | Cross-impl parity test |
| `.github/workflows/ci.yml` | Modify | Add run_convert.sh step |
| `README.md` | Modify | Add convert to pipeline overview |
| `tutorial/README.md` | Modify | Add convert step after filter |

---

### Task C1: Create test fixture

**Files:**
- Create: `tests/fixtures/FilteredReads_small.fna`

- [ ] **Step 1: Write the fixture file**

```
tests/fixtures/FilteredReads_small.fna
```

Content (4 records: 2 for Sample1, 2 for Sample2; last record is 4 nt to test length filtering):

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

Sizes: Sample1 record 1 → 5+4=9; record 2 → 3+2=5; Sample2 record 3 → 10+8=18; record 4 → 1+0=1.
Sequence lengths: records 1–3 are 60 nt; record 4 is 4 nt.

- [ ] **Step 2: Commit**

```bash
git add tests/fixtures/FilteredReads_small.fna
git commit -m "test(fixture): add FilteredReads_small.fna for convert tests"
```

---

### Task C2: Python unit tests (write failing)

**Files:**
- Create: `python/tests/test_convert.py`

- [ ] **Step 1: Write the test file**

```python
# python/tests/test_convert.py
import os
import pytest
from dame.convert import convert, _parse_fasta


SMALL_FNA = """\
>Sample1 Tag1-Tag2.Tag3-Tag4_1 5_4
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
>Sample1 Tag1-Tag2.Tag3-Tag4_2 3_2
GCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCG
>Sample2 Tag5-Tag6.Tag7-Tag8_1 10_8
TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
>Sample2 Tag5-Tag6.Tag7-Tag8_2 1_0
AAAA
"""


def write_fna(tmp_path):
    p = tmp_path / "FilteredReads_small.fna"
    p.write_text(SMALL_FNA)
    return str(p)


def test_parse_fasta(tmp_path):
    fna = write_fna(tmp_path)
    records = list(_parse_fasta(fna))
    assert len(records) == 4
    assert records[0] == (
        "Sample1", 9,
        "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG",
    )
    assert records[1] == (
        "Sample1", 5,
        "GCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCG",
    )
    assert records[2] == ("Sample2", 18, "T" * 60)
    assert records[3] == ("Sample2", 1, "AAAA")


def test_sumaclust_output(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    fna = write_fna(tmp_path)
    out = convert(fna)
    assert out == "FilteredReads.forsumaclust.fna"
    lines = open(out).readlines()
    assert lines[0] == ">Sample1:1 count=9\n"
    assert lines[1] == "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\n"
    assert lines[2] == ">Sample1:2 count=5\n"
    assert lines[3] == "GCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCG\n"
    assert lines[4] == ">Sample2:3 count=18\n"
    assert lines[5] == "T" * 60 + "\n"
    assert lines[6] == ">Sample2:4 count=1\n"
    assert lines[7] == "AAAA\n"
    assert len(lines) == 8


def test_usearch_output(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    fna = write_fna(tmp_path)
    out = convert(fna, usearch=True)
    assert out == "FilteredReads.forusearch.fna"
    lines = open(out).readlines()
    assert lines[0] == ">Sample1;size=9\n"
    assert lines[1] == "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\n"
    assert lines[6] == ">Sample2;size=1\n"
    assert lines[7] == "AAAA\n"
    assert len(lines) == 8


def test_usearch_padding(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    fna = write_fna(tmp_path)
    out = convert(fna, usearch=True, max_length=65)
    lines = open(out).readlines()
    # 60-nt seq padded to 65
    assert lines[1].rstrip("\n") == "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG" + "N" * 5
    assert len(lines[1].rstrip("\n")) == 65
    # 4-nt seq padded to 65
    assert lines[7].rstrip("\n") == "AAAA" + "N" * 61
    assert len(lines[7].rstrip("\n")) == 65


def test_no_padding_in_sumaclust_mode(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    fna = write_fna(tmp_path)
    out = convert(fna, max_length=65)   # max_length set but NOT usearch → no padding
    lines = open(out).readlines()
    assert lines[1].rstrip("\n") == "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"
    assert len(lines[1].rstrip("\n")) == 60


def test_min_length_filter(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    fna = write_fna(tmp_path)
    out = convert(fna, min_length=10)
    lines = open(out).readlines()
    # 4-nt AAAA record dropped; 3 records × 2 lines = 6 lines
    assert len(lines) == 6
    assert all("Sample2:4" not in l and "AAAA" not in l for l in lines)


def test_max_length_filter_sumaclust(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    fna = write_fna(tmp_path)
    out = convert(fna, max_length=10)
    lines = open(out).readlines()
    # Only AAAA (4 nt) passes max_length=10; counter is 1 for first passing record
    assert len(lines) == 2
    assert lines[0] == ">Sample2:1 count=1\n"
    assert lines[1] == "AAAA\n"


def test_sample_fastas_creates_directory(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    fna = write_fna(tmp_path)
    convert(fna, sample_fastas=True)
    assert os.path.isdir("SampleFastas")
    assert os.path.isfile("SampleFastas/Sample1.fixed.fasta")
    assert os.path.isfile("SampleFastas/Sample2.fixed.fasta")


def test_sample_fastas_content(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    fna = write_fna(tmp_path)
    convert(fna, sample_fastas=True)
    s1 = open("SampleFastas/Sample1.fixed.fasta").readlines()
    # Sample1 has 2 records
    assert len(s1) == 4
    assert s1[0] == ">Sample1:1 count=9\n"
    assert s1[2] == ">Sample1:2 count=5\n"
    s2 = open("SampleFastas/Sample2.fixed.fasta").readlines()
    assert len(s2) == 4
    assert s2[0] == ">Sample2:3 count=18\n"


def test_convert_argparser(tmp_path):
    import argparse
    import dame.convert as conv
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers()
    conv.register_subcommand(sub)
    # Canonical flags
    args = parser.parse_args(["convert", "-i", "x.fna", "--min-length", "5",
                               "--max-length", "100", "-u", "-s"])
    assert args.in_fasta == "x.fna"
    assert args.min_length == 5
    assert args.max_length == 100
    assert args.usearch is True
    assert args.sample_fastas is True


def test_convert_argparser_legacy_flags(tmp_path):
    import argparse
    import dame.convert as conv
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers()
    conv.register_subcommand(sub)
    # Legacy v1.0 flag spellings
    args = parser.parse_args(["convert", "--inFasta", "x.fna",
                               "-lmin", "5", "-lmax", "100",
                               "--sampleFastas"])
    assert args.in_fasta == "x.fna"
    assert args.min_length == 5
    assert args.max_length == 100
    assert args.sample_fastas is True
```

- [ ] **Step 2: Run tests to verify they fail (module doesn't exist yet)**

```bash
cd /path/to/DAMe
python -m pytest python/tests/test_convert.py -v 2>&1 | head -20
```

Expected: `ModuleNotFoundError: No module named 'dame.convert'`

---

### Task C3: Python implementation

**Files:**
- Create: `python/dame/convert.py`

- [ ] **Step 1: Write the module**

```python
# python/dame/convert.py
import os


def _parse_fasta(path):
    """Yield (sample, size, sequence) tuples from a FilteredReads.fna file."""
    with open(path) as fh:
        header = None
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                toks = line.split()
                sample = toks[0][1:]
                size = sum(int(x) for x in toks[2].split("_"))
                header = (sample, size)
            elif header is not None:
                yield header[0], header[1], line
                header = None


def convert(in_fasta, min_length=0, max_length=None, usearch=False, sample_fastas=False):
    """
    Convert FilteredReads.fna to USEARCH or sumaclust format.

    Returns the path of the main output file created.
    """
    out_name = "FilteredReads.forusearch.fna" if usearch else "FilteredReads.forsumaclust.fna"

    if sample_fastas:
        os.makedirs("SampleFastas", exist_ok=True)

    sample_handles = {}
    counter = 1

    with open(out_name, "w") as out:
        for sample, size, seq in _parse_fasta(in_fasta):
            if len(seq) < min_length:
                continue
            if max_length is not None and len(seq) > max_length:
                continue

            if usearch:
                hdr = f">{sample};size={size}"
                out_seq = seq.ljust(max_length, "N") if max_length is not None else seq
            else:
                hdr = f">{sample}:{counter} count={size}"
                out_seq = seq
                counter += 1

            out.write(hdr + "\n" + out_seq + "\n")

            if sample_fastas:
                if sample not in sample_handles:
                    sample_handles[sample] = open(
                        f"SampleFastas/{sample}.fixed.fasta", "w"
                    )
                sample_handles[sample].write(hdr + "\n" + out_seq + "\n")

    for fh in sample_handles.values():
        fh.close()

    return out_name


def register_subcommand(subparsers):
    p = subparsers.add_parser(
        "convert",
        description="Convert FilteredReads.fna to USEARCH or sumaclust input format",
    )
    p.add_argument(
        "-i", "--in-fasta", "--inFasta",
        dest="in_fasta", required=True, metavar="FILE",
        help="Input FilteredReads.fna file",
    )
    p.add_argument(
        "--min-length", "-lmin", "--minLength",
        dest="min_length", type=int, default=0, metavar="N",
        help="Drop sequences shorter than N [default 0]",
    )
    p.add_argument(
        "--max-length", "-lmax", "--maxLength",
        dest="max_length", type=int, default=None, metavar="N",
        help="Drop sequences longer than N; pad to N in USEARCH mode",
    )
    p.add_argument(
        "-u", "--usearch",
        dest="usearch", action="store_true",
        help="Write USEARCH format (default: sumaclust)",
    )
    p.add_argument(
        "-s", "--sample-fastas", "--sampleFastas",
        dest="sample_fastas", action="store_true",
        help="Write per-sample fastas to SampleFastas/",
    )
    p.set_defaults(func=run)


def run(args):
    convert(
        in_fasta=args.in_fasta,
        min_length=args.min_length,
        max_length=args.max_length,
        usearch=args.usearch,
        sample_fastas=args.sample_fastas,
    )
```

- [ ] **Step 2: Run tests**

```bash
python -m pytest python/tests/test_convert.py -v
```

Expected: all tests PASS.

- [ ] **Step 3: Commit**

```bash
git add python/dame/convert.py python/tests/test_convert.py
git commit -m "feat(python): add convert subcommand"
```

---

### Task C4: Register Python convert subcommand

**Files:**
- Modify: `python/dame/__main__.py`

- [ ] **Step 1: Add import and registration**

In `python/dame/__main__.py`, add `import dame.convert as convert_mod` alongside the existing imports, and call `convert_mod.register_subcommand(subparsers)` in `main()`.

Current file:
```python
import argparse

from dame import sort, chimera_check, decollapse, rsi
import dame.filter as filter_mod


def main():
    parser = argparse.ArgumentParser(
        prog="dame-py",
        description="DAMe: DNA Metabarcoding pipeline toolkit",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    sort.register_subcommand(subparsers)
    chimera_check.register_subcommand(subparsers)
    filter_mod.register_subcommand(subparsers)
    decollapse.register_subcommand(subparsers)
    rsi.register_subcommand(subparsers)

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
```

Replace with:
```python
import argparse

from dame import sort, chimera_check, decollapse, rsi
import dame.filter as filter_mod
import dame.convert as convert_mod


def main():
    parser = argparse.ArgumentParser(
        prog="dame-py",
        description="DAMe: DNA Metabarcoding pipeline toolkit",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    sort.register_subcommand(subparsers)
    chimera_check.register_subcommand(subparsers)
    filter_mod.register_subcommand(subparsers)
    decollapse.register_subcommand(subparsers)
    rsi.register_subcommand(subparsers)
    convert_mod.register_subcommand(subparsers)

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
```

- [ ] **Step 2: Smoke test**

```bash
dame-py convert --help
```

Expected: prints help with `-i`, `--min-length`, `--max-length`, `-u`, `-s` flags.

- [ ] **Step 3: Run full Python test suite**

```bash
python -m pytest python/tests/ -v
```

Expected: all tests PASS.

- [ ] **Step 4: Commit**

```bash
git add python/dame/__main__.py
git commit -m "feat(python): register convert subcommand in CLI"
```

---

### Task C5: Rust unit tests (write failing)

**Files:**
- Create: `rust/src/convert.rs` (skeleton only, enough to make module exist)

- [ ] **Step 1: Create stub convert.rs with failing tests**

```rust
// rust/src/convert.rs
use anyhow::Result;
use clap::Args;
use std::collections::HashMap;
use std::fs::{self, File};
use std::io::{BufRead, BufReader, BufWriter, Write};

#[derive(Args)]
pub struct ConvertArgs {
    #[arg(short = 'i', long = "in-fasta")]
    pub in_fasta: String,
    #[arg(long = "min-length", default_value_t = 0)]
    pub min_length: usize,
    #[arg(long = "max-length")]
    pub max_length: Option<usize>,
    #[arg(short = 'u', long = "usearch")]
    pub usearch: bool,
    #[arg(short = 's', long = "sample-fastas")]
    pub sample_fastas: bool,
}

pub fn run(args: ConvertArgs) -> Result<()> {
    todo!()
}

fn process<R: BufRead, W: Write>(
    _reader: R,
    _writer: W,
    _min_length: usize,
    _max_length: Option<usize>,
    _usearch: bool,
    _sample_dir: Option<&str>,
) -> Result<()> {
    todo!()
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    const SMALL_FNA: &str = "\
>Sample1 Tag1-Tag2.Tag3-Tag4_1 5_4\n\
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\n\
>Sample1 Tag1-Tag2.Tag3-Tag4_2 3_2\n\
GCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCG\n\
>Sample2 Tag5-Tag6.Tag7-Tag8_1 10_8\n\
TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT\n\
>Sample2 Tag5-Tag6.Tag7-Tag8_2 1_0\n\
AAAA\n";

    #[test]
    fn test_sumaclust_output() {
        let mut out = Vec::new();
        process(Cursor::new(SMALL_FNA), &mut out, 0, None, false, None).unwrap();
        let s = String::from_utf8(out).unwrap();
        let lines: Vec<&str> = s.lines().collect();
        assert_eq!(lines[0], ">Sample1:1 count=9");
        assert_eq!(lines[1], "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG");
        assert_eq!(lines[2], ">Sample1:2 count=5");
        assert_eq!(lines[4], ">Sample2:3 count=18");
        assert_eq!(lines[6], ">Sample2:4 count=1");
        assert_eq!(lines[7], "AAAA");
        assert_eq!(lines.len(), 8);
    }

    #[test]
    fn test_usearch_output() {
        let mut out = Vec::new();
        process(Cursor::new(SMALL_FNA), &mut out, 0, None, true, None).unwrap();
        let s = String::from_utf8(out).unwrap();
        let lines: Vec<&str> = s.lines().collect();
        assert_eq!(lines[0], ">Sample1;size=9");
        assert_eq!(lines[1], "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG");
        assert_eq!(lines[6], ">Sample2;size=1");
        assert_eq!(lines[7], "AAAA");
        assert_eq!(lines.len(), 8);
    }

    #[test]
    fn test_usearch_padding() {
        let mut out = Vec::new();
        process(Cursor::new(SMALL_FNA), &mut out, 0, Some(65), true, None).unwrap();
        let s = String::from_utf8(out).unwrap();
        let lines: Vec<&str> = s.lines().collect();
        assert_eq!(lines[1].len(), 65);
        assert!(lines[1].ends_with("NNNNN"));
        assert_eq!(lines[7].len(), 65);
        assert!(lines[7].starts_with("AAAA"));
        assert!(lines[7].ends_with('N'));
    }

    #[test]
    fn test_no_padding_in_sumaclust_mode() {
        let mut out = Vec::new();
        process(Cursor::new(SMALL_FNA), &mut out, 0, Some(65), false, None).unwrap();
        let s = String::from_utf8(out).unwrap();
        let lines: Vec<&str> = s.lines().collect();
        assert_eq!(lines[1].len(), 60);
    }

    #[test]
    fn test_min_length_filter() {
        let mut out = Vec::new();
        process(Cursor::new(SMALL_FNA), &mut out, 10, None, false, None).unwrap();
        let s = String::from_utf8(out).unwrap();
        let lines: Vec<&str> = s.lines().collect();
        assert_eq!(lines.len(), 6);
        assert!(!s.contains("Sample2:4"));
        assert!(!s.contains("AAAA"));
    }

    #[test]
    fn test_max_length_filter() {
        let mut out = Vec::new();
        process(Cursor::new(SMALL_FNA), &mut out, 0, Some(10), false, None).unwrap();
        let s = String::from_utf8(out).unwrap();
        let lines: Vec<&str> = s.lines().collect();
        assert_eq!(lines.len(), 2);
        assert_eq!(lines[0], ">Sample2:1 count=1");
        assert_eq!(lines[1], "AAAA");
    }
}
```

- [ ] **Step 2: Add convert to lib.rs (so it compiles)**

In `rust/src/lib.rs`, add `pub mod convert;` alongside the existing module declarations:

```rust
pub mod chimera_check;
pub mod convert;
pub mod decollapse;
pub mod filter;
pub mod rsi;
pub mod sort;
```

- [ ] **Step 3: Run tests to verify they fail (todo! panics)**

```bash
cd rust && cargo test convert -- --nocapture 2>&1 | tail -20
```

Expected: tests panic with `not yet implemented`.

---

### Task C6: Rust implementation

**Files:**
- Modify: `rust/src/convert.rs`

- [ ] **Step 1: Implement process and run**

Replace the stub `process` and `run` functions (keep `ConvertArgs`, `#[cfg(test)]` block unchanged):

```rust
pub fn run(args: ConvertArgs) -> Result<()> {
    let reader = BufReader::new(File::open(&args.in_fasta)?);
    let out_name = if args.usearch {
        "FilteredReads.forusearch.fna"
    } else {
        "FilteredReads.forsumaclust.fna"
    };
    let writer = BufWriter::new(File::create(out_name)?);

    let sample_dir = if args.sample_fastas {
        fs::create_dir_all("SampleFastas")?;
        Some("SampleFastas")
    } else {
        None
    };

    process(reader, writer, args.min_length, args.max_length, args.usearch, sample_dir)
}

fn process<R: BufRead, W: Write>(
    reader: R,
    mut writer: W,
    min_length: usize,
    max_length: Option<usize>,
    usearch: bool,
    sample_dir: Option<&str>,
) -> Result<()> {
    let mut sample_handles: HashMap<String, BufWriter<File>> = HashMap::new();
    let mut counter: u64 = 1;

    let mut lines = reader.lines();
    while let Some(hdr_line) = lines.next() {
        let hdr_line = hdr_line?;
        if !hdr_line.starts_with('>') {
            continue;
        }
        let seq_line = match lines.next() {
            Some(l) => l?,
            None => break,
        };
        let seq = seq_line.trim_end();

        // Parse: >Sample TagPair counts_underscore
        let toks: Vec<&str> = hdr_line.splitn(4, ' ').collect();
        if toks.len() < 3 {
            continue;
        }
        let sample = &toks[0][1..];
        let size: u64 = toks[2]
            .split('_')
            .filter_map(|x| x.parse::<u64>().ok())
            .sum();

        if seq.len() < min_length {
            continue;
        }
        if let Some(max_len) = max_length {
            if seq.len() > max_len {
                continue;
            }
        }

        let out_hdr: String;
        let out_seq: String;
        if usearch {
            out_hdr = format!(">{};size={}", sample, size);
            out_seq = if let Some(max_len) = max_length {
                let mut s = seq.to_string();
                while s.len() < max_len {
                    s.push('N');
                }
                s
            } else {
                seq.to_string()
            };
        } else {
            out_hdr = format!(">{}:{} count={}", sample, counter, size);
            out_seq = seq.to_string();
            counter += 1;
        }

        writeln!(writer, "{}", out_hdr)?;
        writeln!(writer, "{}", out_seq)?;

        if let Some(dir) = sample_dir {
            if !sample_handles.contains_key(sample) {
                let path = format!("{}/{}.fixed.fasta", dir, sample);
                let fh = BufWriter::new(File::create(path)?);
                sample_handles.insert(sample.to_string(), fh);
            }
            let fh = sample_handles.get_mut(sample).unwrap();
            writeln!(fh, "{}", out_hdr)?;
            writeln!(fh, "{}", out_seq)?;
        }
    }

    Ok(())
}
```

The full `rust/src/convert.rs` after this step (imports + ConvertArgs + run + process + tests):

```rust
use anyhow::Result;
use clap::Args;
use std::collections::HashMap;
use std::fs::{self, File};
use std::io::{BufRead, BufReader, BufWriter, Write};

#[derive(Args)]
pub struct ConvertArgs {
    #[arg(short = 'i', long = "in-fasta")]
    pub in_fasta: String,
    #[arg(long = "min-length", default_value_t = 0)]
    pub min_length: usize,
    #[arg(long = "max-length")]
    pub max_length: Option<usize>,
    #[arg(short = 'u', long = "usearch")]
    pub usearch: bool,
    #[arg(short = 's', long = "sample-fastas")]
    pub sample_fastas: bool,
}

pub fn run(args: ConvertArgs) -> Result<()> {
    let reader = BufReader::new(File::open(&args.in_fasta)?);
    let out_name = if args.usearch {
        "FilteredReads.forusearch.fna"
    } else {
        "FilteredReads.forsumaclust.fna"
    };
    let writer = BufWriter::new(File::create(out_name)?);

    let sample_dir = if args.sample_fastas {
        fs::create_dir_all("SampleFastas")?;
        Some("SampleFastas")
    } else {
        None
    };

    process(reader, writer, args.min_length, args.max_length, args.usearch, sample_dir)
}

fn process<R: BufRead, W: Write>(
    reader: R,
    mut writer: W,
    min_length: usize,
    max_length: Option<usize>,
    usearch: bool,
    sample_dir: Option<&str>,
) -> Result<()> {
    let mut sample_handles: HashMap<String, BufWriter<File>> = HashMap::new();
    let mut counter: u64 = 1;

    let mut lines = reader.lines();
    while let Some(hdr_line) = lines.next() {
        let hdr_line = hdr_line?;
        if !hdr_line.starts_with('>') {
            continue;
        }
        let seq_line = match lines.next() {
            Some(l) => l?,
            None => break,
        };
        let seq = seq_line.trim_end();

        let toks: Vec<&str> = hdr_line.splitn(4, ' ').collect();
        if toks.len() < 3 {
            continue;
        }
        let sample = &toks[0][1..];
        let size: u64 = toks[2]
            .split('_')
            .filter_map(|x| x.parse::<u64>().ok())
            .sum();

        if seq.len() < min_length {
            continue;
        }
        if let Some(max_len) = max_length {
            if seq.len() > max_len {
                continue;
            }
        }

        let out_hdr: String;
        let out_seq: String;
        if usearch {
            out_hdr = format!(">{};size={}", sample, size);
            out_seq = if let Some(max_len) = max_length {
                let mut s = seq.to_string();
                while s.len() < max_len {
                    s.push('N');
                }
                s
            } else {
                seq.to_string()
            };
        } else {
            out_hdr = format!(">{}:{} count={}", sample, counter, size);
            out_seq = seq.to_string();
            counter += 1;
        }

        writeln!(writer, "{}", out_hdr)?;
        writeln!(writer, "{}", out_seq)?;

        if let Some(dir) = sample_dir {
            if !sample_handles.contains_key(sample) {
                let path = format!("{}/{}.fixed.fasta", dir, sample);
                let fh = BufWriter::new(File::create(path)?);
                sample_handles.insert(sample.to_string(), fh);
            }
            let fh = sample_handles.get_mut(sample).unwrap();
            writeln!(fh, "{}", out_hdr)?;
            writeln!(fh, "{}", out_seq)?;
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    const SMALL_FNA: &str = "\
>Sample1 Tag1-Tag2.Tag3-Tag4_1 5_4\n\
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\n\
>Sample1 Tag1-Tag2.Tag3-Tag4_2 3_2\n\
GCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCG\n\
>Sample2 Tag5-Tag6.Tag7-Tag8_1 10_8\n\
TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT\n\
>Sample2 Tag5-Tag6.Tag7-Tag8_2 1_0\n\
AAAA\n";

    #[test]
    fn test_sumaclust_output() {
        let mut out = Vec::new();
        process(Cursor::new(SMALL_FNA), &mut out, 0, None, false, None).unwrap();
        let s = String::from_utf8(out).unwrap();
        let lines: Vec<&str> = s.lines().collect();
        assert_eq!(lines[0], ">Sample1:1 count=9");
        assert_eq!(lines[1], "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG");
        assert_eq!(lines[2], ">Sample1:2 count=5");
        assert_eq!(lines[4], ">Sample2:3 count=18");
        assert_eq!(lines[6], ">Sample2:4 count=1");
        assert_eq!(lines[7], "AAAA");
        assert_eq!(lines.len(), 8);
    }

    #[test]
    fn test_usearch_output() {
        let mut out = Vec::new();
        process(Cursor::new(SMALL_FNA), &mut out, 0, None, true, None).unwrap();
        let s = String::from_utf8(out).unwrap();
        let lines: Vec<&str> = s.lines().collect();
        assert_eq!(lines[0], ">Sample1;size=9");
        assert_eq!(lines[1], "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG");
        assert_eq!(lines[6], ">Sample2;size=1");
        assert_eq!(lines[7], "AAAA");
        assert_eq!(lines.len(), 8);
    }

    #[test]
    fn test_usearch_padding() {
        let mut out = Vec::new();
        process(Cursor::new(SMALL_FNA), &mut out, 0, Some(65), true, None).unwrap();
        let s = String::from_utf8(out).unwrap();
        let lines: Vec<&str> = s.lines().collect();
        assert_eq!(lines[1].len(), 65);
        assert!(lines[1].ends_with("NNNNN"));
        assert_eq!(lines[7].len(), 65);
        assert!(lines[7].starts_with("AAAA"));
        assert!(lines[7].ends_with('N'));
    }

    #[test]
    fn test_no_padding_in_sumaclust_mode() {
        let mut out = Vec::new();
        process(Cursor::new(SMALL_FNA), &mut out, 0, Some(65), false, None).unwrap();
        let s = String::from_utf8(out).unwrap();
        let lines: Vec<&str> = s.lines().collect();
        assert_eq!(lines[1].len(), 60);
    }

    #[test]
    fn test_min_length_filter() {
        let mut out = Vec::new();
        process(Cursor::new(SMALL_FNA), &mut out, 10, None, false, None).unwrap();
        let s = String::from_utf8(out).unwrap();
        let lines: Vec<&str> = s.lines().collect();
        assert_eq!(lines.len(), 6);
        assert!(!s.contains("Sample2:4"));
        assert!(!s.contains("AAAA"));
    }

    #[test]
    fn test_max_length_filter() {
        let mut out = Vec::new();
        process(Cursor::new(SMALL_FNA), &mut out, 0, Some(10), false, None).unwrap();
        let s = String::from_utf8(out).unwrap();
        let lines: Vec<&str> = s.lines().collect();
        assert_eq!(lines.len(), 2);
        assert_eq!(lines[0], ">Sample2:1 count=1");
        assert_eq!(lines[1], "AAAA");
    }
}
```

- [ ] **Step 2: Run Rust tests**

```bash
cd rust && cargo test convert 2>&1 | tail -20
```

Expected: all 6 convert tests PASS.

- [ ] **Step 3: Commit**

```bash
git add rust/src/convert.rs rust/src/lib.rs
git commit -m "feat(rust): add convert subcommand with unit tests"
```

---

### Task C7: Register Rust convert subcommand

**Files:**
- Modify: `rust/src/main.rs`

- [ ] **Step 1: Add Convert to Commands and dispatch**

Replace `rust/src/main.rs` with:

```rust
use anyhow::Result;
use clap::{Parser, Subcommand};
use dame::{chimera_check, convert, decollapse, filter, rsi, sort};

#[derive(Parser)]
#[command(name = "dame", about = "DNA Metabarcoding toolkit")]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    Sort(sort::SortArgs),
    Chimera(chimera_check::ChimeraArgs),
    Filter(filter::FilterArgs),
    Decollapse(decollapse::DecollapseArgs),
    Rsi(rsi::RsiArgs),
    Convert(convert::ConvertArgs),
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    match cli.command {
        Commands::Sort(args) => sort::run(args),
        Commands::Chimera(args) => chimera_check::run(args),
        Commands::Filter(args) => filter::run(args),
        Commands::Decollapse(args) => decollapse::run(args),
        Commands::Rsi(args) => rsi::run(args),
        Commands::Convert(args) => convert::run(args),
    }
}
```

- [ ] **Step 2: Build and smoke test**

```bash
cd rust && cargo build --release 2>&1 | tail -5
./target/release/dame convert --help
```

Expected: builds cleanly; help shows `-i`, `--min-length`, `--max-length`, `-u`, `-s`.

- [ ] **Step 3: Run full Rust test suite**

```bash
cd rust && cargo test --all 2>&1 | tail -10
```

Expected: all tests PASS.

- [ ] **Step 4: Commit**

```bash
git add rust/src/main.rs
git commit -m "feat(rust): register convert subcommand in CLI"
```

---

### Task C8: Integration test + CI

**Files:**
- Create: `tests/integration/run_convert.sh`
- Modify: `.github/workflows/ci.yml`

- [ ] **Step 1: Write the integration test script**

```bash
#!/usr/bin/env bash
# tests/integration/run_convert.sh
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
FIXTURES="$REPO_ROOT/tests/fixtures"
DAME_BIN="$REPO_ROOT/rust/target/release/dame"
FNA="$FIXTURES/FilteredReads_small.fna"

if [ ! -f "$DAME_BIN" ]; then
    echo "SKIP: dame binary not found at $DAME_BIN (run: cd rust && cargo build --release)"
    exit 0
fi

TMPPY=$(mktemp -d)
TMPRS=$(mktemp -d)
trap "rm -rf '$TMPPY' '$TMPRS'" EXIT

echo "==> Sumaclust mode..."
(cd "$TMPPY" && dame-py convert -i "$FNA")
(cd "$TMPRS" && "$DAME_BIN" convert -i "$FNA")
diff "$TMPPY/FilteredReads.forsumaclust.fna" "$TMPRS/FilteredReads.forsumaclust.fna" \
    || { echo "FAIL: sumaclust output differs"; exit 1; }
echo "PASS: sumaclust"

echo "==> USEARCH mode..."
(cd "$TMPPY" && dame-py convert -i "$FNA" -u)
(cd "$TMPRS" && "$DAME_BIN" convert -i "$FNA" -u)
diff "$TMPPY/FilteredReads.forusearch.fna" "$TMPRS/FilteredReads.forusearch.fna" \
    || { echo "FAIL: usearch output differs"; exit 1; }
echo "PASS: usearch"

echo "==> USEARCH + --max-length 65..."
(cd "$TMPPY" && dame-py convert -i "$FNA" -u --max-length 65)
(cd "$TMPRS" && "$DAME_BIN" convert -i "$FNA" -u --max-length 65)
diff "$TMPPY/FilteredReads.forusearch.fna" "$TMPRS/FilteredReads.forusearch.fna" \
    || { echo "FAIL: usearch padded output differs"; exit 1; }
echo "PASS: usearch + max-length"

echo "==> --min-length 10 filter..."
(cd "$TMPPY" && dame-py convert -i "$FNA" --min-length 10)
(cd "$TMPRS" && "$DAME_BIN" convert -i "$FNA" --min-length 10)
diff "$TMPPY/FilteredReads.forsumaclust.fna" "$TMPRS/FilteredReads.forsumaclust.fna" \
    || { echo "FAIL: min-length output differs"; exit 1; }
echo "PASS: min-length filter"

echo "==> --sample-fastas..."
(cd "$TMPPY" && dame-py convert -i "$FNA" -s)
(cd "$TMPRS" && "$DAME_BIN" convert -i "$FNA" -s)
[ -d "$TMPPY/SampleFastas" ] || { echo "FAIL: dame-py SampleFastas not created"; exit 1; }
[ -d "$TMPRS/SampleFastas" ] || { echo "FAIL: dame SampleFastas not created"; exit 1; }
diff "$TMPPY/SampleFastas/Sample1.fixed.fasta" "$TMPRS/SampleFastas/Sample1.fixed.fasta" \
    || { echo "FAIL: Sample1.fixed.fasta differs"; exit 1; }
diff "$TMPPY/SampleFastas/Sample2.fixed.fasta" "$TMPRS/SampleFastas/Sample2.fixed.fasta" \
    || { echo "FAIL: Sample2.fixed.fasta differs"; exit 1; }
echo "PASS: sample-fastas"

echo "PASS: dame and dame-py convert produce identical output"
```

Make it executable:

```bash
chmod +x tests/integration/run_convert.sh
```

- [ ] **Step 2: Run locally**

```bash
bash tests/integration/run_convert.sh
```

Expected: 5 PASS lines, then `PASS: dame and dame-py convert produce identical output`.

- [ ] **Step 3: Add CI step**

In `.github/workflows/ci.yml`, add after the existing `run_filter.sh` step (before the chimera step):

```yaml
      - name: Run convert integration test
        run: bash tests/integration/run_convert.sh
```

- [ ] **Step 4: Commit**

```bash
git add tests/integration/run_convert.sh .github/workflows/ci.yml
git commit -m "test(integration): add convert cross-impl parity test + CI step"
```

---

### Task C9: Docs + version bump

**Files:**
- Modify: `README.md`
- Modify: `tutorial/README.md`
- Modify: `rust/Cargo.toml`
- Modify: `python/pyproject.toml`

- [ ] **Step 1: Bump versions**

In `rust/Cargo.toml`, change:
```toml
version = "0.6.0"
```
to:
```toml
version = "0.7.0"
```

In `python/pyproject.toml`, change:
```toml
version = "3.0.0"
```
to:
```toml
version = "3.1.0"
```

- [ ] **Step 2: Update README pipeline overview**

In `README.md`, find the pipeline overview block and add the `convert` step after `dame rsi`:

```
dame filter   --ps-info PSinfo.txt --x 2 --y 2 --t 2 --l 50
              → Comparisons_2PCRs.txt (all seqs, all replicates)
              → FilteredReads.fna   (passed all thresholds)

dame convert  -i FilteredReads.fna [-u] [--min-length N] [--max-length N] [-s]
              → FilteredReads.forsumaclust.fna  (sumaclust input, default)
              → FilteredReads.forusearch.fna    (USEARCH input, with -u)
              → SampleFastas/<Sample>.fixed.fasta  (per-sample, with -s)

dame rsi      Comparisons_2PCRs.txt
              → RSI_output.txt
```

(Move `dame convert` between `dame filter` and `dame rsi` in the pipeline block.)

- [ ] **Step 3: Update tutorial/README.md**

In `tutorial/README.md`, add a section after the `dame filter` step (look for "## 4. Filter" or similar heading). Add a new step:

```markdown
## 5. Convert for clustering

Convert the filtered reads for downstream clustering with USEARCH or sumaclust:

```bash
# For sumaclust (default):
dame-py convert -i FilteredReads.fna
# → FilteredReads.forsumaclust.fna

# For USEARCH (adds ;size= tag, optional N-padding to fixed length):
dame-py convert -i FilteredReads.fna -u
dame-py convert -i FilteredReads.fna -u --max-length 313

# Per-sample fastas (adds SampleFastas/ directory):
dame-py convert -i FilteredReads.fna -s
```

Use `--min-length` and `--max-length` to filter by amplicon length (inclusive bounds).

| Flag | Description |
|------|-------------|
| `-i` / `--in-fasta` | Input `FilteredReads.fna` (required) |
| `-u` / `--usearch` | USEARCH output (`>Sample;size=N`); default is sumaclust (`>Sample:N count=N`) |
| `--min-length N` | Drop sequences shorter than N |
| `--max-length N` | Drop sequences longer than N; pad to N in USEARCH mode |
| `-s` / `--sample-fastas` | Write per-sample fastas to `SampleFastas/` |

Original v1.0 flag spellings also accepted in `dame-py`: `--inFasta`, `-lmin`, `-lmax`, `--sampleFastas`.
```

- [ ] **Step 4: Rebuild Rust to confirm version change doesn't break**

```bash
cd rust && cargo build --release 2>&1 | tail -3
```

Expected: builds cleanly.

- [ ] **Step 5: Run full test suites**

```bash
python -m pytest python/tests/ -v 2>&1 | tail -10
cd rust && cargo test --all 2>&1 | tail -10
bash tests/integration/run_convert.sh
```

Expected: all pass.

- [ ] **Step 6: Commit**

```bash
git add rust/Cargo.toml python/pyproject.toml README.md tutorial/README.md
git commit -m "docs: add convert to pipeline overview + tutorial; bump v2.6 (py 3.1.0, crate 0.7.0)"
```

---

## Self-review against spec

**Spec coverage check:**

| Spec requirement | Task |
|-----------------|------|
| Input: parse 3-token header, size = sum of `_`-split counts | C3 (`_parse_fasta`) |
| Sumaclust output: `>Sample:N count=size` | C3, C6 |
| USEARCH output: `>Sample;size=size` | C3, C6 |
| N-padding in USEARCH mode only when `--max-length` set | C3 `ljust`, C6 loop |
| `SampleFastas/` with `makedirs(exist_ok=True)` | C3, C6 |
| `-i`/`--in-fasta` required flag | C3 (`register_subcommand`), C7 |
| `--min-length` inclusive lower bound | C3, C6 (`< min_length`) |
| `--max-length` inclusive upper bound | C3, C6 (`> max_len`) |
| Python backwards-compat aliases (`--inFasta`, `-lmin`, `-lmax`, `--sampleFastas`) | C3 (`register_subcommand`), tested in C2 |
| Rust canonical flags only | C5/C6 `ConvertArgs` |
| Parity test | C8 |
| CI step | C8 |
| Version bump Rust 0.6.0 → 0.7.0 | C9 |
| Version bump Python 3.0.0 → 3.1.0 | C9 |
| README pipeline block | C9 |
| Tutorial section | C9 |
| Fixture file | C1 |
| Commit `DAMe_1.0/bin/convertToUSearch.py` | Done in spec commit |

**Placeholder scan:** No TBD/TODO in plan. All code blocks are complete.

**Type consistency:** `_parse_fasta` → yields `(str, int, str)` tuples; `convert` → uses those tuples. `process` in Rust takes `R: BufRead, W: Write` + `Option<usize>` params; `run` passes `BufReader<File>` and `BufWriter<File>` — consistent throughout.
