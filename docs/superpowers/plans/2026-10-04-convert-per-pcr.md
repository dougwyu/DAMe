# `convert --per-pcr` Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a `--per-pcr` mode to `dame convert` / `dame-py convert` that writes one FASTA record per PCR replicate (labelled for mapping) plus a `PCRinfo.txt` table, so users can build per-PCR OTU tables without the zeros that `filter --y` introduces.

**Architecture:** A new per-PCR module in each implementation (`python/dame/perpcr.py`, `rust/src/perpcr.rs`) does the work; `convert` only parses flags and dispatches. A small read-only PSinfo helper is added next to `filter`'s own PSinfo code so both commands assign PCR numbers identically. Python and Rust must produce byte-identical output and identical error text.

**Tech Stack:** Python 3.11+ (argparse, pytest), Rust stable (clap 4 derive, anyhow, indexmap, tempfile for tests), bash integration scripts, optional vsearch 2.31 and R (dplyr, tidyr) for the end-to-end check.

**Spec:** `docs/superpowers/specs/2026-10-04-convert-per-pcr-design.md` (read it first; this plan implements it).

## Global Constraints

- Without `--per-pcr`, `convert` output must stay byte-identical to v3.1.0 (existing tests must pass unchanged).
- Python and Rust per-PCR outputs (`FilteredReads.perpcr.fna`, `PCRinfo.txt`) must be byte-identical.
- Error and warning text must be identical in both implementations. Errors print to stderr as `Error: <message>` and exit with status 1 (anyhow's default in Rust; `sys.exit("Error: ...")` in Python).
- Output file names: `FilteredReads.perpcr.fna` and `PCRinfo.txt`, written to `<name>.tmp` first and renamed only on success.
- Per-PCR record label: `>{sample}_PCR{k}.{n};size={c};sample={sample}_PCR{k};` with `n` a running counter from 1. No N-padding ever.
- `PCRinfo.txt` header with `--ps-info`: `pcr_id\tsample\tpcr\ttag_pair\tpool\treads_pre_mapping`; without: `pcr_id\tsample\tpcr\ttag_pair\treads_pre_mapping`.
- No em-dashes in any prose (docs, comments, commit messages). Use `--`.
- Version becomes 3.2.0 in `python/pyproject.toml`, `rust/Cargo.toml`, README title, and `rust/tests/cli_version_test.rs`.
- Commit messages end with `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`.

## File Structure

- Create `tests/fixtures/perpcr/`: synthetic data set (generator, `dame sort`-style pool files, `PSinfo.txt`, `filter` outputs) and `expected/` outputs. Shared by Python tests, Rust tests and the integration script.
- Modify `python/dame/modules_filter.py`: add `readPSinfoRows`.
- Create `python/dame/perpcr.py`: per-PCR conversion, messages, flag validation.
- Modify `python/dame/convert.py`: new flags, dispatch.
- Create `python/tests/test_perpcr.py`; modify `python/tests/test_filter.py`, `python/tests/test_convert.py`.
- Modify `rust/src/filter.rs`: add `ps_info_rows`.
- Create `rust/src/perpcr.rs`; modify `rust/src/lib.rs`, `rust/src/convert.rs`.
- Create `rust/tests/perpcr_test.rs`, `rust/tests/convert_perpcr_cli_test.rs`; modify `rust/tests/filter_test.rs`, `rust/tests/cli_version_test.rs`.
- Create `tutorial/perpcr_to_occupancy.R`: the R step of the recipe (used by the tutorial and the integration test).
- Create `tests/integration/run_perpcr.sh`; modify `tests/integration/run_malformed.sh`, `.github/workflows/ci.yml`.
- Modify `README.md`, `tutorial/README.md`, version files.

---

### Task 1: Synthetic per-PCR fixture

**Files:**
- Create: `tests/fixtures/perpcr/make_data.py`
- Create (generated): `tests/fixtures/perpcr/PSinfo.txt`, `tests/fixtures/perpcr/pool1/*.txt`, `tests/fixtures/perpcr/pool2/*.txt`, `tests/fixtures/perpcr/Comparisons_3PCRs.fasta`, `tests/fixtures/perpcr/FilteredReads.fna`
- Create: `tests/fixtures/perpcr/expected/FilteredReads.perpcr.fna`, `expected/PCRinfo.txt`, `expected/PCRinfo_no_psinfo.txt`, `expected/table.tsv`, `expected/survey.tsv`
- Create: `tests/fixtures/perpcr/README.md`

**Interfaces:**
- Produces: the fixture directory `tests/fixtures/perpcr/` with the files above. Later tasks refer to `Comparisons_3PCRs.fasta`, `PSinfo.txt` and `expected/*` by these exact names.

- [ ] **Step 1: Write the generator**

`tests/fixtures/perpcr/make_data.py` (keep the order of random calls exactly as written; the expected files depend on it):

```python
"""Generate the synthetic per-PCR fixture in `dame sort` output form.

Run from this directory: python3 make_data.py
Then: dame filter --ps-info PSinfo.txt --x 3 --y 2 --t 1 --l 100

Sequences (120-bp marker):
  A, B, C  real sequences
  A2       1-bp variant of A that also passes the filter (same OTU at 97%)
  E        2-substitution error copy of B (fails everywhere)
  Bt       B trimmed to 112 bp, inside B (counted by mapping with --mincols 110)
  Bs       90-bp fragment of B (rejected by --mincols 110)
S3 PCR2 gets no reads at all; S4 gets no reads in any PCR.
"""
import os
import random

random.seed(11)


def rs(n):
    return "".join(random.choice("ACGT") for _ in range(n))


def sub(s, i):
    return s[:i] + ("A" if s[i] != "A" else "C") + s[i + 1:]


A = rs(120)
A2 = sub(A, 60)
B = rs(120)
C = rs(120)
E = sub(sub(B, 30), 90)
Bt = B[4:116]
Bs = B[:90]
seqs = dict(A=A, A2=A2, B=B, C=C, E=E, Bt=Bt, Bs=Bs)

# PSinfo: sample, forward tag, reverse tag, pool. 4 samples x 3 PCRs.
ps = []
tag = 1
for s in ["S1", "S2", "S3", "S4"]:
    for k in range(3):
        ps.append((s, f"t{tag}", f"t{tag + 1}", 1 if k < 2 else 2))
        tag += 2
with open("PSinfo.txt", "w") as f:
    f.write("".join("\t".join(map(str, r)) + "\n" for r in ps))

# Reads per sample, per sequence, per PCR
counts = {
    "S1": {"A": (50, 40, 45), "A2": (5, 6, 0), "B": (0, 1, 0), "E": (2, 0, 0)},
    "S2": {"B": (30, 25, 0), "Bt": (0, 0, 3), "C": (0, 0, 1), "A": (0, 1, 0)},
    "S3": {"C": (8, 0, 6), "Bs": (0, 0, 2)},
    "S4": {},
}
for s, ft, rt, pool in ps:
    k = [r for r in ps if r[0] == s].index((s, ft, rt, pool))
    rows = [(n, c[k]) for n, c in counts[s].items() if c[k] > 0]
    if not rows:
        continue  # sort writes no file for a tag pair with no reads
    os.makedirs(f"pool{pool}", exist_ok=True)
    with open(f"pool{pool}/{ft}_{rt}.txt", "w") as f:
        for n, c in rows:
            f.write(f"COI\t{ft}\t{rt}\t{c}\t{seqs[n]}\n")
```

- [ ] **Step 2: Generate the data and run `filter`**

```bash
cd tests/fixtures/perpcr
python3 make_data.py
../../../rust/target/release/dame filter --ps-info PSinfo.txt --x 3 --y 2 --t 1 --l 100
rm PS1_files.txt PS2_files.txt PS3_files.txt Comparisons_3PCRs.txt Comparisons_2outOf3PCRs.txt \
   Comparisons_2outOf3PCRs.countsThreshold1.txt FilteredReads_atLeast2.fasta FilteredReads_atLeast2.threshold.fasta
python3 -c "import hashlib;[print(f, hashlib.md5(open(f,'rb').read()).hexdigest()) for f in ['Comparisons_3PCRs.fasta','FilteredReads.fna']]"
```

(Build the binary first with `cd rust && cargo build --release` if it is missing.)

Expected checksums:
```
Comparisons_3PCRs.fasta d55bf790ae9b23254415bbac01eaa1a3
FilteredReads.fna dc8510b4a94a44f721e5dbe7353037c1
```
If they differ, stop: the generator or `filter` behaves differently from when the expected files below were made.

Expected files present: `pool1/t1_t2.txt t3_t4.txt t7_t8.txt t9_t10.txt t13_t14.txt`, `pool2/t5_t6.txt t11_t12.txt t17_t18.txt`.

- [ ] **Step 3: Write the expected outputs**

`expected/FilteredReads.perpcr.fna` (exact content, 15 records):

```
>S1_PCR1.1;size=2;sample=S1_PCR1;
CCCTTATAAAAGCTGTTGCACCTAGCCAAGATCAACGGCAGCTGCAATGGAAATAGGCAATGACGGATATATATTAAAAAGTGTTTTAAGCTACATTGAGGCCCGTTCGTGCTCCTCGCC
>S1_PCR2.2;size=1;sample=S1_PCR2;
CCCTTATAAAAGCTGTTGCACCTAGCCAAGTTCAACGGCAGCTGCAATGGAAATAGGCAATGACGGATATATATTAAAAAGTGTTTTAAGATACATTGAGGCCCGTTCGTGCTCCTCGCC
>S1_PCR1.3;size=5;sample=S1_PCR1;
TTTCCTCATGCAATTCAAAACCATGTCCGTAATGTAGGCGAAATAGTAAACCATTTTACGAAGGATACCAAATTCCTCCTTATTCAGGACCTAACCTGAGGTAAACCAGGTCTCTCCGCC
>S1_PCR2.4;size=6;sample=S1_PCR2;
TTTCCTCATGCAATTCAAAACCATGTCCGTAATGTAGGCGAAATAGTAAACCATTTTACGAAGGATACCAAATTCCTCCTTATTCAGGACCTAACCTGAGGTAAACCAGGTCTCTCCGCC
>S1_PCR1.5;size=50;sample=S1_PCR1;
TTTCCTCATGCAATTCAAAACCATGTCCGTAATGTAGGCGAAATAGTAAACCATTTTACGGAGGATACCAAATTCCTCCTTATTCAGGACCTAACCTGAGGTAAACCAGGTCTCTCCGCC
>S1_PCR2.6;size=40;sample=S1_PCR2;
TTTCCTCATGCAATTCAAAACCATGTCCGTAATGTAGGCGAAATAGTAAACCATTTTACGGAGGATACCAAATTCCTCCTTATTCAGGACCTAACCTGAGGTAAACCAGGTCTCTCCGCC
>S1_PCR3.7;size=45;sample=S1_PCR3;
TTTCCTCATGCAATTCAAAACCATGTCCGTAATGTAGGCGAAATAGTAAACCATTTTACGGAGGATACCAAATTCCTCCTTATTCAGGACCTAACCTGAGGTAAACCAGGTCTCTCCGCC
>S2_PCR1.8;size=30;sample=S2_PCR1;
CCCTTATAAAAGCTGTTGCACCTAGCCAAGTTCAACGGCAGCTGCAATGGAAATAGGCAATGACGGATATATATTAAAAAGTGTTTTAAGATACATTGAGGCCCGTTCGTGCTCCTCGCC
>S2_PCR2.9;size=25;sample=S2_PCR2;
CCCTTATAAAAGCTGTTGCACCTAGCCAAGTTCAACGGCAGCTGCAATGGAAATAGGCAATGACGGATATATATTAAAAAGTGTTTTAAGATACATTGAGGCCCGTTCGTGCTCCTCGCC
>S2_PCR3.10;size=1;sample=S2_PCR3;
CTGAAGCATTGCTTTGTGAAGAGGGACTTCAGCCAATAGACCTGCATACCGGCTCATTCTTCATGTGCAACCTAGGGAGAATGTGTACATACGCTCTTACTGCGGTCGCGTCTAATAATA
>S2_PCR3.11;size=3;sample=S2_PCR3;
TATAAAAGCTGTTGCACCTAGCCAAGTTCAACGGCAGCTGCAATGGAAATAGGCAATGACGGATATATATTAAAAAGTGTTTTAAGATACATTGAGGCCCGTTCGTGCTCCT
>S2_PCR2.12;size=1;sample=S2_PCR2;
TTTCCTCATGCAATTCAAAACCATGTCCGTAATGTAGGCGAAATAGTAAACCATTTTACGGAGGATACCAAATTCCTCCTTATTCAGGACCTAACCTGAGGTAAACCAGGTCTCTCCGCC
>S3_PCR3.13;size=2;sample=S3_PCR3;
CCCTTATAAAAGCTGTTGCACCTAGCCAAGTTCAACGGCAGCTGCAATGGAAATAGGCAATGACGGATATATATTAAAAAGTGTTTTAAG
>S3_PCR1.14;size=8;sample=S3_PCR1;
CTGAAGCATTGCTTTGTGAAGAGGGACTTCAGCCAATAGACCTGCATACCGGCTCATTCTTCATGTGCAACCTAGGGAGAATGTGTACATACGCTCTTACTGCGGTCGCGTCTAATAATA
>S3_PCR3.15;size=6;sample=S3_PCR3;
CTGAAGCATTGCTTTGTGAAGAGGGACTTCAGCCAATAGACCTGCATACCGGCTCATTCTTCATGTGCAACCTAGGGAGAATGTGTACATACGCTCTTACTGCGGTCGCGTCTAATAATA
```

`expected/PCRinfo.txt` (tab-separated):

```
pcr_id	sample	pcr	tag_pair	pool	reads_pre_mapping
S1_PCR1	S1	1	t1-t2	1	57
S1_PCR2	S1	2	t3-t4	1	47
S1_PCR3	S1	3	t5-t6	2	45
S2_PCR1	S2	1	t7-t8	1	30
S2_PCR2	S2	2	t9-t10	1	26
S2_PCR3	S2	3	t11-t12	2	4
S3_PCR1	S3	1	t13-t14	1	8
S3_PCR2	S3	2	t15-t16	1	0
S3_PCR3	S3	3	t17-t18	2	8
S4_PCR1	S4	1	t19-t20	1	0
S4_PCR2	S4	2	t21-t22	1	0
S4_PCR3	S4	3	t23-t24	2	0
```

`expected/PCRinfo_no_psinfo.txt` (tab-separated):

```
pcr_id	sample	pcr	tag_pair	reads_pre_mapping
S1_PCR1	S1	1	t1-t2	57
S1_PCR2	S1	2	t3-t4	47
S1_PCR3	S1	3	t5-t6	45
S2_PCR1	S2	1	t7-t8	30
S2_PCR2	S2	2	t9-t10	26
S2_PCR3	S2	3	t11-t12	4
S3_PCR1	S3	1	t13-t14	8
S3_PCR2	S3	2	empty	0
S3_PCR3	S3	3	t17-t18	8
```

`expected/table.tsv` (vsearch 2.31 `--otutabout`, tab-separated):

```
#OTU ID	S1_PCR1	S1_PCR2	S1_PCR3	S2_PCR1	S2_PCR2	S2_PCR3	S3_PCR1	S3_PCR3
seq1	50	40	45	0	1	0	0	0
seq2	0	1	0	30	25	3	0	0
seq3	0	0	0	0	0	1	8	6
seq4	5	6	0	0	0	0	0	0
```

(seq1 = A, seq2 = B, seq3 = C, seq4 = A2: `--derep_fulllength --relabel seq` numbers by decreasing total abundance.)

`expected/survey.tsv` (tab-separated):

```
pcr_id	sample	pcr	tag_pair	pool	reads_pre_mapping	OTU_seq1	OTU_seq2	OTU_seq3
S1_PCR1	S1	1	t1-t2	1	57	55	0	0
S1_PCR2	S1	2	t3-t4	1	47	46	1	0
S1_PCR3	S1	3	t5-t6	2	45	45	0	0
S2_PCR1	S2	1	t7-t8	1	30	0	30	0
S2_PCR2	S2	2	t9-t10	1	26	1	25	0
S2_PCR3	S2	3	t11-t12	2	4	0	3	1
S3_PCR1	S3	1	t13-t14	1	8	0	0	8
S3_PCR2	S3	2	t15-t16	1	0	0	0	0
S3_PCR3	S3	3	t17-t18	2	8	0	0	6
S4_PCR1	S4	1	t19-t20	1	0	0	0	0
S4_PCR2	S4	2	t21-t22	1	0	0	0	0
S4_PCR3	S4	3	t23-t24	2	0	0	0	0
```

Make sure the expected files use real tab characters and end with a newline.

- [ ] **Step 4: Write the fixture README**

`tests/fixtures/perpcr/README.md`:

```markdown
# Per-PCR fixture

Synthetic data set for `convert --per-pcr`, in `dame sort` output form. See the docstring in
`make_data.py` for what each sequence and sample exercises.

Regenerate with `python3 make_data.py`, then
`dame filter --ps-info PSinfo.txt --x 3 --y 2 --t 1 --l 100`, keeping only
`Comparisons_3PCRs.fasta` and `FilteredReads.fna` from the filter output.

`expected/` holds the outputs the tests compare against: the per-PCR FASTA and `PCRinfo.txt`
(with and without `--ps-info`), the vsearch 2.31 sequence-level table, and the survey table
written by `tutorial/perpcr_to_occupancy.R`.
```

- [ ] **Step 5: Commit**

```bash
git add tests/fixtures/perpcr
git commit -m "test: add synthetic per-PCR fixture and expected outputs

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 2: Python PSinfo helper

**Files:**
- Modify: `python/dame/modules_filter.py` (add a function after `MakeSampleNameArray`)
- Test: `python/tests/test_filter.py`

**Interfaces:**
- Produces: `readPSinfoRows(PSinfo: str, X: int) -> list[tuple[str, int, str, str]]`, each tuple `(sample, pcr, tag_pair, pool)` with `pcr` 1-based and `tag_pair` = `"<col2>-<col3>"`, in PSinfo line order.

- [ ] **Step 1: Write the failing tests**

Add to `python/tests/test_filter.py` (and add `readPSinfoRows` to its import list from `dame.modules_filter`):

```python
def test_readPSinfoRows_assigns_pcr_numbers(tmp_path):
    psinfo = write_psinfo(tmp_path, [
        "S1\tt1\tt2\t1",
        "S1\tt3\tt4\t1",
        "S1\tt5\tt6\t2",
        "S2\tt7\tt8\t1",
        "S2\tt9\tt10\t1",
        "S2\tt11\tt12\t2",
    ])
    assert readPSinfoRows(psinfo, 3) == [
        ("S1", 1, "t1-t2", "1"), ("S1", 2, "t3-t4", "1"), ("S1", 3, "t5-t6", "2"),
        ("S2", 1, "t7-t8", "1"), ("S2", 2, "t9-t10", "1"), ("S2", 3, "t11-t12", "2"),
    ]


def test_readPSinfoRows_blank_line_consumes_a_slot_like_makePSnumFiles(tmp_path, monkeypatch):
    # A blank line still advances the line number, so the next row lands in the
    # slot after it, exactly as makePSnumFiles assigns replicate files.
    lines = ["S1\tt1\tt2\t1", "", "S1\tt3\tt4\t1"]
    psinfo = write_psinfo(tmp_path, lines)
    assert readPSinfoRows(psinfo, 2) == [("S1", 1, "t1-t2", "1"), ("S1", 1, "t3-t4", "1")]
    monkeypatch.chdir(tmp_path)
    makePSnumFiles(psinfo, X=2, P=1, chimeraChecked=False)
    assert open("PS1_files.txt").read() == "pool1/t1_t2.txt\npool1/t3_t4.txt\n"


def test_readPSinfoRows_skips_short_lines(tmp_path):
    psinfo = write_psinfo(tmp_path, ["S1\tt1\tt2\t1", "S1\tt3", "S1\tt5\tt6\t1"])
    assert readPSinfoRows(psinfo, 3) == [("S1", 1, "t1-t2", "1"), ("S1", 3, "t5-t6", "1")]
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `cd python && pytest tests/test_filter.py -k readPSinfoRows -v`
Expected: FAIL with `ImportError: cannot import name 'readPSinfoRows'`

- [ ] **Step 3: Implement**

Add to `python/dame/modules_filter.py` after `MakeSampleNameArray`:

```python
def readPSinfoRows(PSinfo, X):
    """Return (sample, pcr, tag_pair, pool) for each usable PSinfo line.

    PCR numbers (1..X) follow exactly the line-number rule makePSnumFiles uses
    to assign replicate files: blank and short lines are skipped but still
    advance the line number. Keeping one rule means convert --per-pcr and
    filter cannot disagree about which PCR is which.
    """
    rows = []
    with open(PSinfo) as f:
        for NR, line in enumerate(f, start=1):
            parts = line.split()
            if len(parts) < 4:
                continue
            residue = NR % X
            pcr = residue if residue != 0 else X
            rows.append((parts[0], pcr, "%s-%s" % (parts[1], parts[2]), parts[3]))
    return rows
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `cd python && pytest tests/test_filter.py -v`
Expected: all PASS

- [ ] **Step 5: Commit**

```bash
git add python/dame/modules_filter.py python/tests/test_filter.py
git commit -m "feat(python): add readPSinfoRows PSinfo helper

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 3: Python per-PCR conversion

**Files:**
- Create: `python/dame/perpcr.py`
- Test: `python/tests/test_perpcr.py`

**Interfaces:**
- Consumes: `readPSinfoRows` (Task 2); fixture files (Task 1).
- Produces (in `dame.perpcr`):
  - `FASTA_OUT = "FilteredReads.perpcr.fna"`, `INFO_OUT = "PCRinfo.txt"`
  - `FILTERED_READS_WARNING: str`, `USEARCH_NOTE: str` (exact text below)
  - `class PerPcrError(Exception)`: message is user-facing, without the `Error: ` prefix
  - `flag_error(per_pcr: bool, ps_info: str | None, sample_fastas: bool) -> str | None`
  - `per_pcr_convert(in_fasta, ps_info=None, min_length=0, max_length=None, usearch=False, out_dir=".") -> tuple[str, str]` returning the two output paths; raises `PerPcrError`.

- [ ] **Step 1: Write the failing tests**

`python/tests/test_perpcr.py`:

```python
import os
import shutil
from pathlib import Path

import pytest

from dame.perpcr import (
    FILTERED_READS_WARNING, USEARCH_NOTE, PerPcrError, flag_error, per_pcr_convert,
)

FIXTURE = Path(__file__).resolve().parents[2] / "tests" / "fixtures" / "perpcr"
COMP = str(FIXTURE / "Comparisons_3PCRs.fasta")
PSINFO = str(FIXTURE / "PSinfo.txt")
EXPECTED = FIXTURE / "expected"


def read(path):
    return Path(path).read_text()


def test_fixture_with_psinfo_matches_expected(tmp_path):
    fasta, info = per_pcr_convert(COMP, ps_info=PSINFO, out_dir=str(tmp_path))
    assert read(fasta) == read(EXPECTED / "FilteredReads.perpcr.fna")
    assert read(info) == read(EXPECTED / "PCRinfo.txt")
    assert sorted(os.listdir(tmp_path)) == ["FilteredReads.perpcr.fna", "PCRinfo.txt"]


def test_fixture_without_psinfo_matches_expected(tmp_path):
    fasta, info = per_pcr_convert(COMP, out_dir=str(tmp_path))
    assert read(fasta) == read(EXPECTED / "FilteredReads.perpcr.fna")
    assert read(info) == read(EXPECTED / "PCRinfo_no_psinfo.txt")


def test_one_read_in_one_pcr_becomes_one_record(tmp_path):
    fasta, _ = per_pcr_convert(COMP, out_dir=str(tmp_path))
    assert ">S1_PCR2.2;size=1;sample=S1_PCR2;\n" in read(fasta)


def test_length_filters_drop_records_without_padding(tmp_path):
    fasta, info = per_pcr_convert(COMP, ps_info=PSINFO, min_length=100, max_length=200,
                                  out_dir=str(tmp_path))
    lines = read(fasta).splitlines()
    assert len(lines) == 28                      # 14 records; the 90-bp fragment is gone
    assert all(len(s) in (112, 120) for s in lines[1::2])   # nothing padded to 200
    assert lines[-2] == ">S3_PCR3.14;size=6;sample=S3_PCR3;"
    assert "S3_PCR3\tS3\t3\tt17-t18\t2\t6\n" in read(info)


def test_warning_for_filteredreads_input(tmp_path, capsys):
    src = tmp_path / "FilteredReads_copy.fna"
    shutil.copy(COMP, src)
    out = tmp_path / "out"
    out.mkdir()
    per_pcr_convert(str(src), out_dir=str(out))
    assert FILTERED_READS_WARNING + "\n" in capsys.readouterr().err
    assert (out / "FilteredReads.perpcr.fna").exists()


def test_no_warning_for_comparisons_input(tmp_path, capsys):
    per_pcr_convert(COMP, out_dir=str(tmp_path))
    assert capsys.readouterr().err == ""


def test_usearch_note(tmp_path, capsys):
    fasta, _ = per_pcr_convert(COMP, usearch=True, max_length=200, out_dir=str(tmp_path))
    assert USEARCH_NOTE + "\n" in capsys.readouterr().err
    assert read(fasta) == read(EXPECTED / "FilteredReads.perpcr.fna")


def test_short_header_record_is_skipped(tmp_path):
    src = tmp_path / "in.fasta"
    src.write_text(">S1\nACGT\n>S1\tt1-t2.t3-t4_1\t1_0\nACGT\n")
    out = tmp_path / "out"
    out.mkdir()
    fasta, _ = per_pcr_convert(str(src), out_dir=str(out))
    assert read(fasta) == ">S1_PCR1.1;size=1;sample=S1_PCR1;\nACGT\n"


def test_flag_error():
    assert flag_error(True, None, True) == "--sample-fastas cannot be used with --per-pcr"
    assert flag_error(False, "PSinfo.txt", False) == "--ps-info requires --per-pcr"
    assert flag_error(True, "PSinfo.txt", False) is None
    assert flag_error(False, None, True) is None


ERROR_CASES = [
    ("counts_mismatch", ">S1\tt1-t2.t3-t4_1\t1_2_3\nACGT\n", None,
     "record 1 for sample S1 has 3 counts but 2 tag pairs"),
    ("non_integer", ">S1\tt1-t2.t3-t4_1\t1_x\nACGT\n", None,
     "record 1 for sample S1 has a non-integer count 'x'"),
    ("inconsistent_x", ">S1\tt1-t2.t3-t4_1\t1_0\nACGT\n>S2\tt5-t6.t7-t8.t9-t10_1\t1_0_0\nACGT\n", None,
     "record 2 for sample S2 has 3 PCRs but earlier records have 2"),
    ("inconsistent_tags", ">S1\tt1-t2.t3-t4_1\t1_1\nACGT\n>S1\tt1-t2.empty-empty_2\t1_0\nAACC\n", None,
     "sample S1 has tag pairs t1-t2.empty-empty in record 2 but t1-t2.t3-t4 in an earlier record"),
    ("sample_not_in_psinfo", ">S9\tt1-t2.t3-t4_1\t1_0\nACGT\n", "S1\tt1\tt2\t1\nS1\tt3\tt4\t1\n",
     "sample S9 is in the input but not in {psinfo}"),
    ("tag_mismatch", ">S1\tt1-t2.t3-t4_1\t1_0\nACGT\n", "S1\tt1\tt2\t1\nS1\tt3\tt9\t1\n",
     "sample S1 PCR 2: tag pair t3-t4 in the input but t3-t9 in PSinfo"),
    ("psinfo_rows", ">S1\tt1-t2.t3-t4_1\t1_0\nACGT\n", "S1\tt1\tt2\t1\n",
     "sample S1 does not have exactly one PSinfo row for each PCR 1..2"),
    ("empty_input", "", None, "no records in {input}"),
]


@pytest.mark.parametrize("name,fasta,psinfo,message", ERROR_CASES, ids=[c[0] for c in ERROR_CASES])
def test_errors_leave_no_output(tmp_path, name, fasta, psinfo, message):
    src = tmp_path / "in.fasta"
    src.write_text(fasta)
    ps = None
    if psinfo is not None:
        ps = tmp_path / "PSinfo.txt"
        ps.write_text(psinfo)
    out = tmp_path / "out"
    out.mkdir()
    with pytest.raises(PerPcrError) as exc:
        per_pcr_convert(str(src), ps_info=str(ps) if ps else None, out_dir=str(out))
    assert str(exc.value) == message.format(psinfo=ps, input=src)
    assert os.listdir(out) == []
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `cd python && pytest tests/test_perpcr.py -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'dame.perpcr'`

- [ ] **Step 3: Implement**

`python/dame/perpcr.py`:

```python
"""convert --per-pcr: one FASTA record per PCR replicate, plus PCRinfo.txt.

Input is Comparisons_<X>PCRs.fasta (or FilteredReads.fna) from dame filter.
Output reads keep their PCR as `sample=<sample>_PCR<k>` so a mapper such as
vsearch --usearch_global --otutabout builds a per-PCR table. The Rust
implementation (rust/src/perpcr.rs) must produce byte-identical output and
identical messages.
"""
import os
import sys

from dame.modules_filter import readPSinfoRows

FASTA_OUT = "FilteredReads.perpcr.fna"
INFO_OUT = "PCRinfo.txt"

FILTERED_READS_WARNING = (
    "Warning: --per-pcr input looks like FilteredReads output, which has already been\n"
    "filtered by --y/--t, so sequences that failed in a sample will appear as zeros.\n"
    "For per-PCR OTU tables, use Comparisons_<X>PCRs.fasta instead."
)

USEARCH_NOTE = (
    "Note: -u has no effect with --per-pcr; per-PCR output always uses the\n"
    ";size=N;sample=<pcr_id>; label and is never padded."
)


class PerPcrError(Exception):
    """An input problem that stops convert --per-pcr. The message is shown to the user."""


def flag_error(per_pcr, ps_info, sample_fastas):
    """Return the message for an invalid flag combination, or None."""
    if per_pcr and sample_fastas:
        return "--sample-fastas cannot be used with --per-pcr"
    if ps_info is not None and not per_pcr:
        return "--ps-info requires --per-pcr"
    return None


def _records(path):
    """Yield (ordinal, sample, tag_pairs, count_parts, seq) for each record.

    ordinal counts header lines from 1. A header with fewer than three tokens,
    or one not followed by a sequence line, is skipped, as in per-sample convert.
    """
    ordinal = 0
    pending = None
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                ordinal += 1
                pending = None
                toks = line.split()
                if len(toks) >= 3:
                    tags = toks[1].rsplit("_", 1)[0].split(".")
                    pending = (ordinal, toks[0][1:], tags, toks[2].split("_"))
            elif pending is not None:
                yield pending + (line.rstrip(),)
                pending = None


def _load_ps_info(ps_info, x):
    """Return {sample: [(tag_pair, pool) for PCR 1..x]} in PSinfo first-seen order."""
    by_sample = {}
    for sample, pcr, tag_pair, pool in readPSinfoRows(ps_info, x):
        rows = by_sample.setdefault(sample, {})
        if pcr in rows:
            raise PerPcrError(
                "sample %s does not have exactly one PSinfo row for each PCR 1..%d" % (sample, x))
        rows[pcr] = (tag_pair, pool)
    out = {}
    for sample, rows in by_sample.items():
        if sorted(rows) != list(range(1, x + 1)):
            raise PerPcrError(
                "sample %s does not have exactly one PSinfo row for each PCR 1..%d" % (sample, x))
        out[sample] = [rows[k] for k in range(1, x + 1)]
    return out


def _write(in_fasta, ps_info, min_length, max_length, fasta_tmp, info_tmp):
    x = None
    ps_rows = None
    samples = {}  # sample -> {"tags": [...], "reads": [...]}, in first-seen order
    counter = 0
    with open(fasta_tmp, "w") as out:
        for ordinal, sample, tags, parts, seq in _records(in_fasta):
            counts = []
            for part in parts:
                try:
                    counts.append(int(part))
                except ValueError:
                    raise PerPcrError("record %d for sample %s has a non-integer count '%s'"
                                      % (ordinal, sample, part))
            if len(counts) != len(tags):
                raise PerPcrError("record %d for sample %s has %d counts but %d tag pairs"
                                  % (ordinal, sample, len(counts), len(tags)))
            if x is None:
                x = len(tags)
                if ps_info is not None:
                    ps_rows = _load_ps_info(ps_info, x)
            elif len(tags) != x:
                raise PerPcrError("record %d for sample %s has %d PCRs but earlier records have %d"
                                  % (ordinal, sample, len(tags), x))

            state = samples.get(sample)
            if state is None:
                if ps_rows is not None:
                    if sample not in ps_rows:
                        raise PerPcrError("sample %s is in the input but not in %s"
                                          % (sample, ps_info))
                    for k, (tag_pair, _pool) in enumerate(ps_rows[sample], start=1):
                        if tags[k - 1] != "empty-empty" and tags[k - 1] != tag_pair:
                            raise PerPcrError(
                                "sample %s PCR %d: tag pair %s in the input but %s in PSinfo"
                                % (sample, k, tags[k - 1], tag_pair))
                state = samples[sample] = {"tags": tags, "reads": [0] * x}
            elif tags != state["tags"]:
                raise PerPcrError(
                    "sample %s has tag pairs %s in record %d but %s in an earlier record"
                    % (sample, ".".join(tags), ordinal, ".".join(state["tags"])))

            if len(seq) < min_length or (max_length is not None and len(seq) > max_length):
                continue
            for k, c in enumerate(counts, start=1):
                if c > 0:
                    counter += 1
                    out.write(">%s_PCR%d.%d;size=%d;sample=%s_PCR%d;\n%s\n"
                              % (sample, k, counter, c, sample, k, seq))
                    state["reads"][k - 1] += c

    if x is None:
        raise PerPcrError("no records in %s" % in_fasta)

    with open(info_tmp, "w") as f:
        if ps_rows is None:
            f.write("pcr_id\tsample\tpcr\ttag_pair\treads_pre_mapping\n")
            for sample, state in samples.items():
                for k in range(1, x + 1):
                    tag_pair = state["tags"][k - 1]
                    if tag_pair == "empty-empty":
                        tag_pair = "empty"
                    f.write("%s_PCR%d\t%s\t%d\t%s\t%d\n"
                            % (sample, k, sample, k, tag_pair, state["reads"][k - 1]))
        else:
            f.write("pcr_id\tsample\tpcr\ttag_pair\tpool\treads_pre_mapping\n")
            for sample, rows in ps_rows.items():
                reads = samples[sample]["reads"] if sample in samples else [0] * x
                for k, (tag_pair, pool) in enumerate(rows, start=1):
                    f.write("%s_PCR%d\t%s\t%d\t%s\t%s\t%d\n"
                            % (sample, k, sample, k, tag_pair, pool, reads[k - 1]))


def per_pcr_convert(in_fasta, ps_info=None, min_length=0, max_length=None,
                    usearch=False, out_dir="."):
    """Write FilteredReads.perpcr.fna and PCRinfo.txt into out_dir.

    Outputs are written under .tmp names and renamed only when the whole input
    has been read without error. Returns the two output paths.
    """
    if os.path.basename(in_fasta).startswith("FilteredReads"):
        print(FILTERED_READS_WARNING, file=sys.stderr)
    if usearch:
        print(USEARCH_NOTE, file=sys.stderr)

    fasta_out = os.path.join(out_dir, FASTA_OUT)
    info_out = os.path.join(out_dir, INFO_OUT)
    fasta_tmp = fasta_out + ".tmp"
    info_tmp = info_out + ".tmp"
    try:
        _write(in_fasta, ps_info, min_length, max_length, fasta_tmp, info_tmp)
    except BaseException:
        for path in (fasta_tmp, info_tmp):
            if os.path.exists(path):
                os.remove(path)
        raise
    os.replace(fasta_tmp, fasta_out)
    os.replace(info_tmp, info_out)
    return fasta_out, info_out
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `cd python && pytest tests/test_perpcr.py -v`
Expected: all PASS

- [ ] **Step 5: Commit**

```bash
git add python/dame/perpcr.py python/tests/test_perpcr.py
git commit -m "feat(python): add per-PCR conversion module

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 4: Python `convert` flags and dispatch

**Files:**
- Modify: `python/dame/convert.py` (`register_subcommand`, `run`)
- Test: `python/tests/test_convert.py`

**Interfaces:**
- Consumes: `flag_error`, `per_pcr_convert`, `PerPcrError` from `dame.perpcr` (Task 3).
- Produces: CLI flags `--per-pcr`/`-perPCR` (`args.per_pcr`, bool) and `--ps-info`/`-psInfo` (`args.ps_info`, str or None).

- [ ] **Step 1: Write the failing tests**

Add to `python/tests/test_convert.py`:

```python
from pathlib import Path

PERPCR_FIXTURE = Path(__file__).resolve().parents[2] / "tests" / "fixtures" / "perpcr"


def _parser():
    import argparse
    import dame.convert as conv
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers()
    conv.register_subcommand(sub)
    return parser


def test_convert_argparser_per_pcr_flags():
    args = _parser().parse_args(["convert", "-i", "x.fasta", "--per-pcr", "--ps-info", "PS.txt"])
    assert args.per_pcr is True
    assert args.ps_info == "PS.txt"
    legacy = _parser().parse_args(["convert", "-i", "x.fasta", "-perPCR", "-psInfo", "PS.txt"])
    assert legacy.per_pcr is True
    assert legacy.ps_info == "PS.txt"
    default = _parser().parse_args(["convert", "-i", "x.fasta"])
    assert default.per_pcr is False
    assert default.ps_info is None


def test_run_per_pcr_writes_outputs(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    args = _parser().parse_args([
        "convert", "-i", str(PERPCR_FIXTURE / "Comparisons_3PCRs.fasta"),
        "--per-pcr", "--ps-info", str(PERPCR_FIXTURE / "PSinfo.txt"),
    ])
    args.func(args)
    expected = PERPCR_FIXTURE / "expected"
    assert open("FilteredReads.perpcr.fna").read() == (expected / "FilteredReads.perpcr.fna").read_text()
    assert open("PCRinfo.txt").read() == (expected / "PCRinfo.txt").read_text()
    assert not os.path.exists("FilteredReads.forsumaclust.fna")


@pytest.mark.parametrize("argv,message", [
    (["--per-pcr", "-s"], "Error: --sample-fastas cannot be used with --per-pcr"),
    (["--ps-info", "PS.txt"], "Error: --ps-info requires --per-pcr"),
])
def test_run_rejects_flag_combinations(tmp_path, monkeypatch, argv, message):
    monkeypatch.chdir(tmp_path)
    args = _parser().parse_args(["convert", "-i", "x.fasta"] + argv)
    with pytest.raises(SystemExit) as exc:
        args.func(args)
    assert exc.value.code == message
    assert os.listdir(tmp_path) == []


def test_run_per_pcr_error_exits_with_message(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    src = tmp_path / "in.fasta"
    src.write_text(">S1\tt1-t2.t3-t4_1\t1_2_3\nACGT\n")
    args = _parser().parse_args(["convert", "-i", str(src), "--per-pcr"])
    with pytest.raises(SystemExit) as exc:
        args.func(args)
    assert exc.value.code == "Error: record 1 for sample S1 has 3 counts but 2 tag pairs"
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `cd python && pytest tests/test_convert.py -v`
Expected: the four new tests FAIL (`unrecognized arguments: --per-pcr` / `AttributeError: 'Namespace' object has no attribute 'per_pcr'`); existing tests PASS.

- [ ] **Step 3: Implement**

In `python/dame/convert.py`, add `import sys` at the top and `from dame.perpcr import PerPcrError, flag_error, per_pcr_convert` after the imports. In `register_subcommand`, before `p.set_defaults(func=run)`, add:

```python
    p.add_argument(
        "--per-pcr", "-perPCR",
        dest="per_pcr", action="store_true",
        help="Write one record per PCR replicate (for mapping) and PCRinfo.txt, "
             "instead of summing counts per sample. Use Comparisons_<X>PCRs.fasta as input",
    )
    p.add_argument(
        "--ps-info", "-psInfo",
        dest="ps_info", default=None, metavar="FILE",
        help="With --per-pcr: PSinfo file, adding pool and real tag pairs to PCRinfo.txt "
             "and cross-checking tag pairs",
    )
```

Replace `run` with:

```python
def run(args):
    message = flag_error(args.per_pcr, args.ps_info, args.sample_fastas)
    if message is not None:
        sys.exit("Error: " + message)
    if args.per_pcr:
        try:
            per_pcr_convert(
                in_fasta=args.in_fasta,
                ps_info=args.ps_info,
                min_length=args.min_length,
                max_length=args.max_length,
                usearch=args.usearch,
            )
        except PerPcrError as e:
            sys.exit("Error: %s" % e)
        return
    convert(
        in_fasta=args.in_fasta,
        min_length=args.min_length,
        max_length=args.max_length,
        usearch=args.usearch,
        sample_fastas=args.sample_fastas,
    )
```

Also update the subparser `description` to `"Convert FilteredReads.fna to USEARCH or sumaclust input format, or write per-PCR reads for mapping (--per-pcr)"`.

- [ ] **Step 4: Run all Python tests**

Run: `cd python && pytest -v`
Expected: all PASS (including every pre-existing convert test).

- [ ] **Step 5: Commit**

```bash
git add python/dame/convert.py python/tests/test_convert.py
git commit -m "feat(python): add --per-pcr and --ps-info to convert

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 5: Rust PSinfo helper

**Files:**
- Modify: `rust/src/filter.rs` (add a function after `make_sample_name_array`)
- Test: `rust/tests/filter_test.rs`

**Interfaces:**
- Produces: `pub fn ps_info_rows(ps_info: &str, x: usize) -> anyhow::Result<Vec<(String, usize, String, String)>>`, each `(sample, pcr, tag_pair, pool)`, `pcr` 1-based, `tag_pair` = `"<col2>-<col3>"`, in PSinfo line order. Same behaviour as Python `readPSinfoRows`.

- [ ] **Step 1: Write the failing tests**

Add `ps_info_rows` to the `use dame::filter::{...}` list in `rust/tests/filter_test.rs`, then add:

```rust
#[test]
fn test_ps_info_rows_assigns_pcr_numbers() {
    let dir = tempdir().unwrap();
    let psinfo = write_psinfo(
        dir.path(),
        &[
            "S1\tt1\tt2\t1",
            "S1\tt3\tt4\t1",
            "S1\tt5\tt6\t2",
            "S2\tt7\tt8\t1",
            "S2\tt9\tt10\t1",
            "S2\tt11\tt12\t2",
        ],
    );
    let rows = ps_info_rows(psinfo.to_str().unwrap(), 3).unwrap();
    let expected: Vec<(String, usize, String, String)> = vec![
        ("S1", 1, "t1-t2", "1"),
        ("S1", 2, "t3-t4", "1"),
        ("S1", 3, "t5-t6", "2"),
        ("S2", 1, "t7-t8", "1"),
        ("S2", 2, "t9-t10", "1"),
        ("S2", 3, "t11-t12", "2"),
    ]
    .into_iter()
    .map(|(s, k, t, p)| (s.to_string(), k, t.to_string(), p.to_string()))
    .collect();
    assert_eq!(rows, expected);
}

#[test]
fn test_ps_info_rows_blank_line_consumes_a_slot() {
    let dir = tempdir().unwrap();
    let psinfo = write_psinfo(dir.path(), &["S1\tt1\tt2\t1", "", "S1\tt3\tt4\t1"]);
    let rows = ps_info_rows(psinfo.to_str().unwrap(), 2).unwrap();
    assert_eq!(rows[0].1, 1);
    assert_eq!(rows[1].1, 1);
    assert_eq!(rows[1].2, "t3-t4");
}

#[test]
fn test_ps_info_rows_skips_short_lines() {
    let dir = tempdir().unwrap();
    let psinfo = write_psinfo(dir.path(), &["S1\tt1\tt2\t1", "S1\tt3", "S1\tt5\tt6\t1"]);
    let rows = ps_info_rows(psinfo.to_str().unwrap(), 3).unwrap();
    assert_eq!(rows.len(), 2);
    assert_eq!((rows[1].1, rows[1].2.as_str()), (3, "t5-t6"));
}
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `cd rust && cargo test --test filter_test ps_info_rows`
Expected: compile error, `cannot find function ps_info_rows`.

- [ ] **Step 3: Implement**

Add to `rust/src/filter.rs` after `make_sample_name_array`:

```rust
/// Returns (sample, pcr, tag_pair, pool) for each usable PSinfo line.
///
/// PCR numbers (1..=x) follow exactly the line-number rule `make_ps_num_files`
/// uses to assign replicate files: blank and short lines are skipped but still
/// advance the line number. Keeping one rule means `convert --per-pcr` and
/// `filter` cannot disagree about which PCR is which.
pub fn ps_info_rows(ps_info: &str, x: usize) -> Result<Vec<(String, usize, String, String)>> {
    let reader =
        BufReader::new(File::open(ps_info).with_context(|| format!("opening {}", ps_info))?);
    let mut rows = Vec::new();
    for (nr, line) in reader.lines().enumerate() {
        let line = line?;
        let parts: Vec<&str> = line.split_whitespace().collect();
        if parts.len() < 4 {
            continue;
        }
        let residue = (nr + 1) % x;
        let pcr = if residue != 0 { residue } else { x };
        rows.push((
            parts[0].to_string(),
            pcr,
            format!("{}-{}", parts[1], parts[2]),
            parts[3].to_string(),
        ));
    }
    Ok(rows)
}
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `cd rust && cargo test --test filter_test`
Expected: all PASS

- [ ] **Step 5: Commit**

```bash
git add rust/src/filter.rs rust/tests/filter_test.rs
git commit -m "feat(rust): add ps_info_rows PSinfo helper

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 6: Rust per-PCR conversion

**Files:**
- Create: `rust/src/perpcr.rs`
- Modify: `rust/src/lib.rs` (add `pub mod perpcr;`)
- Test: `rust/tests/perpcr_test.rs`

**Interfaces:**
- Consumes: `dame::filter::ps_info_rows` (Task 5); fixture files (Task 1).
- Produces (in `dame::perpcr`):
  - `pub const FASTA_OUT: &str`, `pub const INFO_OUT: &str`, `pub const FILTERED_READS_WARNING: &str`, `pub const USEARCH_NOTE: &str` (same text as Python)
  - `pub fn flag_error(per_pcr: bool, ps_info: Option<&str>, sample_fastas: bool) -> Option<&'static str>`
  - `pub fn run_per_pcr(in_fasta: &str, ps_info: Option<&str>, min_length: usize, max_length: Option<usize>, usearch: bool, out_dir: &Path) -> anyhow::Result<()>`; errors carry the same message text as Python's `PerPcrError`.

- [ ] **Step 1: Write the failing tests**

`rust/tests/perpcr_test.rs`:

```rust
use dame::perpcr::{flag_error, run_per_pcr};
use std::fs;
use std::path::{Path, PathBuf};
use tempfile::tempdir;

fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../tests/fixtures/perpcr").join(name)
}

fn read(p: &Path) -> String {
    fs::read_to_string(p).unwrap()
}

fn comp() -> String {
    fixture("Comparisons_3PCRs.fasta").to_str().unwrap().to_string()
}

#[test]
fn fixture_with_psinfo_matches_expected() {
    let out = tempdir().unwrap();
    let ps = fixture("PSinfo.txt");
    run_per_pcr(&comp(), Some(ps.to_str().unwrap()), 0, None, false, out.path()).unwrap();
    assert_eq!(
        read(&out.path().join("FilteredReads.perpcr.fna")),
        read(&fixture("expected/FilteredReads.perpcr.fna"))
    );
    assert_eq!(read(&out.path().join("PCRinfo.txt")), read(&fixture("expected/PCRinfo.txt")));
    let mut names: Vec<String> = fs::read_dir(out.path())
        .unwrap()
        .map(|e| e.unwrap().file_name().into_string().unwrap())
        .collect();
    names.sort();
    assert_eq!(names, vec!["FilteredReads.perpcr.fna", "PCRinfo.txt"]);
}

#[test]
fn fixture_without_psinfo_matches_expected() {
    let out = tempdir().unwrap();
    run_per_pcr(&comp(), None, 0, None, false, out.path()).unwrap();
    assert_eq!(
        read(&out.path().join("FilteredReads.perpcr.fna")),
        read(&fixture("expected/FilteredReads.perpcr.fna"))
    );
    assert_eq!(
        read(&out.path().join("PCRinfo.txt")),
        read(&fixture("expected/PCRinfo_no_psinfo.txt"))
    );
}

#[test]
fn length_filters_drop_records_without_padding() {
    let out = tempdir().unwrap();
    let ps = fixture("PSinfo.txt");
    run_per_pcr(&comp(), Some(ps.to_str().unwrap()), 100, Some(200), false, out.path()).unwrap();
    let fasta = read(&out.path().join("FilteredReads.perpcr.fna"));
    let lines: Vec<&str> = fasta.lines().collect();
    assert_eq!(lines.len(), 28);
    assert!(lines.iter().skip(1).step_by(2).all(|s| s.len() == 112 || s.len() == 120));
    assert_eq!(lines[26], ">S3_PCR3.14;size=6;sample=S3_PCR3;");
    assert!(read(&out.path().join("PCRinfo.txt")).contains("S3_PCR3\tS3\t3\tt17-t18\t2\t6\n"));
}

#[test]
fn short_header_record_is_skipped() {
    let dir = tempdir().unwrap();
    let src = dir.path().join("in.fasta");
    fs::write(&src, ">S1\nACGT\n>S1\tt1-t2.t3-t4_1\t1_0\nACGT\n").unwrap();
    let out = dir.path().join("out");
    fs::create_dir(&out).unwrap();
    run_per_pcr(src.to_str().unwrap(), None, 0, None, false, &out).unwrap();
    assert_eq!(
        read(&out.join("FilteredReads.perpcr.fna")),
        ">S1_PCR1.1;size=1;sample=S1_PCR1;\nACGT\n"
    );
}

#[test]
fn flag_errors() {
    assert_eq!(
        flag_error(true, None, true),
        Some("--sample-fastas cannot be used with --per-pcr")
    );
    assert_eq!(flag_error(false, Some("PS.txt"), false), Some("--ps-info requires --per-pcr"));
    assert_eq!(flag_error(true, Some("PS.txt"), false), None);
    assert_eq!(flag_error(false, None, true), None);
}

/// (name, fasta, psinfo, message); {psinfo} and {input} are replaced with paths.
const ERROR_CASES: &[(&str, &str, Option<&str>, &str)] = &[
    ("counts_mismatch", ">S1\tt1-t2.t3-t4_1\t1_2_3\nACGT\n", None,
     "record 1 for sample S1 has 3 counts but 2 tag pairs"),
    ("non_integer", ">S1\tt1-t2.t3-t4_1\t1_x\nACGT\n", None,
     "record 1 for sample S1 has a non-integer count 'x'"),
    ("inconsistent_x", ">S1\tt1-t2.t3-t4_1\t1_0\nACGT\n>S2\tt5-t6.t7-t8.t9-t10_1\t1_0_0\nACGT\n", None,
     "record 2 for sample S2 has 3 PCRs but earlier records have 2"),
    ("inconsistent_tags", ">S1\tt1-t2.t3-t4_1\t1_1\nACGT\n>S1\tt1-t2.empty-empty_2\t1_0\nAACC\n", None,
     "sample S1 has tag pairs t1-t2.empty-empty in record 2 but t1-t2.t3-t4 in an earlier record"),
    ("sample_not_in_psinfo", ">S9\tt1-t2.t3-t4_1\t1_0\nACGT\n", Some("S1\tt1\tt2\t1\nS1\tt3\tt4\t1\n"),
     "sample S9 is in the input but not in {psinfo}"),
    ("tag_mismatch", ">S1\tt1-t2.t3-t4_1\t1_0\nACGT\n", Some("S1\tt1\tt2\t1\nS1\tt3\tt9\t1\n"),
     "sample S1 PCR 2: tag pair t3-t4 in the input but t3-t9 in PSinfo"),
    ("psinfo_rows", ">S1\tt1-t2.t3-t4_1\t1_0\nACGT\n", Some("S1\tt1\tt2\t1\n"),
     "sample S1 does not have exactly one PSinfo row for each PCR 1..2"),
    ("empty_input", "", None, "no records in {input}"),
];

#[test]
fn errors_leave_no_output() {
    for (name, fasta, psinfo, message) in ERROR_CASES {
        let dir = tempdir().unwrap();
        let src = dir.path().join("in.fasta");
        fs::write(&src, fasta).unwrap();
        let ps = psinfo.map(|text| {
            let p = dir.path().join("PSinfo.txt");
            fs::write(&p, text).unwrap();
            p
        });
        let out = dir.path().join("out");
        fs::create_dir(&out).unwrap();
        let err = run_per_pcr(
            src.to_str().unwrap(),
            ps.as_ref().map(|p| p.to_str().unwrap()),
            0,
            None,
            false,
            &out,
        )
        .expect_err(name);
        let expected = message
            .replace("{psinfo}", ps.as_ref().map(|p| p.to_str().unwrap()).unwrap_or(""))
            .replace("{input}", src.to_str().unwrap());
        assert_eq!(err.to_string(), expected, "case {}", name);
        assert_eq!(fs::read_dir(&out).unwrap().count(), 0, "case {} left output", name);
    }
}
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `cd rust && cargo test --test perpcr_test`
Expected: compile error, `could not find perpcr in dame`.

- [ ] **Step 3: Implement**

Add `pub mod perpcr;` to `rust/src/lib.rs` (keep alphabetical order, after `pub mod filter;`).

`rust/src/perpcr.rs`:

```rust
//! `convert --per-pcr`: one FASTA record per PCR replicate, plus PCRinfo.txt.
//!
//! Input is Comparisons_<X>PCRs.fasta (or FilteredReads.fna) from `dame filter`.
//! Output reads keep their PCR as `sample=<sample>_PCR<k>` so a mapper such as
//! `vsearch --usearch_global --otutabout` builds a per-PCR table. The Python
//! implementation (python/dame/perpcr.py) must produce byte-identical output
//! and identical messages.

use anyhow::{bail, Context, Result};
use indexmap::IndexMap;
use std::fs::{self, File};
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

use crate::filter::ps_info_rows;

pub const FASTA_OUT: &str = "FilteredReads.perpcr.fna";
pub const INFO_OUT: &str = "PCRinfo.txt";

pub const FILTERED_READS_WARNING: &str = "Warning: --per-pcr input looks like FilteredReads output, which has already been\nfiltered by --y/--t, so sequences that failed in a sample will appear as zeros.\nFor per-PCR OTU tables, use Comparisons_<X>PCRs.fasta instead.";

pub const USEARCH_NOTE: &str = "Note: -u has no effect with --per-pcr; per-PCR output always uses the\n;size=N;sample=<pcr_id>; label and is never padded.";

/// Returns the message for an invalid flag combination, or None.
pub fn flag_error(per_pcr: bool, ps_info: Option<&str>, sample_fastas: bool) -> Option<&'static str> {
    if per_pcr && sample_fastas {
        return Some("--sample-fastas cannot be used with --per-pcr");
    }
    if ps_info.is_some() && !per_pcr {
        return Some("--ps-info requires --per-pcr");
    }
    None
}

struct SampleState {
    tags: Vec<String>,
    reads: Vec<u64>,
}

/// PSinfo rows per sample in first-seen order; element k-1 is PCR k's (tag pair, pool).
type PsRows = IndexMap<String, Vec<(String, String)>>;

fn load_ps_info(ps_info: &str, x: usize) -> Result<PsRows> {
    let mut by_sample: IndexMap<String, IndexMap<usize, (String, String)>> = IndexMap::new();
    for (sample, pcr, tag_pair, pool) in ps_info_rows(ps_info, x)? {
        let rows = by_sample.entry(sample.clone()).or_default();
        if rows.contains_key(&pcr) {
            bail!("sample {} does not have exactly one PSinfo row for each PCR 1..{}", sample, x);
        }
        rows.insert(pcr, (tag_pair, pool));
    }
    let mut out = PsRows::new();
    for (sample, mut rows) in by_sample {
        if rows.len() != x || !(1..=x).all(|k| rows.contains_key(&k)) {
            bail!("sample {} does not have exactly one PSinfo row for each PCR 1..{}", sample, x);
        }
        let ordered = (1..=x).map(|k| rows.swap_remove(&k).unwrap()).collect();
        out.insert(sample, ordered);
    }
    Ok(out)
}

fn write_outputs(
    in_fasta: &str,
    ps_info: Option<&str>,
    min_length: usize,
    max_length: Option<usize>,
    fasta_tmp: &Path,
    info_tmp: &Path,
) -> Result<()> {
    let reader =
        BufReader::new(File::open(in_fasta).with_context(|| format!("opening {}", in_fasta))?);
    let mut out = BufWriter::new(
        File::create(fasta_tmp).with_context(|| format!("creating {}", fasta_tmp.display()))?,
    );

    let mut x: Option<usize> = None;
    let mut ps_rows: Option<PsRows> = None;
    let mut samples: IndexMap<String, SampleState> = IndexMap::new();
    let mut counter: u64 = 0;
    let mut ordinal: usize = 0;
    // (ordinal, sample, tag pairs, count parts) of a header awaiting its sequence line
    let mut pending: Option<(usize, String, Vec<String>, Vec<String>)> = None;

    for line in reader.lines() {
        let line = line?;
        if line.starts_with('>') {
            ordinal += 1;
            pending = None;
            let toks: Vec<&str> = line.split_whitespace().collect();
            if toks.len() >= 3 {
                let tag_field = toks[1].rsplit_once('_').map_or(toks[1], |(head, _)| head);
                let tags = tag_field.split('.').map(String::from).collect();
                let parts = toks[2].split('_').map(String::from).collect();
                pending = Some((ordinal, toks[0][1..].to_string(), tags, parts));
            }
            continue;
        }
        let Some((ord, sample, tags, parts)) = pending.take() else {
            continue;
        };
        let seq = line.trim_end();

        let mut counts: Vec<i64> = Vec::with_capacity(parts.len());
        for part in &parts {
            match part.parse::<i64>() {
                Ok(c) => counts.push(c),
                Err(_) => bail!(
                    "record {} for sample {} has a non-integer count '{}'",
                    ord, sample, part
                ),
            }
        }
        if counts.len() != tags.len() {
            bail!(
                "record {} for sample {} has {} counts but {} tag pairs",
                ord, sample, counts.len(), tags.len()
            );
        }
        let xv = match x {
            None => {
                x = Some(tags.len());
                if let Some(p) = ps_info {
                    ps_rows = Some(load_ps_info(p, tags.len())?);
                }
                tags.len()
            }
            Some(xv) => {
                if tags.len() != xv {
                    bail!(
                        "record {} for sample {} has {} PCRs but earlier records have {}",
                        ord, sample, tags.len(), xv
                    );
                }
                xv
            }
        };

        if !samples.contains_key(&sample) {
            {
                if let Some(rows) = &ps_rows {
                    let Some(ps) = rows.get(&sample) else {
                        bail!(
                            "sample {} is in the input but not in {}",
                            sample,
                            ps_info.unwrap_or("")
                        );
                    };
                    for (k, (tag_pair, _pool)) in ps.iter().enumerate() {
                        if tags[k] != "empty-empty" && &tags[k] != tag_pair {
                            bail!(
                                "sample {} PCR {}: tag pair {} in the input but {} in PSinfo",
                                sample, k + 1, tags[k], tag_pair
                            );
                        }
                    }
                }
                samples.insert(sample.clone(), SampleState { tags: tags.clone(), reads: vec![0; xv] });
            }
        } else if samples[&sample].tags != tags {
            bail!(
                "sample {} has tag pairs {} in record {} but {} in an earlier record",
                sample, tags.join("."), ord, samples[&sample].tags.join(".")
            );
        }

        if seq.len() < min_length || max_length.map_or(false, |m| seq.len() > m) {
            continue;
        }
        let state = samples.get_mut(&sample).unwrap();
        for (k, &c) in counts.iter().enumerate() {
            if c > 0 {
                counter += 1;
                writeln!(
                    out,
                    ">{s}_PCR{k}.{n};size={c};sample={s}_PCR{k};",
                    s = sample, k = k + 1, n = counter, c = c
                )?;
                writeln!(out, "{}", seq)?;
                state.reads[k] += c as u64;
            }
        }
    }
    out.flush()?;

    let Some(x) = x else {
        bail!("no records in {}", in_fasta);
    };

    let mut info = BufWriter::new(
        File::create(info_tmp).with_context(|| format!("creating {}", info_tmp.display()))?,
    );
    match &ps_rows {
        None => {
            writeln!(info, "pcr_id\tsample\tpcr\ttag_pair\treads_pre_mapping")?;
            for (sample, state) in &samples {
                for k in 0..x {
                    let tag_pair = if state.tags[k] == "empty-empty" { "empty" } else { state.tags[k].as_str() };
                    writeln!(info, "{s}_PCR{k}\t{s}\t{k}\t{t}\t{r}", s = sample, k = k + 1, t = tag_pair, r = state.reads[k])?;
                }
            }
        }
        Some(rows) => {
            writeln!(info, "pcr_id\tsample\tpcr\ttag_pair\tpool\treads_pre_mapping")?;
            for (sample, ps) in rows {
                for (k, (tag_pair, pool)) in ps.iter().enumerate() {
                    let r = samples.get(sample).map_or(0, |s| s.reads[k]);
                    writeln!(info, "{s}_PCR{k}\t{s}\t{k}\t{t}\t{p}\t{r}", s = sample, k = k + 1, t = tag_pair, p = pool, r = r)?;
                }
            }
        }
    }
    info.flush()?;
    Ok(())
}

/// Writes FilteredReads.perpcr.fna and PCRinfo.txt into `out_dir`.
///
/// Outputs are written under `.tmp` names and renamed only when the whole
/// input has been read without error.
pub fn run_per_pcr(
    in_fasta: &str,
    ps_info: Option<&str>,
    min_length: usize,
    max_length: Option<usize>,
    usearch: bool,
    out_dir: &Path,
) -> Result<()> {
    let looks_filtered = Path::new(in_fasta)
        .file_name()
        .map_or(false, |n| n.to_string_lossy().starts_with("FilteredReads"));
    if looks_filtered {
        eprintln!("{}", FILTERED_READS_WARNING);
    }
    if usearch {
        eprintln!("{}", USEARCH_NOTE);
    }

    let fasta_out = out_dir.join(FASTA_OUT);
    let info_out = out_dir.join(INFO_OUT);
    let fasta_tmp = out_dir.join(format!("{}.tmp", FASTA_OUT));
    let info_tmp = out_dir.join(format!("{}.tmp", INFO_OUT));
    if let Err(e) = write_outputs(in_fasta, ps_info, min_length, max_length, &fasta_tmp, &info_tmp) {
        let _ = fs::remove_file(&fasta_tmp);
        let _ = fs::remove_file(&info_tmp);
        return Err(e);
    }
    fs::rename(&fasta_tmp, &fasta_out)?;
    fs::rename(&info_tmp, &info_out)?;
    Ok(())
}
```

Note: Python's `int()` and Rust's `i64::parse` both accept a leading `+` and leading zeros, so they agree on every count `filter` can write.

- [ ] **Step 4: Run the tests to verify they pass**

Run: `cd rust && cargo test --test perpcr_test`
Expected: all PASS. Then `cargo test --all` (all PASS) and `cargo build --release`.

- [ ] **Step 5: Commit**

```bash
git add rust/src/perpcr.rs rust/src/lib.rs rust/tests/perpcr_test.rs
git commit -m "feat(rust): add per-PCR conversion module

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 7: Rust `convert` flags and dispatch

**Files:**
- Modify: `rust/src/convert.rs` (`ConvertArgs`, `run`)
- Test: `rust/tests/convert_perpcr_cli_test.rs`

**Interfaces:**
- Consumes: `dame::perpcr::{flag_error, run_per_pcr, FILTERED_READS_WARNING, USEARCH_NOTE}` (Task 6).
- Produces: CLI flags `--per-pcr` and `--ps-info <FILE>` on `dame convert`.

- [ ] **Step 1: Write the failing tests**

`rust/tests/convert_perpcr_cli_test.rs`:

```rust
use dame::perpcr::{FILTERED_READS_WARNING, USEARCH_NOTE};
use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;
use tempfile::tempdir;

fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../tests/fixtures/perpcr").join(name)
}

fn dame(dir: &Path, args: &[&str]) -> std::process::Output {
    Command::new(env!("CARGO_BIN_EXE_dame"))
        .current_dir(dir)
        .arg("convert")
        .args(args)
        .output()
        .expect("run dame convert")
}

#[test]
fn per_pcr_writes_expected_outputs() {
    let dir = tempdir().unwrap();
    let comp = fixture("Comparisons_3PCRs.fasta");
    let ps = fixture("PSinfo.txt");
    let out = dame(dir.path(), &["-i", comp.to_str().unwrap(), "--per-pcr", "--ps-info", ps.to_str().unwrap()]);
    assert!(out.status.success(), "{}", String::from_utf8_lossy(&out.stderr));
    assert!(out.stderr.is_empty());
    assert_eq!(
        fs::read_to_string(dir.path().join("FilteredReads.perpcr.fna")).unwrap(),
        fs::read_to_string(fixture("expected/FilteredReads.perpcr.fna")).unwrap()
    );
    assert_eq!(
        fs::read_to_string(dir.path().join("PCRinfo.txt")).unwrap(),
        fs::read_to_string(fixture("expected/PCRinfo.txt")).unwrap()
    );
}

#[test]
fn per_pcr_warns_for_filteredreads_and_notes_u() {
    let dir = tempdir().unwrap();
    let src = dir.path().join("FilteredReads.fna");
    fs::copy(fixture("FilteredReads.fna"), &src).unwrap();
    let out = dame(dir.path(), &["-i", "FilteredReads.fna", "--per-pcr", "-u"]);
    assert!(out.status.success());
    let stderr = String::from_utf8(out.stderr).unwrap();
    assert_eq!(stderr, format!("{}\n{}\n", FILTERED_READS_WARNING, USEARCH_NOTE));
}

#[test]
fn rejects_sample_fastas_with_per_pcr() {
    let dir = tempdir().unwrap();
    let out = dame(dir.path(), &["-i", "x.fasta", "--per-pcr", "-s"]);
    assert_eq!(out.status.code(), Some(1));
    assert_eq!(
        String::from_utf8(out.stderr).unwrap(),
        "Error: --sample-fastas cannot be used with --per-pcr\n"
    );
    assert_eq!(fs::read_dir(dir.path()).unwrap().count(), 0);
}

#[test]
fn rejects_ps_info_without_per_pcr() {
    let dir = tempdir().unwrap();
    let out = dame(dir.path(), &["-i", "x.fasta", "--ps-info", "PS.txt"]);
    assert_eq!(out.status.code(), Some(1));
    assert_eq!(String::from_utf8(out.stderr).unwrap(), "Error: --ps-info requires --per-pcr\n");
}

#[test]
fn per_pcr_error_message_and_exit_code() {
    let dir = tempdir().unwrap();
    fs::write(dir.path().join("in.fasta"), ">S1\tt1-t2.t3-t4_1\t1_2_3\nACGT\n").unwrap();
    let out = dame(dir.path(), &["-i", "in.fasta", "--per-pcr"]);
    assert_eq!(out.status.code(), Some(1));
    assert_eq!(
        String::from_utf8(out.stderr).unwrap(),
        "Error: record 1 for sample S1 has 3 counts but 2 tag pairs\n"
    );
}
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `cd rust && cargo test --test convert_perpcr_cli_test`
Expected: FAIL (clap rejects `--per-pcr` as an unexpected argument, exit code 2).

- [ ] **Step 3: Implement**

In `rust/src/convert.rs`, add `use crate::perpcr;` and `use std::path::Path;` to the imports, and `bail` to the anyhow import (`use anyhow::{bail, Result};`). Add to `ConvertArgs` after `sample_fastas`:

```rust
    /// Write one record per PCR replicate (for mapping) and PCRinfo.txt, instead of
    /// summing counts per sample. Use Comparisons_<X>PCRs.fasta as input
    #[arg(long = "per-pcr")]
    pub per_pcr: bool,
    /// With --per-pcr: PSinfo file, adding pool and real tag pairs to PCRinfo.txt
    /// and cross-checking tag pairs
    #[arg(long = "ps-info")]
    pub ps_info: Option<String>,
```

At the start of `run`, before opening the input:

```rust
    if let Some(message) = perpcr::flag_error(args.per_pcr, args.ps_info.as_deref(), args.sample_fastas) {
        bail!("{}", message);
    }
    if args.per_pcr {
        return perpcr::run_per_pcr(
            &args.in_fasta,
            args.ps_info.as_deref(),
            args.min_length,
            args.max_length,
            args.usearch,
            Path::new("."),
        );
    }
```

Update the `in_fasta` doc comment to `/// Input FilteredReads.fna file (Comparisons_<X>PCRs.fasta with --per-pcr)`.

- [ ] **Step 4: Run the tests to verify they pass**

Run: `cd rust && cargo test --all`
Expected: all PASS (including the existing convert unit tests). Then `cargo build --release`.

- [ ] **Step 5: Commit**

```bash
git add rust/src/convert.rs rust/tests/convert_perpcr_cli_test.rs
git commit -m "feat(rust): add --per-pcr and --ps-info to convert

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 8: Integration tests, R step, CI

**Files:**
- Create: `tutorial/perpcr_to_occupancy.R`
- Create: `tests/integration/run_perpcr.sh`
- Create: `tests/fixtures/malformed/Comparisons_perpcr_counts_mismatch.fasta`, `Comparisons_perpcr_non_integer.fasta`, `PSinfo_perpcr_tag_mismatch.txt`
- Modify: `tests/integration/run_malformed.sh` (append a per-PCR section before the final summary line, if any)
- Modify: `.github/workflows/ci.yml`

**Interfaces:**
- Consumes: both CLIs (Tasks 4 and 7), fixture (Task 1).
- Produces: `tutorial/perpcr_to_occupancy.R` with usage `Rscript perpcr_to_occupancy.R table.tsv clusters.uc PCRinfo.txt survey.tsv` (the tutorial in Task 9 shows this script).

- [ ] **Step 1: Add the R script**

`tutorial/perpcr_to_occupancy.R`:

```r
#!/usr/bin/env Rscript
# Turn a per-PCR, sequence-level vsearch table into an occupancy-model survey table.
#
# Usage: Rscript perpcr_to_occupancy.R table.tsv clusters.uc PCRinfo.txt survey.tsv
#
#   table.tsv    vsearch --otutabout output: one row per DAMe-passed sequence,
#                one column per PCR (from convert --per-pcr labels)
#   clusters.uc  vsearch --cluster_* --uc output for the same sequences
#   PCRinfo.txt  from dame convert --per-pcr
#   survey.tsv   written: one row per PCR (PCRinfo columns), one column per OTU
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 4)

# vsearch drops ";size=N;" style annotations from table row names but keeps them
# in .uc files, so compare labels without annotations
strip_annotations <- function(label) sub(";.*$", "", label)

seq_table <- read.delim(args[1], check.names = FALSE) |>
  rename(sequence = 1) |>
  mutate(sequence = strip_annotations(sequence))

# .uc: S rows are centroids (V9 = label), H rows are members (V9 = label, V10 = centroid)
uc <- read.delim(args[2], header = FALSE, stringsAsFactors = FALSE)
membership <- bind_rows(
  uc |> filter(V1 == "S") |> transmute(sequence = V9, centroid = V9),
  uc |> filter(V1 == "H") |> transmute(sequence = V9, centroid = V10)
) |>
  transmute(sequence = strip_annotations(sequence),
            otu = paste0("OTU_", strip_annotations(centroid)))

unassigned <- setdiff(seq_table$sequence, membership$sequence)
if (length(unassigned) > 0) {
  stop("sequences in ", args[1], " missing from ", args[2], ": ",
       paste(head(unassigned, 5), collapse = ", "))
}

pcr_info <- read.delim(args[3], stringsAsFactors = FALSE)

# Sum sequence rows within OTUs
otu_table <- seq_table |>
  inner_join(membership, by = "sequence") |>
  group_by(otu) |>
  summarise(across(-sequence, sum), .groups = "drop")

# One row per PCR, one column per OTU
otu_by_pcr <- otu_table |>
  pivot_longer(-otu, names_to = "pcr_id", values_to = "count") |>
  pivot_wider(names_from = otu, values_from = count)
otu_cols <- setdiff(names(otu_by_pcr), "pcr_id")

# Keep every PCR in PCRinfo; PCRs with no column in table.tsv had no matched reads
survey <- pcr_info |>
  left_join(otu_by_pcr, by = "pcr_id") |>
  mutate(across(all_of(otu_cols), \(x) replace_na(x, 0L)))

write.table(survey, args[4], sep = "\t", quote = FALSE, row.names = FALSE)
```

- [ ] **Step 2: Write `run_perpcr.sh`**

`tests/integration/run_perpcr.sh` (make it executable with `chmod +x`):

```bash
#!/usr/bin/env bash
#
# convert --per-pcr end to end on tests/fixtures/perpcr:
#   filter reproduces the committed Comparisons/FilteredReads files;
#   dame and dame-py write the expected per-PCR FASTA and PCRinfo.txt;
#   optionally, the tutorial recipe (vsearch, then R) gives the expected tables.
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
FIX="$REPO_ROOT/tests/fixtures/perpcr"
EXP="$FIX/expected"
DAME_BIN="$REPO_ROOT/rust/target/release/dame"

if [ ! -f "$DAME_BIN" ]; then
    echo "SKIP: dame binary not found at $DAME_BIN (run: cd rust && cargo build --release)"
    exit 0
fi

WORK=$(mktemp -d)
trap "rm -rf '$WORK'" EXIT
fail() { echo "FAIL: $*"; exit 1; }

echo "==> filter reproduces the committed fixture inputs..."
mkdir "$WORK/filter"
cp -R "$FIX/pool1" "$FIX/pool2" "$FIX/PSinfo.txt" "$WORK/filter/"
(cd "$WORK/filter" && "$DAME_BIN" filter --ps-info PSinfo.txt --x 3 --y 2 --t 1 --l 100 >/dev/null)
diff "$WORK/filter/Comparisons_3PCRs.fasta" "$FIX/Comparisons_3PCRs.fasta" || fail "Comparisons_3PCRs.fasta differs"
diff "$WORK/filter/FilteredReads.fna" "$FIX/FilteredReads.fna" || fail "FilteredReads.fna differs"
echo "PASS: filter"

# run_impl <impl> <outdir> <convert args...>
run_impl() {
    local impl="$1" dir="$2"; shift 2
    mkdir -p "$dir"
    if [ "$impl" = "py" ]; then
        (cd "$dir" && dame-py convert "$@")
    else
        (cd "$dir" && "$DAME_BIN" convert "$@")
    fi
}

echo "==> convert --per-pcr --ps-info..."
for impl in py rs; do
    run_impl "$impl" "$WORK/ps_$impl" -i "$FIX/Comparisons_3PCRs.fasta" --per-pcr --ps-info "$FIX/PSinfo.txt"
    diff "$WORK/ps_$impl/FilteredReads.perpcr.fna" "$EXP/FilteredReads.perpcr.fna" || fail "$impl FASTA differs from expected"
    diff "$WORK/ps_$impl/PCRinfo.txt" "$EXP/PCRinfo.txt" || fail "$impl PCRinfo.txt differs from expected"
done
echo "PASS: per-PCR with PSinfo (dame and dame-py)"

echo "==> convert --per-pcr without PSinfo..."
for impl in py rs; do
    run_impl "$impl" "$WORK/nops_$impl" -i "$FIX/Comparisons_3PCRs.fasta" --per-pcr
    diff "$WORK/nops_$impl/PCRinfo.txt" "$EXP/PCRinfo_no_psinfo.txt" || fail "$impl PCRinfo.txt (no PSinfo) differs"
done
echo "PASS: per-PCR without PSinfo"

echo "==> convert --per-pcr with length filters (parity)..."
for impl in py rs; do
    run_impl "$impl" "$WORK/len_$impl" -i "$FIX/Comparisons_3PCRs.fasta" --per-pcr --ps-info "$FIX/PSinfo.txt" --min-length 100 --max-length 200
done
diff "$WORK/len_py/FilteredReads.perpcr.fna" "$WORK/len_rs/FilteredReads.perpcr.fna" || fail "length-filtered FASTA differs"
diff "$WORK/len_py/PCRinfo.txt" "$WORK/len_rs/PCRinfo.txt" || fail "length-filtered PCRinfo.txt differs"
echo "PASS: length filters"

echo "==> warning text parity..."
for impl in py rs; do
    mkdir -p "$WORK/warn_$impl"
    cp "$FIX/FilteredReads.fna" "$WORK/warn_$impl/"
    run_impl "$impl" "$WORK/warn_$impl" -i FilteredReads.fna --per-pcr -u 2>"$WORK/warn_$impl.err"
done
diff "$WORK/warn_py.err" "$WORK/warn_rs.err" || fail "warning/note text differs"
grep -q "^Warning: --per-pcr input looks like FilteredReads output" "$WORK/warn_rs.err" || fail "warning missing"
echo "PASS: warning and note"

if ! command -v vsearch >/dev/null 2>&1; then
    echo "SKIP: vsearch not on PATH; recipe check skipped"
    echo "PASS: per-PCR integration"
    exit 0
fi

echo "==> recipe: dereplicate, map, cluster..."
R="$WORK/recipe"
mkdir "$R"
(cd "$R" && "$DAME_BIN" convert -i "$FIX/FilteredReads.fna" -u >/dev/null)
(cd "$R" && vsearch --derep_fulllength FilteredReads.forusearch.fna --sizein --sizeout --relabel seq --output passed.fna --quiet)
(cd "$R" && vsearch --usearch_global "$WORK/ps_rs/FilteredReads.perpcr.fna" --db passed.fna \
    --id 1.0 --mincols 110 --query_cov 1.0 --otutabout table.tsv --quiet)
diff "$R/table.tsv" "$EXP/table.tsv" || fail "vsearch table differs from expected"
(cd "$R" && vsearch --cluster_size passed.fna --sizein --id 0.97 --uc clusters.uc --quiet)
echo "PASS: vsearch recipe"

if ! command -v Rscript >/dev/null 2>&1 || ! Rscript -e 'library(dplyr); library(tidyr)' >/dev/null 2>&1; then
    echo "SKIP: Rscript with dplyr and tidyr not available; R step skipped"
    echo "PASS: per-PCR integration"
    exit 0
fi

echo "==> recipe: R join..."
Rscript "$REPO_ROOT/tutorial/perpcr_to_occupancy.R" "$R/table.tsv" "$R/clusters.uc" \
    "$WORK/ps_rs/PCRinfo.txt" "$R/survey.tsv" 2>/dev/null
diff "$R/survey.tsv" "$EXP/survey.tsv" || fail "survey table differs from expected"
echo "PASS: R join"

echo "PASS: per-PCR integration"
```

- [ ] **Step 3: Run it**

Run: `bash tests/integration/run_perpcr.sh`
Expected: every `PASS:` line, ending `PASS: per-PCR integration` (vsearch and R steps run on a machine that has them; otherwise `SKIP:` lines).

- [ ] **Step 4: Add malformed-input cases**

Create the fixtures:

`tests/fixtures/malformed/Comparisons_perpcr_counts_mismatch.fasta`:
```
>S1	t1-t2.t3-t4_1	1_2_3
ACGT
```

`tests/fixtures/malformed/Comparisons_perpcr_non_integer.fasta`:
```
>S1	t1-t2.t3-t4_1	1_x
ACGT
```

`tests/fixtures/malformed/PSinfo_perpcr_tag_mismatch.txt` (pair with the fixture's `Comparisons_3PCRs.fasta`; the second line's reverse tag is wrong):
```
S1	t1	t2	1
S1	t3	t9	1
S1	t5	t6	2
S2	t7	t8	1
S2	t9	t10	1
S2	t11	t12	2
S3	t13	t14	1
S3	t15	t16	1
S3	t17	t18	2
S4	t19	t20	1
S4	t21	t22	1
S4	t23	t24	2
```

Append to `tests/integration/run_malformed.sh`, before its final line (keep the existing final summary line last):

```bash
# check_perpcr_error <name> <expected message> <convert args...>
# Both implementations must fail with exit 1, the same message, and no output files.
check_perpcr_error() {
    local name="$1" message="$2"; shift 2
    local d
    for impl in py rs; do
        d="$WORK/perpcr_${name}_$impl"
        mkdir -p "$d"
        if [ "$impl" = "py" ]; then
            (cd "$d" && $TIMEOUT dame-py convert "$@" 2>stderr.txt) && fail "$name: dame-py succeeded"
        else
            (cd "$d" && $TIMEOUT "$DAME_BIN" convert "$@" 2>stderr.txt) && fail "$name: dame succeeded"
        fi
        [ "$(cat "$d/stderr.txt")" = "Error: $message" ] || fail "$name ($impl): got '$(cat "$d/stderr.txt")'"
        [ ! -e "$d/FilteredReads.perpcr.fna" ] && [ ! -e "$d/PCRinfo.txt" ] || fail "$name ($impl): left output"
    done
    echo "PASS: convert --per-pcr, $name"
}

echo "==> convert --per-pcr with damaged inputs..."
check_perpcr_error counts_mismatch "record 1 for sample S1 has 3 counts but 2 tag pairs" \
    -i "$MALFORMED/Comparisons_perpcr_counts_mismatch.fasta" --per-pcr
check_perpcr_error non_integer "record 1 for sample S1 has a non-integer count 'x'" \
    -i "$MALFORMED/Comparisons_perpcr_non_integer.fasta" --per-pcr
check_perpcr_error tag_mismatch "sample S1 PCR 2: tag pair t3-t4 in the input but t3-t9 in PSinfo" \
    -i "$FIXTURES/perpcr/Comparisons_3PCRs.fasta" --per-pcr --ps-info "$MALFORMED/PSinfo_perpcr_tag_mismatch.txt"
```

Also add one line to the header comment list at the top of `run_malformed.sh`:
```
#   per-PCR convert errors                      ->  new in v3.2.0; both must refuse identically
```

Run: `bash tests/integration/run_malformed.sh`
Expected: all existing `PASS:` lines plus three `PASS: convert --per-pcr, ...` lines.

- [ ] **Step 5: Add the CI step**

In `.github/workflows/ci.yml`, after the `Run convert integration test` step, add:

```yaml
      - name: Run per-PCR convert integration test (vsearch/R steps skip if absent)
        run: bash tests/integration/run_perpcr.sh
```

- [ ] **Step 6: Commit**

```bash
git add tutorial/perpcr_to_occupancy.R tests/integration/run_perpcr.sh tests/integration/run_malformed.sh \
    tests/fixtures/malformed/Comparisons_perpcr_counts_mismatch.fasta \
    tests/fixtures/malformed/Comparisons_perpcr_non_integer.fasta \
    tests/fixtures/malformed/PSinfo_perpcr_tag_mismatch.txt .github/workflows/ci.yml
git commit -m "test: per-PCR integration, malformed cases, and the tutorial R step

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 9: Documentation and version 3.2.0

**Files:**
- Modify: `README.md` (title line 1; Pipeline overview around line 157; changelog after entry 15 around line 343; Testing around line 417; Repository layout around lines 436 and 449)
- Modify: `tutorial/README.md` (new section after the convert section, around line 330)
- Modify: `python/pyproject.toml`, `rust/Cargo.toml`, `rust/Cargo.lock` (via cargo), `rust/tests/cli_version_test.rs`

**Interfaces:**
- Consumes: everything above. No code interfaces.

- [ ] **Step 1: Bump the version, test first**

In `rust/tests/cli_version_test.rs`, rename the test to `cli_reports_version_3_2_0` and change the expected string to `"dame 3.2.0\n"`.

Run: `cd rust && cargo test --test cli_version_test`
Expected: FAIL (`dame 3.1.0`).

Set `version = "3.2.0"` in `rust/Cargo.toml` and `python/pyproject.toml`. Run `cd rust && cargo build --release` (updates `Cargo.lock`) and `pip install -e python/` (or `uv tool install --editable ./python --force`).

Run: `cd rust && cargo test --test cli_version_test && dame-py --version`
Expected: PASS, and `dame-py 3.2.0`.

- [ ] **Step 2: Update README.md**

1. Line 1: `# DAMe v3.2.0: DNA Metabarcoding toolkit`.
2. In the opening paragraph, change `converts filtered reads to USEARCH or sumaclust input (**convert**)` to `converts filtered reads to USEARCH or sumaclust input, or writes per-PCR reads for building occupancy-model tables (**convert**)`.
3. In "Pipeline overview", after the existing `dame convert` lines, add:
```
dame convert  -i Comparisons_2PCRs.fasta --per-pcr [--ps-info PSinfo.txt] [--min-length N] [--max-length N]
              → FilteredReads.perpcr.fna  (one record per PCR replicate, for mapping)
              → PCRinfo.txt               (one row per PCR)
```
4. After changelog entry 15, add:
```
16. **DAMe v3.2.0 -- Per-PCR output for occupancy models.**  `convert --per-pcr`
    writes one record per PCR replicate, labelled `;size=N;sample=<sample>_PCR<k>;`,
    plus `PCRinfo.txt` (one row per PCR, with `--ps-info` adding pools and
    cross-checking tag pairs). Run on the unfiltered `Comparisons_<X>PCRs.fasta`
    and mapped onto the DAMe-passed sequences, it gives a per-PCR table without
    the zeros that `filter --y` writes into samples where a real sequence was
    seen in too few PCRs. See the tutorial section "Per-PCR OTU tables for
    occupancy and detection models".
```
5. In "Testing", after `bash tests/integration/run_convert.sh`, add `bash tests/integration/run_perpcr.sh`.
6. In "Repository layout", add `perpcr.py                    convert --per-pcr: per-PCR reads and PCRinfo.txt` under the Python list and `perpcr.rs                    convert --per-pcr: per-PCR reads and PCRinfo.txt` under the Rust list, aligned with their neighbours.

Use `--` rather than an em-dash in entry 16 even though older entries use `—`.

- [ ] **Step 3: Add the tutorial section**

Insert into `tutorial/README.md` after the convert section (before the rsi section):

````markdown
## Per-PCR OTU tables for occupancy and detection models

Occupancy and detection models such as occPlus and occJSDM need one row per PCR reaction, and
they estimate false positives and false negatives themselves. The usual DAMe route gets in their
way twice:

1. **Filtered cells become zeros.** `filter --y 2` keeps a sequence in a sample only if it is in
   at least 2 of that sample's PCRs. A real sequence seen once, in one PCR, is dropped from that
   sample even though it passes elsewhere, so its OTU column survives but the cell is 0.
2. **PCRs are summed.** `convert` adds the replicate counts into one record per sample.

The fix is to keep DAMe's decision about *which* sequences are real, and rebuild the *counts* by
mapping every read, PCR by PCR, onto those sequences. Ji et al. (2025) built their occPlus table
this way.

### Example

`tests/fixtures/perpcr/` holds a small synthetic data set (4 samples x 3 PCRs, a 120-bp marker).
After `filter --x 3 --y 2 --t 1 --l 100`, the usual route gives:

```
sample  OTU_A  OTU_B  OTU_C
S1        146      0      0     B was seen in S1 PCR2 (1 read)
S2          0     55      0     C and A were each seen once in S2
S3          0      0     14
```

The per-PCR route below gives one row per PCR, keeps those detections, and lists every PCR,
including S3 PCR2 (no reads) and sample S4 (no reads in any PCR):

```
pcr_id   sample pcr tag_pair pool reads_pre_mapping OTU_seq1 OTU_seq2 OTU_seq3
S1_PCR1  S1     1   t1-t2    1    57                55       0        0
S1_PCR2  S1     2   t3-t4    1    47                46       1        0
S1_PCR3  S1     3   t5-t6    2    45                45       0        0
S2_PCR1  S2     1   t7-t8    1    30                0        30       0
S2_PCR2  S2     2   t9-t10   1    26                1        25       0
S2_PCR3  S2     3   t11-t12  2    4                 0        3        1
S3_PCR1  S3     1   t13-t14  1    8                 0        0        8
S3_PCR2  S3     2   t15-t16  1    0                 0        0        0
S3_PCR3  S3     3   t17-t18  2    8                 0        0        6
S4_PCR1  S4     1   t19-t20  1    0                 0        0        0
S4_PCR2  S4     2   t21-t22  1    0                 0        0        0
S4_PCR3  S4     3   t23-t24  2    0                 0        0        0
```

(OTU_seq1 is A plus its 1-bp variant A2, OTU_seq2 is B, OTU_seq3 is C.)

### Recipe

```bash
# 1. Per-PCR reads from the UNFILTERED comparisons file, plus PCRinfo.txt
dame convert -i Comparisons_3PCRs.fasta --per-pcr --ps-info PSinfo.txt

# 2. Reference: the DAMe-passed sequences, dereplicated
#    (no --max-length here: with -u it pads with N, and padded sequences cannot match exactly;
#     use vsearch --minseqlength/--maxseqlength for a length range)
dame convert -i FilteredReads.fna -u
vsearch --derep_fulllength FilteredReads.forusearch.fna --sizein --sizeout --relabel seq --output passed.fna

# 3. Map: exact, contained matches of at least N bases (N a little below the amplicon length)
vsearch --usearch_global FilteredReads.perpcr.fna --db passed.fna \
    --id 1.0 --mincols 110 --query_cov 1.0 --otutabout table.tsv

# 4. Cluster the passed sequences into OTUs
vsearch --cluster_size passed.fna --sizein --id 0.97 --uc clusters.uc

# 5. Sum sequences within OTUs, one row per PCR, all PCRs listed
Rscript perpcr_to_occupancy.R table.tsv clusters.uc PCRinfo.txt survey.tsv
```

`table.tsv` has one row per passed sequence, so every step can be inspected and other
clusterings (or LULU-style curation) can be tried without remapping.
`perpcr_to_occupancy.R` is in this directory. Its output's first columns (`pcr_id` to
`reads_pre_mapping`) are the per-PCR covariates; the `OTU_*` columns are the read counts.

### Why these mapping options

Checked with vsearch 2.31:

- `--id 1.0` on its own is not an exact-match test: vsearch ignores end gaps when computing
  identity, so any fragment of a reference scores 100% (a 30-bp fragment was counted as the full
  sequence). `--mincols N` sets the shortest match accepted.
- `--query_cov 1.0` makes every base of a read align, so a read must lie inside its reference
  and cannot be longer than it. Without it, a read with extra bases at one end is accepted.
  `--maxqt 1.0` is not a substitute: a read shifted off one end of its reference passes it.
- `--minseqlength` does not filter reads in a search; it removes reference sequences. To drop
  short reads, use `--min-length` on `convert --per-pcr`.
- For markers whose length varies a lot, `--target_cov F` (fraction of the reference aligned)
  is the length-relative alternative to `--mincols`.
- A read can match two passed sequences equally if they differ only beyond its ends; vsearch
  then picks one. This matters only when the two are in different OTUs.

With these options a read counts only if it matches a DAMe-passed sequence exactly. Sequences
that failed in every sample (errors, chimeras) match nothing. Mapping at a lower identity
(`--id 0.97 --mincols N --query_cov 1.0`) lets error copies add their reads to their parent
sequence, at some risk of also absorbing chimeras or reads of rare relatives.

### Notes

- **Use `Comparisons_<X>PCRs.fasta`, not `FilteredReads.fna`.** `FilteredReads.fna` has already
  been filtered by `--y`/`--t`, so converting it per PCR brings the zeros back. `convert` warns
  when the input file name starts with `FilteredReads`.
- **`reads_pre_mapping`** counts every read written for that PCR, including reads that later
  match nothing. It measures sequencing depth; it is not the row total of the mapped table.
- **PCRs with no reads** have no column in `table.tsv`; `PCRinfo.txt` lists them with
  `reads_pre_mapping = 0` and the R script adds them as all-zero rows. DAMe cannot tell a PCR
  that was run and gave no reads from one that was never run; that is for you to decide.
- **Without `--ps-info`**, samples with no reads in any PCR do not appear at all, empty PCRs have
  `tag_pair = empty`, and there is no `pool` column.
- Tag names containing `.`, `-` or `_` break the header format, as for `convert` and `rsi`.
````

Then copy the fixture-based example numbers exactly from `tests/fixtures/perpcr/expected/survey.tsv` if any differ from the table above (they should not).

- [ ] **Step 4: Check prose rules and run everything**

Run:
```bash
awk '/^## Per-PCR OTU tables/{on=1;next} /^## /{on=0} on' tutorial/README.md | grep -c "—"   # expect 0
awk '/^16\. \*\*DAMe v3.2.0/{on=1} /^17\./{on=0} on' README.md | grep -c "—"                  # expect 0
cd python && pytest -q && cd ../rust && cargo test --all -q && cd .. && \
  for t in run_convert run_perpcr run_malformed run_filter run_pipeline; do bash tests/integration/$t.sh; done
```
Expected: no em-dash in the new tutorial section or changelog entry 16; all tests PASS.

- [ ] **Step 5: Commit**

```bash
git add README.md tutorial/README.md python/pyproject.toml rust/Cargo.toml rust/Cargo.lock rust/tests/cli_version_test.rs
git commit -m "docs: document convert --per-pcr; release 3.2.0

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

## After the plan

- Open a PR from `convert-per-pcr` to `master` on `dougwyu/DAMe` only when Doug asks.
- occJSDM follow-up (separate repo, separate PR, later): a one-line pointer from Lesson 0's input-data section to the DAMe tutorial section.
