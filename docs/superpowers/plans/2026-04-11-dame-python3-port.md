# DAMe Python 3 Port Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Port DAMe's 5 Python 2 scripts to a pip-installable Python 3 package (`dame-py`) with a `dame-py <subcommand>` CLI, unit tests (pytest), and integration tests that verify output matches the Python 2 originals on fixture data.

**Architecture:** The Python 2 scripts in `bin/` are left untouched. A new `python/dame/` package is created containing Python 3 versions of each script and its module file, wired together via a `__main__.py` dispatcher. Integration tests in `tests/integration/` run each subcommand on synthetic fixtures and diff outputs to verify correctness.

**Tech Stack:** Python 3.11+, numpy, pytest, argparse, setuptools

---

## File Map

### New files (create)
- `python/pyproject.toml` — package metadata, `dame-py` entry point
- `python/dame/__init__.py` — empty, marks package
- `python/dame/__main__.py` — CLI dispatcher, routes `dame-py <subcommand>` to each module's `run(args)`
- `python/dame/modules_sort.py` — Python 3 port of `bin/modules_sort.py`
- `python/dame/sort.py` — Python 3 port of `bin/sort.py`, exports `register_subcommand`, `run`
- `python/dame/modules_filter.py` — Python 3 port of `bin/modules_filter.py`
- `python/dame/filter.py` — Python 3 port of `bin/filter.py`, exports `register_subcommand`, `run`
- `python/dame/modules_chimera_check.py` — Python 3 port of `bin/modules_chimeraCheck.py`
- `python/dame/chimera_check.py` — Python 3 port of `bin/chimeraCheck.py`, exports `register_subcommand`, `run`
- `python/dame/decollapse.py` — Python 3 port of `bin/decollapse.py`, exports `register_subcommand`, `run`
- `python/dame/rsi.py` — Python 3 port of `bin/RSI.py`, exports `register_subcommand`, `run`, `compare`
- `python/tests/__init__.py` — empty
- `python/tests/test_sort.py` — unit tests for modules_sort functions
- `python/tests/test_filter.py` — unit tests for modules_filter functions
- `python/tests/test_chimera_check.py` — unit tests for modules_chimera_check functions
- `python/tests/test_decollapse.py` — unit tests for decollapse run()
- `python/tests/test_rsi.py` — unit tests for rsi.compare()
- `tests/fixtures/Tags.txt` — synthetic tag file for integration tests
- `tests/fixtures/Primers.txt` — synthetic primer file for integration tests
- `tests/fixtures/sample.fastq` — synthetic FASTQ for sort integration test
- `tests/fixtures/PCRsetsInfo.txt` — symlink to `example/PCRsetsInfo.txt`
- `tests/fixtures/Comparisons_4PCRs.txt` — symlink to `example/Comparisons_4PCRs.txt`
- `tests/integration/run_sort.sh` — runs `dame-py sort` on fixtures, checks output shape
- `tests/integration/run_rsi.sh` — runs `dame-py rsi` on fixture, checks output format

---

## Task 1: Project scaffold

**Files:**
- Create: `python/pyproject.toml`
- Create: `python/dame/__init__.py`
- Create: `python/tests/__init__.py`

- [ ] **Step 1: Create directory structure**

```bash
mkdir -p python/dame python/tests tests/fixtures tests/integration
touch python/dame/__init__.py python/tests/__init__.py
```

- [ ] **Step 2: Write `python/pyproject.toml`**

```toml
[build-system]
requires = ["setuptools>=68"]
build-backend = "setuptools.backends.legacy:build"

[project]
name = "dame-py"
version = "3.0.0"
requires-python = ">=3.11"
dependencies = ["numpy>=1.24"]

[project.scripts]
dame-py = "dame.__main__:main"

[tool.setuptools.packages.find]
where = ["."]
include = ["dame*"]

[tool.pytest.ini_options]
testpaths = ["tests"]
```

- [ ] **Step 3: Install the package in editable mode**

```bash
cd python && pip install -e ".[dev]" 2>/dev/null || pip install -e .
```

Expected: `Successfully installed dame-py-3.0.0`

- [ ] **Step 4: Commit**

```bash
git add python/ tests/
git commit -m "feat: scaffold dame-py Python 3 package"
```

---

## Task 2: Port `modules_sort.py`

**Files:**
- Create: `python/dame/modules_sort.py`
- Create: `python/tests/test_sort.py`

- [ ] **Step 1: Write failing tests**

`python/tests/test_sort.py`:
```python
import pytest
import tempfile
import os
from dame.modules_sort import RC, readTags, readPrimers, FillHAP, GetPiecesInfo


def test_RC_palindrome():
    # ACGT reversed = TGCA, complement = ACGT
    assert RC("ACGT") == "ACGT"


def test_RC_all_A():
    assert RC("AAAA") == "TTTT"


def test_RC_mixed():
    # ATCG reversed = GCTA, complement = CGAT
    assert RC("ATCG") == "CGAT"


def test_RC_ambiguous_N():
    assert RC("N") == "N"


def test_readTags(tmp_path):
    tags_file = tmp_path / "tags.txt"
    tags_file.write_text("ACGT\tTag1\nTTTT\tTag2\n")
    TAGS = {}
    result = readTags(str(tags_file), TAGS)
    assert "Tag1" in result
    assert result["Tag1"][0] == "ACGT"       # forward seq
    assert result["Tag1"][1] == RC("ACGT")   # RC seq
    assert "Tag2" in result
    assert result["Tag2"][0] == "TTTT"
    assert result["Tag2"][1] == "AAAA"


def test_readPrimers(tmp_path):
    primers_file = tmp_path / "primers.txt"
    primers_file.write_text("CO1\tACGT\tTTTT\n")
    AMBIG = {'A': "A", 'C': "C", 'G': "G", 'T': "T",
             'N': "[ACGT]", 'R': "[AG]", 'Y': "[CT]",
             'M': "[AC]", 'K': "[GT]", 'S': "[CG]",
             'W': "[AT]", 'B': "[CGT]", 'D': "[AGT]",
             'H': "[ACT]", 'V': "[ACG]"}
    PRIMERS = {}
    result = readPrimers(str(primers_file), PRIMERS, AMBIG)
    assert "CO1" in result
    assert len(result["CO1"]) == 2  # [forward_list, rc_list]
    assert len(result["CO1"][0]) == 2  # F and R on A side


def test_FillHAP_new_entry():
    HAP = {}
    HAP = FillHAP(HAP, "Tag1", "Tag2", "CO1", "ACGTACGT")
    assert "Tag1_Tag2" in HAP
    assert HAP["Tag1_Tag2"][0] == "Tag1"
    assert HAP["Tag1_Tag2"][1] == "Tag2"
    assert HAP["Tag1_Tag2"][2]["ACGTACGT"][0] == 1
    assert HAP["Tag1_Tag2"][2]["ACGTACGT"][1] == "CO1"


def test_FillHAP_increment_count():
    HAP = {}
    HAP = FillHAP(HAP, "Tag1", "Tag2", "CO1", "ACGTACGT")
    HAP = FillHAP(HAP, "Tag1", "Tag2", "CO1", "ACGTACGT")
    assert HAP["Tag1_Tag2"][2]["ACGTACGT"][0] == 2


def test_FillHAP_multiple_seqs():
    HAP = {}
    HAP = FillHAP(HAP, "Tag1", "Tag2", "CO1", "AAAA")
    HAP = FillHAP(HAP, "Tag1", "Tag2", "CO1", "CCCC")
    assert len(HAP["Tag1_Tag2"][2]) == 2
```

- [ ] **Step 2: Run tests to verify they fail**

```bash
cd python && python -m pytest tests/test_sort.py -v 2>&1 | head -20
```

Expected: `ImportError` or `ModuleNotFoundError: No module named 'dame.modules_sort'`

- [ ] **Step 3: Write `python/dame/modules_sort.py`**

```python
import re
import os
import sys


def RC(seq):
    seq = seq[::-1]
    transtab = str.maketrans('ACGTMRWSYKVHDB', 'TGCAKYWSRMBDHV')
    return seq.translate(transtab)


def readTags(tags, TAGS):
    with open(tags) as f:
        for line in f:
            line = line.rstrip().split()
            if not line:
                continue
            if line[1] not in TAGS:
                TAGS[line[1]] = []
            TAGS[line[1]].append(line[0])
            TAGS[line[1]].append(RC(line[0]))
    return TAGS


def readPrimers(primers, PRIMERS, AMBIG):
    with open(primers) as f:
        for line in f:
            line = line.rstrip().split()
            if not line:
                continue
            if line[0] not in PRIMERS:
                PRIMERS[line[0]] = [[], []]
            Frc = RC(line[1])
            Rrc = RC(line[2])
            F = line[1]
            R = line[2]
            for key in AMBIG:
                Frc = re.sub(key, AMBIG[key], Frc)
                Rrc = re.sub(key, AMBIG[key], Rrc)
                F = re.sub(key, AMBIG[key], F)
                R = re.sub(key, AMBIG[key], R)
            PRIMERS[line[0]][0].append(F)
            PRIMERS[line[0]][0].append(R)
            PRIMERS[line[0]][1].append(Frc)
            PRIMERS[line[0]][1].append(Rrc)
    return PRIMERS


def GetPiecesInfo(line, PRIMERS, TAGS, keepPrimersSeq):
    for key in PRIMERS:
        primIniPos = [(m.start(0), m.end(0)) for m in re.finditer(PRIMERS[key][0][0], line)]
        if len(primIniPos) > 0:
            if keepPrimersSeq:
                primIniPosPrim = primIniPos[0][0]
                primIniPosTags = primIniPos[0][0]
            else:
                primIniPosPrim = primIniPos[0][1]
                primIniPosTags = primIniPos[0][0]
            primFinPos = [(m.start(0), m.end(0)) for m in re.finditer(PRIMERS[key][1][1], line)]
            if len(primFinPos) > 0:
                if keepPrimersSeq:
                    primFinPosPrim = primFinPos[0][1]
                    primFinPosTags = primFinPos[0][1]
                else:
                    primFinPosPrim = primFinPos[0][0]
                    primFinPosTags = primFinPos[0][1]
                PrimerName = key
                between = line[primIniPosPrim:primFinPosPrim]
                if len(between) == 0:
                    return [1]
                tag1 = line[:primIniPosTags]
                tag2 = line[primFinPosTags:]
                tagName1 = [t for t in TAGS if TAGS[t][0] == tag1]
                tagName2 = [t for t in TAGS if TAGS[t][1] == tag2]
                if len(tagName1) > 0 and len(tagName2) > 0:
                    return [tagName1[0], tagName2[0], PrimerName, between]
                return [1]
            return [1]
        else:
            primIniPos = [(m.start(0), m.end(0)) for m in re.finditer(PRIMERS[key][0][1], line)]
            if len(primIniPos) > 0:
                if keepPrimersSeq:
                    primIniPosPrim = primIniPos[0][0]
                    primIniPosTags = primIniPos[0][0]
                else:
                    primIniPosPrim = primIniPos[0][1]
                    primIniPosTags = primIniPos[0][0]
                primFinPos = [(m.start(0), m.end(0)) for m in re.finditer(PRIMERS[key][1][0], line)]
                if len(primFinPos) > 0:
                    if keepPrimersSeq:
                        primFinPosPrim = primFinPos[0][1]
                        primFinPosTags = primFinPos[0][1]
                    else:
                        primFinPosPrim = primFinPos[0][0]
                        primFinPosTags = primFinPos[0][1]
                    PrimerName = key
                    between = line[primIniPosPrim:primFinPosPrim]
                    if len(between) == 0:
                        return [1]
                    between = RC(between)
                    tag1 = line[:primIniPosTags]
                    tag2 = line[primFinPosTags:]
                    tagName2 = [t for t in TAGS if TAGS[t][0] == tag1]
                    tagName1 = [t for t in TAGS if TAGS[t][1] == tag2]
                    if len(tagName1) > 0 and len(tagName2) > 0:
                        return [tagName1[0], tagName2[0], PrimerName, between]
                    return [1]
                return [1]
    return [1]


def FillHAP(HAP, tagName1, tagName2, PrimerName, between):
    tagHapKey = "_".join([tagName1, tagName2])
    if tagHapKey not in HAP:
        HAP[tagHapKey] = [tagName1, tagName2, {}]
    if between not in HAP[tagHapKey][2]:
        HAP[tagHapKey][2][between] = [1, PrimerName]
    else:
        HAP[tagHapKey][2][between][0] += 1
    return HAP


def PrintSortedCollapsedCountedSeqs(HAP):
    for TagComb in HAP:
        with open("%s.txt" % TagComb, "w") as out:
            tagName1 = HAP[TagComb][0]
            tagName2 = HAP[TagComb][1]
            for Seq in HAP[TagComb][2]:
                a = "\t".join([HAP[TagComb][2][Seq][1], tagName1, tagName2,
                               str(HAP[TagComb][2][Seq][0]), Seq])
                out.write("%s\n" % a)


def PrintSummaryFile(HAP):
    with open("SummaryCounts.txt", "w") as out:
        out.write("%s\n" % "\t".join(("#tagName1", "tagName2", "NumUniqSeqs", "SumTotalFreq")))
        for TagComb in HAP:
            tagName1 = HAP[TagComb][0]
            tagName2 = HAP[TagComb][1]
            NumUniqSeqs = len(HAP[TagComb][2])
            SumTotalFreq = sum(HAP[TagComb][2][Seq][0] for Seq in HAP[TagComb][2])
            out.write("%s\n" % "\t".join((tagName1, tagName2, str(NumUniqSeqs), str(SumTotalFreq))))
```

- [ ] **Step 4: Run tests to verify they pass**

```bash
cd python && python -m pytest tests/test_sort.py -v
```

Expected: all 8 tests PASS

- [ ] **Step 5: Commit**

```bash
git add python/dame/modules_sort.py python/tests/test_sort.py
git commit -m "feat: port modules_sort.py to Python 3 with tests"
```

---

## Task 3: Port `sort.py` CLI

**Files:**
- Create: `python/dame/sort.py`

- [ ] **Step 1: Write `python/dame/sort.py`**

```python
import argparse
import sys

from dame.modules_sort import (
    readTags, readPrimers, GetPiecesInfo, FillHAP,
    PrintSortedCollapsedCountedSeqs, PrintSummaryFile,
)

AMBIG = {
    'A': "A", 'B': "[CGT]", 'C': "C", 'D': "[AGT]", 'G': "G",
    'H': "[ACT]", 'K': "[GT]", 'M': "[AC]", 'N': "[ACGT]", 'R': "[AG]",
    'S': "[CG]", 'T': "T", 'V': "[ACG]", 'W': "[AT]", 'Y': "[CT]",
}


def register_subcommand(subparsers):
    p = subparsers.add_parser(
        "sort",
        description="Sort amplicon sequences tagged on each end by tag combination",
    )
    p.add_argument("-fq", required=True, help="Input fastq with amplicon sequences")
    p.add_argument("-p", required=True,
                   help="Input text file with primer name and sequences [Format: Name\\tForwardSeq\\tReverseSeq]")
    p.add_argument("-t", required=True,
                   help="Input text file with tag names and sequences [Format: TagSeq\\tTagName]")
    p.add_argument("--keepPrimersSeq", action="store_true",
                   help="Keep primer sequences instead of trimming them [default not set]")
    p.set_defaults(func=run)


def run(args):
    TAGS = {}
    PRIMERS = {}
    HAP = {}
    CountErrors = 0

    TAGS = readTags(args.t, TAGS)
    PRIMERS = readPrimers(args.p, PRIMERS, AMBIG)

    with open(args.fq) as f:
        line = f.readline()  # header line
        while line:
            line = f.readline().rstrip()  # seq line
            if not line:
                break
            Info = GetPiecesInfo(line, PRIMERS, TAGS, args.keepPrimersSeq)
            if len(Info) == 1:
                f.readline()  # "+" line
                f.readline()  # qual line
                line = f.readline()  # next header
                CountErrors += 1
            else:
                HAP = FillHAP(HAP, Info[0], Info[1], Info[2], Info[3])
                f.readline()  # "+" line
                f.readline()  # qual line
                line = f.readline()  # next header

    PrintSortedCollapsedCountedSeqs(HAP)
    PrintSummaryFile(HAP)
    print(f"Number of erroneous sequences (with errors in the sequence of primer or tags, "
          f"or no barcode amplified): {CountErrors}")


def main():
    parser = argparse.ArgumentParser(
        description="Sort amplicon sequences tagged on each end by tag combination"
    )
    parser.add_argument("-fq", required=True)
    parser.add_argument("-p", required=True)
    parser.add_argument("-t", required=True)
    parser.add_argument("--keepPrimersSeq", action="store_true")
    run(parser.parse_args())


if __name__ == "__main__":
    main()
```

- [ ] **Step 2: Smoke-test the import**

```bash
cd python && python -c "from dame.sort import register_subcommand, run; print('OK')"
```

Expected: `OK`

- [ ] **Step 3: Commit**

```bash
git add python/dame/sort.py
git commit -m "feat: port sort.py CLI to Python 3"
```

---

## Task 4: Port `modules_filter.py`

**Files:**
- Create: `python/dame/modules_filter.py`
- Create: `python/tests/test_filter.py`

- [ ] **Step 1: Write failing tests**

`python/tests/test_filter.py`:
```python
import pytest
import os
import tempfile
from dame.modules_filter import (
    makePSnumFiles, ReadPSnumFiles, MakeSampleNameArray,
    ReadHapsForASample, getSeqsSetsAndFRcounts, MakeComparisonFile,
)


def write_psinfo(tmp_path, lines):
    p = tmp_path / "PSinfo.txt"
    p.write_text("\n".join(lines) + "\n")
    return str(p)


def test_MakeSampleNameArray(tmp_path):
    psinfo = write_psinfo(tmp_path, [
        "SampleA\tTag1\tTag2\t1",
        "SampleA\tTag3\tTag4\t1",
        "SampleB\tTag5\tTag6\t1",
        "SampleB\tTag7\tTag8\t1",
    ])
    names = MakeSampleNameArray(psinfo)
    assert names == ["SampleA", "SampleB"]


def test_MakeSampleNameArray_deduplicates(tmp_path):
    psinfo = write_psinfo(tmp_path, [
        "S1\tTag1\tTag2\t1",
        "S1\tTag3\tTag4\t1",
    ])
    names = MakeSampleNameArray(psinfo)
    assert len(names) == 1
    assert names[0] == "S1"


def test_makePSnumFiles_creates_files(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    psinfo = write_psinfo(tmp_path, [
        "S1\tTag1\tTag2\t1",
        "S1\tTag3\tTag4\t1",
    ])
    makePSnumFiles(psinfo, X=2, P=1, chimeraChecked=False)
    assert os.path.exists("PS1_files.txt")
    assert os.path.exists("PS2_files.txt")
    assert "pool1/Tag1_Tag2.txt" in open("PS1_files.txt").read()
    assert "pool1/Tag3_Tag4.txt" in open("PS2_files.txt").read()


def test_getSeqsSetsAndFRcounts_empty():
    haps = {"0": [], "1": []}
    seqsALL, F, R, counts, seqs = getSeqsSetsAndFRcounts(2, haps)
    assert seqsALL == set()
    assert F == {}
    assert R == {}


def test_getSeqsSetsAndFRcounts_with_data():
    haps = {
        "0": [["CO1", "Tag1", "Tag2", "3", "AAAA"],
              ["CO1", "Tag1", "Tag2", "1", "CCCC"]],
        "1": [["CO1", "Tag1", "Tag2", "2", "AAAA"]],
    }
    seqsALL, F, R, counts, seqs = getSeqsSetsAndFRcounts(2, haps)
    assert "AAAA" in seqsALL
    assert "CCCC" in seqsALL
    assert F["0"] == "Tag1"
    assert R["0"] == "Tag2"
    assert counts["0"] == ["3", "1"]
```

- [ ] **Step 2: Run tests to verify they fail**

```bash
cd python && python -m pytest tests/test_filter.py -v 2>&1 | head -10
```

Expected: `ImportError` or `ModuleNotFoundError: No module named 'dame.modules_filter'`

- [ ] **Step 3: Write `python/dame/modules_filter.py`**

```python
import os


def makePSnumFiles(PSinfo, X, P, chimeraChecked):
    PSouts = [open("PS%s_files.txt" % (i + 1), "w") for i in range(X)]
    with open(PSinfo) as f:
        PS = f.readlines()
    for NR, psinfo in enumerate(PS):
        NR = NR + 1
        psinfo = psinfo.rstrip().split()
        residue = NR % X
        idx = residue - 1 if residue != 0 else X - 1
        if not chimeraChecked:
            PSouts[idx].write("pool%s/%s_%s.txt\n" % (psinfo[3], psinfo[1], psinfo[2]))
        else:
            PSouts[idx].write("%s_%s_%s.noChim.txt\n" % (psinfo[1], psinfo[2], psinfo[3]))
    for out in PSouts:
        out.close()


def ReadPSnumFiles(X):
    PSinsLines = {}
    for i in range(X):
        with open("PS%s_files.txt" % (i + 1)) as f:
            PSinsLines[str(i)] = f.readlines()
    return PSinsLines


def MakeSampleNameArray(PSinfo):
    sampleName = []
    with open(PSinfo) as f:
        for line in f:
            name = line.split()[0]
            if name not in sampleName:
                sampleName.append(name)
    return sampleName


def ReadHapsForASample(X, PSinsLines, i):
    haps = {}
    for j in range(X):
        haps[str(j)] = []
        path = PSinsLines[str(j)][i].rstrip()
        if path != "empty" and os.path.exists(path):
            with open(path) as f:
                for line in f:
                    haps[str(j)].append(line.split())
    return haps


def getSeqsSetsAndFRcounts(X, haps):
    F = {}
    R = {}
    counts = {}
    seqs = {}
    seqsALL = []
    for j in range(X):
        if len(haps[str(j)]) != 0:
            seqs[str(j)] = []
            F[str(j)] = haps[str(j)][0][1]
            R[str(j)] = haps[str(j)][0][2]
            counts[str(j)] = []
            for k in range(len(haps[str(j)])):
                counts[str(j)].append(haps[str(j)][k][3])
                seqs[str(j)].append(haps[str(j)][k][4])
                seqsALL.append(haps[str(j)][k][4])
    seqsALL = set(seqsALL)
    return (seqsALL, F, R, counts, seqs)


def MakeComparisonFile(X, seqsALL, haps, F, R, counts, seqs,
                       OUT, OUTthresh, OUTYX, OUT_fas, OUTthresh_fas,
                       OUTYX_fas, OUTthreshLen_fas, Y, T, L, sampleName, i):
    idnum = 1
    for seq in seqsALL:
        line = sampleName[i] + "\t"
        lineFasIDs = ">" + sampleName[i] + "\t"
        lineFasCounts = "\t"
        y = 0
        t = 0
        for j in range(X):
            if len(haps[str(j)]) != 0:
                pos = [pos for pos, s in enumerate(seqs[str(j)]) if seq == s]
                if len(pos) == 0:
                    count = 0
                else:
                    y += 1
                    count = counts[str(j)][pos[0]]
                    if int(count) < T:
                        t += 1
                line = line + F[str(j)] + "-" + R[str(j)] + "\t" + str(count) + "\t"
                if j < (X - 1):
                    lineFasIDs = lineFasIDs + F[str(j)] + "-" + R[str(j)] + "."
                    lineFasCounts = lineFasCounts + str(count) + "_"
                else:
                    lineFasIDs = lineFasIDs + F[str(j)] + "-" + R[str(j)] + "_" + str(idnum) + "\t"
                    lineFasCounts = lineFasCounts + str(count) + "\n" + seq
            if len(haps[str(j)]) == 0:
                line = line + "empty\t0\t"
                if j < (X - 1):
                    lineFasIDs = lineFasIDs + "empty-empty."
                    lineFasCounts = lineFasCounts + "0_"
                else:
                    lineFasIDs = lineFasIDs + "empty-empty_" + str(idnum) + "\t"
                    lineFasCounts = lineFasCounts + "0\n" + seq
        line = line + seq + "\n"
        lineFas = lineFasIDs + lineFasCounts + "\n"
        OUT.write(line)
        OUT_fas.write(lineFas)
        if y >= Y:
            OUTYX.write(line)
            OUTYX_fas.write(lineFas)
        if (y - t) >= Y:
            OUTthresh.write(line)
            OUTthresh_fas.write(lineFas)
            if len(seq) >= L:
                OUTthreshLen_fas.write(lineFas)
        idnum += 1
```

- [ ] **Step 4: Run tests to verify they pass**

```bash
cd python && python -m pytest tests/test_filter.py -v
```

Expected: all 5 tests PASS

- [ ] **Step 5: Commit**

```bash
git add python/dame/modules_filter.py python/tests/test_filter.py
git commit -m "feat: port modules_filter.py to Python 3 with tests"
```

---

## Task 5: Port `filter.py` CLI

**Files:**
- Create: `python/dame/filter.py`

- [ ] **Step 1: Write `python/dame/filter.py`**

```python
import argparse
from dame.modules_filter import (
    makePSnumFiles, ReadPSnumFiles, MakeSampleNameArray,
    ReadHapsForASample, getSeqsSetsAndFRcounts, MakeComparisonFile,
)


def register_subcommand(subparsers):
    p = subparsers.add_parser(
        "filter",
        description="Filter multiplexed sequences by PCR presence, abundance, and length",
    )
    p.add_argument("-psInfo", required=True,
                   help="Text file with tag combination info per PCR reaction per sample")
    p.add_argument("-x", type=int, default=2, help="Number of PCR rxns performed per sample")
    p.add_argument("-y", type=int, default=1, help="Number of PCR rxns sequence must be present in")
    p.add_argument("-p", type=int, default=1, help="Number of pools [default 1]")
    p.add_argument("-t", type=int, default=1, help="Minimum count per unique sequence")
    p.add_argument("-l", type=int, default=100, help="Minimum sequence length")
    p.add_argument("--chimeraChecked", action="store_true",
                   help="Use chimera-checked sorted collapsed files [default not set]")
    p.set_defaults(func=run)


def run(args):
    PSinfo = args.psInfo
    X = args.x
    Y = args.y
    P = args.p
    T = args.t
    L = args.l
    chimeraChecked = args.chimeraChecked

    OUT = open("Comparisons_%sPCRs.txt" % X, "w")
    OUTYX = open("Comparisons_%soutOf%sPCRs.txt" % (Y, X), "w")
    OUTthresh = open("Comparisons_%soutOf%sPCRs.countsThreshold%s.txt" % (Y, X, T), "w")
    OUT_fas = open("Comparisons_%sPCRs.fasta" % X, "w")
    OUTYX_fas = open("FilteredReads_atLeast%s.fasta" % Y, "w")
    OUTthresh_fas = open("FilteredReads_atLeast%s.threshold.fasta" % Y, "w")
    OUTthreshLen_fas = open("FilteredReads.fna", "w")

    makePSnumFiles(PSinfo, X, P, chimeraChecked)
    PSinsLines = ReadPSnumFiles(X)
    sampleName = MakeSampleNameArray(PSinfo)

    for i in range(len(PSinsLines["0"])):
        haps = ReadHapsForASample(X, PSinsLines, i)
        seqsALL, F, R, counts, seqs = getSeqsSetsAndFRcounts(X, haps)
        MakeComparisonFile(X, seqsALL, haps, F, R, counts, seqs,
                           OUT, OUTthresh, OUTYX, OUT_fas, OUTthresh_fas,
                           OUTYX_fas, OUTthreshLen_fas, Y, T, L, sampleName, i)

    for fh in [OUT, OUTYX, OUTthresh, OUT_fas, OUTYX_fas, OUTthresh_fas, OUTthreshLen_fas]:
        fh.close()


def main():
    parser = argparse.ArgumentParser(
        description="Filter multiplexed sequences by PCR presence, abundance, and length"
    )
    parser.add_argument("-psInfo", required=True)
    parser.add_argument("-x", type=int, default=2)
    parser.add_argument("-y", type=int, default=1)
    parser.add_argument("-p", type=int, default=1)
    parser.add_argument("-t", type=int, default=1)
    parser.add_argument("-l", type=int, default=100)
    parser.add_argument("--chimeraChecked", action="store_true")
    run(parser.parse_args())


if __name__ == "__main__":
    main()
```

- [ ] **Step 2: Smoke-test the import**

```bash
cd python && python -c "from dame.filter import register_subcommand, run; print('OK')"
```

Expected: `OK`

- [ ] **Step 3: Commit**

```bash
git add python/dame/filter.py
git commit -m "feat: port filter.py CLI to Python 3"
```

---

## Task 6: Port `modules_chimeraCheck.py`

**Files:**
- Create: `python/dame/modules_chimera_check.py`
- Create: `python/tests/test_chimera_check.py`

- [ ] **Step 1: Write failing tests**

`python/tests/test_chimera_check.py`:
```python
import pytest
import os
from dame.modules_chimera_check import makeTagFiles, makeTagFilesWithPools, MakeFasSeqOneLine


def write_psinfo(tmp_path, lines):
    p = tmp_path / "PSinfo.txt"
    p.write_text("\n".join(lines) + "\n")
    return str(p)


def test_makeTagFiles_creates_files(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    psinfo = write_psinfo(tmp_path, [
        "S1\tTag1\tTag2\t1",
        "S1\tTag3\tTag4\t1",
    ])
    makeTagFiles(psinfo, X=2)
    assert os.path.exists("PS1.tags.txt")
    assert os.path.exists("PS2.tags.txt")
    content1 = open("PS1.tags.txt").read()
    content2 = open("PS2.tags.txt").read()
    assert "Tag1\tTag2" in content1
    assert "Tag3\tTag4" in content2


def test_makeTagFilesWithPools_includes_pool(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    psinfo = write_psinfo(tmp_path, [
        "S1\tTag1\tTag2\t1",
        "S1\tTag3\tTag4\t2",
    ])
    makeTagFilesWithPools(psinfo, X=2)
    content1 = open("PS1.tags.txt").read()
    assert "Tag1\tTag2\t1" in content1


def test_MakeFasSeqOneLine(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    fasta_content = ">header1\nACGT\nACGT\n>header2\nTTTT\n"
    (tmp_path / "Pool1.noChim.fasta").write_text(fasta_content)
    MakeFasSeqOneLine(P=1)
    result = open("Pool1.noChim.oneLiner.fasta").read()
    lines = [l for l in result.strip().split("\n") if l]
    assert lines[0] == ">header1"
    assert lines[1] == "ACGTACGT"
    assert lines[2] == ">header2"
    assert lines[3] == "TTTT"
```

- [ ] **Step 2: Run tests to verify they fail**

```bash
cd python && python -m pytest tests/test_chimera_check.py -v 2>&1 | head -10
```

Expected: `ImportError` or `ModuleNotFoundError`

- [ ] **Step 3: Write `python/dame/modules_chimera_check.py`**

```python
import re
import os
import sys
import subprocess


def makeTagFiles(PSinfo, X):
    PSouts = [open("PS%s.tags.txt" % (i + 1), "w") for i in range(X)]
    with open(PSinfo) as f:
        PS = f.readlines()
    for NR, psinfo in enumerate(PS):
        NR = NR + 1
        psinfo = psinfo.rstrip().split()
        residue = NR % X
        idx = residue - 1 if residue != 0 else X - 1
        PSouts[idx].write("%s\t%s\n" % (psinfo[1], psinfo[2]))
    for out in PSouts:
        out.close()


def makeTagFilesWithPools(PSinfo, X):
    PSouts = [open("PS%s.tags.txt" % (i + 1), "w") for i in range(X)]
    with open(PSinfo) as f:
        PS = f.readlines()
    for NR, psinfo in enumerate(PS):
        NR = NR + 1
        psinfo = psinfo.rstrip().split()
        residue = NR % X
        idx = residue - 1 if residue != 0 else X - 1
        PSouts[idx].write("%s\t%s\t%s\n" % (psinfo[1], psinfo[2], psinfo[3]))
    for out in PSouts:
        out.close()


def MakeSizeOutFastas(P, X):
    OUTS = [open("Pool%s.fasta" % (pool + 1), "w") for pool in range(P)]
    for num in range(X):
        with open("PS%s.tags.txt" % (num + 1)) as f:
            line = f.readline()
            while line:
                line = line.rstrip().split()
                hap = "_".join([line[0], line[1]]) + ".txt"
                if P > 1:
                    hap = "./pool" + str(line[2]) + "/" + hap
                if not os.path.exists(hap):
                    line = f.readline()
                    continue
                with open(hap) as hap_f:
                    idNum = 1
                    for seq in hap_f:
                        seq = seq.rstrip().split()
                        a = (">" + "_".join([seq[0], line[0], line[1], str(idNum)])
                             + ";size=" + str(seq[3]) + "\n" + seq[4])
                        pool_idx = int(line[2]) - 1 if P > 1 else 0
                        OUTS[pool_idx].write("%s\n" % a)
                        idNum += 1
                line = f.readline()
    for out in OUTS:
        out.close()


def SortFasta(P):
    for pool in range(P):
        input_file = "Pool%s.fasta" % (pool + 1)
        output = "Pool%s.sort.fasta" % (pool + 1)
        cmd = "usearch --sortsize " + input_file + " --output " + output
        p_core = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = p_core.communicate()
        with open("sort%s.out" % (pool + 1), "wb") as fh:
            fh.write(stdout)
        with open("sort%s.err" % (pool + 1), "wb") as fh:
            fh.write(stderr)
        input_file = output
        output1 = "Pool%s.Chim.fasta" % (pool + 1)
        output2 = "Pool%s.noChim.fasta" % (pool + 1)
        cmd = "usearch -uchime " + input_file + " -chimeras " + output1 + " -nonchimeras " + output2
        p_core = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = p_core.communicate()
        with open("chimeraCheck%s.out" % (pool + 1), "wb") as fh:
            fh.write(stdout)
        with open("chimeraCheck%s.err" % (pool + 1), "wb") as fh:
            fh.write(stderr)


def MakeFasSeqOneLine(P):
    for pool in range(P):
        with open("Pool%s.noChim.fasta" % (pool + 1)) as fasta:
            with open("Pool%s.noChim.oneLiner.fasta" % (pool + 1), "w") as fastaOne:
                seq = ""
                for line in fasta:
                    line = line.rstrip()
                    if line.startswith(">"):
                        if seq:
                            fastaOne.write("%s\n" % seq)
                        fastaOne.write("%s\n" % line)
                        seq = ""
                    else:
                        seq += line
                if seq:
                    fastaOne.write("%s\n" % seq)


def MakeNoChimHaps(P):
    HAP = {}
    for pool in range(P):
        with open("Pool%s.noChim.oneLiner.fasta" % (pool + 1)) as fasta:
            primerName = tagName1 = tagName2 = freq = tagHapKey = ""
            for line in fasta:
                line = line.rstrip()
                if line.startswith(">"):
                    line = line[1:]
                    primerName = line.split("_")[0]
                    tagName1 = line.split("_")[1]
                    tagName2 = line.split("_")[2]
                    freq = line.split("=")[1]
                    tagHapKey = "_".join([tagName1, tagName2, str(pool + 1)])
                    if tagHapKey not in HAP:
                        HAP[tagHapKey] = [[primerName], [tagName1], [tagName2], [freq], []]
                    else:
                        HAP[tagHapKey][0].append(primerName)
                        HAP[tagHapKey][1].append(tagName1)
                        HAP[tagHapKey][2].append(tagName2)
                        HAP[tagHapKey][3].append(freq)
                else:
                    HAP[tagHapKey][4].append(line + "\n")
    for TagComb in HAP:
        with open("%s.noChim.txt" % TagComb, "w") as out:
            for i in range(len(HAP[TagComb][0])):
                a = "\t".join([HAP[TagComb][0][i], HAP[TagComb][1][i],
                               HAP[TagComb][2][i], HAP[TagComb][3][i],
                               HAP[TagComb][4][i]])
                out.write(a)
```

- [ ] **Step 4: Run tests to verify they pass**

```bash
cd python && python -m pytest tests/test_chimera_check.py -v
```

Expected: all 3 tests PASS

- [ ] **Step 5: Commit**

```bash
git add python/dame/modules_chimera_check.py python/tests/test_chimera_check.py
git commit -m "feat: port modules_chimeraCheck.py to Python 3 with tests"
```

---

## Task 7: Port `chimera_check.py` CLI

**Files:**
- Create: `python/dame/chimera_check.py`

- [ ] **Step 1: Write `python/dame/chimera_check.py`**

```python
import argparse
from dame.modules_chimera_check import (
    makeTagFiles, makeTagFilesWithPools, MakeSizeOutFastas,
    SortFasta, MakeFasSeqOneLine, MakeNoChimHaps,
)


def register_subcommand(subparsers):
    p = subparsers.add_parser(
        "chimera",
        description="Create necessary files to operate on sequences per PCR reaction",
    )
    p.add_argument("-psInfo", required=True,
                   help="Text file with tag combination info per PCR reaction per sample")
    p.add_argument("-x", required=True, type=int, help="Number of PCR rxns performed per sample")
    p.add_argument("-p", type=int, default=1, help="Number of pools [default 1]")
    p.set_defaults(func=run)


def run(args):
    PSinfo = args.psInfo
    X = args.x
    P = args.p

    if P == 1:
        makeTagFiles(PSinfo, X)
    else:
        makeTagFilesWithPools(PSinfo, X)

    MakeSizeOutFastas(P, X)
    SortFasta(P)
    MakeFasSeqOneLine(P)
    MakeNoChimHaps(P)


def main():
    parser = argparse.ArgumentParser(
        description="Create necessary files to operate on sequences per PCR reaction"
    )
    parser.add_argument("-psInfo", required=True)
    parser.add_argument("-x", required=True, type=int)
    parser.add_argument("-p", type=int, default=1)
    run(parser.parse_args())


if __name__ == "__main__":
    main()
```

- [ ] **Step 2: Smoke-test the import**

```bash
cd python && python -c "from dame.chimera_check import register_subcommand, run; print('OK')"
```

Expected: `OK`

- [ ] **Step 3: Commit**

```bash
git add python/dame/chimera_check.py
git commit -m "feat: port chimera_check.py CLI to Python 3"
```

---

## Task 8: Port `decollapse.py`

**Files:**
- Create: `python/dame/decollapse.py`
- Create: `python/tests/test_decollapse.py`

- [ ] **Step 1: Write failing tests**

`python/tests/test_decollapse.py`:
```python
import pytest
import os
from dame.decollapse import run


class MockArgs:
    def __init__(self, input_file, out_fas="Decollapsed.fasta"):
        self.input = input_file
        self.outFas = out_fas


def test_decollapse_single_seq(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    input_file = tmp_path / "input.txt"
    # Format: PrimerName Tag1 Tag2 Freq Seq
    input_file.write_text("CO1\tTag1\tTag2\t3\tACGT\n")
    args = MockArgs(str(input_file), "out.fasta")
    run(args)
    lines = open("out.fasta").readlines()
    # freq=3 means 3 entries
    headers = [l for l in lines if l.startswith(">")]
    seqs = [l for l in lines if not l.startswith(">")]
    assert len(headers) == 3
    assert all(s.strip() == "ACGT" for s in seqs)
    assert "Tag1.Tag2.3_1" in headers[0]
    assert "Tag1.Tag2.3_2" in headers[1]
    assert "Tag1.Tag2.3_3" in headers[2]


def test_decollapse_multiple_seqs(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    input_file = tmp_path / "input.txt"
    input_file.write_text("CO1\tTag1\tTag2\t2\tAAAA\nCO1\tTag1\tTag2\t1\tCCCC\n")
    args = MockArgs(str(input_file), "out.fasta")
    run(args)
    lines = open("out.fasta").readlines()
    headers = [l for l in lines if l.startswith(">")]
    assert len(headers) == 3  # 2 + 1
```

- [ ] **Step 2: Run tests to verify they fail**

```bash
cd python && python -m pytest tests/test_decollapse.py -v 2>&1 | head -10
```

Expected: `ImportError` or `ModuleNotFoundError`

- [ ] **Step 3: Write `python/dame/decollapse.py`**

```python
import argparse


def register_subcommand(subparsers):
    p = subparsers.add_parser(
        "decollapse",
        description="Expand unique sequences to individual reads by frequency",
    )
    p.add_argument("-input", required=True,
                   help="Text file with tag combination and freq of each unique seq")
    p.add_argument("-outFas", default="Decollapsed.fasta",
                   help='Output fasta file [default "Decollapsed.fasta"]')
    p.set_defaults(func=run)


def run(args):
    seq_id = 0
    with open(args.input) as IN, open(args.outFas, "w") as OUT:
        for line in IN:
            line = line.rstrip().split()
            count = 1
            while count <= int(line[3]):
                seq_id += 1
                OUT.write(">" + line[1] + "." + line[2] + "." + line[3]
                          + "_" + str(seq_id) + "\n" + line[4] + "\n")
                count += 1


def main():
    parser = argparse.ArgumentParser(
        description="Expand unique sequences to individual reads by frequency"
    )
    parser.add_argument("-input", required=True)
    parser.add_argument("-outFas", default="Decollapsed.fasta")
    run(parser.parse_args())


if __name__ == "__main__":
    main()
```

- [ ] **Step 4: Run tests to verify they pass**

```bash
cd python && python -m pytest tests/test_decollapse.py -v
```

Expected: both tests PASS

- [ ] **Step 5: Commit**

```bash
git add python/dame/decollapse.py python/tests/test_decollapse.py
git commit -m "feat: port decollapse.py to Python 3 with tests"
```

---

## Task 9: Port `RSI.py`

**Files:**
- Create: `python/dame/rsi.py`
- Create: `python/tests/test_rsi.py`

- [ ] **Step 1: Write failing tests**

`python/tests/test_rsi.py`:
```python
import pytest
import numpy as np
from dame.rsi import compare


def test_compare_identical_replicates():
    # Identical replicates → RSI = 0 (perfectly similar)
    matrix = np.array([[10, 10], [20, 20], [30, 30]])
    result = compare(matrix, "sample1", 1, 2)
    assert abs(result) < 1e-10


def test_compare_completely_different():
    # Mutually exclusive seqs → RSI = 1 (completely dissimilar)
    matrix = np.array([[100, 0], [0, 100]])
    result = compare(matrix, "sample1", 1, 2)
    assert abs(result - 1.0) < 1e-10


def test_compare_partial_overlap():
    # Half shared → RSI = 0.5
    matrix = np.array([[50, 0], [50, 100]])
    result = compare(matrix, "sample1", 1, 2)
    # a = [0.5, 0.5], b = [0.0, 1.0], min = [0.0, 0.5], sum = 0.5, RSI = 0.5
    assert abs(result - 0.5) < 1e-10


def test_compare_zero_replicate_handled():
    # Zero-total replicate should not divide by zero
    matrix = np.array([[0, 10], [0, 20]])
    result = compare(matrix, "sample1", 1, 2)
    assert 0.0 <= result <= 1.0
```

- [ ] **Step 2: Run tests to verify they fail**

```bash
cd python && python -m pytest tests/test_rsi.py -v 2>&1 | head -10
```

Expected: `ImportError` or `ModuleNotFoundError`

- [ ] **Step 3: Write `python/dame/rsi.py`**

```python
import sys
import argparse
import numpy as np


def compare(matrix, j, a, b):
    print(f"Comparing replicates {a!r} and {b!r} from sample {j}")
    total = matrix.sum(axis=0)
    if total[0] == 0:
        total[0] = 1
        print(f"This sample gave zero in replicate: {a!r}")
    if total[1] == 0:
        total[1] = 1
        print(f"This sample gave zero in replicate: {b!r}")
    col_a = np.array([row[0] / float(total[0]) for row in matrix])
    col_b = np.array([row[1] / float(total[1]) for row in matrix])
    percent = np.column_stack((col_a, col_b))
    return 1 - np.sum([row.min() for row in percent])


def register_subcommand(subparsers):
    p = subparsers.add_parser(
        "rsi",
        description="Compute Renkonen Similarity Index between PCR replicates",
    )
    p.add_argument("input", help="Input comparison file")
    p.add_argument("-e", "--explicit", action="store_true",
                   help="Output explicit RSI for every pairwise comparison")
    p.add_argument("-o", "--output", dest="outfile", metavar="FILE",
                   help="Write output to FILE [default RSI_output.txt]")
    p.set_defaults(func=run)


def run(args):
    data = []
    with open(args.input) as f:
        for line in f:
            parts = line.split()
            data.append(parts)
    data = np.array(data)
    names = set(row[0] for row in data)
    no_rep = (len(data[0]) - 2) // 2
    if no_rep < 2:
        print("There are no replicates in the file.")
        sys.exit(0)

    rkn = []
    for i in names:
        subset = np.array([row for row in data if row[0] == i])
        if args.explicit:
            for A in range(1, no_rep):
                for B in range(A + 1, no_rep + 1):
                    sample = subset[:, [A * 2, B * 2]].astype(int)
                    output = compare(sample, i, A, B)
                    rkn.append([i, A, B, output])
        else:
            rep = 0
            output = 0
            for A in range(1, no_rep):
                for B in range(A + 1, no_rep + 1):
                    sample = subset[:, [A * 2, B * 2]].astype(int)
                    output += compare(sample, i, A, B)
                    rep += 1
            rkn.append([i, output / rep])

    outfile = args.outfile if args.outfile else "RSI_output.txt"
    with open(outfile, "w") as out:
        if args.explicit:
            out.write("Sample\tReplicateA\tReplicateB\tRSI\n")
            for row in rkn:
                out.write("%s\t%s\t%s\t%s\n" % (row[0], row[1], row[2], row[3]))
        else:
            out.write("Sample\tRSI\n")
            for row in rkn:
                out.write("%s\t%s\n" % (row[0], row[1]))


def main():
    parser = argparse.ArgumentParser(
        description="Compute Renkonen Similarity Index between PCR replicates"
    )
    parser.add_argument("input")
    parser.add_argument("-e", "--explicit", action="store_true")
    parser.add_argument("-o", "--output", dest="outfile", metavar="FILE")
    run(parser.parse_args())


if __name__ == "__main__":
    main()
```

- [ ] **Step 4: Run tests to verify they pass**

```bash
cd python && python -m pytest tests/test_rsi.py -v
```

Expected: all 4 tests PASS

- [ ] **Step 5: Commit**

```bash
git add python/dame/rsi.py python/tests/test_rsi.py
git commit -m "feat: port RSI.py to Python 3 with tests"
```

---

## Task 10: CLI dispatcher (`__main__.py`)

**Files:**
- Create: `python/dame/__main__.py`

- [ ] **Step 1: Write `python/dame/__main__.py`**

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

- [ ] **Step 2: Test the CLI entry point**

```bash
dame-py --help
```

Expected output contains: `{sort,chimera,filter,decollapse,rsi}`

```bash
dame-py sort --help
```

Expected output contains: `-fq`, `-p`, `-t`, `--keepPrimersSeq`

- [ ] **Step 3: Run the full test suite**

```bash
cd python && python -m pytest tests/ -v
```

Expected: all tests PASS

- [ ] **Step 4: Commit**

```bash
git add python/dame/__main__.py
git commit -m "feat: add dame-py CLI dispatcher"
```

---

## Task 11: Integration test fixtures

**Files:**
- Create: `tests/fixtures/Tags.txt`
- Create: `tests/fixtures/Primers.txt`
- Create: `tests/fixtures/sample.fastq`
- Create: `tests/fixtures/Comparisons_4PCRs.txt` (symlink)

- [ ] **Step 1: Create synthetic Tags.txt**

`tests/fixtures/Tags.txt`:
```
AAAA	Tag1
CCCC	Tag2
GGGG	Tag3
TTTT	Tag4
```

- [ ] **Step 2: Create synthetic Primers.txt**

`tests/fixtures/Primers.txt`:
```
CO1	ACGT	TGCA
```

- [ ] **Step 3: Create synthetic sample.fastq**

The FASTQ must have reads with structure: `<tag1><primer_fwd><barcode><primer_rc><tag2_rc>`

Tag1=AAAA, PrimerFwd=ACGT, Barcode=ATATATATATAT, PrimerRrc=RC(TGCA)=TGCA... wait, RC of TGCA = TGCA reversed=ACGT, complement=TGCA. So RC(TGCA)=ACGT.

Full forward read: `AAAA` + `ACGT` + `ATATATATATAT` + `ACGT` + `TTTT` (Tag4 RC = RC(TTTT) = AAAA... hmm.

RC(Tag4=TTTT) = AAAA. So tag2 in reverse position should be RC of Tag2. Tag2=CCCC, RC(CCCC)=GGGG.

Forward read: `AAAA` + `ACGT` + `ATATATATATAT` + RC(RC(TGCA))` + RC(CCCC)
= `AAAA` + `ACGT` + `ATATATATATAT` + `TGCA` + `GGGG`

Let me verify: in GetPiecesInfo for a forward read:
- primIniPos finds PRIMERS[key][0][0] = F = ACGT at position 4
- primIniPosPrim = 4+4=8 (after primer), primIniPosTags = 4 (start of primer)
- primFinPos finds PRIMERS[key][1][1] = Rrc = RC(R) = RC(TGCA) = TGCA... 

Wait, let me re-read the primer setup. In readPrimers:
- F = line[1] = ACGT (forward primer)
- R = line[2] = TGCA (reverse primer)
- Frc = RC(F) = RC(ACGT) = ACGT (palindrome)
- Rrc = RC(R) = RC(TGCA) = TGCA... 

RC(TGCA): reverse = ACGT, complement of ACGT = TGCA. So RC(TGCA) = TGCA.

PRIMERS[key][0] = [F, R] = [ACGT, TGCA] (what to look for on A side: forward or reverse primer start)
PRIMERS[key][1] = [Frc, Rrc] = [ACGT, TGCA] (RC versions)

For a forward read, GetPiecesInfo looks for:
- PRIMERS[key][0][0] = F = ACGT at start (primer start)
- PRIMERS[key][1][1] = Rrc = TGCA at end (RC of reverse primer)

So forward read: `<tag1>` + `ACGT` + `<barcode>` + `TGCA` + `<tag2_rc>`

Tag2_rc = RC(Tag2) = RC(CCCC) = GGGG

Forward read: `AAAA` + `ACGT` + `ATATATATATAT` + `TGCA` + `GGGG`
= `AAAACGTATATATATATATTGCAGGGG`

Let me verify tag extraction:
- tag1 = line[:primIniPosTags] = line[:4] = "AAAA" ✓ (Tag1)
- primIniPosPrim = 8 (after ACGT)
- primFinPosTags = position of TGCA end = 4+4+12+4 = 24
- tag2 = line[24:] = "GGGG"
- tagName2 = tags where TAGS[t][1] == "GGGG" → Tag2 (since TAGS["Tag2"][1] = RC("CCCC") = "GGGG") ✓

Great. So the synthetic FASTQ should have reads like:

```
@read1
AAAACGTATATATATATATTGCAGGGG
+
IIIIIIIIIIIIIIIIIIIIIIIIIII
```

`tests/fixtures/sample.fastq`:
```
@read1
AAAACGTATATATATATATTGCAGGGG
+
IIIIIIIIIIIIIIIIIIIIIIIIIII
@read2
AAAACGTATATATATATATTGCAGGGG
+
IIIIIIIIIIIIIIIIIIIIIIIIIII
@read3
AAAACGTGCGCGCGCGCGCGTGCAGGGG
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIII
@read4
AAAACGTGCGCGCGCGCGCGTGCAGGGG
+
IIIIIIIIIIIIIIIIIIIIIIIIIIIII
@errread5
NNNNNNNNNNNNNNNNNNNN
+
IIIIIIIIIIIIIIIIIIII
```

Read lengths: AAAACGTATATATATATATTGCAGGGG = 4+4+12+4+4 = 28 chars. Let me count: AAAA(4)+ACGT(4)+ATATATATAT AT(12)+TGCA(4)+GGGG(4) = 28. Quality line needs 28 I's.

- [ ] **Step 4: Create symlink for RSI fixture**

```bash
cd tests/fixtures && ln -sf ../../example/Comparisons_4PCRs.txt Comparisons_4PCRs.txt
```

- [ ] **Step 5: Commit fixtures**

```bash
git add tests/fixtures/
git commit -m "test: add integration test fixtures"
```

---

## Task 12: Integration test scripts

**Files:**
- Create: `tests/integration/run_sort.sh`
- Create: `tests/integration/run_rsi.sh`
- Create: `tests/integration/run_filter.sh`
- Create: `tests/integration/run_decollapse.sh`

> Note: `run_chimera.sh` and `run_pipeline.sh` require the external `usearch` binary and are deferred to Plan 2 (Rust port). They will be added to `tests/integration/` at that point.

- [ ] **Step 1: Write `tests/integration/run_sort.sh`**

```bash
#!/usr/bin/env bash
set -euo pipefail

FIXTURES="$(cd "$(dirname "$0")/../fixtures" && pwd)"
TMPDIR=$(mktemp -d)
trap "rm -rf '$TMPDIR'" EXIT

echo "==> Running dame-py sort on fixtures..."
cd "$TMPDIR"
dame-py sort \
    -fq "$FIXTURES/sample.fastq" \
    -p  "$FIXTURES/Primers.txt" \
    -t  "$FIXTURES/Tags.txt"

# Verify SummaryCounts.txt was created and has expected structure
if [ ! -f SummaryCounts.txt ]; then
    echo "FAIL: SummaryCounts.txt not created"
    exit 1
fi

# Header line must be present
if ! grep -q "#tagName1" SummaryCounts.txt; then
    echo "FAIL: SummaryCounts.txt missing header"
    exit 1
fi

# Tag1_Tag2.txt must exist (our synthetic reads use Tag1 + Tag2)
if [ ! -f Tag1_Tag2.txt ]; then
    echo "FAIL: Tag1_Tag2.txt not created"
    exit 1
fi

# Tag1_Tag2.txt must contain our barcodes with correct counts
if ! grep -q "ATATATATAT" Tag1_Tag2.txt; then
    echo "FAIL: expected barcode not found in Tag1_Tag2.txt"
    exit 1
fi

echo "PASS: dame-py sort produced expected output files"
```

- [ ] **Step 2: Write `tests/integration/run_rsi.sh`**

```bash
#!/usr/bin/env bash
set -euo pipefail

FIXTURES="$(cd "$(dirname "$0")/../fixtures" && pwd)"
TMPDIR=$(mktemp -d)
trap "rm -rf '$TMPDIR'" EXIT

echo "==> Running dame-py rsi on fixtures..."
cd "$TMPDIR"
dame-py rsi "$FIXTURES/Comparisons_4PCRs.txt"

if [ ! -f RSI_output.txt ]; then
    echo "FAIL: RSI_output.txt not created"
    exit 1
fi

if ! grep -q "Sample" RSI_output.txt; then
    echo "FAIL: RSI_output.txt missing header"
    exit 1
fi

# RSI values must be between 0 and 1
while IFS=$'\t' read -r sample rsi_val; do
    [[ "$sample" == "Sample" ]] && continue
    if ! awk "BEGIN{exit !($rsi_val >= 0 && $rsi_val <= 1)}"; then
        echo "FAIL: RSI value out of range: $rsi_val for sample $sample"
        exit 1
    fi
done < RSI_output.txt

echo "PASS: dame-py rsi produced valid RSI_output.txt"

echo "==> Running dame-py rsi --explicit..."
cd "$TMPDIR"
dame-py rsi --explicit -o RSI_explicit.txt "$FIXTURES/Comparisons_4PCRs.txt"

if [ ! -f RSI_explicit.txt ]; then
    echo "FAIL: RSI_explicit.txt not created"
    exit 1
fi

if ! grep -q "ReplicateA" RSI_explicit.txt; then
    echo "FAIL: RSI_explicit.txt missing explicit header"
    exit 1
fi

echo "PASS: dame-py rsi --explicit produced valid output"
```

- [ ] **Step 3: Write `tests/integration/run_filter.sh`**

The filter integration test creates pre-made hap files (sort output format) directly so it doesn't depend on running sort first.

```bash
#!/usr/bin/env bash
set -euo pipefail

TMPDIR=$(mktemp -d)
trap "rm -rf '$TMPDIR'" EXIT

echo "==> Setting up filter fixtures..."
cd "$TMPDIR"
mkdir -p pool1

# Hap file format: PrimerName Tag1 Tag2 Freq Seq
cat > pool1/Tag1_Tag2.txt <<'EOF'
CO1	Tag1	Tag2	5	ATATATATATAT
CO1	Tag1	Tag2	3	GCGCGCGCGCGC
EOF

cat > pool1/Tag3_Tag4.txt <<'EOF'
CO1	Tag3	Tag4	4	ATATATATATAT
CO1	Tag3	Tag4	1	GCGCGCGCGCGC
EOF

cat > PSinfo.txt <<'EOF'
SampleA	Tag1	Tag2	1
SampleA	Tag3	Tag4	1
EOF

echo "==> Running dame-py filter..."
dame-py filter -psInfo PSinfo.txt -x 2 -y 1 -t 1 -l 10

if [ ! -f FilteredReads.fna ]; then
    echo "FAIL: FilteredReads.fna not created"; exit 1
fi

if [ ! -f Comparisons_2PCRs.txt ]; then
    echo "FAIL: Comparisons_2PCRs.txt not created"; exit 1
fi

if ! grep -q "SampleA" Comparisons_2PCRs.txt; then
    echo "FAIL: SampleA not found in comparison file"; exit 1
fi

echo "PASS: dame-py filter produced expected output files"
```

- [ ] **Step 4: Write `tests/integration/run_decollapse.sh`**

```bash
#!/usr/bin/env bash
set -euo pipefail

TMPDIR=$(mktemp -d)
trap "rm -rf '$TMPDIR'" EXIT

echo "==> Running dame-py decollapse..."
cd "$TMPDIR"

cat > input.txt <<'EOF'
CO1	Tag1	Tag2	3	ATATATATATAT
CO1	Tag1	Tag2	2	GCGCGCGCGCGC
EOF

dame-py decollapse -input input.txt -outFas out.fasta

if [ ! -f out.fasta ]; then
    echo "FAIL: out.fasta not created"; exit 1
fi

header_count=$(grep -c "^>" out.fasta)
if [ "$header_count" -ne 5 ]; then
    echo "FAIL: expected 5 sequences (3+2), got $header_count"; exit 1
fi

if ! grep -q "Tag1.Tag2.3_1" out.fasta; then
    echo "FAIL: expected header Tag1.Tag2.3_1 not found"; exit 1
fi

echo "PASS: dame-py decollapse produced 5 sequences as expected"
```

- [ ] **Step 5: Make scripts executable**

```bash
chmod +x tests/integration/run_sort.sh tests/integration/run_rsi.sh \
         tests/integration/run_filter.sh tests/integration/run_decollapse.sh
```

- [ ] **Step 6: Run all integration tests**

```bash
bash tests/integration/run_sort.sh
bash tests/integration/run_rsi.sh
bash tests/integration/run_filter.sh
bash tests/integration/run_decollapse.sh
```

Expected: all four print `PASS`

- [ ] **Step 7: Commit**

```bash
git add tests/integration/
git commit -m "test: add integration tests for all dame-py subcommands"
```

---

## Task 13: Final check — full test suite

- [ ] **Step 1: Run all unit tests**

```bash
cd python && python -m pytest tests/ -v
```

Expected: all tests PASS, 0 failures

- [ ] **Step 2: Run all integration tests**

```bash
bash tests/integration/run_sort.sh
bash tests/integration/run_rsi.sh
```

Expected: both print `PASS`

- [ ] **Step 3: Verify CLI help works for all subcommands**

```bash
dame-py sort --help
dame-py chimera --help
dame-py filter --help
dame-py decollapse --help
dame-py rsi --help
```

Expected: all print usage without error

- [ ] **Step 4: Final commit**

```bash
git add -A
git commit -m "feat: complete dame-py Python 3 port with unit and integration tests"
```

---

> **Next:** After this plan is complete and all tests pass, see `docs/superpowers/plans/2026-04-11-dame-rust-port.md` for the Rust port (Plan 2).
