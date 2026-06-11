# Sort Mismatch Tolerance Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add `-m` (primer mismatch) and `-mt` (tag mismatch) options to `dame-py sort` using a tag-anchored search strategy.

**Architecture:** When `m=0` and `mt=0` (defaults), the existing `GetPiecesInfo` regex path runs unchanged. When either is non-zero, the new `GetPiecesInfoMismatch` function is called instead. It anchors to known tag positions and checks primers at the exact expected offset, avoiding false matches inside the amplicon. Best-match (fewest total mismatches) wins; ties are discarded as ambiguous.

**Tech Stack:** Python 3.11+, pytest, no new dependencies.

---

## File Map

| File | Change |
|------|--------|
| `python/dame/modules_sort.py` | Add `IUPAC_SETS`, `iupac_match`, `hamming_iupac`; extend `readPrimers` to store raw sequences at indices [2][3]; add `GetPiecesInfoMismatch` |
| `python/dame/sort.py` | Add `-m` / `-mt` CLI args; dispatch to `GetPiecesInfoMismatch` when non-zero |
| `python/tests/test_sort.py` | Add tests for new helpers and `GetPiecesInfoMismatch`; update `test_readPrimers` for new structure |

---

## Task 1: IUPAC helper functions

**Files:**
- Modify: `python/dame/modules_sort.py` (add after the `RC` function)
- Test: `python/tests/test_sort.py`

- [ ] **Step 1: Write the failing tests**

Add to `python/tests/test_sort.py`:

```python
from dame.modules_sort import RC, readTags, readPrimers, FillHAP, GetPiecesInfo, \
    iupac_match, hamming_iupac, GetPiecesInfoMismatch


def test_iupac_match_exact_bases():
    assert iupac_match('A', 'A') is True
    assert iupac_match('C', 'C') is True
    assert iupac_match('A', 'C') is False


def test_iupac_match_N_matches_all():
    for b in ('A', 'C', 'G', 'T'):
        assert iupac_match('N', b) is True


def test_iupac_match_R_matches_purines():
    assert iupac_match('R', 'A') is True
    assert iupac_match('R', 'G') is True
    assert iupac_match('R', 'C') is False
    assert iupac_match('R', 'T') is False


def test_iupac_match_case_insensitive():
    assert iupac_match('a', 'A') is True
    assert iupac_match('N', 't') is True


def test_hamming_iupac_no_mismatch():
    assert hamming_iupac('ACGT', 'ACGT') == 0


def test_hamming_iupac_one_mismatch():
    assert hamming_iupac('ACGT', 'ACGC') == 1


def test_hamming_iupac_iupac_not_counted():
    assert hamming_iupac('ACGN', 'ACGT') == 0
    assert hamming_iupac('ACGR', 'ACGA') == 0
    assert hamming_iupac('ACGR', 'ACGC') == 1
```

- [ ] **Step 2: Run tests to confirm they fail**

```bash
cd /Users/jiyinqiu/src/cowork/DAMe-master-v2.3/python
python -m pytest tests/test_sort.py::test_iupac_match_exact_bases -v
```

Expected: `ImportError: cannot import name 'iupac_match'`

- [ ] **Step 3: Add IUPAC helpers to `modules_sort.py`**

Insert after the `RC` function (after line 11), before `readTags`:

```python
IUPAC_SETS = {
    'A': {'A'}, 'C': {'C'}, 'G': {'G'}, 'T': {'T'},
    'M': {'A', 'C'}, 'R': {'A', 'G'}, 'W': {'A', 'T'},
    'S': {'C', 'G'}, 'Y': {'C', 'T'}, 'K': {'G', 'T'},
    'V': {'A', 'C', 'G'}, 'H': {'A', 'C', 'T'},
    'D': {'A', 'G', 'T'}, 'B': {'C', 'G', 'T'},
    'N': {'A', 'C', 'G', 'T'},
}


def iupac_match(p, r):
    """Return True if read base r satisfies IUPAC primer base p."""
    return r.upper() in IUPAC_SETS.get(p.upper(), set())


def hamming_iupac(primer, region):
    """Count positions where region does not satisfy the IUPAC primer pattern."""
    return sum(0 if iupac_match(p, r) else 1 for p, r in zip(primer, region))
```

- [ ] **Step 4: Run tests to confirm they pass**

```bash
cd /Users/jiyinqiu/src/cowork/DAMe-master-v2.3/python
python -m pytest tests/test_sort.py -k "iupac or hamming" -v
```

Expected: 8 tests PASSED

- [ ] **Step 5: Commit**

```bash
cd /Users/jiyinqiu/src/cowork/DAMe-master-v2.3
git add python/dame/modules_sort.py python/tests/test_sort.py
git commit -m "feat: add iupac_match and hamming_iupac helpers to modules_sort"
```

---

## Task 2: Extend `readPrimers` to store raw sequences

**Files:**
- Modify: `python/dame/modules_sort.py` (`readPrimers` function)
- Test: `python/tests/test_sort.py` (update `test_readPrimers`)

- [ ] **Step 1: Update the existing test first**

Replace `test_readPrimers` in `python/tests/test_sort.py`:

```python
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
    assert len(result["CO1"]) == 4          # [regex_fwd, regex_rc, raw_fwd, raw_rc]
    assert len(result["CO1"][0]) == 2       # F and R regex patterns
    assert result["CO1"][2][0] == "ACGT"    # F raw (unchanged, no IUPAC here)
    assert result["CO1"][2][1] == "TTTT"    # R raw
    assert result["CO1"][3][0] == RC("ACGT")  # RC(F) raw
    assert result["CO1"][3][1] == RC("TTTT")  # RC(R) raw
```

- [ ] **Step 2: Run the updated test to confirm it fails**

```bash
cd /Users/jiyinqiu/src/cowork/DAMe-master-v2.3/python
python -m pytest tests/test_sort.py::test_readPrimers -v
```

Expected: `AssertionError: assert 2 == 4`

- [ ] **Step 3: Modify `readPrimers` in `modules_sort.py`**

Replace the entire `readPrimers` function:

```python
def readPrimers(primers, PRIMERS, AMBIG):
    with smart_open(primers) as f:
        for line in f:
            line = line.rstrip().split()
            if not line:
                continue
            if line[0] not in PRIMERS:
                PRIMERS[line[0]] = [[], [], [], []]
            Frc = RC(line[1])
            Rrc = RC(line[2])
            F = line[1]
            R = line[2]
            # Store raw IUPAC sequences before regex substitution (indices 2, 3)
            PRIMERS[line[0]][2].append(F)
            PRIMERS[line[0]][2].append(R)
            PRIMERS[line[0]][3].append(Frc)
            PRIMERS[line[0]][3].append(Rrc)
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
```

- [ ] **Step 4: Run all sort tests to confirm they pass**

```bash
cd /Users/jiyinqiu/src/cowork/DAMe-master-v2.3/python
python -m pytest tests/test_sort.py -v
```

Expected: all tests PASSED (the existing `GetPiecesInfo` tests still pass because indices [0] and [1] are unchanged)

- [ ] **Step 5: Commit**

```bash
cd /Users/jiyinqiu/src/cowork/DAMe-master-v2.3
git add python/dame/modules_sort.py python/tests/test_sort.py
git commit -m "feat: extend readPrimers to store raw IUPAC sequences at indices [2][3]"
```

---

## Task 3: Implement `GetPiecesInfoMismatch`

**Files:**
- Modify: `python/dame/modules_sort.py` (add after `GetPiecesInfo`)
- Test: `python/tests/test_sort.py`

- [ ] **Step 1: Write the failing tests**

Add to `python/tests/test_sort.py`:

```python
# ── shared fixtures for GetPiecesInfoMismatch tests ──────────────────────────

def _make_tags():
    return {
        "Tag1": ["AAAA", RC("AAAA")],   # RC = TTTT
        "Tag2": ["CCCC", RC("CCCC")],   # RC = GGGG
        "Tag3": ["GGGG", RC("GGGG")],   # RC = CCCC
        "Tag4": ["TTTT", RC("TTTT")],   # RC = AAAA
    }


def _make_primers():
    # F=ACGT, R=TGCA, RC(F)=RC("ACGT")="ACGT", RC(R)=RC("TGCA")="TGCA"
    # (both palindromes — no IUPAC ambiguity in this fixture)
    raw = ["ACGT", "TGCA"]
    raw_rc = [RC("ACGT"), RC("TGCA")]
    return {"CO1": [raw[:], raw_rc[:], raw[:], raw_rc[:]]}


# Forward read: AAAA + ACGT + ATATATATATAT + TGCA + GGGG
# tag2_rc = TAGS["Tag2"][1] = RC("CCCC") = "GGGG"
_FWD = "AAAA" + "ACGT" + "ATATATATATAT" + "TGCA" + "GGGG"
_AMP = "ATATATATATAT"


def test_mismatch_exact_match():
    result = GetPiecesInfoMismatch(_FWD, _make_primers(), _make_tags(), False, 0, 0)
    assert result == ["Tag1", "Tag2", "CO1", _AMP]


def test_mismatch_exact_keep_primers():
    result = GetPiecesInfoMismatch(_FWD, _make_primers(), _make_tags(), True, 0, 0)
    assert result == ["Tag1", "Tag2", "CO1", "ACGT" + _AMP + "TGCA"]


def test_primer_1mm_rejected_at_m0():
    # F primer ACGT → ACGG (1 mismatch at pos 3: T→G)
    read = "AAAA" + "ACGG" + _AMP + "TGCA" + "GGGG"
    result = GetPiecesInfoMismatch(read, _make_primers(), _make_tags(), False, 0, 0)
    assert result == [1]


def test_primer_1mm_accepted_at_m1():
    read = "AAAA" + "ACGG" + _AMP + "TGCA" + "GGGG"
    result = GetPiecesInfoMismatch(read, _make_primers(), _make_tags(), False, 1, 0)
    assert result == ["Tag1", "Tag2", "CO1", _AMP]


def test_primer_2mm_rejected_at_m1():
    # F primer ACGT → ACCC (2 mismatches: G→C, T→C)
    read = "AAAA" + "ACCC" + _AMP + "TGCA" + "GGGG"
    result = GetPiecesInfoMismatch(read, _make_primers(), _make_tags(), False, 1, 0)
    assert result == [1]


def test_primer_2mm_accepted_at_m2():
    read = "AAAA" + "ACCC" + _AMP + "TGCA" + "GGGG"
    result = GetPiecesInfoMismatch(read, _make_primers(), _make_tags(), False, 2, 0)
    assert result == ["Tag1", "Tag2", "CO1", _AMP]


def test_tag1_1mm_rejected_at_mt0():
    # tag1 AAAA → AAAT (1 mismatch at pos 3)
    read = "AAAT" + "ACGT" + _AMP + "TGCA" + "GGGG"
    result = GetPiecesInfoMismatch(read, _make_primers(), _make_tags(), False, 0, 0)
    assert result == [1]


def test_tag1_1mm_accepted_at_mt1():
    read = "AAAT" + "ACGT" + _AMP + "TGCA" + "GGGG"
    result = GetPiecesInfoMismatch(read, _make_primers(), _make_tags(), False, 0, 1)
    assert result == ["Tag1", "Tag2", "CO1", _AMP]


def test_empty_amplicon_rejected():
    # Primer and tag back-to-back with no amplicon
    read = "AAAA" + "ACGT" + "TGCA" + "GGGG"
    result = GetPiecesInfoMismatch(read, _make_primers(), _make_tags(), False, 0, 0)
    assert result == [1]


def test_no_primer_match_returns_1():
    result = GetPiecesInfoMismatch("NNNNNNNNNNNN", _make_primers(), _make_tags(), False, 0, 0)
    assert result == [1]
```

- [ ] **Step 2: Run tests to confirm they fail**

```bash
cd /Users/jiyinqiu/src/cowork/DAMe-master-v2.3/python
python -m pytest tests/test_sort.py::test_mismatch_exact_match -v
```

Expected: `ImportError: cannot import name 'GetPiecesInfoMismatch'`

- [ ] **Step 3: Add `GetPiecesInfoMismatch` to `modules_sort.py`**

Add after the `GetPiecesInfo` function:

```python
def GetPiecesInfoMismatch(line, PRIMERS, TAGS, keepPrimersSeq, max_primer_mm, max_tag_mm):
    """Tag-anchored primer search allowing up to max_primer_mm / max_tag_mm mismatches.

    Checks each known tag as a prefix/suffix (Hamming ≤ max_tag_mm), then verifies
    the adjacent primer at the anchored position (IUPAC-aware Hamming ≤ max_primer_mm).
    Returns the combination with the fewest total mismatches; ties are discarded.
    """
    read_len = len(line)
    best_mm = None
    best = None
    ambiguous = False

    for key in PRIMERS:
        F_raw, R_raw = PRIMERS[key][2][0], PRIMERS[key][2][1]
        Frc_raw, Rrc_raw = PRIMERS[key][3][0], PRIMERS[key][3][1]

        # ── forward orientation: [tag1_fwd][F][amplicon][RC(R)][tag2_rc] ──
        for t1 in TAGS:
            t1_seq = TAGS[t1][0]
            t1_len = len(t1_seq)
            if t1_len > read_len:
                continue
            t1_mm = hamming_iupac(t1_seq, line[:t1_len])
            if t1_mm > max_tag_mm:
                continue
            f_end = t1_len + len(F_raw)
            if f_end > read_len:
                continue
            f_mm = hamming_iupac(F_raw, line[t1_len:f_end])
            if f_mm > max_primer_mm:
                continue

            for t2 in TAGS:
                t2_seq = TAGS[t2][1]
                t2_len = len(t2_seq)
                rrc_end = read_len - t2_len
                rrc_start = rrc_end - len(Rrc_raw)
                if rrc_start < f_end or rrc_start < 0:
                    continue
                t2_mm = hamming_iupac(t2_seq, line[rrc_end:])
                if t2_mm > max_tag_mm:
                    continue
                rrc_mm = hamming_iupac(Rrc_raw, line[rrc_start:rrc_end])
                if rrc_mm > max_primer_mm:
                    continue
                between = line[t1_len:rrc_end] if keepPrimersSeq else line[f_end:rrc_start]
                if not between:
                    continue
                total = t1_mm + f_mm + rrc_mm + t2_mm
                if best_mm is None or total < best_mm:
                    best_mm, best, ambiguous = total, [t1, t2, key, between], False
                elif total == best_mm:
                    ambiguous = True

        # ── reverse orientation: [tag1_fwd][R][RC(amplicon)][RC(F)][tag2_rc] ──
        for t1 in TAGS:
            t1_seq = TAGS[t1][0]
            t1_len = len(t1_seq)
            if t1_len > read_len:
                continue
            t1_mm = hamming_iupac(t1_seq, line[:t1_len])
            if t1_mm > max_tag_mm:
                continue
            r_end = t1_len + len(R_raw)
            if r_end > read_len:
                continue
            r_mm = hamming_iupac(R_raw, line[t1_len:r_end])
            if r_mm > max_primer_mm:
                continue

            for t2 in TAGS:
                t2_seq = TAGS[t2][1]
                t2_len = len(t2_seq)
                frc_end = read_len - t2_len
                frc_start = frc_end - len(Frc_raw)
                if frc_start < r_end or frc_start < 0:
                    continue
                t2_mm = hamming_iupac(t2_seq, line[frc_end:])
                if t2_mm > max_tag_mm:
                    continue
                frc_mm = hamming_iupac(Frc_raw, line[frc_start:frc_end])
                if frc_mm > max_primer_mm:
                    continue
                between = RC(line[t1_len:frc_end] if keepPrimersSeq else line[r_end:frc_start])
                if not between:
                    continue
                total = t1_mm + r_mm + frc_mm + t2_mm
                if best_mm is None or total < best_mm:
                    best_mm, best, ambiguous = total, [t1, t2, key, between], False
                elif total == best_mm:
                    ambiguous = True

    if best is None or ambiguous:
        return [1]
    return best
```

- [ ] **Step 4: Run all new tests**

```bash
cd /Users/jiyinqiu/src/cowork/DAMe-master-v2.3/python
python -m pytest tests/test_sort.py -k "mismatch or primer_1mm or primer_2mm or tag1_1mm or empty_amplicon or no_primer" -v
```

Expected: 10 tests PASSED

- [ ] **Step 5: Run the full test suite to confirm nothing regressed**

```bash
cd /Users/jiyinqiu/src/cowork/DAMe-master-v2.3/python
python -m pytest tests/ -v
```

Expected: all tests PASSED

- [ ] **Step 6: Commit**

```bash
cd /Users/jiyinqiu/src/cowork/DAMe-master-v2.3
git add python/dame/modules_sort.py python/tests/test_sort.py
git commit -m "feat: add GetPiecesInfoMismatch with tag-anchored mismatch search"
```

---

## Task 4: Add `-m` / `-mt` CLI arguments to `sort.py`

**Files:**
- Modify: `python/dame/sort.py`
- Test: `python/tests/test_sort.py`

- [ ] **Step 1: Write the failing test**

Add to `python/tests/test_sort.py`:

```python
def test_sort_cli_m_mt_args():
    import argparse
    from dame.sort import register_subcommand
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers()
    register_subcommand(sub)
    args = parser.parse_args(["sort", "-fq", "x.fq", "-p", "p.txt", "-t", "t.txt",
                               "-m", "2", "-mt", "1"])
    assert args.m == 2
    assert args.mt == 1


def test_sort_cli_m_mt_defaults():
    import argparse
    from dame.sort import register_subcommand
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers()
    register_subcommand(sub)
    args = parser.parse_args(["sort", "-fq", "x.fq", "-p", "p.txt", "-t", "t.txt"])
    assert args.m == 0
    assert args.mt == 0
```

- [ ] **Step 2: Run tests to confirm they fail**

```bash
cd /Users/jiyinqiu/src/cowork/DAMe-master-v2.3/python
python -m pytest tests/test_sort.py::test_sort_cli_m_mt_args -v
```

Expected: `error: unrecognized arguments: -m 2`

- [ ] **Step 3: Update `sort.py` — add arguments and dispatch**

In `register_subcommand`, add two new arguments before `p.set_defaults(func=run)`:

```python
    p.add_argument("-m", dest="m", type=int, default=0,
                   help="Max mismatches allowed per primer (F and R each) [default 0]")
    p.add_argument("-mt", dest="mt", type=int, default=0,
                   help="Max mismatches allowed per tag (tag1 and tag2 each) [default 0]")
```

In `run`, update the import at the top of `sort.py` to include `GetPiecesInfoMismatch`:

```python
from dame.modules_sort import (
    readTags, readPrimers, GetPiecesInfo, GetPiecesInfoMismatch, FillHAP,
    PrintSortedCollapsedCountedSeqs, PrintSummaryFile,
    PrintSplitSummaryFile, read_valid_pairs,
)
```

In `run`, replace the single `GetPiecesInfo` call inside the read loop:

```python
            if args.m == 0 and args.mt == 0:
                Info = GetPiecesInfo(line, PRIMERS, TAGS, args.keepPrimersSeq)
            else:
                Info = GetPiecesInfoMismatch(line, PRIMERS, TAGS, args.keepPrimersSeq,
                                             args.m, args.mt)
```

In `main`, add the same two arguments before `run(parser.parse_args())`:

```python
    parser.add_argument("-m", dest="m", type=int, default=0)
    parser.add_argument("-mt", dest="mt", type=int, default=0)
```

- [ ] **Step 4: Run all tests**

```bash
cd /Users/jiyinqiu/src/cowork/DAMe-master-v2.3/python
python -m pytest tests/ -v
```

Expected: all tests PASSED

- [ ] **Step 5: Smoke-test the CLI**

```bash
cd /Users/jiyinqiu/src/cowork/DAMe-master-v2.3/python
dame-py sort --help
```

Expected: output includes `-m` and `-mt` in the help text.

- [ ] **Step 6: Final commit**

```bash
cd /Users/jiyinqiu/src/cowork/DAMe-master-v2.3
git add python/dame/sort.py python/tests/test_sort.py
git commit -m "feat: add -m/-mt mismatch options to dame-py sort"
```
