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


# ── malformed PSinfo ──────────────────────────────────────────────────────────
# filter.rs skips blank and short lines. These pin the Python side to the same
# behaviour, including that a skipped line still consumes a replicate slot so
# the PS file assignment stays in step between implementations.

def test_blank_psinfo_line_is_skipped(tmp_path, monkeypatch):
    psinfo = tmp_path / "PSinfo.txt"
    psinfo.write_text("Sample1\ttag1\ttag2\t1\nSample1\ttag3\ttag4\t2\n\n")
    monkeypatch.chdir(tmp_path)
    makePSnumFiles(str(psinfo), 2, 1, False)
    assert (tmp_path / "PS1_files.txt").read_text() == "pool1/tag1_tag2.txt\n"
    assert (tmp_path / "PS2_files.txt").read_text() == "pool2/tag3_tag4.txt\n"


def test_blank_psinfo_line_skipped_by_sample_name_array(tmp_path):
    psinfo = tmp_path / "PSinfo.txt"
    psinfo.write_text("Sample1\ttag1\ttag2\t1\n\nSample2\ttag3\ttag4\t2\n")
    assert MakeSampleNameArray(str(psinfo)) == ["Sample1", "Sample2"]


def test_short_psinfo_line_is_skipped(tmp_path, monkeypatch):
    # A three-column line has no pool number and is dropped, but still consumes
    # its replicate slot, so line 3 lands in PS1 rather than shifting up.
    psinfo = tmp_path / "PSinfo.txt"
    psinfo.write_text("Sample1\ttag1\ttag2\t1\nSample1\ttag3\ttag4\nSample2\ttag5\ttag6\t1\n")
    monkeypatch.chdir(tmp_path)
    makePSnumFiles(str(psinfo), 2, 1, False)
    assert (tmp_path / "PS1_files.txt").read_text() == "pool1/tag1_tag2.txt\npool1/tag5_tag6.txt\n"
    assert (tmp_path / "PS2_files.txt").read_text() == ""
