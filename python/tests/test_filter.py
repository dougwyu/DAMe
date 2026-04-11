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
