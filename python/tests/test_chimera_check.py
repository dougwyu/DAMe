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
