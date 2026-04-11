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
