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
    out = convert(fna, max_length=65)   # max_length set but NOT usearch -> no padding
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

    # Also check long-form legacy spellings
    args2 = parser.parse_args(["convert", "--inFasta", "y.fna",
                                "--minLength", "3", "--maxLength", "200"])
    assert args2.in_fasta == "y.fna"
    assert args2.min_length == 3
    assert args2.max_length == 200


def test_sample_fastas_usearch_padded(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    fna = write_fna(tmp_path)
    convert(fna, usearch=True, max_length=65, sample_fastas=True)
    assert os.path.isfile("SampleFastas/Sample1.fixed.fasta")
    s1 = open("SampleFastas/Sample1.fixed.fasta").readlines()
    # Sample1 record 1: 60-nt seq padded to 65
    assert s1[0] == ">Sample1;size=9\n"
    assert s1[1].rstrip("\n") == "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG" + "N" * 5
    assert len(s1[1].rstrip("\n")) == 65
