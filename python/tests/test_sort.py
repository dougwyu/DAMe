import pytest
import tempfile
import os
from dame.modules_sort import RC, readTags, readPrimers, FillHAP, GetPiecesInfo, iupac_matches, find_primer


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


def test_iupac_matches_exact():
    assert iupac_matches("A", "A")
    assert not iupac_matches("A", "C")


def test_iupac_matches_ambiguous():
    # R = A or G
    assert iupac_matches("R", "A")
    assert iupac_matches("R", "G")
    assert not iupac_matches("R", "C")
    # N = anything
    assert iupac_matches("N", "A")
    assert iupac_matches("N", "T")


def test_iupac_matches_non_acgt_read():
    # A read base of N never satisfies any primer code (matches Rust).
    assert not iupac_matches("N", "N")
    assert not iupac_matches("A", "N")


def test_find_primer_exact():
    assert find_primer("ACGT", "XXXXACGTXXXX", 0) == (4, 8)
    assert find_primer("ACGT", "AAAAAAAA", 0) is None


def test_find_primer_leftmost():
    assert find_primer("ACGT", "ACGTXXXXACGT", 0) == (0, 4)


def test_find_primer_one_mismatch():
    # ACGA differs from ACGT at the last base.
    assert find_primer("ACGT", "XXXXACGAXXXX", 0) is None
    assert find_primer("ACGT", "XXXXACGAXXXX", 1) == (4, 8)


def test_find_primer_two_mismatches_rejected_at_one():
    assert find_primer("ACGT", "XXXXAAGAXXXX", 1) is None
    assert find_primer("ACGT", "XXXXAAGAXXXX", 2) == (4, 8)


def test_find_primer_leftmost_within_budget():
    # GCTTGC (1 mismatch) before exact GCATGC — leftmost-within-budget wins.
    assert find_primer("GCATGC", "TTGCTTGCATGC", 1) == (2, 8)


def test_find_primer_longer_than_seq():
    assert find_primer("ACGT", "AC", 0) is None
