import pytest
import tempfile
import os
from dame.modules_sort import (
    RC, readTags, readPrimers, FillHAP, GetPiecesInfo, iupac_matches, find_primer,
    hamming_iupac, GetPiecesInfoMismatch, min_equal_length_tag_distance,
)


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
    PRIMERS = {}
    result = readPrimers(str(primers_file), PRIMERS)
    assert "CO1" in result
    assert len(result["CO1"]) == 2          # [forward_list, rc_list]
    assert result["CO1"][0] == ["ACGT", "TTTT"]      # F, R (raw)
    assert result["CO1"][1] == ["ACGT", "AAAA"]      # RC(F), RC(R)


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


def _mismatch_fixtures(tmp_path):
    tags_file = tmp_path / "tags.txt"
    tags_file.write_text("AAAA\tTag1\nCCCC\tTag2\nGGGG\tTag3\nTTTT\tTag4\n")
    primers_file = tmp_path / "primers.txt"
    primers_file.write_text("CO1\tACGT\tTGCA\n")
    TAGS = readTags(str(tags_file), {})
    PRIMERS = readPrimers(str(primers_file), {})
    return PRIMERS, TAGS


def test_GetPiecesInfo_one_mismatch(tmp_path):
    PRIMERS, TAGS = _mismatch_fixtures(tmp_path)
    # Forward primer ACGT miscalled as ACGA.
    line = "AAAAACGAATATATTGCAGGGG"
    # N=0 -> rejected (error sentinel [1])
    assert GetPiecesInfo(line, PRIMERS, TAGS, False, 0) == [1]
    # N=1 -> sorts to Tag1_Tag2 with barcode ATATAT
    info = GetPiecesInfo(line, PRIMERS, TAGS, False, 1)
    assert info[0] == "Tag1"
    assert info[1] == "Tag2"
    assert info[2] == "CO1"
    assert info[3] == "ATATAT"


def test_GetPiecesInfo_end_primer_mismatch(tmp_path):
    PRIMERS, TAGS = _mismatch_fixtures(tmp_path)
    # Forward primer ACGT exact; END primer RC(R)=TGCA miscalled as TGCT.
    line = "AAAAACGTATATATTGCTGGGG"
    assert GetPiecesInfo(line, PRIMERS, TAGS, False, 0) == [1]
    info = GetPiecesInfo(line, PRIMERS, TAGS, False, 1)
    assert info[0] == "Tag1"
    assert info[1] == "Tag2"
    assert info[2] == "CO1"
    assert info[3] == "ATATAT"


def test_GetPiecesInfo_reverse_orientation_mismatch(tmp_path):
    tags_file = tmp_path / "tags.txt"
    tags_file.write_text("AAAA\tTag1\nCCCC\tTag2\nGGGG\tTag3\nTTTT\tTag4\n")
    primers_file = tmp_path / "primers.txt"
    primers_file.write_text("CO1\tAAAG\tCCGT\n")
    TAGS = readTags(str(tags_file), {})
    PRIMERS = readPrimers(str(primers_file), {})
    # Reverse-orientation read; start primer R=CCGT miscalled CCGA.
    line = "CCCCCCGAGGGGGGCTTTTTTT"
    assert GetPiecesInfo(line, PRIMERS, TAGS, False, 0) == [1]
    info = GetPiecesInfo(line, PRIMERS, TAGS, False, 1)
    assert info[0] == "Tag1"
    assert info[1] == "Tag2"
    assert info[2] == "CO1"
    assert info[3] == "CCCCCC"


def test_GetPiecesInfo_lowercase_read(tmp_path):
    PRIMERS, TAGS = _mismatch_fixtures(tmp_path)
    # Soft-masked / lowercase reads are uppercased before matching.
    line = "aaaaacgtatatattgcagggg"
    info = GetPiecesInfo(line, PRIMERS, TAGS, False, 0)
    assert info[0] == "Tag1"
    assert info[1] == "Tag2"
    assert info[2] == "CO1"
    assert info[3] == "ATATAT"


def test_hamming_iupac():
    assert hamming_iupac("ACGT", "ACGT") == 0
    assert hamming_iupac("ACGT", "ACGA") == 1
    assert hamming_iupac("GCRTGC", "GCATGC") == 0   # R matches A
    assert hamming_iupac("GCRTGC", "GCCTGC") == 1
    # length mismatch -> sentinel greater than any real budget
    assert hamming_iupac("ACGT", "ACG") > 3


def test_min_equal_length_tag_distance():
    assert min_equal_length_tag_distance({"T1": ["AAAA", "TTTT"],
                                          "T2": ["AATT", "AATT"],
                                          "T3": ["AAAT", "ATTT"]}) == 1
    assert min_equal_length_tag_distance({"T1": ["AAAA", "TTTT"]}) is None


def test_anchored_forward_tag_mismatch(tmp_path):
    PRIMERS, TAGS = _mismatch_fixtures(tmp_path)
    line = "AAAGACGTATATATTGCAGGGG"  # tag1 AAAA->AAAG
    assert GetPiecesInfoMismatch(line, PRIMERS, TAGS, False, 0, 0) == [1]
    info = GetPiecesInfoMismatch(line, PRIMERS, TAGS, False, 0, 1)
    assert info == ["Tag1", "Tag2", "CO1", "ATATAT"]


def test_anchored_primer_only_parity(tmp_path):
    PRIMERS, TAGS = _mismatch_fixtures(tmp_path)
    line = "AAAAACGAATATATTGCAGGGG"  # primer ACGT->ACGA
    assert GetPiecesInfoMismatch(line, PRIMERS, TAGS, False, 0, 0) == [1]
    info = GetPiecesInfoMismatch(line, PRIMERS, TAGS, False, 1, 0)
    assert info == ["Tag1", "Tag2", "CO1", "ATATAT"]


def test_anchored_combined(tmp_path):
    PRIMERS, TAGS = _mismatch_fixtures(tmp_path)
    line = "AAAGACGAATATATTGCAGGGG"  # tag1 and primer each 1 off
    assert GetPiecesInfoMismatch(line, PRIMERS, TAGS, False, 0, 1) == [1]
    assert GetPiecesInfoMismatch(line, PRIMERS, TAGS, False, 1, 0) == [1]
    assert GetPiecesInfoMismatch(line, PRIMERS, TAGS, False, 1, 1) == ["Tag1", "Tag2", "CO1", "ATATAT"]


def test_anchored_reverse_orientation_tag_mismatch(tmp_path):
    tags_file = tmp_path / "t.txt"
    tags_file.write_text("AAAA\tTag1\nCCCC\tTag2\nGGGG\tTag3\nTTTT\tTag4\n")
    primers_file = tmp_path / "p.txt"
    primers_file.write_text("CO1\tAAAG\tCCGT\n")
    TAGS = readTags(str(tags_file), {})
    PRIMERS = readPrimers(str(primers_file), {})
    line = "CCCCCCGAGGGGGGCTTTTTTA"  # trailing tag TTTT->TTTA
    assert GetPiecesInfoMismatch(line, PRIMERS, TAGS, False, 0, 0) == [1]
    info = GetPiecesInfoMismatch(line, PRIMERS, TAGS, False, 1, 1)
    assert info == ["Tag1", "Tag2", "CO1", "CCCCCC"]


def test_anchored_ambiguous_tag_is_discarded(tmp_path):
    tags_file = tmp_path / "t.txt"
    tags_file.write_text("AAAA\tTagA\nAATT\tTagB\nGGGG\tTagR\n")
    primers_file = tmp_path / "p.txt"
    primers_file.write_text("CO1\tACGT\tTGCA\n")
    TAGS = readTags(str(tags_file), {})
    PRIMERS = readPrimers(str(primers_file), {})
    line = "AAATACGTATATATTGCACCCC"  # AAAT is Hamming-1 from both AAAA and AATT
    assert GetPiecesInfoMismatch(line, PRIMERS, TAGS, False, 0, 1) == [1]


def test_sort_argparser_accepts_tag_mismatches():
    import argparse
    import dame.sort as sortmod
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers()
    sortmod.register_subcommand(sub)
    args = parser.parse_args(["sort", "-fq", "x", "-p", "y", "-t", "z", "-mt", "1"])
    assert args.tag_mismatches == 1
    args2 = parser.parse_args(["sort", "--fq", "x", "--primers", "y", "--tags", "z",
                               "--tag-mismatches", "2"])
    assert args2.tag_mismatches == 2
