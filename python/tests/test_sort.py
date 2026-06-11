import pytest
import tempfile
import os
from dame.modules_sort import RC, readTags, readPrimers, FillHAP, GetPiecesInfo, \
    iupac_match, hamming_iupac, GetPiecesInfoMismatch


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
    assert len(result["CO1"]) == 4          # [regex_fwd, regex_rc, raw_fwd, raw_rc]
    assert len(result["CO1"][0]) == 2       # F and R regex patterns
    assert result["CO1"][2][0] == "ACGT"    # F raw (unchanged, no IUPAC in this fixture)
    assert result["CO1"][2][1] == "TTTT"    # R raw
    assert result["CO1"][3][0] == RC("ACGT")  # RC(F) raw
    assert result["CO1"][3][1] == RC("TTTT")  # RC(R) raw


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


# ── shared fixtures for GetPiecesInfoMismatch tests ──────────────────────────

def _make_tags():
    return {
        "Tag1": ["AAAA", RC("AAAA")],   # RC("AAAA") = "TTTT"
        "Tag2": ["CCCC", RC("CCCC")],   # RC("CCCC") = "GGGG"
        "Tag3": ["GGGG", RC("GGGG")],   # RC("GGGG") = "CCCC"
        "Tag4": ["TTTT", RC("TTTT")],   # RC("TTTT") = "AAAA"
    }


def _make_primers():
    # F=ACGT, R=TGCA, RC(F)=RC("ACGT")="ACGT", RC(R)=RC("TGCA")="TGCA" (both palindromes)
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
    # Primers back-to-back with no amplicon between them
    read = "AAAA" + "ACGT" + "TGCA" + "GGGG"
    result = GetPiecesInfoMismatch(read, _make_primers(), _make_tags(), False, 0, 0)
    assert result == [1]


def test_no_match_returns_1():
    result = GetPiecesInfoMismatch("NNNNNNNNNNNN", _make_primers(), _make_tags(), False, 0, 0)
    assert result == [1]


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
