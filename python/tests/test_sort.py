import argparse
import pytest
import tempfile
import os
from dame.modules_sort import (
    RC, readTags, readPrimers, FillHAP, GetPiecesInfo, iupac_matches, find_primer,
    hamming_iupac, GetPiecesInfoMismatch, min_equal_length_tag_distance,
    compile_primer_regexes,
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


def test_compile_primer_regexes_expands_iupac():
    raw = {"CO1": [["GCRTGC", "CTGACT"], ["GCAYGC", "AGTCAG"]]}
    rx = compile_primer_regexes(raw)
    assert rx["CO1"][0][0] == "GC[AG]TGC"   # R -> [AG]
    assert rx["CO1"][0][1] == "CTGACT"      # no ambiguity codes


def test_GetPiecesInfo_exact_forward(tmp_path):
    # Fast exact (zero-mismatch) regex matcher on a clean forward read.
    PRIMERS, TAGS = _mismatch_fixtures(tmp_path)
    PRIMERS_RX = compile_primer_regexes(PRIMERS)
    line = "AAAAACGTATATATTGCAGGGG"
    assert GetPiecesInfo(line, PRIMERS_RX, TAGS, False) == ["Tag1", "Tag2", "CO1", "ATATAT"]
    assert GetPiecesInfo("NNNNNNNNNNNNNNNNNNNN", PRIMERS_RX, TAGS, False) == [1]


def test_GetPiecesInfo_lowercase_read(tmp_path):
    # Soft-masked / lowercase reads are uppercased before matching (exact path).
    PRIMERS, TAGS = _mismatch_fixtures(tmp_path)
    PRIMERS_RX = compile_primer_regexes(PRIMERS)
    line = "aaaaacgtatatattgcagggg"
    assert GetPiecesInfo(line, PRIMERS_RX, TAGS, False) == ["Tag1", "Tag2", "CO1", "ATATAT"]


def test_anchored_end_primer_mismatch(tmp_path):
    PRIMERS, TAGS = _mismatch_fixtures(tmp_path)
    # Forward primer ACGT exact; END primer RC(R)=TGCA miscalled as TGCT.
    line = "AAAAACGTATATATTGCTGGGG"
    assert GetPiecesInfoMismatch(line, PRIMERS, TAGS, False, 0, 0) == [1]
    assert GetPiecesInfoMismatch(line, PRIMERS, TAGS, False, 1, 0) == ["Tag1", "Tag2", "CO1", "ATATAT"]


def test_anchored_reverse_primer_mismatch(tmp_path):
    tags_file = tmp_path / "tags.txt"
    tags_file.write_text("AAAA\tTag1\nCCCC\tTag2\nGGGG\tTag3\nTTTT\tTag4\n")
    primers_file = tmp_path / "primers.txt"
    primers_file.write_text("CO1\tAAAG\tCCGT\n")
    TAGS = readTags(str(tags_file), {})
    PRIMERS = readPrimers(str(primers_file), {})
    # Reverse-orientation read; start primer R=CCGT miscalled CCGA.
    line = "CCCCCCGAGGGGGGCTTTTTTT"
    assert GetPiecesInfoMismatch(line, PRIMERS, TAGS, False, 0, 0) == [1]
    assert GetPiecesInfoMismatch(line, PRIMERS, TAGS, False, 1, 0) == ["Tag1", "Tag2", "CO1", "CCCCCC"]


def test_hamming_iupac():
    assert hamming_iupac("ACGT", "ACGT") == 0
    assert hamming_iupac("ACGT", "ACGA") == 1
    assert hamming_iupac("GCRTGC", "GCATGC") == 0   # R matches A
    assert hamming_iupac("GCRTGC", "GCCTGC") == 1
    # length mismatch -> sentinel greater than any real budget
    assert hamming_iupac("ACGT", "ACG") > 3


@pytest.mark.parametrize("primer_base", "ACGTRYSWKMBDHVN")
@pytest.mark.parametrize("read_base", "ACGTNX")
def test_hamming_iupac_matches_single_base_predicate(primer_base, read_base):
    expected = 0 if iupac_matches(primer_base, read_base) else 1
    assert hamming_iupac(primer_base, read_base) == expected


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


def _write_fastq(path, records):
    """Write records as 4-line FASTQ. Each record is (header, seq)."""
    with open(path, "w") as fh:
        for hdr, seq in records:
            fh.write("@%s\n%s\n+\n%s\n" % (hdr, seq, "I" * len(seq)))


def test_sort_empty_sequence_line_does_not_truncate(tmp_path, monkeypatch):
    """A blank sequence line is a malformed record, not end of input.

    DAMe v1.0 counted it as an error and carried on; an early `break` here
    silently discarded every read after it.
    """
    tags_file = tmp_path / "tags.txt"
    tags_file.write_text("AACCGGT\ttag1\nTTGGCCA\ttag2\n")
    primers_file = tmp_path / "primers.txt"
    primers_file.write_text("CO1\tGCRTGC\tCTGACT\n")

    good = "AACCGGT" + "GCATGC" + "TTTTTTTTTT" + "AGTCAG" + "TGGCCAA"
    fq = tmp_path / "in.fastq"
    # A blank sequence sits between two otherwise sortable reads.
    _write_fastq(fq, [("r1", good), ("r2", ""), ("r3", good)])

    monkeypatch.chdir(tmp_path)
    from dame import sort as sort_mod

    args = argparse.Namespace(
        fq=str(fq), p=str(primers_file), t=str(tags_file),
        keepPrimersSeq=False, primer_mismatches=0, tag_mismatches=0,
    )
    sort_mod.run(args)

    # r1 and r3 must both survive; only the blank record is an error.
    summary = (tmp_path / "SummaryCounts.txt").read_text().splitlines()
    assert len(summary) == 2, summary          # header + one tag pair
    assert summary[1].split("\t")[3] == "2"    # SumTotalFreq counts r1 and r3


def test_readTags_rejects_duplicate_sequence(tmp_path):
    """A repeated tag sequence is refused rather than silently resolved.

    The name a read is assigned to would otherwise depend on lookup order, and
    that name becomes an output filename.
    """
    tags_file = tmp_path / "tags.txt"
    tags_file.write_text("AAAA\tTag1\nCCCC\tTag2\nAAAA\tTagX\n")
    with pytest.raises(ValueError) as exc:
        readTags(str(tags_file), {})
    msg = str(exc.value)
    assert "duplicate tag sequence 'AAAA'" in msg
    assert "'Tag1'" in msg and "'TagX'" in msg
    assert "line 3" in msg


def test_readTags_rejects_duplicate_name(tmp_path):
    """A repeated tag name is refused: names and sequences are one to one.

    Only TAGS[name][0] is ever consulted, so a second sequence under the same
    name was unreachable here while Rust matched it on the exact path but not
    the anchored one.
    """
    tags_file = tmp_path / "tags.txt"
    tags_file.write_text("AAAA\tTag1\nCCCC\tTag1\n")
    with pytest.raises(ValueError) as exc:
        readTags(str(tags_file), {})
    msg = str(exc.value)
    assert "duplicate tag name 'Tag1'" in msg
    assert "'AAAA'" in msg and "'CCCC'" in msg
    assert "line 2" in msg


def test_readTags_accepts_a_well_formed_panel(tmp_path):
    """The checks must not fire on a normal one-to-one tags file."""
    tags_file = tmp_path / "tags.txt"
    tags_file.write_text("AAAA\tTag1\nCCCC\tTag2\nGGGG\tTag3\n")
    result = readTags(str(tags_file), {})
    assert sorted(result) == ["Tag1", "Tag2", "Tag3"]
    assert result["Tag1"][0] == "AAAA"


def test_sort_exits_nonzero_on_duplicate_tag_sequence(tmp_path, monkeypatch):
    tags_file = tmp_path / "tags.txt"
    tags_file.write_text("AAAA\tTag1\nAAAA\tTagX\n")
    primers_file = tmp_path / "primers.txt"
    primers_file.write_text("CO1\tACGT\tTGCA\n")
    fq = tmp_path / "in.fastq"
    _write_fastq(fq, [("r1", "AAAAACGTATATTGCATTTT")])

    monkeypatch.chdir(tmp_path)
    from dame import sort as sort_mod
    args = argparse.Namespace(
        fq=str(fq), p=str(primers_file), t=str(tags_file),
        keepPrimersSeq=False, primer_mismatches=0, tag_mismatches=0,
    )
    with pytest.raises(SystemExit) as exc:
        sort_mod.run(args)
    assert exc.value.code == 1
    # Refused before writing anything.
    assert not (tmp_path / "SummaryCounts.txt").exists()


def test_readTags_rejects_incomplete_entry(tmp_path):
    """A line with a sequence but no name is refused, not skipped.

    Skipping would drop the tag, and every read carrying it would land in the
    erroneous-sequence count, indistinguishable from a real tag error.
    """
    tags_file = tmp_path / "tags.txt"
    tags_file.write_text("AAAA\tTag1\nCCCC\tTag2\nGGGG\n")
    with pytest.raises(ValueError) as exc:
        readTags(str(tags_file), {})
    msg = str(exc.value)
    assert "incomplete tag entry" in msg
    assert "line 3" in msg
    assert "'GGGG'" in msg


def test_readPrimers_rejects_incomplete_entry(tmp_path):
    """A primer line missing the reverse sequence is refused, not skipped.

    Skipping would drop the primer set, leaving nothing for reads to match.
    """
    f = tmp_path / "primers.txt"
    f.write_text("CO1\tACGT\tTGCA\nCO2\tGGGG\n")
    with pytest.raises(ValueError) as exc:
        readPrimers(str(f), {})
    msg = str(exc.value)
    assert "incomplete primer entry" in msg
    assert "line 2" in msg


def test_readPrimers_rejects_duplicate_name(tmp_path):
    """A repeated primer set name is refused.

    Only PRIMERS[name][0][0:2] is consulted here so the first definition won,
    while Rust's IndexMap insert replaced the entry and the last won.
    """
    f = tmp_path / "primers.txt"
    f.write_text("CO1\tACGT\tTGCA\nCO1\tGGGG\tCCCC\n")
    with pytest.raises(ValueError) as exc:
        readPrimers(str(f), {})
    msg = str(exc.value)
    assert "duplicate primer set name 'CO1'" in msg
    assert "'ACGT/TGCA'" in msg and "'GGGG/CCCC'" in msg


def test_readPrimers_allows_shared_sequence_across_sets(tmp_path):
    """Two sets sharing a forward primer is a legitimate multiplex design."""
    f = tmp_path / "primers.txt"
    f.write_text("CO1\tACGT\tTGCA\nCO2\tACGT\tGGGG\n")
    result = readPrimers(str(f), {})
    assert sorted(result) == ["CO1", "CO2"]
