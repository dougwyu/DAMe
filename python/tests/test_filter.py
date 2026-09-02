import pytest
import os
from io import StringIO
from dame.modules_filter import (
    makePSnumFiles, ReadPSnumFiles, MakeSampleNameArray,
    ReadHapsForASample, buildReplicateIndexes, allSequences, MakeComparisonFile,
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


def test_buildReplicateIndexes_preserves_raw_first_count_and_union():
    haps = {
        "0": [["CO1", "Tag1", "Tag2", "007", "AAAA"],
              ["CO1", "Tag1", "Tag2", "9", "AAAA"],
              ["CO1", "Tag1", "Tag2", "bad", "CCCC"]],
        "1": [["CO1", "Tag3", "Tag4", "2", "GGGG"]],
    }
    indexes = buildReplicateIndexes(2, haps)

    assert indexes["0"]["forward_tag"] == "Tag1"
    assert indexes["0"]["reverse_tag"] == "Tag2"
    assert indexes["0"]["counts_by_sequence"] == {"AAAA": "007", "CCCC": "bad"}
    assert allSequences(indexes) == {"AAAA", "CCCC", "GGGG"}


def comparison_outputs(haps):
    indexes = buildReplicateIndexes(1, haps)
    outputs = [StringIO() for _ in range(7)]
    MakeComparisonFile(
        1, allSequences(indexes), haps, indexes, *outputs,
        1, 1, 0, ["SampleA"], 0,
    )
    return outputs


def test_MakeComparisonFile_round_trips_noncanonical_count_text():
    haps = {"0": [["CO1", "Tag1", "Tag2", "007", "AAAA"]]}

    outputs = comparison_outputs(haps)

    assert outputs[0].getvalue() == "SampleA\tTag1-Tag2\t007\tAAAA\n"


def test_MakeComparisonFile_invalid_count_still_raises_before_writing_row():
    haps = {"0": [["CO1", "Tag1", "Tag2", "invalid", "AAAA"]]}
    indexes = buildReplicateIndexes(1, haps)
    outputs = [StringIO() for _ in range(7)]

    with pytest.raises(ValueError, match="invalid literal"):
        MakeComparisonFile(
            1, allSequences(indexes), haps, indexes, *outputs,
            1, 1, 0, ["SampleA"], 0,
        )

    assert outputs[0].getvalue() == ""


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
