import os
import shutil
from pathlib import Path

import pytest

from dame.perpcr import (
    FILTERED_READS_WARNING, USEARCH_NOTE, PerPcrError, flag_error, per_pcr_convert,
)

FIXTURE = Path(__file__).resolve().parents[2] / "tests" / "fixtures" / "perpcr"
COMP = str(FIXTURE / "Comparisons_3PCRs.fasta")
PSINFO = str(FIXTURE / "PSinfo.txt")
EXPECTED = FIXTURE / "expected"


def read(path):
    return Path(path).read_text()


def test_fixture_with_psinfo_matches_expected(tmp_path):
    fasta, info = per_pcr_convert(COMP, ps_info=PSINFO, out_dir=str(tmp_path))
    assert read(fasta) == read(EXPECTED / "FilteredReads.perpcr.fna")
    assert read(info) == read(EXPECTED / "PCRinfo.txt")
    assert sorted(os.listdir(tmp_path)) == ["FilteredReads.perpcr.fna", "PCRinfo.txt"]


def test_fixture_without_psinfo_matches_expected(tmp_path):
    fasta, info = per_pcr_convert(COMP, out_dir=str(tmp_path))
    assert read(fasta) == read(EXPECTED / "FilteredReads.perpcr.fna")
    assert read(info) == read(EXPECTED / "PCRinfo_no_psinfo.txt")


def test_one_read_in_one_pcr_becomes_one_record(tmp_path):
    fasta, _ = per_pcr_convert(COMP, out_dir=str(tmp_path))
    assert ">S1_PCR2.2;size=1;sample=S1_PCR2;\n" in read(fasta)


def test_length_filters_drop_records_without_padding(tmp_path):
    fasta, info = per_pcr_convert(COMP, ps_info=PSINFO, min_length=100, max_length=200,
                                  out_dir=str(tmp_path))
    lines = read(fasta).splitlines()
    assert len(lines) == 28                      # 14 records; the 90-bp fragment is gone
    assert all(len(s) in (112, 120) for s in lines[1::2])   # nothing padded to 200
    assert lines[-2] == ">S3_PCR3.14;size=6;sample=S3_PCR3;"
    assert "S3_PCR3\tS3\t3\tt17-t18\t2\t6\n" in read(info)


def test_warning_for_filteredreads_input(tmp_path, capsys):
    src = tmp_path / "FilteredReads_copy.fna"
    shutil.copy(COMP, src)
    out = tmp_path / "out"
    out.mkdir()
    per_pcr_convert(str(src), out_dir=str(out))
    assert FILTERED_READS_WARNING + "\n" in capsys.readouterr().err
    assert (out / "FilteredReads.perpcr.fna").exists()


def test_no_warning_for_comparisons_input(tmp_path, capsys):
    per_pcr_convert(COMP, out_dir=str(tmp_path))
    assert capsys.readouterr().err == ""


def test_usearch_note(tmp_path, capsys):
    fasta, _ = per_pcr_convert(COMP, usearch=True, max_length=200, out_dir=str(tmp_path))
    assert USEARCH_NOTE + "\n" in capsys.readouterr().err
    assert read(fasta) == read(EXPECTED / "FilteredReads.perpcr.fna")


def test_short_header_record_is_skipped(tmp_path):
    src = tmp_path / "in.fasta"
    src.write_text(">S1\nACGT\n>S1\tt1-t2.t3-t4_1\t1_0\nACGT\n")
    out = tmp_path / "out"
    out.mkdir()
    fasta, _ = per_pcr_convert(str(src), out_dir=str(out))
    assert read(fasta) == ">S1_PCR1.1;size=1;sample=S1_PCR1;\nACGT\n"


def test_flag_error():
    assert flag_error(True, None, True) == "--sample-fastas cannot be used with --per-pcr"
    assert flag_error(False, "PSinfo.txt", False) == "--ps-info requires --per-pcr"
    assert flag_error(True, "PSinfo.txt", False) is None
    assert flag_error(False, None, True) is None


ERROR_CASES = [
    ("counts_mismatch", ">S1\tt1-t2.t3-t4_1\t1_2_3\nACGT\n", None,
     "record 1 for sample S1 has 3 counts but 2 tag pairs"),
    ("non_integer", ">S1\tt1-t2.t3-t4_1\t1_x\nACGT\n", None,
     "record 1 for sample S1 has a non-integer count 'x'"),
    ("inconsistent_x", ">S1\tt1-t2.t3-t4_1\t1_0\nACGT\n>S2\tt5-t6.t7-t8.t9-t10_1\t1_0_0\nACGT\n", None,
     "record 2 for sample S2 has 3 PCRs but earlier records have 2"),
    ("inconsistent_tags", ">S1\tt1-t2.t3-t4_1\t1_1\nACGT\n>S1\tt1-t2.empty-empty_2\t1_0\nAACC\n", None,
     "sample S1 has tag pairs t1-t2.empty-empty in record 2 but t1-t2.t3-t4 in an earlier record"),
    ("sample_not_in_psinfo", ">S9\tt1-t2.t3-t4_1\t1_0\nACGT\n", "S1\tt1\tt2\t1\nS1\tt3\tt4\t1\n",
     "sample S9 is in the input but not in {psinfo}"),
    ("tag_mismatch", ">S1\tt1-t2.t3-t4_1\t1_0\nACGT\n", "S1\tt1\tt2\t1\nS1\tt3\tt9\t1\n",
     "sample S1 PCR 2: tag pair t3-t4 in the input but t3-t9 in PSinfo"),
    ("psinfo_rows", ">S1\tt1-t2.t3-t4_1\t1_0\nACGT\n", "S1\tt1\tt2\t1\n",
     "sample S1 does not have exactly one PSinfo row for each PCR 1..2"),
    ("empty_input", "", None, "no records in {input}"),
]


@pytest.mark.parametrize("name,fasta,psinfo,message", ERROR_CASES, ids=[c[0] for c in ERROR_CASES])
def test_errors_leave_no_output(tmp_path, name, fasta, psinfo, message):
    src = tmp_path / "in.fasta"
    src.write_text(fasta)
    ps = None
    if psinfo is not None:
        ps = tmp_path / "PSinfo.txt"
        ps.write_text(psinfo)
    out = tmp_path / "out"
    out.mkdir()
    with pytest.raises(PerPcrError) as exc:
        per_pcr_convert(str(src), ps_info=str(ps) if ps else None, out_dir=str(out))
    assert str(exc.value) == message.format(psinfo=ps, input=src)
    assert os.listdir(out) == []
