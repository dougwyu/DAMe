import pytest
import numpy as np
import subprocess
import sys
from dame.rsi import compare


def test_importing_cli_does_not_import_numpy():
    code = "import sys; import dame.__main__; assert 'numpy' not in sys.modules"
    completed = subprocess.run([sys.executable, "-c", code], check=False)
    assert completed.returncode == 0


def test_compare_identical_replicates():
    # Identical replicates → RSI = 0 (perfectly similar)
    matrix = np.array([[10, 10], [20, 20], [30, 30]])
    result = compare(matrix, "sample1", 1, 2)
    assert abs(result) < 1e-10


def test_compare_completely_different():
    # Mutually exclusive seqs → RSI = 1 (completely dissimilar)
    matrix = np.array([[100, 0], [0, 100]])
    result = compare(matrix, "sample1", 1, 2)
    assert abs(result - 1.0) < 1e-10


def test_compare_partial_overlap():
    # Half shared → RSI = 0.5
    matrix = np.array([[50, 0], [50, 100]])
    result = compare(matrix, "sample1", 1, 2)
    # a = [0.5, 0.5], b = [0.0, 1.0], min = [0.0, 0.5], sum = 0.5, RSI = 0.5
    assert abs(result - 0.5) < 1e-10


def test_compare_zero_replicate_handled():
    # Zero-total replicate should not divide by zero
    matrix = np.array([[0, 10], [0, 20]])
    result = compare(matrix, "sample1", 1, 2)
    assert 0.0 <= result <= 1.0


# ── deterministic sample ordering ─────────────────────────────────────────────
# A Python set iterates in hash-seed-dependent order, so RSI rows used to come
# out in a different order on every run and never matched the Rust output, which
# sorts sample names.

def test_rsi_sample_order_is_sorted(tmp_path, monkeypatch):
    import argparse
    from dame import rsi as rsi_mod

    inp = tmp_path / "cmp.txt"
    inp.write_text(
        "sB\tt1-t2\t5\tt3-t4\t7\tTTTT\n"
        "sD\tt1-t2\t6\tt3-t4\t2\tCCCC\n"
        "sA\tt1-t2\t10\tt3-t4\t8\tACGT\n"
        "sC\tt1-t2\t3\tt3-t4\t9\tGGGG\n"
    )
    out = tmp_path / "RSI_output.txt"
    monkeypatch.chdir(tmp_path)
    rsi_mod.run(argparse.Namespace(
        input=str(inp), explicit=False, outfile=str(out)
    ))
    names = [ln.split("\t")[0] for ln in out.read_text().splitlines()[1:]]
    assert names == ["sA", "sB", "sC", "sD"], names
