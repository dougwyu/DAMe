import pytest
import numpy as np
from dame.rsi import compare


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
