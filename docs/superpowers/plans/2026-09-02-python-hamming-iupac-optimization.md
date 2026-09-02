# Python IUPAC Hamming Optimization Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Reduce Python v3 mismatch-sort time by removing per-base generator and function-call overhead while preserving every existing IUPAC and ambiguity result.

**Architecture:** Keep the public `hamming_iupac(pattern, region)` interface and its length-mismatch sentinel unchanged. Replace only its generator expression with a direct loop that reads the existing `_IUPAC` table; do not change matcher traversal, tie-breaking, tag/primer data structures, or output ownership.

**Tech Stack:** Python 3.11, pytest, cProfile, the pinned Docker v3 benchmark harness.

**Spec:** `docs/superpowers/plans/2026-09-02-dame-v3-performance-optimization.md`, Task 9. Its final profile measured 74,045,111 calls for each mismatch workload and about 7.0 of 8.2 profiled seconds in `hamming_iupac`, including 20,874,000 generator iterations.

## Global Constraints

- Preserve byte-for-byte Python/Rust output parity for exact, primer-mismatch, and tag-mismatch sort paths.
- Preserve the public `iupac_matches` and `hamming_iupac` behavior, including unknown read bases and the `len(pattern) + 1` unequal-length sentinel.
- Add no dependency and make no matcher traversal, candidate-selection, or tie-breaking change.
- Retain the change only if both Python mismatch benchmark cases improve by at least 15% and every other common case stays within the 2% regression guard.

---

### Task 1: Inline the hot IUPAC comparison loop

**Files:**
- Modify: `python/dame/modules_sort.py:55-60`
- Test: `python/tests/test_sort.py`

**Interfaces:**
- Consumes: `_IUPAC: dict[str, frozenset[str]]` and `hamming_iupac(pattern: str, region: str) -> int`.
- Produces: the same `hamming_iupac` result without calling `iupac_matches` once per base.

- [ ] **Step 1: Add exhaustive equivalence tests**

Add a parametrized test covering every supported IUPAC code against `A`, `C`, `G`, `T`, `N`, and `X`, and compare the function against the public single-base predicate:

```python
import pytest


@pytest.mark.parametrize("primer_base", "ACGTRYSWKMBDHVN")
@pytest.mark.parametrize("read_base", "ACGTNX")
def test_hamming_iupac_matches_single_base_predicate(primer_base, read_base):
    expected = 0 if iupac_matches(primer_base, read_base) else 1
    assert hamming_iupac(primer_base, read_base) == expected
```

Keep the existing unequal-length assertion. Run:

```bash
python3 -m pytest python/tests/test_sort.py -q
```

Expected: PASS before and after the implementation; this is a characterization test for an optimized internal path.

- [ ] **Step 2: Replace the generator with a direct loop**

Use the existing lookup table directly while retaining the length guard:

```python
def hamming_iupac(pattern, region):
    """Count positions where region fails pattern's IUPAC constraint.
    Returns a sentinel greater than len(pattern) when lengths differ."""
    if len(pattern) != len(region):
        return len(pattern) + 1
    mismatches = 0
    for primer_base, read_base in zip(pattern, region):
        allowed = _IUPAC.get(primer_base)
        if allowed is None or read_base not in allowed:
            mismatches += 1
    return mismatches
```

- [ ] **Step 3: Run focused correctness tests**

```bash
python3 -m pytest python/tests/test_sort.py -q
bash tests/integration/run_sort_mismatch.sh
bash tests/integration/run_sort_tag_mismatch.sh
```

Expected: all tests pass and Python/Rust canonical outputs remain identical.

### Task 2: Prove the speedup and stop if it is too small

**Files:**
- Modify only if the gate passes: `python/dame/modules_sort.py`
- Test: `benchmark/run_v3_benchmarks.sh`

**Interfaces:**
- Consumes: a fresh benchmark result from the unchanged predecessor and the Task 1 candidate.
- Produces: an accept/reject decision based on both mismatch targets and all-case guards.

- [ ] **Step 1: Run an adjacent predecessor/candidate pair**

Create a detached worktree at the Task 1 predecessor. Run the pinned benchmark in the predecessor first and candidate second, with no unrelated workload between them:

```bash
./benchmark/run_v3_benchmarks.sh python-hamming-baseline
./benchmark/run_v3_benchmarks.sh python-hamming-candidate
```

- [ ] **Step 2: Apply the two-target gate**

```bash
python3 benchmark/compare_results.py \
  benchmark/results/python-hamming-baseline.json \
  benchmark/results/python-hamming-candidate.json \
  --target python_sort_primer_mm1 \
  --target python_sort_tag_mm1 \
  --min-improvement 0.15 \
  --guard-all --max-regression 0.02
```

Expected: both mismatch cases improve by at least 15%, all canonical output hashes match, and no guarded case regresses by more than 2%. Revert the implementation and test if the gate fails after one reverse-order confirmation pair.

- [ ] **Step 3: Re-profile the retained candidate**

Run cProfile on both 98,000-read mismatch cases in the pinned candidate image. Confirm that `hamming_iupac` cumulative time falls materially and record the new call/time totals in the commit message body or review notes. Do not expand scope to candidate hoisting in this plan.

### Task 3: Verify and commit the retained change

**Files:**
- Modify: `python/dame/modules_sort.py`
- Modify: `python/tests/test_sort.py`

**Interfaces:**
- Consumes: the gate-passing implementation from Task 2.
- Produces: one independently revertible performance commit.

- [ ] **Step 1: Run final verification**

```bash
python3 -m pytest python/tests -q
cargo test --manifest-path rust/Cargo.toml
cargo build --locked --release --manifest-path rust/Cargo.toml
for script in tests/integration/run_*.sh; do bash "$script"; done
git diff --check
```

Expected: every command passes.

- [ ] **Step 2: Review the exact diff**

```bash
git diff -- python/dame/modules_sort.py python/tests/test_sort.py
git status --short
```

Expected: only the direct Hamming loop and its characterization test are changed.

- [ ] **Step 3: Commit**

```bash
git add python/dame/modules_sort.py python/tests/test_sort.py
git commit -m "perf(python): inline IUPAC hamming comparisons"
```

Stop after this commit. Any precompiled-mask representation or matcher-loop restructuring requires a new profile and a separate plan.
