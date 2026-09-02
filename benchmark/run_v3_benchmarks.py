#!/usr/bin/env python3
"""Run interleaved DAMe v3 benchmarks and emit machine-readable results."""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
import shutil
import statistics
import subprocess
import sys
import tempfile
import time
from collections.abc import Sequence
from dataclasses import dataclass
from datetime import UTC, datetime
from pathlib import Path


@dataclass(frozen=True)
class Case:
    name: str
    command: Sequence[str]
    cwd: Path
    units: int
    unit_name: str
    cleanup_globs: tuple[str, ...] = ()


def rotated(names: list[str], round_number: int) -> list[str]:
    """Rotate case names so a different implementation runs first."""
    offset = round_number % len(names)
    return names[offset:] + names[:offset]


def summarize(
    samples_ms: list[float], units: int, unit_name: str
) -> dict[str, object]:
    """Summarize raw wall-time samples with median throughput."""
    median_ms = statistics.median(samples_ms)
    return {
        "samples_ms": samples_ms,
        "median_ms": median_ms,
        "units": units,
        "unit_name": unit_name,
        "units_per_second": units / (median_ms / 1_000),
    }


def _clean_case_outputs(case: Case) -> None:
    for pattern in case.cleanup_globs:
        for path in case.cwd.glob(pattern):
            if path.is_file() or path.is_symlink():
                path.unlink()


def _run_case(case: Case) -> float:
    _clean_case_outputs(case)
    started = time.perf_counter_ns()
    subprocess.run(
        case.command,
        cwd=case.cwd,
        check=True,
        capture_output=True,
        text=True,
    )
    return (time.perf_counter_ns() - started) / 1_000_000


def run_interleaved(
    cases: list[Case], warmups: int, rounds: int
) -> dict[str, dict[str, object]]:
    """Run warmups, then measured rounds with rotating case order."""
    by_name = {case.name: case for case in cases}
    names = [case.name for case in cases]
    samples = {name: [] for name in names}

    for _ in range(warmups):
        for case in cases:
            _run_case(case)

    for round_number in range(rounds):
        for name in rotated(names, round_number):
            samples[name].append(_run_case(by_name[name]))

    return {
        case.name: summarize(samples[case.name], case.units, case.unit_name)
        for case in cases
    }


def canonical_hashes(root: Path, patterns: Sequence[str]) -> dict[str, str]:
    """Return SHA-256 hashes keyed by normalized relative output path."""
    selected = {
        path
        for pattern in patterns
        for path in root.glob(pattern)
        if path.is_file()
    }
    return {
        path.relative_to(root).as_posix(): hashlib.sha256(path.read_bytes()).hexdigest()
        for path in sorted(selected)
    }


def verify_pair(
    left_root: Path, right_root: Path, patterns: Sequence[str]
) -> dict[str, str]:
    """Return canonical hashes when paired outputs are byte-for-byte equal."""
    left = canonical_hashes(left_root, patterns)
    right = canonical_hashes(right_root, patterns)
    if not left and not right:
        raise RuntimeError(
            f"{left_root.name} and {right_root.name} produced no canonical outputs"
        )
    if left != right:
        missing_left = sorted(set(right) - set(left))
        missing_right = sorted(set(left) - set(right))
        changed = sorted(
            path for path in set(left) & set(right) if left[path] != right[path]
        )
        raise RuntimeError(
            f"output mismatch between {left_root.name} and {right_root.name}; "
            f"missing from left={missing_left}; missing from right={missing_right}; "
            f"different hashes={changed}"
        )
    return left


def execute_benchmark(
    cases: list[Case],
    pairs: dict[str, tuple[str, str]],
    patterns_by_pair: dict[str, Sequence[str]],
    workload_hashes: dict[str, str],
    warmups: int,
    rounds: int,
) -> dict[str, object]:
    """Time cases and reject the document unless every pair has equal output."""
    results = run_interleaved(cases, warmups=warmups, rounds=rounds)
    roots = {case.name: case.cwd for case in cases}
    pair_results = {}
    for pair_name, (left_name, right_name) in pairs.items():
        hashes = verify_pair(
            roots[left_name], roots[right_name], patterns_by_pair[pair_name]
        )
        pair_results[pair_name] = {
            "cases": [left_name, right_name],
            "outputs": sorted(hashes),
            "hashes": hashes,
        }
    return {
        "schema_version": 2,
        "parity_ok": True,
        "cases": results,
        "pairs": pair_results,
        "workloads": workload_hashes,
    }


def workload_hashes(sort_input: Path, filter_input: Path) -> dict[str, str]:
    """Fingerprint every generated input using stable fixture-relative paths."""
    return {
        **{
            f"sort/{path}": digest
            for path, digest in canonical_hashes(sort_input, ("**/*",)).items()
        },
        **{
            f"filter/{path}": digest
            for path, digest in canonical_hashes(filter_input, ("**/*",)).items()
        },
    }


SORT_OUTPUTS = ("SummaryCounts.txt", "tag*.txt")
FILTER_OUTPUTS = ("Comparisons_*.txt", "Comparisons_*.fasta", "FilteredReads*")
FILTER_CLEANUP = ("PS*_files.txt", "Comparisons_*", "FilteredReads*")


def build_suite(
    work_root: Path,
    sort_input: Path,
    filter_input: Path,
    python_bin: Path,
    rust_bin: Path,
) -> tuple[
    list[Case],
    dict[str, tuple[str, str]],
    dict[str, tuple[str, ...]],
]:
    """Create isolated case directories and all Python/Rust v3 commands."""
    work_root.mkdir(parents=True, exist_ok=True)
    cases = []
    binaries = {"python": str(python_bin), "rust": str(rust_bin)}
    sort_base = [
        "sort",
        "--fq", str(sort_input / "Pool1.fastq"),
        "--primers", str(sort_input / "Primers.txt"),
        "--tags", str(sort_input / "Tags.txt"),
    ]
    sort_variants = {
        "default": [],
        "primer_mm1": ["--primer-mismatches", "1"],
        "tag_mm1": ["--tag-mismatches", "1"],
    }
    for variant, extra_args in sort_variants.items():
        for implementation in ("python", "rust"):
            name = f"{implementation}_sort_{variant}"
            cwd = work_root / name
            cwd.mkdir()
            cases.append(
                Case(
                    name,
                    [binaries[implementation], *sort_base, *extra_args],
                    cwd,
                    98_000,
                    "reads",
                    cleanup_globs=SORT_OUTPUTS,
                )
            )

    filter_units = {"tiny": 4, "scaled": 40_000}
    for size, units in filter_units.items():
        for implementation in ("python", "rust"):
            name = f"{implementation}_filter_{size}"
            cwd = work_root / name
            shutil.copytree(filter_input / size, cwd)
            cases.append(
                Case(
                    name,
                    [
                        binaries[implementation],
                        "filter",
                        "--ps-info", "PSinfo.txt",
                        "--x", "2",
                        "--y", "1",
                        "--t", "1",
                        "--l", "100",
                    ],
                    cwd,
                    units,
                    "records",
                    cleanup_globs=FILTER_CLEANUP,
                )
            )

    pairs = {
        "sort_default": ("python_sort_default", "rust_sort_default"),
        "sort_primer_mm1": ("python_sort_primer_mm1", "rust_sort_primer_mm1"),
        "sort_tag_mm1": ("python_sort_tag_mm1", "rust_sort_tag_mm1"),
        "filter_tiny": ("python_filter_tiny", "rust_filter_tiny"),
        "filter_scaled": ("python_filter_scaled", "rust_filter_scaled"),
    }
    patterns = {
        "sort_default": SORT_OUTPUTS,
        "sort_primer_mm1": SORT_OUTPUTS,
        "sort_tag_mm1": SORT_OUTPUTS,
        "filter_tiny": FILTER_OUTPUTS,
        "filter_scaled": FILTER_OUTPUTS,
    }
    return cases, pairs, patterns


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--warmups", type=int, default=2)
    parser.add_argument("--rounds", type=int, default=10)
    parser.add_argument(
        "--repo-root", type=Path, default=Path(__file__).resolve().parents[1]
    )
    parser.add_argument("--python-bin", default="dame-py")
    parser.add_argument("--rust-bin", type=Path)
    return parser.parse_args(argv)


def _version(command: Sequence[str]) -> str:
    completed = subprocess.run(command, check=True, capture_output=True, text=True)
    return (completed.stdout or completed.stderr).strip()


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    repo_root = args.repo_root.resolve()
    python_path = shutil.which(args.python_bin)
    if python_path is None:
        raise RuntimeError(f"Python executable not found: {args.python_bin}")
    rust_bin = (args.rust_bin or repo_root / "rust/target/release/dame").resolve()
    if not rust_bin.is_file():
        raise RuntimeError(f"Rust executable not found: {rust_bin}")

    with tempfile.TemporaryDirectory(prefix="dame-v3-benchmark-") as temp:
        temp_root = Path(temp)
        sort_input = temp_root / "input/sort"
        sort_input.mkdir(parents=True)
        subprocess.run(
            [sys.executable, str(repo_root / "benchmark/generate_benchmark_data.py")],
            cwd=sort_input,
            check=True,
        )
        filter_input = temp_root / "input/filter"
        subprocess.run(
            [
                sys.executable,
                str(repo_root / "benchmark/generate_filter_benchmark_data.py"),
                str(filter_input),
            ],
            check=True,
        )
        cases, pairs, patterns = build_suite(
            temp_root / "runs",
            sort_input,
            filter_input,
            Path(python_path),
            rust_bin,
        )
        document = execute_benchmark(
            cases,
            pairs,
            patterns,
            workload_hashes(sort_input, filter_input),
            warmups=args.warmups,
            rounds=args.rounds,
        )

    document["environment"] = {
        "python_runtime": platform.python_version(),
        "numpy_version": _version(
            [sys.executable, "-c", "import numpy; print(numpy.__version__)"]
        ),
        "rustc_version": _version(["rustc", "--version"]),
    }
    document["metadata"] = {
        "generated_at": datetime.now(UTC).isoformat(),
        "platform": platform.platform(),
        "warmups": args.warmups,
        "rounds": args.rounds,
        "python_version": _version([python_path, "--version"]),
        "rust_version": _version([str(rust_bin), "--version"]),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(document, indent=2, sort_keys=True) + "\n")
    print(f"wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
