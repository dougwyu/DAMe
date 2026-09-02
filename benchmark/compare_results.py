#!/usr/bin/env python3
"""Compare two DAMe v3 benchmark JSON results against performance gates."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any


def relative_change(baseline_ms: float, candidate_ms: float) -> float:
    """Return a positive fraction when the candidate is faster."""
    return (baseline_ms - candidate_ms) / baseline_ms


def assert_gate(
    baseline: dict[str, object],
    candidate: dict[str, object],
    target_keys: set[str],
    min_improvement: float,
    guard_keys: set[str],
    max_regression: float,
) -> None:
    """Raise one assertion containing every missed target and guard."""
    failures = []
    for key in sorted(target_keys):
        baseline_ms = baseline[key]["median_ms"]
        candidate_ms = candidate[key]["median_ms"]
        gain = relative_change(baseline_ms, candidate_ms)
        if gain < min_improvement:
            failures.append(
                f"{key}: gain {gain:.1%} is below {min_improvement:.1%}"
            )
    for key in sorted(guard_keys):
        baseline_ms = baseline[key]["median_ms"]
        candidate_ms = candidate[key]["median_ms"]
        gain = relative_change(baseline_ms, candidate_ms)
        if gain < -max_regression:
            failures.append(
                f"{key}: regression {-gain:.1%} exceeds {max_regression:.1%}"
            )
    if failures:
        raise AssertionError("\n".join(failures))


def validate_compatible_documents(
    baseline: dict[str, Any], candidate: dict[str, Any]
) -> None:
    """Reject benchmark documents that do not describe the same work."""
    if baseline.get("schema_version") != candidate.get("schema_version"):
        raise ValueError("benchmark schema_version values differ")
    if baseline.get("schema_version") != 2:
        raise ValueError("unsupported benchmark schema_version (expected 2)")

    baseline_cases = baseline.get("cases", {})
    candidate_cases = candidate.get("cases", {})
    if set(baseline_cases) != set(candidate_cases):
        raise ValueError("benchmark case names differ")
    for name in sorted(baseline_cases):
        for field in ("units", "unit_name"):
            if baseline_cases[name].get(field) != candidate_cases[name].get(field):
                raise ValueError(f"{name}: {field} differs between results")

    if baseline.get("workloads") != candidate.get("workloads"):
        raise ValueError("workload fingerprints differ between results")
    if baseline.get("pairs") != candidate.get("pairs"):
        raise ValueError("paired output fingerprints differ between results")
    if baseline.get("environment") != candidate.get("environment"):
        raise ValueError("benchmark environments differ between results")


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("baseline", type=Path)
    parser.add_argument("candidate", type=Path)
    parser.add_argument("--target", action="append", default=[])
    parser.add_argument("--guard", action="append", default=[])
    parser.add_argument("--guard-all", action="store_true")
    parser.add_argument("--min-improvement", type=float, default=0.05)
    parser.add_argument("--max-regression", type=float, default=0.02)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    baseline_document = json.loads(args.baseline.read_text())
    candidate_document = json.loads(args.candidate.read_text())
    if not baseline_document.get("parity_ok", False):
        print(f"baseline parity failed: {args.baseline}", file=sys.stderr)
        return 1
    if not candidate_document.get("parity_ok", False):
        print(f"candidate parity failed: {args.candidate}", file=sys.stderr)
        return 1
    try:
        validate_compatible_documents(baseline_document, candidate_document)
    except ValueError as error:
        print(error, file=sys.stderr)
        return 1

    baseline = baseline_document["cases"]
    candidate = candidate_document["cases"]
    targets = set(args.target)
    guards = set(args.guard)
    unknown = (targets | guards) - (set(baseline) & set(candidate))
    if unknown:
        print(
            "unknown benchmark case: " + ", ".join(sorted(unknown)),
            file=sys.stderr,
        )
        print("available cases: " + ", ".join(sorted(baseline)), file=sys.stderr)
        return 1
    if args.guard_all:
        guards.update(set(baseline) & set(candidate) - targets)

    try:
        assert_gate(
            baseline,
            candidate,
            target_keys=targets,
            min_improvement=args.min_improvement,
            guard_keys=guards,
            max_regression=args.max_regression,
        )
    except AssertionError as error:
        print(error, file=sys.stderr)
        return 1

    for key in sorted(targets | guards):
        change = relative_change(
            baseline[key]["median_ms"], candidate[key]["median_ms"]
        )
        print(f"{key}: {change:+.1%}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
