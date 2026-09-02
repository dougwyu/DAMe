import json
import subprocess
import sys
from pathlib import Path

import pytest

from benchmark.compare_results import (
    assert_gate,
    relative_change,
    validate_compatible_documents,
)
from benchmark.run_v3_benchmarks import (
    Case,
    build_suite,
    canonical_hashes,
    execute_benchmark,
    rotated,
    run_interleaved,
    summarize,
    verify_pair,
)
from benchmark.generate_filter_benchmark_data import generate_filter_fixture


def result(median_ms):
    return {
        "median_ms": median_ms,
        "units": 100,
        "unit_name": "records",
    }


def benchmark_document(*, schema_version=2, workload_digest="same"):
    return {
        "schema_version": schema_version,
        "parity_ok": True,
        "cases": {"target": result(100.0), "guard": result(100.0)},
        "pairs": {
            "pair": {
                "cases": ["target", "guard"],
                "hashes": {"result.txt": "same-output"},
            }
        },
        "workloads": {"fixture.txt": workload_digest},
        "environment": {
            "python_runtime": "3.11.2",
            "numpy_version": "2.4.6",
            "rustc_version": "rustc 1.85.1",
        },
    }


def test_rotated_changes_first_case_each_round():
    names = ["python", "rust", "reference"]

    assert rotated(names, 0) == ["python", "rust", "reference"]
    assert rotated(names, 1) == ["rust", "reference", "python"]
    assert rotated(names, 3) == ["python", "rust", "reference"]


def test_summarize_reports_median_and_throughput():
    assert summarize([30.0, 10.0, 20.0], 100, "reads") == {
        "samples_ms": [30.0, 10.0, 20.0],
        "median_ms": 20.0,
        "units": 100,
        "unit_name": "reads",
        "units_per_second": 5_000.0,
    }


def test_relative_change_is_positive_for_a_speedup():
    assert relative_change(100.0, 75.0) == pytest.approx(0.25)


def test_assert_gate_accepts_target_gain_and_guard_noise():
    baseline = {"target": result(100.0), "guard": result(100.0)}
    candidate = {"target": result(90.0), "guard": result(101.0)}

    assert_gate(
        baseline,
        candidate,
        target_keys={"target"},
        min_improvement=0.05,
        guard_keys={"guard"},
        max_regression=0.02,
    )


def test_assert_gate_rejects_a_missed_target():
    baseline = {"target": result(100.0)}
    candidate = {"target": result(97.0)}

    with pytest.raises(AssertionError, match="target: gain 3.0% is below 5.0%"):
        assert_gate(
            baseline,
            candidate,
            target_keys={"target"},
            min_improvement=0.05,
            guard_keys=set(),
            max_regression=0.02,
        )


def test_assert_gate_rejects_a_guard_regression():
    baseline = {"guard": result(100.0)}
    candidate = {"guard": result(103.0)}

    with pytest.raises(AssertionError, match="guard: regression 3.0% exceeds 2.0%"):
        assert_gate(
            baseline,
            candidate,
            target_keys=set(),
            min_improvement=0.05,
            guard_keys={"guard"},
            max_regression=0.02,
        )


def test_compatible_documents_reject_changed_schema():
    baseline = benchmark_document(schema_version=2)
    candidate = benchmark_document(schema_version=3)

    with pytest.raises(ValueError, match="schema_version"):
        validate_compatible_documents(baseline, candidate)


def test_compatible_documents_reject_changed_case_units():
    baseline = benchmark_document()
    candidate = benchmark_document()
    candidate["cases"]["target"]["units"] = 99

    with pytest.raises(ValueError, match="target.*units"):
        validate_compatible_documents(baseline, candidate)


def test_compatible_documents_reject_changed_workload_fingerprint():
    baseline = benchmark_document(workload_digest="baseline")
    candidate = benchmark_document(workload_digest="candidate")

    with pytest.raises(ValueError, match="workload fingerprints"):
        validate_compatible_documents(baseline, candidate)


def test_compatible_documents_reject_changed_outputs():
    baseline = benchmark_document()
    candidate = benchmark_document()
    candidate["pairs"]["pair"]["hashes"]["result.txt"] = "changed-output"

    with pytest.raises(ValueError, match="paired output fingerprints"):
        validate_compatible_documents(baseline, candidate)


def test_run_interleaved_rotates_measured_case_order(tmp_path):
    log = tmp_path / "order.log"
    script = (
        "from pathlib import Path; import sys; "
        "Path(sys.argv[1]).open('a').write(sys.argv[2] + '\\n')"
    )
    cases = [
        Case("python", [sys.executable, "-c", script, str(log), "python"], tmp_path, 5, "reads"),
        Case("rust", [sys.executable, "-c", script, str(log), "rust"], tmp_path, 5, "reads"),
    ]

    results = run_interleaved(cases, warmups=1, rounds=3)

    assert log.read_text().splitlines() == [
        "python", "rust",  # warmup
        "python", "rust",  # round 0
        "rust", "python",  # round 1
        "python", "rust",  # round 2
    ]
    assert set(results) == {"python", "rust"}
    assert len(results["python"]["samples_ms"]) == 3
    assert results["python"]["units"] == 5
    assert results["python"]["unit_name"] == "reads"


def test_run_interleaved_cleans_declared_outputs_before_each_run(tmp_path):
    script = (
        "from pathlib import Path; "
        "p = Path('result.txt'); "
        "assert not p.exists(); "
        "p.write_text('same\\n')"
    )
    case = Case(
        "only",
        [sys.executable, "-c", script],
        tmp_path,
        1,
        "record",
        cleanup_globs=("result.txt",),
    )

    run_interleaved([case], warmups=1, rounds=2)

    assert (tmp_path / "result.txt").read_text() == "same\n"


def test_canonical_hashes_selects_relative_output_paths(tmp_path):
    (tmp_path / "SummaryCounts.txt").write_text("summary\n")
    (tmp_path / "nested").mkdir()
    (tmp_path / "nested" / "counts.tsv").write_text("counts\n")
    (tmp_path / "ignored.log").write_text("ignored\n")

    hashes = canonical_hashes(tmp_path, ["*.txt", "**/*.tsv"])

    assert set(hashes) == {"SummaryCounts.txt", "nested/counts.tsv"}
    assert all(len(digest) == 64 for digest in hashes.values())


def test_verify_pair_rejects_different_output_content(tmp_path):
    python_root = tmp_path / "python"
    rust_root = tmp_path / "rust"
    python_root.mkdir()
    rust_root.mkdir()
    (python_root / "SummaryCounts.txt").write_text("python\n")
    (rust_root / "SummaryCounts.txt").write_text("rust\n")

    with pytest.raises(RuntimeError, match="output mismatch"):
        verify_pair(python_root, rust_root, ["SummaryCounts.txt"])


def test_verify_pair_rejects_when_both_sides_produce_no_outputs(tmp_path):
    python_root = tmp_path / "python"
    rust_root = tmp_path / "rust"
    python_root.mkdir()
    rust_root.mkdir()

    with pytest.raises(RuntimeError, match="produced no canonical outputs"):
        verify_pair(python_root, rust_root, ["SummaryCounts.txt"])


def test_execute_benchmark_records_successful_pair_parity(tmp_path):
    cases = []
    for name in ("python", "rust"):
        case_root = tmp_path / name
        case_root.mkdir()
        cases.append(
            Case(
                name,
                [
                    sys.executable,
                    "-c",
                    "from pathlib import Path; Path('result.txt').write_text('same\\n')",
                ],
                case_root,
                2,
                "records",
                cleanup_globs=("result.txt",),
            )
        )

    document = execute_benchmark(
        cases,
        pairs={"filter_tiny": ("python", "rust")},
        patterns_by_pair={"filter_tiny": ("result.txt",)},
        workload_hashes={"fixture.txt": "fixture-hash"},
        warmups=0,
        rounds=1,
    )

    assert document["parity_ok"] is True
    assert set(document["cases"]) == {"python", "rust"}
    assert document["pairs"]["filter_tiny"]["outputs"] == ["result.txt"]
    assert document["workloads"] == {"fixture.txt": "fixture-hash"}


def test_generate_filter_fixture_has_requested_records_and_union(tmp_path):
    records = generate_filter_fixture(
        tmp_path,
        samples=1,
        records_per_replicate=3,
        union_per_sample=4,
    )

    assert records == 6
    assert len((tmp_path / "PSinfo.txt").read_text().splitlines()) == 2
    hap_files = sorted((tmp_path / "pool1").glob("*.txt"))
    assert len(hap_files) == 2
    assert [len(path.read_text().splitlines()) for path in hap_files] == [3, 3]
    sequences = {
        line.split()[4]
        for path in hap_files
        for line in path.read_text().splitlines()
    }
    assert len(sequences) == 4


def test_build_suite_defines_all_v3_cases_and_copies_filter_inputs(tmp_path):
    sort_input = tmp_path / "sort-input"
    sort_input.mkdir()
    for name in ("Pool1.fastq", "Primers.txt", "Tags.txt"):
        (sort_input / name).write_text("fixture\n")
    filter_input = tmp_path / "filter-input"
    for size in ("tiny", "scaled"):
        fixture = filter_input / size
        (fixture / "pool1").mkdir(parents=True)
        (fixture / "PSinfo.txt").write_text("Sample\tF\tR\t1\n")
        (fixture / "pool1" / "F_R.txt").write_text("CO1\tF\tR\t1\tAAAA\n")

    cases, pairs, patterns = build_suite(
        tmp_path / "runs",
        sort_input,
        filter_input,
        python_bin=Path("/tools/dame-py"),
        rust_bin=Path("/tools/dame"),
    )

    assert {case.name for case in cases} == {
        "python_sort_default", "rust_sort_default",
        "python_sort_primer_mm1", "rust_sort_primer_mm1",
        "python_sort_tag_mm1", "rust_sort_tag_mm1",
        "python_filter_tiny", "rust_filter_tiny",
        "python_filter_scaled", "rust_filter_scaled",
    }
    primer_case = next(case for case in cases if case.name == "python_sort_primer_mm1")
    assert primer_case.command[-2:] == ["--primer-mismatches", "1"]
    scaled_case = next(case for case in cases if case.name == "rust_filter_scaled")
    assert scaled_case.units == 40_000
    assert (scaled_case.cwd / "PSinfo.txt").is_file()
    assert set(pairs) == {"sort_default", "sort_primer_mm1", "sort_tag_mm1", "filter_tiny", "filter_scaled"}
    assert patterns["sort_default"] == ("SummaryCounts.txt", "tag*.txt")


def test_compare_cli_guard_all_checks_non_targets(tmp_path):
    baseline = benchmark_document()
    candidate = benchmark_document()
    candidate["cases"]["target"]["median_ms"] = 90.0
    candidate["cases"]["guard"]["median_ms"] = 103.0
    baseline_path = tmp_path / "baseline.json"
    candidate_path = tmp_path / "candidate.json"
    baseline_path.write_text(json.dumps(baseline))
    candidate_path.write_text(json.dumps(candidate))
    script = Path(__file__).parents[1] / "compare_results.py"

    completed = subprocess.run(
        [
            sys.executable,
            str(script),
            str(baseline_path),
            str(candidate_path),
            "--target", "target",
            "--min-improvement", "0.05",
            "--guard-all",
            "--max-regression", "0.02",
        ],
        check=False,
        capture_output=True,
        text=True,
    )

    assert completed.returncode == 1
    assert "guard: regression 3.0% exceeds 2.0%" in completed.stderr


def test_compare_cli_reports_unknown_case_without_traceback(tmp_path):
    baseline_path = tmp_path / "baseline.json"
    candidate_path = tmp_path / "candidate.json"
    baseline_path.write_text(json.dumps(benchmark_document()))
    candidate_path.write_text(json.dumps(benchmark_document()))
    script = Path(__file__).parents[1] / "compare_results.py"

    completed = subprocess.run(
        [
            sys.executable,
            str(script),
            str(baseline_path),
            str(candidate_path),
            "--target", "missing",
        ],
        check=False,
        capture_output=True,
        text=True,
    )

    assert completed.returncode == 1
    assert "unknown benchmark case: missing" in completed.stderr
    assert "Traceback" not in completed.stderr


def test_benchmark_wrapper_rejects_unsafe_result_label():
    script = Path(__file__).parents[1] / "run_v3_benchmarks.sh"

    completed = subprocess.run(
        ["bash", str(script), "../escape"],
        check=False,
        capture_output=True,
        text=True,
    )

    assert completed.returncode == 2
    assert "invalid benchmark label" in completed.stderr
