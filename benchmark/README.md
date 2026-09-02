# DAMe v3 benchmarks

The benchmark gate compares the Python and Rust v3 implementations for exact
and mismatch-enabled `sort`, plus tiny and scaled `filter` workloads. It checks
output parity before accepting timing data.

## Prerequisites

```bash
docker info
```

For local integration tests, install the Python CLI from this working tree in
editable mode:

```bash
python3 -m pip install -e ./python
```

This prevents the integration scripts from silently exercising an older
installed `dame-py`. Build the current Rust release binary before running
those scripts with `cargo build --release --manifest-path rust/Cargo.toml`.

## Run

```bash
./benchmark/run_v3_benchmarks.sh baseline
./benchmark/run_v3_benchmarks.sh candidate
python3 benchmark/compare_results.py \
  benchmark/results/baseline.json \
  benchmark/results/candidate.json \
  --guard-all --max-regression 0.02
```

The container builds both v3 implementations once. Each case receives two
warmups and ten measured runs in rotating interleaved order. Results are saved
as ignored JSON files under `benchmark/results/` and include median wall time,
all raw samples, and reads or records per second.

The sort fixture contains 98,000 reads and eight tags. The scaled filter
fixture contains 100 samples, two replicates per sample, 200 collapsed records
per replicate, and 250 union sequences per sample. The tiny filter fixture
separates command startup from sustained filtering work.

Use repeatable `--target NAME` arguments with `--min-improvement` to require a
speedup. Use repeatable `--guard NAME` arguments or `--guard-all` with
`--max-regression` to reject slowdowns. A result whose Python/Rust outputs do
not match is rejected before gates are evaluated.

## Caveats

Absolute timings vary substantially with machine and background load, and this
synthetic set has high amplicon duplication. Compare baseline and candidate
results from the same host and session; do not compare absolute values across
machines.
