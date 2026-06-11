# Sort benchmark

Reproducible throughput benchmark for the `sort` step (the hot path), comparing
the Python (`dame-py`) and Rust (`dame`) implementations on a large synthetic
dataset.

## Prerequisites

```bash
pip install -e python/                                   # dame-py on PATH
cargo build --release --manifest-path rust/Cargo.toml    # release binary
```

## Run

```bash
bash benchmark/run_sort_benchmark.sh
```

This generates 196,000 synthetic reads (2 pools × 98k, the 8-tag CO1 tutorial
panel) in a temp dir via `generate_benchmark_data.py`, then times `sort` for
both implementations at the default, `--primer-mismatches 1`, and
`--tag-mismatches 1` settings (best of 3 each).

## Caveats

Absolute timings vary substantially with machine and background load, and this
synthetic set has high amplicon duplication, so the numbers are **not**
comparable across machines or to the historical figures in the top-level
README. The meaningful signal is the *within-session ratios* (Rust vs Python,
and default vs anchored path).
