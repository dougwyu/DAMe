# Filter performance: DAMe v1 and v3

Benchmarked 3 September 2026. DAMe v3 provides Python and Rust
implementations developed in parallel; the results below compare each with
upstream DAMe v1 and compare the two v3 implementations directly.

## Results

The primary workload contained 100 samples, two replicates per sample and
40,000 collapsed input records:

| Implementation | Median | Throughput | Comparison |
|---|---:|---:|---:|
| Upstream DAMe v1 (Python 2.7) | 291.781 ms | 137,089 records/s | baseline |
| DAMe Python v3 | 81.895 ms | 488,428 records/s | **3.6× v1; 1.3× Rust v3** |
| DAMe Rust v3 | 104.316 ms | 383,451 records/s | **2.8× v1** |

## Where the speedup comes from

Both v3 implementations index counts by sequence, replacing repeated linear
searches through each replicate. Python v3 also loads NumPy only when RSI is
run, which removes a large startup cost from `filter`; on this workload those
changes make Python v3 1.3× faster than Rust. Rust retains compiled loops and
low per-record overhead, but full-command time also includes input discovery,
deterministic output and file writing.

## Method

All implementations ran in the same Linux/aarch64 Docker image. Each received
two warmups and ten measured runs, interleaved with rotating run order. Times
include startup, input discovery and reading, filtering, and output writing.
All implementations produced semantically equivalent output.

The v1 source was upstream commit
`f8aa1fc7b64f9c1bfaeb4d44378f9667a7e5a689`; both v3 implementations were
benchmarked from the v3 performance branch at commit
`bf91630e24525c10f1bb9ab8a39762cd067f5dfb`. A tiny fixture was excluded
because CLI startup dominated it. Absolute
times are host- and workload-specific; the within-session ratios are the useful
comparison.
