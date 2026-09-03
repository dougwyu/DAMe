# Filter performance: DAMe v1 and v3

Benchmarked 3 September 2026. DAMe v3 provides Python and Rust
implementations developed in parallel; the results below compare each with
upstream DAMe v1 and compare the two v3 implementations directly.

## Results

The primary workload contained 100 samples, two replicates per sample and
40,000 collapsed input records:

| Implementation | Median | Throughput | Comparison |
|---|---:|---:|---:|
| Upstream DAMe v1 (Python 2.7) | 296.564 ms | 134,878 records/s | baseline |
| DAMe Python v3 | 86.354 ms | 463,208 records/s | **3.4× v1; 1.2× Rust v3** |
| DAMe Rust v3 | 105.994 ms | 377,380 records/s | **2.8× v1** |

## Where the speedup comes from

Both v3 implementations index counts by sequence, replacing repeated linear
searches through each replicate. Python v3 also loads NumPy only when RSI is
run, which removes a large startup cost from `filter`; on this workload those
changes make Python v3 1.2× faster than Rust. Rust retains compiled loops and
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
`6daf673087a1352acd895059c13fa3441738bc42`. A tiny fixture was excluded
because CLI startup dominated it. Absolute
times are host- and workload-specific; the within-session ratios are the useful
comparison.
