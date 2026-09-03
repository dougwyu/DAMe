# Sort performance: DAMe v1 and v3

Benchmarked 3 September 2026. DAMe v3 provides Python and Rust
implementations developed in parallel; the results below compare each with
upstream DAMe v1, then compare the two v3 implementations directly.

## Results

The shared default path processed the same 98,000-read synthetic FASTQ:

| Implementation | Median | Throughput | Comparison |
|---|---:|---:|---:|
| Upstream DAMe v1 (Python 2.7) | 334.361 ms | 293,096 reads/s | baseline |
| DAMe Python v3 | 241.093 ms | 406,483 reads/s | **1.4× v1** |
| DAMe Rust v3 | 50.342 ms | 1,946,668 reads/s | **6.6× v1; 4.8× Python v3** |

Only v3 supports mismatch-tolerant sorting:

| v3 sort mode | Python v3 | Rust v3 | Rust speedup |
|---|---:|---:|---:|
| Default / exact | 241.093 ms | 50.342 ms | **4.8×** |
| `--primer-mismatches 1` | 1,844.793 ms | 43.736 ms | **42.2×** |
| `--tag-mismatches 1` | 1,853.773 ms | 43.501 ms | **42.6×** |

## Where the speedup comes from

The Rust implementation parses FASTQ with `needletail`, performs IUPAC
matching on bytes, uses prebuilt hash maps for exact tag lookup, and avoids
many per-read allocations. Python v3 now loads NumPy only for RSI, reducing
`sort` startup time. The mismatch gap remains larger because Python's tolerant
matcher performs nested loops in the interpreter, while Rust keeps the
fixed-offset matching work in compiled code.

## Method

All implementations ran in the same Linux/aarch64 Docker image. Each received
two warmups and ten measured runs, interleaved with rotating run order. Times
cover the complete CLI process, including startup, parsing, matching and output
writing. Default outputs matched across all three implementations; mismatch
outputs matched between Python v3 and Rust v3.

The v1 source was upstream commit
`f8aa1fc7b64f9c1bfaeb4d44378f9667a7e5a689`; both v3 implementations were
benchmarked from the v3 performance branch at commit
`bf91630e24525c10f1bb9ab8a39762cd067f5dfb`. Absolute times are specific to
this synthetic workload and benchmark host; the within-session ratios are the
useful comparison.
