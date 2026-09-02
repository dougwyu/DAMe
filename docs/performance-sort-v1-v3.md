# Sort performance: DAMe v1 and v3

Benchmarked 2 September 2026. DAMe v3 provides Python and Rust
implementations developed in parallel; the results below compare each with
upstream DAMe v1, then compare the two v3 implementations directly.

## Results

The shared default path processed the same 98,000-read synthetic FASTQ:

| Implementation | Median | Throughput | Comparison |
|---|---:|---:|---:|
| Upstream DAMe v1 (Python 2.7) | 328.882 ms | 297,979 reads/s | baseline |
| DAMe Python v3 | 302.139 ms | 324,354 reads/s | **1.09× v1** |
| DAMe Rust v3 | 44.836 ms | 2,185,754 reads/s | **7.34× v1; 6.74× Python v3** |

Only v3 supports mismatch-tolerant sorting:

| v3 sort mode | Python v3 | Rust v3 | Rust speedup |
|---|---:|---:|---:|
| Default / exact | 302.139 ms | 44.836 ms | **6.74×** |
| `--primer-mismatches 1` | 1,939.159 ms | 37.301 ms | **51.99×** |
| `--tag-mismatches 1` | 1,967.515 ms | 37.557 ms | **52.39×** |

## Where the speedup comes from

The Rust implementation parses FASTQ with `needletail`, performs IUPAC
matching on bytes, uses prebuilt hash maps for exact tag lookup, and avoids
many per-read allocations. The mismatch gap is larger because Python's
tolerant matcher performs nested loops in the interpreter, while Rust keeps
the fixed-offset matching work in compiled code.

## Method

All implementations ran in the same Linux/aarch64 Docker image. Each received
two warmups and ten measured runs, interleaved with rotating run order. Times
cover the complete CLI process, including startup, parsing, matching and output
writing. Default outputs matched across all three implementations; mismatch
outputs matched between Python v3 and Rust v3.

The v1 source was upstream commit
`f8aa1fc7b64f9c1bfaeb4d44378f9667a7e5a689`; both v3 implementations were
benchmarked from repository tag `v3.0.0` at commit
`5f215f67a92470b5d118703ef1a224c21b178405`. Absolute times are specific to
this synthetic workload and benchmark host; the within-session ratios are the
useful comparison.
