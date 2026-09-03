# Sort performance: DAMe v1 and v3

Benchmarked 3 September 2026. DAMe v3 provides Python and Rust
implementations developed in parallel; the results below compare each with
upstream DAMe v1, then compare the two v3 implementations directly.

## Results

The shared default path processed the same 98,000-read synthetic FASTQ:

| Implementation | Median | Throughput | Comparison |
|---|---:|---:|---:|
| Upstream DAMe v1 (Python 2.7) | 335.299 ms | 292,276 reads/s | baseline |
| DAMe Python v3 | 242.219 ms | 404,592 reads/s | **1.4× v1** |
| DAMe Rust v3 | 50.552 ms | 1,938,589 reads/s | **6.6× v1; 4.8× Python v3** |

Only v3 supports mismatch-tolerant sorting:

| v3 sort mode | Python v3 | Rust v3 | Rust speedup |
|---|---:|---:|---:|
| Default / exact | 242.219 ms | 50.552 ms | **4.8×** |
| `--primer-mismatches 1` | 1,228.854 ms | 44.087 ms | **27.9×** |
| `--tag-mismatches 1` | 1,222.898 ms | 44.155 ms | **27.7×** |

## Where the speedup comes from

The Rust implementation parses FASTQ with `needletail`, performs IUPAC
matching on bytes, uses prebuilt hash maps for exact tag lookup, and avoids
many per-read allocations. Python v3 now loads NumPy only for RSI, reducing
`sort` startup time. Its tolerant matcher also compares each IUPAC window in a
single direct loop, avoiding a generator and one Python function call per DNA
letter. The mismatch gap remains larger because Python still tests candidates
in interpreted loops, while Rust keeps that fixed-offset work in compiled code.

## Method

All implementations ran in the same Linux/aarch64 Docker image. Each received
two warmups and ten measured runs, interleaved with rotating run order. Times
cover the complete CLI process, including startup, parsing, matching and output
writing. Default outputs matched across all three implementations; mismatch
outputs matched between Python v3 and Rust v3.

The v1 source was upstream commit
`f8aa1fc7b64f9c1bfaeb4d44378f9667a7e5a689`; both v3 implementations were
benchmarked from the v3 performance branch at commit
`6daf673087a1352acd895059c13fa3441738bc42`. Absolute times are specific to
this synthetic workload and benchmark host; the within-session ratios are the
useful comparison.
