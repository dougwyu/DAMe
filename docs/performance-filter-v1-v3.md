# Filter performance: DAMe v1 and v3

Benchmarked 2 September 2026. DAMe v3 provides Python and Rust
implementations developed in parallel; the results below compare each with
upstream DAMe v1 and compare the two v3 implementations directly.

## Results

The primary workload contained 100 samples, two replicates per sample and
40,000 collapsed input records:

| Implementation | Median | Throughput | Comparison |
|---|---:|---:|---:|
| Upstream DAMe v1 (Python 2.7) | 258.881 ms | 154,511 records/s | baseline |
| DAMe Python v3 | 306.381 ms | 130,556 records/s | v1 was **1.18× faster** |
| DAMe Rust v3 | 86.611 ms | 461,835 records/s | **2.99× v1; 3.54× Python v3** |

## Where the speedup comes from

Rust reduces per-record overhead through compiled loops, buffered I/O and
lower-cost hash and object handling. Python v3 also sorts sequence sets for
deterministic output and starts through a central dispatcher that imports all
subcommands; those costs contribute to its result relative to the direct v1
script. The benchmark measures the full command rather than an isolated inner
loop.

## Method

All implementations ran in the same Linux/aarch64 Docker image. Each received
two warmups and ten measured runs, interleaved with rotating run order. Times
include startup, input discovery and reading, filtering, and output writing.
All implementations produced semantically equivalent output.

The v1 source was upstream commit
`f8aa1fc7b64f9c1bfaeb4d44378f9667a7e5a689`; both v3 implementations were
benchmarked from repository tag `v3.0.0` at commit
`5f215f67a92470b5d118703ef1a224c21b178405`. A second fixture with only 40
collapsed records was excluded because CLI startup dominated it. Absolute
times are host- and workload-specific; the within-session ratios are the useful
comparison.
