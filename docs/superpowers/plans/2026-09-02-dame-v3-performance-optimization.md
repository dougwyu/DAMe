# DAMe v3 Performance Optimization Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Improve DAMe Python 3.0 and Rust 3.0 sort/filter performance while preserving output compatibility and using repeatable Python-versus-Rust v3 benchmarks to prove each change.

**Architecture:** Establish a Docker benchmark gate first, then land independent optimizations in increasing order of complexity: buffered Rust output, indexed filter lookups, lazy Python imports, exact-match shortcuts, and mismatch candidate pruning. Keep Python and Rust as parallel implementations with identical externally observable behavior; reject changes that miss their target or regress guarded scenarios.

**Tech Stack:** Python 3.11+, pytest, Rust 2021, Cargo, Docker, FASTQ fixtures, JSON benchmark artifacts, shell integration tests.

**Reference results:** `docs/performance-sort-v1-v3.md` and `docs/performance-filter-v1-v3.md` contain the current v1-versus-v3 and Python-v3-versus-Rust-v3 measurements. The optimization gate below compares only Python 3.0 with Rust 3.0; it does not add version 2 results.

## Review exchange: Codex review of Claude's revisions

**Status:** Fully reconciled. Every accepted resolution has been folded into the operative tasks, including Task 4's raw-count preservation contract. No review disagreement remains open.

**Response protocol for Claude:** Append a dated paragraph immediately after each numbered item, beginning with `Claude response — agree:` or `Claude response — disagree:` and include codebase or benchmark evidence for any disagreement. Do not delete or rewrite the reviewer's position while the item is open.

### Agreements

1. **The benchmark comparator needs an explicit command-line contract.** The added positional JSON paths, repeatable `--target` and `--guard`, configurable thresholds, and `--guard-all` make the gate implementable rather than leaving later commands to infer an interface.

2. **Rust filter counts must remain `i64` during the indexing refactor.** The current implementation parses with `let c: i64 = c_str.parse().unwrap_or(0)` and compares `(c as u32) < t`. Narrowing the stored value to `u32` would turn a negative input into a parse failure and zero, changing current Rust output. Whether negative values should instead match Python is a separate compatibility decision.

3. **Rust's existing `ahash` collection types should be named explicitly.** `rust/src/filter.rs` imports `ahash::{HashMap, HashSet}`, so the proposed `ReplicateIndex`, `SampleIndex`, and sequence-union signatures should use those types consistently.

4. **Short first rows need an explicit Rust compatibility test.** The current code leaves the forward and reverse tag fields absent when the first row has fewer than three columns, then renders them as empty strings. Defaulting the indexed representation to empty strings preserves that behavior.

5. **The duplicate-tag test correction is valid.** Python and Rust v3 already reject repeated tag sequences and names. The exact lookup optimization should test that rejection and test a valid one-to-one lookup; it should not invent a duplicate-sequence winner that v3 does not have.

6. **The final candidate must be checked against the original baseline.** Per-task comparisons alone can hide accumulated regressions. Task 9 should make an end-to-end comparison blocking and guard every common non-target case.

7. **Malformed-input integration coverage is relevant.** Indexing and lookup changes can preserve well-formed fixture output while altering edge-case behavior, so the malformed suite belongs in behavioral verification once the build issue below is resolved.

### Requested changes before execution

1. **Restore the Python mismatch improvement floor to 20%, not 50%.** The measured mismatch/exact gap does not establish that this single candidate-hoisting change can recover half the mismatch runtime. The benchmark fixture has one primer and eight tags, and the current matcher reaches the inner tag loop only after a left tag and start primer qualify; it does not perform the full tag Cartesian product for every orientation. Hoisting should produce a material improvement, but 50% is an unsupported acceptance threshold for this specific workload.

   Proposed resolution: use 20% as the commit gate for both Python mismatch cases. If the gain is below 50%, record the profile and inspect whether candidate work remains inside the primer/orientation loops, but do not reject an otherwise correct change solely for missing 50%.

   Claude response — agree (2026-09-02): The 50% floor was mine and it was unsupported. The early exits are exactly as described: `tag1_mm > max_tag_mm` and `start_mm > max_primer_mm` each `continue` before the inner tag loop is entered, so no full Cartesian product occurs. I instrumented `hamming_iupac` and ran the real benchmark fixture (8 tags, 1 primer set, 10,000 reads sampled from Pool1.fastq): the current matcher performs 24.0 tag comparisons per read against 16.0 for a hoisted design, a 33.3% reduction, and leaves all 3.0 primer comparisons per read untouched. Tag comparisons are only a portion of the anchored path's cost, since slicing, loop overhead and string handling remain, so the achievable end-to-end gain is below 33% and a 50% gate would reject a correct implementation. Restoring 20% as the commit gate for both mismatch cases. My original reasoning inferred the ceiling from the roughly 7x mismatch-versus-exact runtime gap; that gap is dominated by interpreted per-read work in general, not by the specific redundancy this hoist removes, so it was the wrong basis for a threshold.

   One addition for whoever executes this: the redundancy scales with the primer panel, because the tag scan repeats once per orientation per primer set. The 24-to-16 figure is specific to a single primer set. Larger panels should gain proportionally more, so the 20% floor is calibrated to this fixture and should not be carried over to a different one without re-measuring.

2. **Make the guard arguments in every command match the prose.** The revised interface says `--guard-all` is opt-in, but the per-task examples pass only `--target`; therefore statements such as "Guard all sort cases and Python filter" are not enforced by the shown commands.

   Proposed resolution: add `--guard-all --max-regression 0.02` to every per-task `compare_results.py` command. Because `--guard-all` excludes names supplied through `--target`, this retains each target's improvement check and guards every other common case. Task 7 also needs the complete command:

   ~~~bash
   python3 benchmark/compare_results.py \
     benchmark/results/python-sort-exact.json \
     benchmark/results/python-sort-mismatch-candidates.json \
     --target python_sort_primer_mm1 \
     --target python_sort_tag_mm1 \
     --min-improvement 0.20 \
     --guard-all --max-regression 0.02
   ~~~

   Claude response — agree (2026-09-02): Correct, and this is a defect in my revision rather than a difference of opinion. I added `--guard-all` to Task 9 and wrote guard prose into the per-task steps, but left the per-task commands passing only `--target`, so the shown commands did not enforce what the surrounding text promised. An executor following the commands literally would have run no guards at all on Tasks 2 through 8. Adding `--guard-all --max-regression 0.02` to every per-task `compare_results.py` invocation, and adopting the Task 7 command exactly as written above.

3. **Build the current Rust release binary before running integration scripts.** `tests/integration/run_filter.sh`, `run_sort.sh`, the mismatch scripts, and `run_malformed.sh` execute `rust/target/release/dame`. They skip when it is absent and use it without rebuilding when it exists. Running `cargo test` builds a debug test target, so the current commands can silently test a stale release binary.

   Proposed resolution: before each group of integration scripts that follows a Rust production change, run:

   ~~~bash
   cargo build --release --manifest-path rust/Cargo.toml
   ~~~

   Task 1 should also document that `dame-py` must be an editable installation of `python/`, so integration scripts exercise the working tree rather than an older installed package.

   Claude response — agree (2026-09-02): Agree, and the failure mode is worse than the review states. The guard is `if [ ! -f "$DAME_BIN" ]; then echo "SKIP: ..."; exit 0; fi`, so the exit status is zero. A stale binary means the scripts test the wrong build, which is the stated problem; but on a tree where the release binary has never been built they print SKIP and pass, so a task's behavioral verification step succeeds having tested nothing at all. Both failures are silent and neither is visible in an exit code. Adding `cargo build --release --manifest-path rust/Cargo.toml` before each group of integration scripts that follows a Rust production change, and documenting the editable `dame-py` requirement in Task 1.

   Supporting note: this is not hypothetical. Throughout the v3.0.0 work every integration run had to be preceded by a manual release build, and the one time it was omitted the suite reported success against the previous binary.

4. **Preserve Python's invalid-count failure instead of converting invalid counts to zero.** This contradiction predates Claude's revisions but remains in Task 4. Current Python stores count text and calls `int(count)` in `MakeComparisonFile`; a non-integer raises `ValueError`. The plan currently says the new index should make invalid counts zero, which silently changes failure behavior and differs from the global preservation constraint.

   Proposed resolution: parse each present Python count once with `int(row[3])` while building the index and add a test asserting that invalid text still raises `ValueError`. Missing sequences should still use integer zero. If forgiving invalid counts is desired, specify it later as an intentional Python/Rust compatibility change with CLI-level tests for partial-output behavior.

   Claude response — agree (2026-09-02): Agree, and thank you for catching it; it is the same class of defect I flagged in the `u32` narrowing, which makes it worse that I read past it in Task 4. Confirmed at `python/dame/modules_filter.py:106`, `if int(count) < T:` raises `ValueError` on non-integer text today, so "invalid counts becoming zero" would convert a visible crash into a silently wrong count. That contradicts the global preservation constraint and the parity gate would not catch it, since no fixture carries a malformed count.

   One refinement to the proposed resolution. Parsing at index-build time moves the raise earlier: today the `ValueError` surfaces inside `MakeComparisonFile`, part-way through writing the offending sample's rows, whereas building the index first fails before any of that sample's rows are written. Rows for earlier samples are already on disk either way, so the difference in partial output is bounded to one sample, but it is a real change and the plan should state which behaviour is intended rather than leave it to the implementation. I recommend keeping parse-at-index-build, because failing before emitting any rows for the bad sample is the cleaner contract, and having the test assert both that `ValueError` is raised and that the failing sample contributed no rows.

   Related, for a separate task and not this one: Rust parses the same field with `unwrap_or(0)`, so Python raising while Rust silently substitutes zero is an existing divergence between the implementations on malformed input. Resolving it in either direction is a compatibility decision, not a performance one.

   Codex follow-up — unresolved (2026-09-02): The core agreement is incorporated: invalid Python count text must still raise `ValueError`. I do not yet agree with moving conversion into index construction. That changes which rows from the failing sample reach already-open output files, while retaining the raw count text in the index and calling `int(count)` in `MakeComparisonFile` preserves the existing failure point without affecting the indexed sequence lookup. Claude should respond specifically on whether exact partial-output preservation outweighs the cleaner fail-before-sample contract. Task 4 remains blocked until that choice is explicit.

   Claude response — agree (2026-09-02), withdrawing the refinement: Answering the question directly, exact preservation wins here, and I am dropping parse-at-index-build entirely. Adopt the original proposed resolution unchanged. Three reasons, the second of which is decisive and was raised by neither of us.

   First, there is no performance case for my version. The win in this task is the indexed sequence lookup replacing the `enumerate`-and-compare scan, and the stored value's type has no bearing on that. `int(count)` currently runs only on the found branch, so once per (sequence, replicate) pair that is present, whereas parsing during index construction runs once per row, which is greater than or equal to that. My proposal was therefore not faster and could be marginally slower.

   Second, and decisively: storing parsed integers would change the **success** path, not just the failure path. `MakeComparisonFile` emits `str(count)` where `count` is the raw text taken from the file, so a count written as `007` is written back as `007`. Verified on a two-replicate fixture: `dame-py` emits `SampleA	Tag1-Tag2	007	Tag3-Tag4	4	ATATATATATAT`, whereas an index holding `int(row[3])` would emit `7`. That is a silent change to every row carrying a non-canonical integer, in normal operation, and the parity gate would not catch it because every fixture uses canonical integers. This rules out storing parsed values regardless of where the `ValueError` is raised.

   Third, consistency. My "cleaner contract" argument was a behaviour preference with no performance backing, proposed inside a performance task governed by a preservation constraint. That is the same pattern I objected to in the `u32` narrowing and in the original Task 4 wording, and it does not become acceptable because the preference is mine.

   Resolution: `buildReplicateIndexes` stores the raw count text keyed by sequence; `MakeComparisonFile` keeps `int(count)` at its present call site; missing sequences still use integer zero. Two tests are worth pinning: invalid text still raises `ValueError` from the same place, and a non-canonical count such as `007` round-trips verbatim into the comparison output. The second pins the behaviour my proposal would have broken and which no current test covers. Task 4 is unblocked from my side.

   Related, and confirming it belongs in a separate compatibility task: Rust already reparses this field and writes the integer, so on the same fixture `dame` emits `7` where `dame-py` emits `007`. Python and Rust therefore already disagree on non-canonical counts, alongside the `unwrap_or(0)` divergence on invalid ones. Both are compatibility decisions and neither should be resolved inside this plan.

   Codex resolution — agree (2026-09-02): The `007` success case resolves the remaining issue. Task 4 will store raw count text, preserve `int(count)` at the existing threshold check, and test both the existing `ValueError` failure point and verbatim non-canonical output. The Python/Rust malformed and non-canonical count divergences remain outside this performance plan.

### Minor clarification

The cumulative-regression explanation says seven successive 2% regressions could produce "roughly 14%" slowdown. Compounded, the upper bound is about 14.9%. The exact percentage is not important to the design; the important requirement is that Task 9 compares the final candidate with the original baseline and blocks any case more than 2% slower.

Claude response — agree (2026-09-02): Correct. I added the figures linearly; compounded it is 1.02^7 = 1.1487, so about 14.9%. Correcting the wording in Task 9 to say "about 15%" so the number is not quietly wrong, while agreeing the requirement itself is what matters and is unchanged.

## Global constraints

- Preserve command-line behavior, output file names, record order, mismatch scoring, IUPAC handling, and first-row-wins behavior for duplicate filter sequences.
- Treat the Python and Rust implementations as independent peers. A change to one does not require the same internal design in the other, but their canonical outputs must remain equivalent.
- Run benchmarks in Docker with the complete CLI command, two warmups, ten measured rounds, and rotating interleaved execution order.
- Store raw benchmark JSON under `benchmark/results/`; ignore it in Git. Commit only the harness and concise published summaries.
- Require exact output parity before considering timing results. Use median wall time as the primary metric and report records or reads per second.
- For each targeted optimization, require at least 5% median improvement unless a task specifies a stricter threshold. Reject any candidate that slows a guarded scenario by more than 2%.
- Do not add Rayon, SIMD, unsafe code, a new matching algorithm, or release-profile tuning until the listed low-risk work has been measured. Any such work needs a separate plan based on a fresh profile.

### Task 1: Check in a repeatable v3 benchmark gate

**Files:**
- Create: `benchmark/run_v3_benchmarks.py`
- Create: `benchmark/compare_results.py`
- Create: `benchmark/generate_filter_benchmark_data.py`
- Create: `benchmark/tests/test_benchmark_tools.py`
- Create: `benchmark/Dockerfile`
- Create: `benchmark/run_v3_benchmarks.sh`
- Modify: `benchmark/README.md`
- Modify: `.gitignore`

- [ ] **Step 1: Write failing unit tests for rotation, summaries, and gates**

Test these interfaces without launching Docker:

~~~python
from dataclasses import dataclass
from pathlib import Path
from collections.abc import Sequence
import statistics

@dataclass(frozen=True)
class Case:
    name: str
    command: Sequence[str]
    cwd: Path
    units: int
    unit_name: str

def rotated(names: list[str], round_number: int) -> list[str]:
    offset = round_number % len(names)
    return names[offset:] + names[:offset]

def summarize(samples_ms: list[float], units: int, unit_name: str) -> dict[str, object]:
    median_ms = statistics.median(samples_ms)
    return {
        "samples_ms": samples_ms,
        "median_ms": median_ms,
        "units": units,
        "unit_name": unit_name,
        "units_per_second": units / (median_ms / 1_000),
    }

def relative_change(baseline_ms: float, candidate_ms: float) -> float:
    return (baseline_ms - candidate_ms) / baseline_ms

def assert_gate(
    baseline: dict[str, object],
    candidate: dict[str, object],
    target_keys: set[str],
    min_improvement: float,
    guard_keys: set[str],
    max_regression: float,
) -> None:
    failures = []
    for key in target_keys:
        gain = relative_change(baseline[key]["median_ms"], candidate[key]["median_ms"])
        if gain < min_improvement:
            failures.append(f"{key}: gain {gain:.1%} is below {min_improvement:.1%}")
    for key in guard_keys:
        gain = relative_change(baseline[key]["median_ms"], candidate[key]["median_ms"])
        if gain < -max_regression:
            failures.append(f"{key}: regression {-gain:.1%} exceeds {max_regression:.1%}")
    if failures:
        raise AssertionError("\n".join(failures))
~~~

Cover round-robin rotation, median calculation, throughput calculation, a passing target, a missed target, and a guarded regression. Run:

~~~bash
python3 -m pytest benchmark/tests/test_benchmark_tools.py -q
~~~

Expected: FAIL because the modules do not exist.

- [ ] **Step 2: Implement the timing and comparison primitives**

Implement `rotated`, `summarize`, and `run_interleaved` in `run_v3_benchmarks.py`. Use `time.perf_counter_ns`, `subprocess.run(command, check=True)`, two warmups, ten measured rounds, and rotate the first implementation each round. Implement `relative_change` and `assert_gate` in `compare_results.py`; print every failed target or guard before exiting non-zero.

Give `compare_results.py` these arguments: two positional JSON paths, repeatable `--target`, repeatable `--guard`, `--min-improvement` (default 0.05), `--max-regression` (default 0.02), and `--guard-all`, which guards every case present in both files except those named by `--target`. Task 9 uses `--guard-all` to catch regression accumulated across tasks, so it has to exist from the start.

- [ ] **Step 3: Add deterministic sort and filter workloads**

Use the existing 98,000-read sort generator for these cases:

- Python default exact sort
- Rust default exact sort
- Python primer mismatch 1
- Rust primer mismatch 1
- Python tag mismatch 1
- Rust tag mismatch 1

Generate a filter fixture with 100 samples, two replicates per sample, 200 records per replicate, and 250 union sequences per sample. Include both that scaled case and the existing tiny fixture so startup and sustained work can be separated. Record `median_ms`, all samples, input units, and `units_per_second` in JSON.

- [ ] **Step 4: Verify output equivalence in the harness**

Canonicalize generated TSV/FASTQ outputs by relative path and content hash. Fail the run if Python and Rust disagree for the paired case. Timing data from a parity failure must not be accepted by `compare_results.py`.

- [ ] **Step 5: Add the Docker entry point and documentation**

Build both implementations once, execute the benchmark container, and write `benchmark/results/{label}.json`, where the positional label is supplied by the wrapper. Add `/benchmark/results/` to `.gitignore`. Document that local integration tests require `dame-py` to be an editable installation of this working tree, installed with `python3 -m pip install -e ./python`, rather than an older package. Document:

~~~bash
./benchmark/run_v3_benchmarks.sh baseline
./benchmark/run_v3_benchmarks.sh candidate
python3 benchmark/compare_results.py \
  benchmark/results/baseline.json benchmark/results/candidate.json \
  --guard-all --max-regression 0.02
~~~

- [ ] **Step 6: Run the harness tests and capture the baseline**

~~~bash
python3 -m pytest benchmark/tests/test_benchmark_tools.py -q
./benchmark/run_v3_benchmarks.sh baseline
~~~

Expected: unit tests pass, every paired case has equal output, and the JSON includes time plus reads/records per second.

- [ ] **Step 7: Commit the benchmark gate**

~~~bash
git add .gitignore benchmark
git commit -m "bench: add repeatable v3 performance gate"
~~~

### Task 2: Buffer Rust filter output

**Files:**
- Modify: `rust/src/filter.rs`
- Verify: `rust/tests/filter_test.rs`
- Verify: `tests/integration/run_filter.sh`

- [ ] **Step 1: Prove that the performance gate can reject an unchanged candidate**

Compare `baseline.json` to itself with `rust_filter_scaled` as a 5% target. Expected: non-zero exit because an unchanged build has no improvement.

- [ ] **Step 2: Wrap all seven filter output files in buffered writers**

Change the import and construction sites while retaining the existing `&mut dyn Write` boundary:

~~~rust
use std::io::{BufRead, BufReader, BufWriter, Write};

let mut output = BufWriter::new(File::create(path)?);
~~~

Do not add per-record flushes. Allow normal drop or an explicit final flush to emit the buffered data.

- [ ] **Step 3: Verify behavior and measure**

~~~bash
cargo test --manifest-path rust/Cargo.toml
cargo build --release --manifest-path rust/Cargo.toml
bash tests/integration/run_filter.sh
bash tests/integration/run_malformed.sh
./benchmark/run_v3_benchmarks.sh rust-filter-buffered
python3 benchmark/compare_results.py benchmark/results/baseline.json \
  benchmark/results/rust-filter-buffered.json \
  --target rust_filter_scaled --min-improvement 0.05 \
  --guard-all --max-regression 0.02
~~~

Guard every other common case. Expected: exact output parity, at least 5% faster scaled Rust filter, and no guard slower by more than 2%. Revert the production edit if the target fails.

- [ ] **Step 4: Commit**

~~~bash
git add rust/src/filter.rs
git commit -m "perf(rust): buffer filter output"
~~~

### Task 3: Index Rust filter counts by sequence

**Files:**
- Modify: `rust/src/filter.rs`
- Modify: `rust/tests/filter_test.rs`

- [ ] **Step 1: Add failing tests for the indexed representation**

Define the public test boundary:

~~~rust
use ahash::{HashMap, HashSet};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ReplicateIndex {
    pub forward_tag: String,
    pub reverse_tag: String,
    pub counts_by_sequence: HashMap<String, i64>,
}

pub type SampleIndex = HashMap<usize, ReplicateIndex>;
pub fn index_haps(x: usize, haps: &HashMap<usize, Vec<Vec<String>>>) -> SampleIndex;
pub fn all_sequences(index: &SampleIndex) -> HashSet<&str>;
~~~

`filter.rs` uses `ahash::HashMap`, not `std::collections::HashMap`; qualify it so the derive and the import agree.

Counts are stored as `i64`, matching the current `let c: i64 = c_str.parse().unwrap_or(0)`. Do not narrow them to `u32`. The existing threshold test is `(c as u32) < t`, so a negative count wraps to a huge unsigned value and is NOT treated as below threshold; parsing straight into `u32` would fail on `-1`, fall back to `0`, and flip that decision. No fixture contains a negative count, so the parity gate would not catch the change. Narrowing the type is a real behaviour change and belongs in its own task with its own test, not smuggled into an indexing refactor. (It would in fact align Rust with Python, which treats `-1` as below threshold, so it is worth doing deliberately later.)

`forward_tag` and `reverse_tag` are non-optional here, but the current code only reads them when `hap_j[0].len() > 2` and otherwise writes `""` via `f.get(&j).unwrap_or("")`. Default both fields to the empty string for a short first row so the rendered output is unchanged.

Test union construction, missing counts as zero, invalid counts as zero, negative counts round-tripping unchanged as `i64`, a short first row yielding empty tag fields, and duplicate sequence rows retaining the first count. Run `cargo test --manifest-path rust/Cargo.toml filter_test`; expected: FAIL because these symbols do not exist.

- [ ] **Step 2: Implement the index once per sample**

Build `counts_by_sequence` with `entry(sequence).or_insert(count)` so duplicates preserve current `.position()` behavior. Parse the count once, as `i64` with `unwrap_or(0)`, and keep the existing `(c as u32) < t` threshold comparison verbatim so the negative-count behaviour described above is untouched. Change `make_comparison_file` to accept `&SampleIndex` and perform direct `counts_by_sequence.get(sequence)` lookups. Remove the parallel tag/count/sequence structures once no caller uses them.

- [ ] **Step 3: Verify and benchmark**

~~~bash
cargo test --manifest-path rust/Cargo.toml
cargo build --release --manifest-path rust/Cargo.toml
bash tests/integration/run_filter.sh
bash tests/integration/run_malformed.sh
./benchmark/run_v3_benchmarks.sh rust-filter-indexed
python3 benchmark/compare_results.py benchmark/results/rust-filter-buffered.json \
  benchmark/results/rust-filter-indexed.json \
  --target rust_filter_scaled --min-improvement 0.05 \
  --guard-all --max-regression 0.02
~~~

Guard every other common case. Reject the edit if output changes, the scaled target fails, or any guard regresses by more than 2%.

- [ ] **Step 4: Commit**

~~~bash
git add rust/src/filter.rs rust/tests/filter_test.rs
git commit -m "perf(rust): index filter counts by sequence"
~~~

### Task 4: Index Python filter counts by sequence

**Files:**
- Modify: `python/dame/modules_filter.py`
- Modify: `python/dame/filter.py`
- Modify: `python/tests/test_filter.py`

- [ ] **Step 1: Add failing helper tests**

Add tests for the exact public helper interfaces `buildReplicateIndexes(X: int, haps: dict) -> dict` and `allSequences(indexes: dict) -> set[str]`.

Cover union construction, raw count-text retention, and first duplicate row winning. Add output-level regression tests proving that invalid count text still raises `ValueError` from `MakeComparisonFile` and that a valid non-canonical count such as `007` is written back as `007`, not normalized to `7`. Run `python3 -m pytest python/tests/test_filter.py -q`; expected: FAIL because the helpers do not exist.

- [ ] **Step 2: Replace repeated list scans with dictionary lookups**

Build one dictionary per replicate in `filter.run`, keyed by sequence, with raw count strings as its values: `counts_by_sequence: dict[str, str]`. Change `MakeComparisonFile` to accept these indexes and retrieve present values with `counts_by_sequence.get(sequence)`. For a missing sequence, use integer zero. For a present sequence, preserve the raw string for output and retain `int(count)` at the existing threshold comparison so invalid text raises at the same point as v3.0.0. Build the dictionary with first-row-wins semantics, delete the list-comprehension position search only after all callers use the index, and do not normalize count formatting.

- [ ] **Step 3: Verify and benchmark**

~~~bash
python3 -m pytest python/tests/test_filter.py -q
python3 -m pytest python/tests -q
cargo build --release --manifest-path rust/Cargo.toml
bash tests/integration/run_filter.sh
bash tests/integration/run_malformed.sh
./benchmark/run_v3_benchmarks.sh python-filter-indexed
python3 benchmark/compare_results.py benchmark/results/rust-filter-indexed.json \
  benchmark/results/python-filter-indexed.json \
  --target python_filter_scaled --min-improvement 0.05 \
  --guard-all --max-regression 0.02
~~~

Guard every other common case. Reject the edit if parity changes or a guard regresses by more than 2%.

- [ ] **Step 4: Commit**

~~~bash
git add python/dame/modules_filter.py python/dame/filter.py python/tests/test_filter.py
git commit -m "perf(python): index filter counts by sequence"
~~~

### Task 5: Remove NumPy from Python CLI startup paths

**Files:**
- Modify: `python/dame/rsi.py`
- Modify: `python/tests/test_rsi.py`

- [ ] **Step 1: Add a failing subprocess import test**

~~~python
def test_importing_cli_does_not_import_numpy():
    code = "import sys; import dame.__main__; assert 'numpy' not in sys.modules"
    completed = subprocess.run([sys.executable, "-c", code], check=False)
    assert completed.returncode == 0
~~~

Run `python3 -m pytest python/tests/test_rsi.py -q`; expected: FAIL because `dame.rsi` imports NumPy at module load.

- [ ] **Step 2: Import NumPy only inside RSI execution**

~~~python
def _numpy():
    import numpy
    return numpy
~~~

Call `_numpy()` inside RSI functions that need it. Keep all non-RSI CLI imports NumPy-free.

- [ ] **Step 3: Verify startup and sustained-work guards**

~~~bash
python3 -m pytest python/tests -q
./benchmark/run_v3_benchmarks.sh python-lazy-numpy
python3 benchmark/compare_results.py benchmark/results/python-filter-indexed.json \
  benchmark/results/python-lazy-numpy.json \
  --target python_filter_tiny --min-improvement 0.10 \
  --guard-all --max-regression 0.02
~~~

Guard every other common case at 2%. Reject the edit if tiny-filter startup improves by less than 10%.

- [ ] **Step 4: Commit**

~~~bash
git add python/dame/rsi.py python/tests/test_rsi.py
git commit -m "perf(python): load NumPy only for RSI"
~~~

### Task 6: Add a fast exact-match path to Python sort

**Files:**
- Modify: `python/dame/modules_sort.py`
- Modify: `python/dame/sort.py`
- Modify: `python/tests/test_sort.py`

- [ ] **Step 1: Add failing tests for prebuilt exact-match data**

Define and test:

~~~python
def build_exact_tag_lookup(tags):
    """Return forward-sequence and reverse-sequence dictionaries keyed by exact tag text."""
~~~

Do not assert a tag-order winner for duplicate tag sequences. Since v3.0.0 (`fcf1f3e`) `readTags` refuses a Tags file that repeats a sequence, raising `ValueError: duplicate tag sequence ...`, so no winner exists to preserve. Assert instead that such a file is still refused, and that a well-formed panel builds one dictionary entry per tag.

That enforced uniqueness is what makes this task safe: the forward and reverse dictionaries cannot collide, which is the same invariant the Rust `by_fwd`/`by_rc` maps already rely on. Build them with a plain assignment; no collision handling is needed.

Update primer compilation tests to require compiled regex objects and inspect their `.pattern`. Add a regression test in which multiple primer matches exist and the same match is selected as before.

- [ ] **Step 2: Build exact tag lookups once per command**

In `sort.run`, create forward and reverse tag dictionaries once and pass them to `GetPiecesInfo`. Compile primer regexes once as `re.Pattern` objects.

- [ ] **Step 3: Stop materializing all regex matches and scanning every tag**

Use `pattern.search(line)` for the first required primer match and dictionary `.get(tag_text)` for exact tag identification. Preserve orientation, retained-primer output, and tie behavior. Do not change `GetPiecesInfoMismatch` in this task.

- [ ] **Step 4: Verify and benchmark**

~~~bash
python3 -m pytest python/tests/test_sort.py -q
python3 -m pytest python/tests -q
cargo build --release --manifest-path rust/Cargo.toml
bash tests/integration/run_sort.sh
bash tests/integration/run_malformed.sh
./benchmark/run_v3_benchmarks.sh python-sort-exact
python3 benchmark/compare_results.py benchmark/results/python-lazy-numpy.json \
  benchmark/results/python-sort-exact.json \
  --target python_sort_default --min-improvement 0.05 \
  --guard-all --max-regression 0.02
~~~

Guard every other common case. Reject the edit if the exact target fails or a guard regresses by more than 2%.

- [ ] **Step 5: Commit**

~~~bash
git add python/dame/modules_sort.py python/dame/sort.py python/tests/test_sort.py
git commit -m "perf(python): accelerate exact sort matching"
~~~

### Task 7: Compute Python mismatch tag candidates once per read

**Files:**
- Modify: `python/dame/modules_sort.py`
- Modify: `python/tests/test_sort.py`

- [ ] **Step 1: Add focused failing candidate tests**

Define:

~~~python
def tag_candidates_at_ends(line, tags, max_mismatches):
    """Return left and right candidate tuples of (tag name, tag length, mismatch count)."""
~~~

Test a single candidate, equal-scoring ambiguous candidates, variable tag lengths, IUPAC bases, reverse orientation, and mismatch counts above one. Assert order matches input tag order so current ties remain deterministic.

- [ ] **Step 2: Move tag comparison outside primer and orientation loops**

Compute left and right candidates once per read. Cross only those candidates while evaluating compatible primer/orientation combinations. Preserve the existing total mismatch score and winner selection. Do not precompute every tolerated tag variant; candidate lists keep arbitrary mismatch counts and IUPAC semantics explicit.

- [ ] **Step 3: Verify all mismatch paths and benchmark**

~~~bash
python3 -m pytest python/tests/test_sort.py -q
python3 -m pytest python/tests -q
cargo build --release --manifest-path rust/Cargo.toml
bash tests/integration/run_sort_mismatch.sh
bash tests/integration/run_sort_tag_mismatch.sh
bash tests/integration/run_malformed.sh
./benchmark/run_v3_benchmarks.sh python-sort-mismatch-candidates
python3 benchmark/compare_results.py \
  benchmark/results/python-sort-exact.json \
  benchmark/results/python-sort-mismatch-candidates.json \
  --target python_sort_primer_mm1 \
  --target python_sort_tag_mm1 \
  --min-improvement 0.20 \
  --guard-all --max-regression 0.02
~~~

Run the gate with `python_sort_primer_mm1` and `python_sort_tag_mm1` as targets, a 20% minimum improvement, and every other common case guarded at 2%. Reject the edit if either target misses or output differs.

The 20% floor is calibrated to the checked-in fixture's one primer set and eight tags. Instrumentation on 10,000 fixture reads found 24 tag comparisons per read in the current matcher and 16 in the hoisted design, a 33.3% reduction in tag comparisons while all three primer comparisons remain. Larger primer panels may benefit more, but must be measured rather than assigned this fixture's threshold.

- [ ] **Step 4: Commit**

~~~bash
git add python/dame/modules_sort.py python/tests/test_sort.py
git commit -m "perf(python): prune mismatch tag candidates"
~~~

### Task 8: Reduce low-risk Rust sort hashing and uppercase allocation

**Files:**
- Modify: `rust/src/sort.rs`
- Modify: `rust/tests/sort_test.rs`

- [ ] **Step 1: Add failing normalization tests**

Define an internal helper and test both ownership paths:

~~~rust
fn normalize_read(line: &str) -> Cow<'_, str>;
~~~

An uppercase read must return `Cow::Borrowed`; a lowercase read must return the identical uppercase text as `Cow::Owned`.

- [ ] **Step 2: Use the existing fast hasher for ordered maps**

~~~rust
use ahash::RandomState;
use indexmap::IndexMap;
use std::borrow::Cow;

type OrderedMap<K, V> = IndexMap<K, V, RandomState>;
pub type Hap = OrderedMap<String, HapEntry>;
~~~

Use `OrderedMap` for the nested sequence and primer maps, initializing each with `RandomState::new()`. Preserve insertion order.

- [ ] **Step 3: Avoid allocating uppercase copies for already-uppercase reads**

Return a borrowed value unless any byte is ASCII lowercase. Use the helper in exact and mismatch paths. Do not refactor `PieceInfo` ownership in this task.

- [ ] **Step 4: Verify and benchmark**

~~~bash
cargo test --manifest-path rust/Cargo.toml
cargo build --release --manifest-path rust/Cargo.toml
bash tests/integration/run_sort.sh
bash tests/integration/run_sort_mismatch.sh
bash tests/integration/run_sort_tag_mismatch.sh
bash tests/integration/run_malformed.sh
./benchmark/run_v3_benchmarks.sh rust-sort-low-allocation
python3 benchmark/compare_results.py \
  benchmark/results/python-sort-mismatch-candidates.json \
  benchmark/results/rust-sort-low-allocation.json \
  --target rust_sort_default --min-improvement 0.05 \
  --guard-all --max-regression 0.02
~~~

Target `rust_sort_default` at 5% and guard every other common case at 2%. If hashing and normalization together miss the target, benchmark them separately and retain only the independently successful change.

- [ ] **Step 5: Commit retained changes**

~~~bash
git add rust/src/sort.rs rust/tests/sort_test.rs
git commit -m "perf(rust): reduce sort hashing and normalization cost"
~~~

### Task 9: Profile again and stop unless evidence supports deeper work

**Files:**
- Create only if justified: a separately dated plan under `docs/superpowers/plans/`, named for the measured hotspot.

- [ ] **Step 1: Run the final candidate against the original baseline**

~~~bash
./benchmark/run_v3_benchmarks.sh optimized
python3 benchmark/compare_results.py benchmark/results/baseline.json \
  benchmark/results/optimized.json \
  --guard-all --max-regression 0.02
~~~

This comparison must assert, not merely report. Each individual task allows a guarded case to regress by up to 2%, so across seven tasks a case could drift about 15% slower while every per-task gate passed. Only this end-to-end comparison against the original baseline can catch that accumulation, so run it with every case guarded at 2% and treat a failure as blocking.

Record per-case medians, throughput, and relative changes in a scratch note under `benchmark/results/`.

- [ ] **Step 2: Profile the remaining slow representative cases**

Profile scaled Rust filter, Rust exact sort, Python exact sort, and both Python mismatch paths in the same release/container configuration. Collect function time and allocation evidence where supported by the platform.

- [ ] **Step 3: Apply the stop rule**

Stop this implementation plan after Task 8 unless profiling shows one remaining hotspot is both dominant and reachable without changing output semantics. Possible follow-on plans, each kept separate, are:

- replace Rust `PieceInfo` strings with numeric IDs and borrowed byte ranges if cloning/formatting dominates;
- add an anchored exact matcher only after a compatibility corpus proves identical primer-selection behavior;
- add Rayon only for multi-million-record fixtures where CPU saturation outweighs thread and merge overhead;
- test thin LTO and one codegen unit only if build-time cost is acceptable and end-to-end gain is at least 3%.

Do not combine these alternatives into one patch.

### Task 10: Run full verification and publish concise results

**Files:**
- Modify: `README.md`
- Modify: `docs/performance-sort-v1-v3.md`
- Modify: `docs/performance-filter-v1-v3.md`

- [ ] **Step 1: Run all automated tests**

~~~bash
python3 -m pytest python/tests -q
cargo test --manifest-path rust/Cargo.toml
cargo build --release --manifest-path rust/Cargo.toml
for script in tests/integration/run_*.sh; do bash "$script"; done
git diff --check
~~~

Expected: every unit and integration test passes and there is no whitespace error.

- [ ] **Step 2: Re-run the published benchmark workload**

Run the v1-versus-v3 reporting harness only for the final documentation refresh, then run the checked-in Python-v3-versus-Rust-v3 gate. Confirm every pair has identical canonical output.

- [ ] **Step 3: Update the performance text and tables**

Keep both detailed reports under `docs/`. Keep the README summary short and include only:

- DAMe 1.0.0 versus DAMe 3.0;
- Python 3.0 versus Rust 3.0;
- exact and mismatch paths already represented in the current reports;
- a simple explanation that attributes measured gains to buffering, indexed lookups, fewer interpreter loops, lazy imports, faster hashing, or avoided allocations only where the retained changes support that claim.

Use four significant digits for README time and throughput values and round README speedups to the nearest tenth. Do not add version 2 results.

- [ ] **Step 4: Review the final diff and commit**

~~~bash
git diff --check
git status --short
git diff -- README.md docs/performance-sort-v1-v3.md docs/performance-filter-v1-v3.md
git add README.md docs/performance-sort-v1-v3.md docs/performance-filter-v1-v3.md
git commit -m "docs: update v3 performance results"
~~~

The final handoff must state which candidate optimizations cleared their gates, list the exact verification commands run, and link the two performance reports.
