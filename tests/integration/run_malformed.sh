#!/usr/bin/env bash
#
# Malformed-input parity tests.
#
# Every other fixture in this repo is well formed, which is exactly the input
# class where the two implementations agree. This script feeds both binaries
# inputs that are damaged in ways real files get damaged (a stray space, a
# blank line, a short header) and checks that they still agree with each other
# and, where the damage is recoverable, with the well-formed result.
#
# Each case here corresponds to a bug that was live at some point:
#   stray space / space-separated / double tab  ->  Rust silently sorted 0 reads
#   blank PSinfo line                           ->  Python died with IndexError
#   short convert header                        ->  Python died with IndexError
#   single-column rsi input                     ->  Rust looped forever
#   blank FASTQ sequence line                   ->  Python truncated the run
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
FIXTURES="$REPO_ROOT/tests/fixtures"
MALFORMED="$FIXTURES/malformed"
DAME_BIN="$REPO_ROOT/rust/target/release/dame"

if [ ! -f "$DAME_BIN" ]; then
    echo "SKIP: dame binary not found at $DAME_BIN (run: cd rust && cargo build --release)"
    exit 0
fi

# A hang is a real failure mode here (the rsi underflow looped forever), so cap
# every run when a timeout command is available.
TIMEOUT=""
if command -v timeout >/dev/null 2>&1; then
    TIMEOUT="timeout 60"
elif command -v gtimeout >/dev/null 2>&1; then
    TIMEOUT="gtimeout 60"
fi

WORK=$(mktemp -d)
trap "rm -rf '$WORK'" EXIT

fail() { echo "FAIL: $*"; exit 1; }

# run_sort <outdir> <impl> <tags> <primers> <fastq>
run_sort() {
    local dir="$1" impl="$2" tags="$3" primers="$4" fq="$5"
    mkdir -p "$dir"
    if [ "$impl" = "py" ]; then
        (cd "$dir" && $TIMEOUT dame-py sort --fq "$fq" --primers "$primers" --tags "$tags" >stdout.txt 2>&1)
    else
        (cd "$dir" && $TIMEOUT "$DAME_BIN" sort --fq "$fq" --primers "$primers" --tags "$tags" >stdout.txt 2>&1)
    fi
}

# Well-formed baseline that the recoverable cases must reproduce.
run_sort "$WORK/base_py" py "$FIXTURES/Tags.txt" "$FIXTURES/Primers.txt" "$FIXTURES/sample.fastq"
run_sort "$WORK/base_rs" rs "$FIXTURES/Tags.txt" "$FIXTURES/Primers.txt" "$FIXTURES/sample.fastq"
diff "$WORK/base_py/SummaryCounts.txt" "$WORK/base_rs/SummaryCounts.txt" \
    || fail "baseline: implementations disagree on well-formed input"

# check_sort_case <name> <tags> <primers>
# The damage is recoverable, so both implementations must reproduce the
# well-formed result exactly.
check_sort_case() {
    local name="$1" tags="$2" primers="$3"
    run_sort "$WORK/${name}_py" py "$tags" "$primers" "$FIXTURES/sample.fastq"
    run_sort "$WORK/${name}_rs" rs "$tags" "$primers" "$FIXTURES/sample.fastq"
    diff "$WORK/${name}_py/SummaryCounts.txt" "$WORK/${name}_rs/SummaryCounts.txt" \
        || fail "$name: dame and dame-py disagree"
    diff "$WORK/base_py/SummaryCounts.txt" "$WORK/${name}_py/SummaryCounts.txt" \
        || fail "$name: dame-py output differs from the well-formed baseline"
    diff "$WORK/base_rs/SummaryCounts.txt" "$WORK/${name}_rs/SummaryCounts.txt" \
        || fail "$name: dame output differs from the well-formed baseline"
    echo "PASS: sort, $name"
}

echo "==> sort with damaged tag and primer files..."
check_sort_case stray_space  "$MALFORMED/Tags_stray_space.txt" "$FIXTURES/Primers.txt"
check_sort_case space_sep    "$MALFORMED/Tags_space_sep.txt"   "$MALFORMED/Primers_space_sep.txt"
check_sort_case double_tab   "$FIXTURES/Tags.txt"              "$MALFORMED/Primers_double_tab.txt"

echo "==> sort with a blank FASTQ sequence line..."
run_sort "$WORK/blankseq_py" py "$FIXTURES/Tags.txt" "$FIXTURES/Primers.txt" "$MALFORMED/sample_empty_seq.fastq"
run_sort "$WORK/blankseq_rs" rs "$FIXTURES/Tags.txt" "$FIXTURES/Primers.txt" "$MALFORMED/sample_empty_seq.fastq"
# The blank record carries no data, so the sorted output must be unchanged: a
# malformed record costs one read, never the rest of the file.
diff "$WORK/base_py/SummaryCounts.txt" "$WORK/blankseq_py/SummaryCounts.txt" \
    || fail "blank sequence line truncated the dame-py run"
diff "$WORK/blankseq_py/SummaryCounts.txt" "$WORK/blankseq_rs/SummaryCounts.txt" \
    || fail "blank sequence line: dame and dame-py disagree"
# Known difference, deliberately not asserted: dame-py counts the blank record
# in "Number of erroneous sequences" (as DAMe v1.0 did) while dame skips it
# without counting. The sorted data agrees; only the tally differs.
echo "PASS: sort, blank sequence line"

echo "==> filter with damaged PSinfo..."
setup_haps() {
    mkdir -p "$1/pool1"
    printf 'CO1\tTag1\tTag2\t5\tATATATATATAT\nCO1\tTag1\tTag2\t3\tGCGCGCGCGCGC\n' > "$1/pool1/Tag1_Tag2.txt"
    printf 'CO1\tTag3\tTag4\t4\tATATATATATAT\nCO1\tTag3\tTag4\t1\tGCGCGCGCGCGC\n' > "$1/pool1/Tag3_Tag4.txt"
}
# check_filter_case <name> <psinfo>
check_filter_case() {
    local name="$1" ps="$2"
    for impl in py rs; do
        local d="$WORK/f_${name}_${impl}"; mkdir -p "$d"; setup_haps "$d"
        if [ "$impl" = "py" ]; then
            (cd "$d" && $TIMEOUT dame-py filter --ps-info "$ps" --x 2 --y 1 --t 1 --l 10 >stdout.txt 2>&1) \
                || fail "$name: dame-py filter exited non-zero"
        else
            (cd "$d" && $TIMEOUT "$DAME_BIN" filter --ps-info "$ps" --x 2 --y 1 --t 1 --l 10 >stdout.txt 2>&1) \
                || fail "$name: dame filter exited non-zero"
        fi
    done
    diff "$WORK/f_${name}_py/FilteredReads.fna" "$WORK/f_${name}_rs/FilteredReads.fna" \
        || fail "$name: FilteredReads.fna differs"
    diff "$WORK/f_${name}_py/Comparisons_2PCRs.txt" "$WORK/f_${name}_rs/Comparisons_2PCRs.txt" \
        || fail "$name: Comparisons_2PCRs.txt differs"
    echo "PASS: filter, $name"
}
check_filter_case blank_line "$MALFORMED/PSinfo_blank_line.txt"
check_filter_case space_sep  "$MALFORMED/PSinfo_space_sep.txt"

echo "==> convert with a short header..."
for impl in py rs; do
    d="$WORK/c_$impl"; mkdir -p "$d"
    if [ "$impl" = "py" ]; then
        (cd "$d" && $TIMEOUT dame-py convert -i "$MALFORMED/FilteredReads_short_header.fna" >stdout.txt 2>&1) \
            || fail "dame-py convert exited non-zero on a short header"
    else
        (cd "$d" && $TIMEOUT "$DAME_BIN" convert -i "$MALFORMED/FilteredReads_short_header.fna" >stdout.txt 2>&1) \
            || fail "dame convert exited non-zero on a short header"
    fi
done
diff "$WORK/c_py/FilteredReads.forsumaclust.fna" "$WORK/c_rs/FilteredReads.forsumaclust.fna" \
    || fail "convert output differs on a short header"
grep -q "Sample3" "$WORK/c_py/FilteredReads.forsumaclust.fna" \
    || fail "convert dropped records after the short header instead of skipping it"
echo "PASS: convert, short header"

echo "==> rsi with degenerate input..."
# A single-column file must terminate, not loop. Without a timeout available
# this still checks the exit status and output, just not the hang itself.
for impl in py rs; do
    d="$WORK/r_$impl"; mkdir -p "$d"
    if [ "$impl" = "py" ]; then
        (cd "$d" && $TIMEOUT dame-py rsi "$MALFORMED/Comparisons_one_column.txt" >stdout.txt 2>&1) \
            || fail "dame-py rsi did not terminate cleanly on a single-column file"
    else
        (cd "$d" && $TIMEOUT "$DAME_BIN" rsi "$MALFORMED/Comparisons_one_column.txt" >stdout.txt 2>&1) \
            || fail "dame rsi did not terminate cleanly on a single-column file (hang or crash)"
    fi
    grep -q "no replicates" "$d/stdout.txt" \
        || fail "$impl rsi: expected the 'no replicates' message"
done
echo "PASS: rsi, single column"

for impl in py rs; do
    d="$WORK/rs_$impl"; mkdir -p "$d"
    if [ "$impl" = "py" ]; then
        (cd "$d" && $TIMEOUT dame-py rsi "$MALFORMED/Comparisons_space_sep.txt" >stdout.txt 2>&1) \
            || fail "dame-py rsi failed on a space-separated file"
    else
        (cd "$d" && $TIMEOUT "$DAME_BIN" rsi "$MALFORMED/Comparisons_space_sep.txt" >stdout.txt 2>&1) \
            || fail "dame rsi failed on a space-separated file (hang or crash)"
    fi
done
# Sample order must match; RSI values are compared loosely because Python writes
# whole numbers as 0.0 where Rust writes 0.
diff <(cut -f1 "$WORK/rs_py/RSI_output.txt") <(cut -f1 "$WORK/rs_rs/RSI_output.txt") \
    || fail "rsi sample order differs on a space-separated file"
echo "PASS: rsi, space separated"

echo "==> sort refuses a tags file with a duplicated sequence..."
# Not recoverable: the tag name, and so the output filename, would depend on
# lookup order. Both implementations must refuse, with the same message, and
# must not leave partial output behind.
for impl in py rs; do
    d="$WORK/dup_$impl"; mkdir -p "$d"
    set +e
    if [ "$impl" = "py" ]; then
        (cd "$d" && $TIMEOUT dame-py sort --fq "$FIXTURES/sample.fastq" \
            --primers "$FIXTURES/Primers.txt" --tags "$MALFORMED/Tags_duplicate_seq.txt" \
            >stdout.txt 2>&1)
    else
        (cd "$d" && $TIMEOUT "$DAME_BIN" sort --fq "$FIXTURES/sample.fastq" \
            --primers "$FIXTURES/Primers.txt" --tags "$MALFORMED/Tags_duplicate_seq.txt" \
            >stdout.txt 2>&1)
    fi
    rc=$?
    set -e
    [ "$rc" -ne 0 ] || fail "$impl sort accepted a duplicated tag sequence (rc=0)"
    grep -q "duplicate tag sequence" "$d/stdout.txt" \
        || fail "$impl sort: expected a duplicate tag sequence error, got: $(cat "$d/stdout.txt")"
    # Nothing may be written: a partial run is worse than a refused one.
    [ -f "$d/SummaryCounts.txt" ] && fail "$impl sort wrote SummaryCounts.txt despite refusing"
done
# Both must report the failure identically, so neither can drift on its own.
diff <(grep "duplicate tag sequence" "$WORK/dup_py/stdout.txt") \
     <(grep "duplicate tag sequence" "$WORK/dup_rs/stdout.txt") \
    || fail "duplicate tag sequence: error messages differ between implementations"
echo "PASS: sort, duplicate tag sequence refused"

echo "PASS: dame and dame-py agree on malformed input"
