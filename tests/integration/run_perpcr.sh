#!/usr/bin/env bash
#
# convert --per-pcr end to end on tests/fixtures/perpcr:
#   filter reproduces the committed Comparisons/FilteredReads files;
#   dame and dame-py write the expected per-PCR FASTA and PCRinfo.txt;
#   optionally, the tutorial recipe (vsearch, then R) gives the expected tables.
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
FIX="$REPO_ROOT/tests/fixtures/perpcr"
EXP="$FIX/expected"
DAME_BIN="$REPO_ROOT/rust/target/release/dame"

if [ ! -f "$DAME_BIN" ]; then
    echo "SKIP: dame binary not found at $DAME_BIN (run: cd rust && cargo build --release)"
    exit 0
fi

WORK=$(mktemp -d)
trap "rm -rf '$WORK'" EXIT
fail() { echo "FAIL: $*"; exit 1; }

echo "==> filter reproduces the committed fixture inputs..."
mkdir "$WORK/filter"
cp -R "$FIX/pool1" "$FIX/pool2" "$FIX/PSinfo.txt" "$WORK/filter/"
(cd "$WORK/filter" && "$DAME_BIN" filter --ps-info PSinfo.txt --x 3 --y 2 --t 1 --l 100 >/dev/null)
diff "$WORK/filter/Comparisons_3PCRs.fasta" "$FIX/Comparisons_3PCRs.fasta" || fail "Comparisons_3PCRs.fasta differs"
diff "$WORK/filter/FilteredReads.fna" "$FIX/FilteredReads.fna" || fail "FilteredReads.fna differs"
echo "PASS: filter"

# run_impl <impl> <outdir> <convert args...>
run_impl() {
    local impl="$1" dir="$2"; shift 2
    mkdir -p "$dir"
    if [ "$impl" = "py" ]; then
        (cd "$dir" && dame-py convert "$@")
    else
        (cd "$dir" && "$DAME_BIN" convert "$@")
    fi
}

echo "==> convert --per-pcr --ps-info..."
for impl in py rs; do
    run_impl "$impl" "$WORK/ps_$impl" -i "$FIX/Comparisons_3PCRs.fasta" --per-pcr --ps-info "$FIX/PSinfo.txt"
    diff "$WORK/ps_$impl/FilteredReads.perpcr.fna" "$EXP/FilteredReads.perpcr.fna" || fail "$impl FASTA differs from expected"
    diff "$WORK/ps_$impl/PCRinfo.txt" "$EXP/PCRinfo.txt" || fail "$impl PCRinfo.txt differs from expected"
done
echo "PASS: per-PCR with PSinfo (dame and dame-py)"

echo "==> convert --per-pcr without PSinfo..."
for impl in py rs; do
    run_impl "$impl" "$WORK/nops_$impl" -i "$FIX/Comparisons_3PCRs.fasta" --per-pcr
    diff "$WORK/nops_$impl/PCRinfo.txt" "$EXP/PCRinfo_no_psinfo.txt" || fail "$impl PCRinfo.txt (no PSinfo) differs"
done
echo "PASS: per-PCR without PSinfo"

echo "==> convert --per-pcr with length filters (parity)..."
for impl in py rs; do
    run_impl "$impl" "$WORK/len_$impl" -i "$FIX/Comparisons_3PCRs.fasta" --per-pcr --ps-info "$FIX/PSinfo.txt" --min-length 100 --max-length 200
done
diff "$WORK/len_py/FilteredReads.perpcr.fna" "$WORK/len_rs/FilteredReads.perpcr.fna" || fail "length-filtered FASTA differs"
diff "$WORK/len_py/PCRinfo.txt" "$WORK/len_rs/PCRinfo.txt" || fail "length-filtered PCRinfo.txt differs"
echo "PASS: length filters"

echo "==> warning text parity..."
for impl in py rs; do
    mkdir -p "$WORK/warn_$impl"
    cp "$FIX/FilteredReads.fna" "$WORK/warn_$impl/"
    run_impl "$impl" "$WORK/warn_$impl" -i FilteredReads.fna --per-pcr -u 2>"$WORK/warn_$impl.err"
done
diff "$WORK/warn_py.err" "$WORK/warn_rs.err" || fail "warning/note text differs"
grep -q "^Warning: --per-pcr input looks like FilteredReads output" "$WORK/warn_rs.err" || fail "warning missing"
echo "PASS: warning and note"

if ! command -v vsearch >/dev/null 2>&1; then
    echo "SKIP: vsearch not on PATH; recipe check skipped"
    echo "PASS: per-PCR integration"
    exit 0
fi

echo "==> recipe: dereplicate, map, cluster..."
R="$WORK/recipe"
mkdir "$R"
(cd "$R" && "$DAME_BIN" convert -i "$FIX/FilteredReads.fna" -u >/dev/null)
(cd "$R" && vsearch --derep_fulllength FilteredReads.forusearch.fna --sizein --sizeout --relabel seq --output passed.fna --quiet)
(cd "$R" && vsearch --usearch_global "$WORK/ps_rs/FilteredReads.perpcr.fna" --db passed.fna \
    --id 1.0 --mincols 110 --query_cov 1.0 --otutabout table.tsv --quiet)
diff "$R/table.tsv" "$EXP/table.tsv" || fail "vsearch table differs from expected"
(cd "$R" && vsearch --cluster_size passed.fna --sizein --id 0.97 --uc clusters.uc --quiet)
echo "PASS: vsearch recipe"

if ! command -v Rscript >/dev/null 2>&1 || ! Rscript -e 'library(dplyr); library(tidyr)' >/dev/null 2>&1; then
    echo "SKIP: Rscript with dplyr and tidyr not available; R step skipped"
    echo "PASS: per-PCR integration"
    exit 0
fi

echo "==> recipe: R join..."
Rscript "$REPO_ROOT/tutorial/perpcr_to_occupancy.R" "$R/table.tsv" "$R/clusters.uc" \
    "$WORK/ps_rs/PCRinfo.txt" "$R/survey.tsv" 2>/dev/null
diff "$R/survey.tsv" "$EXP/survey.tsv" || fail "survey table differs from expected"
echo "PASS: R join"

echo "PASS: per-PCR integration"
