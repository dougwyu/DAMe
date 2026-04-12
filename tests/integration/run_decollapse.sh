#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
DAME_BIN="$REPO_ROOT/rust/target/release/dame"

if [ ! -f "$DAME_BIN" ]; then
    echo "SKIP: dame binary not found at $DAME_BIN (run: cd rust && cargo build --release)"
    exit 0
fi

TMPPY=$(mktemp -d)
TMPRS=$(mktemp -d)
trap "rm -rf '$TMPPY' '$TMPRS'" EXIT

# Write the same input to both tempdirs
INPUT_DATA='CO1	Tag1	Tag2	3	ATATATATATAT
CO1	Tag1	Tag2	2	GCGCGCGCGCGC'

echo "$INPUT_DATA" > "$TMPPY/input.txt"
echo "$INPUT_DATA" > "$TMPRS/input.txt"

echo "==> Running dame-py decollapse..."
cd "$TMPPY"
dame-py decollapse -input input.txt -outFas out.fasta

echo "==> Running dame decollapse..."
cd "$TMPRS"
"$DAME_BIN" decollapse --input input.txt --out-fas out.fasta

# Verify output files exist
if [ ! -f "$TMPPY/out.fasta" ]; then
    echo "FAIL: dame-py did not produce out.fasta"; exit 1
fi
if [ ! -f "$TMPRS/out.fasta" ]; then
    echo "FAIL: dame did not produce out.fasta"; exit 1
fi

# Verify header count == 5 in both
py_headers=$(grep -c "^>" "$TMPPY/out.fasta")
rs_headers=$(grep -c "^>" "$TMPRS/out.fasta")

if [ "$py_headers" -ne 5 ]; then
    echo "FAIL: dame-py out.fasta has $py_headers headers, expected 5"; exit 1
fi
if [ "$rs_headers" -ne 5 ]; then
    echo "FAIL: dame out.fasta has $rs_headers headers, expected 5"; exit 1
fi

echo "==> Comparing out.fasta outputs..."
if ! diff <(sort "$TMPPY/out.fasta") <(sort "$TMPRS/out.fasta"); then
    echo "FAIL: out.fasta differs between dame-py and dame"
    exit 1
fi

echo "PASS: dame and dame-py decollapse produce identical output (5 sequences)"
