#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
FIXTURES="$REPO_ROOT/tests/fixtures"
DAME_BIN="$REPO_ROOT/rust/target/release/dame"

if [ ! -f "$DAME_BIN" ]; then
    echo "SKIP: dame binary not found at $DAME_BIN (run: cd rust && cargo build --release)"
    exit 0
fi
if ! command -v dame-py >/dev/null 2>&1; then
    echo "SKIP: dame-py not found on PATH (run: pip install -e python/)"
    exit 0
fi

TMPPY=$(mktemp -d)
TMPRS=$(mktemp -d)
TMPN0=$(mktemp -d)
trap "rm -rf '$TMPPY' '$TMPRS' '$TMPN0'" EXIT

echo "==> dame-py sort -mt 1..."
cd "$TMPPY"
dame-py sort -fq "$FIXTURES/sample_tag_err.fastq" -p "$FIXTURES/Primers.txt" -t "$FIXTURES/Tags.txt" -mt 1

echo "==> dame sort --tag-mismatches 1..."
cd "$TMPRS"
"$DAME_BIN" sort --fq "$FIXTURES/sample_tag_err.fastq" --primers "$FIXTURES/Primers.txt" --tags "$FIXTURES/Tags.txt" --tag-mismatches 1

for D in "$TMPPY" "$TMPRS"; do
    if [ ! -f "$D/Tag1_Tag2.txt" ]; then
        echo "FAIL: $D did not produce Tag1_Tag2.txt at tag-mismatches=1"; exit 1
    fi
done

echo "==> Comparing Tag1_Tag2.txt (py vs rust)..."
if ! diff "$TMPPY/Tag1_Tag2.txt" "$TMPRS/Tag1_Tag2.txt"; then
    echo "FAIL: Tag1_Tag2.txt differs between dame-py and dame at tag-mismatches=1"; exit 1
fi
if ! grep -q "ATATATATATAT" "$TMPRS/Tag1_Tag2.txt"; then
    echo "FAIL: expected recovered barcode ATATATATATAT not found"; exit 1
fi

echo "==> Confirming rejection at tag-mismatches=0 (both tools)..."
cd "$TMPN0"
"$DAME_BIN" sort --fq "$FIXTURES/sample_tag_err.fastq" --primers "$FIXTURES/Primers.txt" --tags "$FIXTURES/Tags.txt"
if [ -f "$TMPN0/Tag1_Tag2.txt" ]; then
    echo "FAIL: rust read should NOT be recovered at tag-mismatches=0"; exit 1
fi
TMPN0PY=$(mktemp -d)
cd "$TMPN0PY"
dame-py sort -fq "$FIXTURES/sample_tag_err.fastq" -p "$FIXTURES/Primers.txt" -t "$FIXTURES/Tags.txt"
if [ -f "$TMPN0PY/Tag1_Tag2.txt" ]; then
    echo "FAIL: dame-py read should NOT be recovered at tag-mismatches=0"; rm -rf "$TMPN0PY"; exit 1
fi
rm -rf "$TMPN0PY"

echo "PASS: dame and dame-py sort agree at tag-mismatches=1 and reject at 0"
