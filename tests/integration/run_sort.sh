#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
FIXTURES="$REPO_ROOT/tests/fixtures"
DAME_BIN="$REPO_ROOT/rust/target/release/dame"

if [ ! -f "$DAME_BIN" ]; then
    echo "SKIP: dame binary not found at $DAME_BIN (run: cd rust && cargo build --release)"
    exit 0
fi

TMPPY=$(mktemp -d)
TMPRS=$(mktemp -d)
trap "rm -rf '$TMPPY' '$TMPRS'" EXIT

echo "==> Running dame-py sort..."
cd "$TMPPY"
dame-py sort \
    -fq "$FIXTURES/sample.fastq" \
    -p  "$FIXTURES/Primers.txt" \
    -t  "$FIXTURES/Tags.txt"

echo "==> Running dame sort..."
cd "$TMPRS"
"$DAME_BIN" sort \
    --fq "$FIXTURES/sample.fastq" \
    --primers "$FIXTURES/Primers.txt" \
    --tags "$FIXTURES/Tags.txt"

# Verify required files exist in both outputs
for f in SummaryCounts.txt Tag1_Tag2.txt; do
    if [ ! -f "$TMPPY/$f" ]; then
        echo "FAIL: dame-py did not produce $f"; exit 1
    fi
    if [ ! -f "$TMPRS/$f" ]; then
        echo "FAIL: dame did not produce $f"; exit 1
    fi
done

echo "==> Comparing SummaryCounts.txt..."
if ! diff <(sort "$TMPPY/SummaryCounts.txt") <(sort "$TMPRS/SummaryCounts.txt"); then
    echo "FAIL: SummaryCounts.txt differs between dame-py and dame"
    exit 1
fi

echo "==> Comparing Tag1_Tag2.txt..."
if ! diff <(sort "$TMPPY/Tag1_Tag2.txt") <(sort "$TMPRS/Tag1_Tag2.txt"); then
    echo "FAIL: Tag1_Tag2.txt differs between dame-py and dame"
    exit 1
fi

echo "PASS: dame and dame-py sort produce identical output"
