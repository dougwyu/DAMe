#!/usr/bin/env bash
set -euo pipefail

FIXTURES="$(cd "$(dirname "$0")/../fixtures" && pwd)"
TMPDIR=$(mktemp -d)
trap "rm -rf '$TMPDIR'" EXIT

echo "==> Running dame-py sort on fixtures..."
cd "$TMPDIR"
dame-py sort \
    -fq "$FIXTURES/sample.fastq" \
    -p  "$FIXTURES/Primers.txt" \
    -t  "$FIXTURES/Tags.txt"

if [ ! -f SummaryCounts.txt ]; then
    echo "FAIL: SummaryCounts.txt not created"; exit 1
fi

if ! grep -q "#tagName1" SummaryCounts.txt; then
    echo "FAIL: SummaryCounts.txt missing header"; exit 1
fi

if [ ! -f Tag1_Tag2.txt ]; then
    echo "FAIL: Tag1_Tag2.txt not created"; exit 1
fi

if ! grep -q "ATATATATAT" Tag1_Tag2.txt; then
    echo "FAIL: expected barcode not found in Tag1_Tag2.txt"; exit 1
fi

echo "PASS: dame-py sort produced expected output files"
