#!/usr/bin/env bash
set -euo pipefail

FIXTURES="$(cd "$(dirname "$0")/../fixtures" && pwd)"
TMPDIR=$(mktemp -d)
trap "rm -rf '$TMPDIR'" EXIT

echo "==> Running dame-py rsi on fixtures..."
cd "$TMPDIR"
dame-py rsi "$FIXTURES/Comparisons_4PCRs.txt"

if [ ! -f RSI_output.txt ]; then
    echo "FAIL: RSI_output.txt not created"; exit 1
fi

if ! grep -q "Sample" RSI_output.txt; then
    echo "FAIL: RSI_output.txt missing header"; exit 1
fi

# RSI values must be between 0 and 1
while IFS=$'\t' read -r sample rsi_val; do
    [[ "$sample" == "Sample" ]] && continue
    if ! awk "BEGIN{exit !($rsi_val >= 0 && $rsi_val <= 1)}"; then
        echo "FAIL: RSI value out of range: $rsi_val for sample $sample"; exit 1
    fi
done < RSI_output.txt

echo "PASS: dame-py rsi produced valid RSI_output.txt"

echo "==> Running dame-py rsi --explicit..."
dame-py rsi --explicit -o RSI_explicit.txt "$FIXTURES/Comparisons_4PCRs.txt"

if [ ! -f RSI_explicit.txt ]; then
    echo "FAIL: RSI_explicit.txt not created"; exit 1
fi

if ! grep -q "ReplicateA" RSI_explicit.txt; then
    echo "FAIL: RSI_explicit.txt missing explicit header"; exit 1
fi

echo "PASS: dame-py rsi --explicit produced valid output"
