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

echo "==> Running dame-py rsi..."
cd "$TMPPY"
dame-py rsi "$FIXTURES/Comparisons_4PCRs.txt"

echo "==> Running dame rsi..."
cd "$TMPRS"
"$DAME_BIN" rsi "$FIXTURES/Comparisons_4PCRs.txt"

# Verify output files exist
if [ ! -f "$TMPPY/RSI_output.txt" ]; then
    echo "FAIL: dame-py did not produce RSI_output.txt"; exit 1
fi
if [ ! -f "$TMPRS/RSI_output.txt" ]; then
    echo "FAIL: dame did not produce RSI_output.txt"; exit 1
fi

echo "==> Comparing RSI_output.txt with numeric tolerance..."
python3 - "$TMPPY/RSI_output.txt" "$TMPRS/RSI_output.txt" <<'PYEOF'
import sys

def parse_rsi(path):
    header = None
    rows = []
    with open(path) as f:
        for line in f:
            line = line.rstrip('\n')
            if not line:
                continue
            parts = line.split('\t')
            # Detect header row (non-numeric first column like "Sample")
            if header is None:
                try:
                    float(parts[1]) if len(parts) > 1 else None
                    rows.append(parts)
                except (ValueError, IndexError):
                    header = parts
            else:
                rows.append(parts)
    # Sort data rows by first column (sample name) for order-independent comparison
    rows.sort(key=lambda r: r[0])
    return header, rows

py_header, py_rows = parse_rsi(sys.argv[1])
rs_header, rs_rows = parse_rsi(sys.argv[2])

# Compare headers if both have them
if py_header is not None and rs_header is not None:
    if py_header != rs_header:
        print(f"FAIL: headers differ: dame-py={py_header} dame={rs_header}")
        sys.exit(1)

if len(py_rows) != len(rs_rows):
    print(f"FAIL: RSI_output.txt row count differs: dame-py={len(py_rows)} dame={len(rs_rows)}")
    sys.exit(1)

for i, (py_row, rs_row) in enumerate(zip(py_rows, rs_rows)):
    if len(py_row) != len(rs_row):
        print(f"FAIL: row {i} column count differs: {py_row} vs {rs_row}")
        sys.exit(1)
    for j, (py_val, rs_val) in enumerate(zip(py_row, rs_row)):
        # Try numeric comparison with tolerance
        try:
            py_num = float(py_val)
            rs_num = float(rs_val)
            if abs(py_num - rs_num) > 1e-9:
                print(f"FAIL: row {i} col {j}: dame-py={py_val} dame={rs_val} (diff > 1e-9)")
                sys.exit(1)
        except ValueError:
            # Non-numeric: exact string comparison
            if py_val != rs_val:
                print(f"FAIL: row {i} col {j}: dame-py={py_val!r} dame={rs_val!r}")
                sys.exit(1)

print("OK: RSI_output.txt values match within tolerance")
PYEOF

echo "==> Testing dame-py rsi --explicit..."
cd "$TMPPY"
dame-py rsi --explicit -o RSI_explicit.txt "$FIXTURES/Comparisons_4PCRs.txt"

echo "==> Testing dame rsi --explicit..."
cd "$TMPRS"
"$DAME_BIN" rsi --explicit -o RSI_explicit.txt "$FIXTURES/Comparisons_4PCRs.txt"

if [ ! -f "$TMPPY/RSI_explicit.txt" ]; then
    echo "FAIL: dame-py did not produce RSI_explicit.txt"; exit 1
fi
if [ ! -f "$TMPRS/RSI_explicit.txt" ]; then
    echo "FAIL: dame did not produce RSI_explicit.txt"; exit 1
fi

echo "==> Comparing RSI_explicit.txt with numeric tolerance..."
python3 - "$TMPPY/RSI_explicit.txt" "$TMPRS/RSI_explicit.txt" <<'PYEOF'
import sys

def parse_rsi(path):
    header = None
    rows = []
    with open(path) as f:
        for line in f:
            line = line.rstrip('\n')
            if not line:
                continue
            parts = line.split('\t')
            if header is None:
                try:
                    float(parts[1]) if len(parts) > 1 else None
                    rows.append(parts)
                except (ValueError, IndexError):
                    header = parts
            else:
                rows.append(parts)
    rows.sort(key=lambda r: r[0])
    return header, rows

py_header, py_rows = parse_rsi(sys.argv[1])
rs_header, rs_rows = parse_rsi(sys.argv[2])

if py_header is not None and rs_header is not None:
    if py_header != rs_header:
        print(f"FAIL: headers differ: dame-py={py_header} dame={rs_header}")
        sys.exit(1)

if len(py_rows) != len(rs_rows):
    print(f"FAIL: RSI_explicit.txt row count differs: dame-py={len(py_rows)} dame={len(rs_rows)}")
    sys.exit(1)

for i, (py_row, rs_row) in enumerate(zip(py_rows, rs_rows)):
    if len(py_row) != len(rs_row):
        print(f"FAIL: row {i} column count differs: {py_row} vs {rs_row}")
        sys.exit(1)
    for j, (py_val, rs_val) in enumerate(zip(py_row, rs_row)):
        try:
            py_num = float(py_val)
            rs_num = float(rs_val)
            if abs(py_num - rs_num) > 1e-9:
                print(f"FAIL: row {i} col {j}: dame-py={py_val} dame={rs_val} (diff > 1e-9)")
                sys.exit(1)
        except ValueError:
            if py_val != rs_val:
                print(f"FAIL: row {i} col {j}: dame-py={py_val!r} dame={rs_val!r}")
                sys.exit(1)

print("OK: RSI_explicit.txt values match within tolerance")
PYEOF

echo "PASS: dame and dame-py rsi produce equivalent output"
