#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
FIXTURES="$REPO_ROOT/tests/fixtures"
DAME_BIN="$REPO_ROOT/rust/target/release/dame"
FNA="$FIXTURES/FilteredReads_small.fna"

if [ ! -f "$DAME_BIN" ]; then
    echo "SKIP: dame binary not found at $DAME_BIN (run: cd rust && cargo build --release)"
    exit 0
fi

TMPPY=$(mktemp -d)
TMPRS=$(mktemp -d)
trap "rm -rf '$TMPPY' '$TMPRS'" EXIT

echo "==> Sumaclust mode..."
(cd "$TMPPY" && dame-py convert -i "$FNA")
(cd "$TMPRS" && "$DAME_BIN" convert -i "$FNA")
diff "$TMPPY/FilteredReads.forsumaclust.fna" "$TMPRS/FilteredReads.forsumaclust.fna" \
    || { echo "FAIL: sumaclust output differs"; exit 1; }
echo "PASS: sumaclust"

echo "==> USEARCH mode..."
(cd "$TMPPY" && dame-py convert -i "$FNA" -u)
(cd "$TMPRS" && "$DAME_BIN" convert -i "$FNA" -u)
diff "$TMPPY/FilteredReads.forusearch.fna" "$TMPRS/FilteredReads.forusearch.fna" \
    || { echo "FAIL: usearch output differs"; exit 1; }
echo "PASS: usearch"

echo "==> USEARCH + --max-length 65..."
(cd "$TMPPY" && dame-py convert -i "$FNA" -u --max-length 65)
(cd "$TMPRS" && "$DAME_BIN" convert -i "$FNA" -u --max-length 65)
diff "$TMPPY/FilteredReads.forusearch.fna" "$TMPRS/FilteredReads.forusearch.fna" \
    || { echo "FAIL: usearch padded output differs"; exit 1; }
echo "PASS: usearch + max-length"

echo "==> --min-length 10 filter..."
(cd "$TMPPY" && dame-py convert -i "$FNA" --min-length 10)
(cd "$TMPRS" && "$DAME_BIN" convert -i "$FNA" --min-length 10)
diff "$TMPPY/FilteredReads.forsumaclust.fna" "$TMPRS/FilteredReads.forsumaclust.fna" \
    || { echo "FAIL: min-length output differs"; exit 1; }
echo "PASS: min-length filter"

echo "==> --sample-fastas..."
(cd "$TMPPY" && dame-py convert -i "$FNA" -s)
(cd "$TMPRS" && "$DAME_BIN" convert -i "$FNA" -s)
[ -d "$TMPPY/SampleFastas" ] || { echo "FAIL: dame-py SampleFastas not created"; exit 1; }
[ -d "$TMPRS/SampleFastas" ] || { echo "FAIL: dame SampleFastas not created"; exit 1; }
diff "$TMPPY/SampleFastas/Sample1.fixed.fasta" "$TMPRS/SampleFastas/Sample1.fixed.fasta" \
    || { echo "FAIL: Sample1.fixed.fasta differs"; exit 1; }
diff "$TMPPY/SampleFastas/Sample2.fixed.fasta" "$TMPRS/SampleFastas/Sample2.fixed.fasta" \
    || { echo "FAIL: Sample2.fixed.fasta differs"; exit 1; }
echo "PASS: sample-fastas"

echo "PASS: dame and dame-py convert produce identical output"
