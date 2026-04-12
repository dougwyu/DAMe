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

# Helper: set up pool/PSinfo fixtures in a directory
setup_fixtures() {
    local dir="$1"
    mkdir -p "$dir/pool1"

    cat > "$dir/pool1/Tag1_Tag2.txt" <<'EOF'
CO1	Tag1	Tag2	5	ATATATATATAT
CO1	Tag1	Tag2	3	GCGCGCGCGCGC
EOF

    cat > "$dir/pool1/Tag3_Tag4.txt" <<'EOF'
CO1	Tag3	Tag4	4	ATATATATATAT
CO1	Tag3	Tag4	1	GCGCGCGCGCGC
EOF

    cat > "$dir/PSinfo.txt" <<'EOF'
SampleA	Tag1	Tag2	1
SampleA	Tag3	Tag4	1
EOF
}

echo "==> Setting up filter fixtures..."
setup_fixtures "$TMPPY"
setup_fixtures "$TMPRS"

echo "==> Running dame-py filter..."
cd "$TMPPY"
dame-py filter -psInfo PSinfo.txt -x 2 -y 1 -t 1 -l 10

echo "==> Running dame filter..."
cd "$TMPRS"
"$DAME_BIN" filter --ps-info PSinfo.txt --x 2 --y 1 --t 1 --l 10

# Verify required output files exist in both
for f in FilteredReads.fna Comparisons_2PCRs.txt; do
    if [ ! -f "$TMPPY/$f" ]; then
        echo "FAIL: dame-py did not produce $f"; exit 1
    fi
    if [ ! -f "$TMPRS/$f" ]; then
        echo "FAIL: dame did not produce $f"; exit 1
    fi
done

echo "==> Comparing FilteredReads.fna..."
# Sort FASTA by sequence line for order-independent comparison.
# dame-py and dame may assign different _N indices but must contain the same sequences.
fasta_sequences() {
    # Extract sequence lines only (non-header lines), sort them
    grep -v "^>" "$1" | sort
}
if ! diff <(fasta_sequences "$TMPPY/FilteredReads.fna") <(fasta_sequences "$TMPRS/FilteredReads.fna"); then
    echo "FAIL: FilteredReads.fna sequences differ between dame-py and dame"
    exit 1
fi

echo "==> Comparing Comparisons_2PCRs.txt..."
if ! diff <(sort "$TMPPY/Comparisons_2PCRs.txt") <(sort "$TMPRS/Comparisons_2PCRs.txt"); then
    echo "FAIL: Comparisons_2PCRs.txt differs between dame-py and dame"
    exit 1
fi

echo "PASS: dame and dame-py filter produce identical output"
