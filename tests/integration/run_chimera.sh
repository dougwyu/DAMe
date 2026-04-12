#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
DAME_BIN="$REPO_ROOT/rust/target/release/dame"

if [ ! -f "$DAME_BIN" ]; then
    echo "SKIP: dame binary not found at $DAME_BIN"
    exit 0
fi

if ! command -v usearch &> /dev/null; then
    echo "SKIP: usearch not found on PATH — skipping chimera integration test"
    exit 0
fi

TMPPY=$(mktemp -d)
TMPRS=$(mktemp -d)
trap "rm -rf '$TMPPY' '$TMPRS'" EXIT

setup_dir() {
    local d="$1"
    cat > "$d/Tag1_Tag2.txt" <<'EOF'
CO1	Tag1	Tag2	5	ATATATATATAT
CO1	Tag1	Tag2	3	GCGCGCGCGCGC
EOF
    cat > "$d/PSinfo.txt" <<'EOF'
SampleA	Tag1	Tag2	1
SampleA	Tag1	Tag2	1
EOF
}

echo "==> Running dame-py chimera..."
setup_dir "$TMPPY"
cd "$TMPPY" && dame-py chimera -psInfo PSinfo.txt -x 2 && cd "$REPO_ROOT"

echo "==> Running dame chimera..."
setup_dir "$TMPRS"
cd "$TMPRS" && "$DAME_BIN" chimera --ps-info PSinfo.txt --x 2 && cd "$REPO_ROOT"

# Compare tag files (don't require usearch output — just the file creation)
for i in 1 2; do
    [ -f "$TMPPY/PS${i}.tags.txt" ] || { echo "FAIL: dame-py PS${i}.tags.txt missing"; exit 1; }
    [ -f "$TMPRS/PS${i}.tags.txt" ] || { echo "FAIL: dame PS${i}.tags.txt missing"; exit 1; }
    if ! diff <(sort "$TMPPY/PS${i}.tags.txt") <(sort "$TMPRS/PS${i}.tags.txt") > /dev/null 2>&1; then
        echo "FAIL: PS${i}.tags.txt differs"
        diff <(sort "$TMPPY/PS${i}.tags.txt") <(sort "$TMPRS/PS${i}.tags.txt")
        exit 1
    fi
done

echo "PASS: dame chimera and dame-py chimera produce identical PS tag files"
