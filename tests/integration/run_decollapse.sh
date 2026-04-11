#!/usr/bin/env bash
set -euo pipefail

TMPDIR=$(mktemp -d)
trap "rm -rf '$TMPDIR'" EXIT

echo "==> Running dame-py decollapse..."
cd "$TMPDIR"

cat > input.txt <<'EOF'
CO1	Tag1	Tag2	3	ATATATATATAT
CO1	Tag1	Tag2	2	GCGCGCGCGCGC
EOF

dame-py decollapse -input input.txt -outFas out.fasta

if [ ! -f out.fasta ]; then
    echo "FAIL: out.fasta not created"; exit 1
fi

header_count=$(grep -c "^>" out.fasta)
if [ "$header_count" -ne 5 ]; then
    echo "FAIL: expected 5 sequences (3+2), got $header_count"; exit 1
fi

if ! grep -q "Tag1.Tag2.3_1" out.fasta; then
    echo "FAIL: expected header Tag1.Tag2.3_1 not found"; exit 1
fi

echo "PASS: dame-py decollapse produced 5 sequences as expected"
