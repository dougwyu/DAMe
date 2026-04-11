#!/usr/bin/env bash
set -euo pipefail

TMPDIR=$(mktemp -d)
trap "rm -rf '$TMPDIR'" EXIT

echo "==> Setting up filter fixtures..."
cd "$TMPDIR"
mkdir -p pool1

cat > pool1/Tag1_Tag2.txt <<'EOF'
CO1	Tag1	Tag2	5	ATATATATATAT
CO1	Tag1	Tag2	3	GCGCGCGCGCGC
EOF

cat > pool1/Tag3_Tag4.txt <<'EOF'
CO1	Tag3	Tag4	4	ATATATATATAT
CO1	Tag3	Tag4	1	GCGCGCGCGCGC
EOF

cat > PSinfo.txt <<'EOF'
SampleA	Tag1	Tag2	1
SampleA	Tag3	Tag4	1
EOF

echo "==> Running dame-py filter..."
dame-py filter -psInfo PSinfo.txt -x 2 -y 1 -t 1 -l 10

if [ ! -f FilteredReads.fna ]; then
    echo "FAIL: FilteredReads.fna not created"; exit 1
fi

if [ ! -f Comparisons_2PCRs.txt ]; then
    echo "FAIL: Comparisons_2PCRs.txt not created"; exit 1
fi

if ! grep -q "SampleA" Comparisons_2PCRs.txt; then
    echo "FAIL: SampleA not found in comparison file"; exit 1
fi

echo "PASS: dame-py filter produced expected output files"
