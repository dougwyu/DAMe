#!/usr/bin/env bash
# Benchmark `sort` throughput (default + anchored paths) for both implementations.
# Requires: dame-py installed (pip install -e python/) and the release binary
# built (cargo build --release --manifest-path rust/Cargo.toml).
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
DAME_BIN="$REPO_ROOT/rust/target/release/dame"
GEN="$REPO_ROOT/benchmark/generate_benchmark_data.py"

if [ ! -f "$DAME_BIN" ]; then
    echo "SKIP: dame binary not found (run: cargo build --release --manifest-path rust/Cargo.toml)"; exit 0
fi
if ! command -v dame-py >/dev/null 2>&1; then
    echo "SKIP: dame-py not found on PATH (run: pip install -e python/)"; exit 0
fi

WORK=$(mktemp -d)
trap "rm -rf '$WORK'" EXIT
cd "$WORK"
echo "==> Generating benchmark dataset in $WORK ..."
python3 "$GEN"

ts() { python3 -c "import time;print(int(time.time()*1000))"; }
cleanup() { rm -f SummaryCounts.txt tag1_tag2.txt tag5_tag6.txt tag3_tag4.txt tag7_tag8.txt; }

bench() { # $1=label  $2=binary  $3=args
    best=99999
    for _ in 1 2 3; do
        cleanup
        T0=$(ts); $2 sort $3 >/dev/null 2>&1; T1=$(ts)
        d=$((T1 - T0)); [ "$d" -lt "$best" ] && best=$d
    done
    printf "  %-26s %5d ms (best of 3)\n" "$1" "$best"
}

echo "==> sort, Pool1 (98k reads), best of 3:"
bench "Rust  default"            "$DAME_BIN" "--fq Pool1.fastq --primers Primers.txt --tags Tags.txt"
bench "Rust  --primer-mismatches 1" "$DAME_BIN" "--fq Pool1.fastq --primers Primers.txt --tags Tags.txt --primer-mismatches 1"
bench "Rust  --tag-mismatches 1"  "$DAME_BIN" "--fq Pool1.fastq --primers Primers.txt --tags Tags.txt --tag-mismatches 1"
bench "Py    default"            "dame-py" "-fq Pool1.fastq -p Primers.txt -t Tags.txt"
bench "Py    --primer-mismatches 1" "dame-py" "-fq Pool1.fastq -p Primers.txt -t Tags.txt -m 1"
bench "Py    --tag-mismatches 1"  "dame-py" "-fq Pool1.fastq -p Primers.txt -t Tags.txt -mt 1"
cleanup
echo "Done."
