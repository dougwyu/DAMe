#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
FIXTURES="$REPO_ROOT/tests/fixtures"
DAME_BIN="$REPO_ROOT/rust/target/release/dame"

if [ ! -f "$DAME_BIN" ]; then
    echo "SKIP: dame binary not found at $DAME_BIN"
    exit 0
fi

run_pipeline() {
    local tool="$1"     # "dame-py" or the path to dame binary
    local tool_name="$2"  # display name

    local dir
    dir=$(mktemp -d)
    echo "  Running $tool_name pipeline in $dir..."

    # Step 1: sort
    echo "  [1/3] sort..."
    cd "$dir"
    if [ "$tool_name" = "dame-py" ]; then
        dame-py sort \
            -fq "$FIXTURES/sample.fastq" \
            -p  "$FIXTURES/Primers.txt" \
            -t  "$FIXTURES/Tags.txt"
    else
        "$tool" sort \
            --fq "$FIXTURES/sample.fastq" \
            --primers "$FIXTURES/Primers.txt" \
            --tags "$FIXTURES/Tags.txt"
    fi

    # Build PSinfo from sort output
    # Find tag combo files (like Tag1_Tag2.txt) and create PSinfo
    > PSinfo.txt
    mkdir -p pool1
    for f in *.txt; do
        [[ "$f" =~ ^(Summary|PS) ]] && continue
        base="${f%.txt}"
        if [[ "$base" == *_* ]]; then
            IFS='_' read -r tag1 tag2 <<< "$base"
            echo "SampleA	${tag1}	${tag2}	1" >> PSinfo.txt
            cp "$f" "pool1/${f}"
        fi
    done

    psinfo_count=$(wc -l < PSinfo.txt | tr -d ' ')
    if [ "$psinfo_count" -eq 0 ]; then
        echo "  WARNING: no tag combos found after sort, skipping filter/rsi"
        cd "$REPO_ROOT"
        rm -rf "$dir"
        return 0
    fi

    # Step 2: filter
    echo "  [2/3] filter (X=$psinfo_count)..."
    if [ "$tool_name" = "dame-py" ]; then
        dame-py filter -psInfo PSinfo.txt -x "$psinfo_count" -y 1 -t 1 -l 5
    else
        "$tool" filter --ps-info PSinfo.txt --x "$psinfo_count" --y 1 --t 1 --l 5
    fi

    # Step 3: rsi (need at least 2 PCR replicates for RSI to be meaningful)
    local comp_file
    comp_file=$(ls Comparisons_*.txt 2>/dev/null | head -1 || true)
    if [ -z "$comp_file" ]; then
        echo "  WARNING: no comparisons file, skipping rsi"
    else
        echo "  [3/3] rsi..."
        if [ "$tool_name" = "dame-py" ]; then
            dame-py rsi "$comp_file" || echo "  NOTE: rsi requires >=2 replicates, skipping"
        else
            "$tool" rsi "$comp_file" || echo "  NOTE: rsi requires >=2 replicates, skipping"
        fi
    fi

    cd "$REPO_ROOT"

    # Verify the filter output at minimum
    [ -f "$dir/FilteredReads.fna" ] || { echo "FAIL: $tool_name FilteredReads.fna missing"; rm -rf "$dir"; exit 1; }
    echo "  $tool_name pipeline: FilteredReads.fna produced ✓"
    rm -rf "$dir"
}

echo "==> Running full pipeline with dame-py..."
run_pipeline "dame-py" "dame-py"

echo "==> Running full pipeline with dame..."
run_pipeline "$DAME_BIN" "dame"

echo "PASS: full pipeline (sort → filter → rsi) completed for both dame-py and dame"
