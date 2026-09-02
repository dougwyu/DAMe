#!/usr/bin/env bash
set -euo pipefail

if [ "$#" -ne 1 ]; then
    echo "usage: $0 LABEL" >&2
    exit 2
fi

label="$1"
case "$label" in
    *[!A-Za-z0-9._-]*|"")
        echo "invalid benchmark label: $label" >&2
        exit 2
        ;;
esac

repo_root="$(cd "$(dirname "$0")/.." && pwd)"
results="$repo_root/benchmark/results"
image="dame-v3-benchmark"
mkdir -p "$results"

docker build --file "$repo_root/benchmark/Dockerfile" --tag "$image" "$repo_root"
docker run --rm \
    --volume "$results:/results" \
    "$image" \
    --output "/results/$label.json"
