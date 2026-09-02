#!/usr/bin/env python3
"""Generate deterministic tiny and scaled filter benchmark fixtures."""

from __future__ import annotations

import argparse
from pathlib import Path


def _sequence(number: int, length: int = 120) -> str:
    bases = "ACGT"
    chars = ["A"] * length
    position = length - 1
    while number:
        number, remainder = divmod(number, 4)
        chars[position] = bases[remainder]
        position -= 1
    return "".join(chars)


def generate_filter_fixture(
    root: Path,
    samples: int,
    records_per_replicate: int,
    union_per_sample: int,
) -> int:
    """Write a two-replicate filter fixture and return collapsed record count."""
    if not records_per_replicate <= union_per_sample <= 2 * records_per_replicate:
        raise ValueError(
            "union_per_sample must be between one and two replicate record counts"
        )

    pool = root / "pool1"
    pool.mkdir(parents=True, exist_ok=True)
    psinfo_lines = []
    for sample_number in range(samples):
        sample = f"Sample{sample_number:03d}"
        starts = (0, union_per_sample - records_per_replicate)
        for replicate, start in enumerate(starts, 1):
            forward = f"F{sample_number:03d}{replicate}"
            reverse = f"R{sample_number:03d}{replicate}"
            psinfo_lines.append(f"{sample}\t{forward}\t{reverse}\t1")
            rows = []
            sequence_base = sample_number * union_per_sample
            for sequence_offset in range(start, start + records_per_replicate):
                count = sequence_offset % 9 + 1
                sequence = _sequence(sequence_base + sequence_offset)
                rows.append(
                    f"CO1\t{forward}\t{reverse}\t{count}\t{sequence}"
                )
            (pool / f"{forward}_{reverse}.txt").write_text("\n".join(rows) + "\n")

    (root / "PSinfo.txt").write_text("\n".join(psinfo_lines) + "\n")
    return samples * 2 * records_per_replicate


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    generate_filter_fixture(
        args.output / "tiny",
        samples=1,
        records_per_replicate=2,
        union_per_sample=2,
    )
    generate_filter_fixture(
        args.output / "scaled",
        samples=100,
        records_per_replicate=200,
        union_per_sample=250,
    )


if __name__ == "__main__":
    main()
