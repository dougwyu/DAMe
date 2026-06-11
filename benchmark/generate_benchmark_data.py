#!/usr/bin/env python3
"""Generate the large synthetic DAMe sort benchmark dataset.

Writes Primers.txt, Tags.txt, Pool1.fastq, Pool2.fastq into the current
directory: 2 pools of ~98k reads each (196k total), the 8-tag CO1 tutorial
panel, ~10 unique 100 nt amplicons per tag pair. All reads are valid
(exact tags + primers), so they sort cleanly in both implementations.

Usage:
    mkdir -p /tmp/dame_bench && cd /tmp/dame_bench
    python3 /path/to/benchmark/generate_benchmark_data.py
"""
import random

random.seed(42)

def RC(seq):
    return seq[::-1].translate(str.maketrans("ACGTMRWSYKVHDB", "TGCAKYWSRMBDHV"))

TAGS = {
    "tag1": "AACCGGT", "tag2": "TTGGCCA", "tag3": "CCGGAAT", "tag4": "GGCCTTA",
    "tag5": "AATCCGG", "tag6": "TTAAGGC", "tag7": "GGTAACC", "tag8": "CCAATTG",
}
FWD_PRIMER = "GCATGC"          # concrete instance of IUPAC GCRTGC (R=A)
RC_REV_PRIMER = RC("CTGACT")   # rc of the reverse primer -> AGTCAG

BASES = "ACGT"
AMPLICONS = ["".join(random.choice(BASES) for _ in range(100)) for _ in range(20)]

POOLS = {
    "Pool1.fastq": [("tag1", "tag2"), ("tag5", "tag6")],
    "Pool2.fastq": [("tag3", "tag4"), ("tag7", "tag8")],
}
READS_PER_POOL = 98000


def make_read(fwd_tag, rev_tag, amplicon, rid):
    seq = TAGS[fwd_tag] + FWD_PRIMER + amplicon + RC_REV_PRIMER + RC(TAGS[rev_tag])
    return f"@{rid}\n{seq}\n+\n{'I' * len(seq)}\n"


def main():
    with open("Primers.txt", "w") as fh:
        fh.write("CO1\tGCRTGC\tCTGACT\n")
    with open("Tags.txt", "w") as fh:
        for name, seq in TAGS.items():
            fh.write(f"{seq}\t{name}\n")
    for fname, pairs in POOLS.items():
        with open(fname, "w") as fh:
            for n in range(READS_PER_POOL):
                fwd, rev = pairs[n % len(pairs)]
                amp = AMPLICONS[n % len(AMPLICONS)]
                fh.write(make_read(fwd, rev, amp, f"{fname}_{n}"))
        print(f"{fname}: {READS_PER_POOL} reads")


if __name__ == "__main__":
    main()
