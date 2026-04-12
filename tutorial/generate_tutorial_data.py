#!/usr/bin/env python3
"""
generate_tutorial_data.py
=========================
Generate synthetic FASTQ files (Pool1.fastq and Pool2.fastq) for the DAMe tutorial.

Design
------
Two samples, each with two PCR replicates split across two pools:

  Sample1 Rep1: tag1 (fwd) + tag2 (rev) -> Pool1
  Sample1 Rep2: tag3 (fwd) + tag4 (rev) -> Pool2
  Sample2 Rep1: tag5 (fwd) + tag6 (rev) -> Pool1
  Sample2 Rep2: tag7 (fwd) + tag8 (rev) -> Pool2

Four amplicon types per sample designed to exercise each filter criterion:

  AMP_PASSES_ALL  (60 nt)  high count in both replicates   -> passes all filters
  AMP_FAILS_PROP  (60 nt)  present only in one replicate   -> fails --y 2
  AMP_FAILS_COUNT (60 nt)  count=1 in both replicates      -> fails --t 2
  AMP_FAILS_LENGTH (20 nt) short, high count               -> fails --l 50

Usage
-----
  python generate_tutorial_data.py

Outputs
-------
  Pool1.fastq   (~440 reads)
  Pool2.fastq   (~380 reads)

Requirements
------------
  Python 3.6+, stdlib only.
"""

import random

SEED = 42
random.seed(SEED)

# ---------------------------------------------------------------------------
# IUPAC-aware reverse complement
# ---------------------------------------------------------------------------
_COMP = str.maketrans(
    "ACGTacgtMRWSYKVHDBNmrwsykvhdbn",
    "TGCAtgcaKYWSRMBDHVNkywsrmbdhvn",
)


def rc(seq: str) -> str:
    """Return the reverse complement of a DNA sequence (IUPAC-aware)."""
    return seq[::-1].translate(_COMP)


# ---------------------------------------------------------------------------
# Tags: (fwd_seq, tag_name)
# DAMe Tags.txt format: TagSeq TAB TagName
# sort.rs stores tags[name] = [fwd_seq, rc(fwd_seq)]
# tag1 match (forward orientation): tag1_str == fwd_seq
# tag2 match (forward orientation): tag2_str == rc(fwd_seq)
# ---------------------------------------------------------------------------
TAGS = {
    "tag1": "AACCGGT",
    "tag2": "TTGGCCA",
    "tag3": "CCGGAAT",
    "tag4": "GGCCTTA",
    "tag5": "AATCCGG",
    "tag6": "TTAAGGC",
    "tag7": "GGTAACC",
    "tag8": "CCAATTG",
}

# ---------------------------------------------------------------------------
# Primers
# DAMe Primers.txt format: Name TAB FwdSeq TAB RevSeq
# GCRTGC has ambiguity R=[AG]: we alternate GCATGC / GCGTGC across reads
# ---------------------------------------------------------------------------
FWD_PRIMER_A = "GCATGC"   # R resolved as A
FWD_PRIMER_B = "GCGTGC"   # R resolved as G
REV_PRIMER   = "CTGACT"

# rc(REV_PRIMER) is placed at the 3' end of a forward-orientation read
RC_REV_PRIMER = rc(REV_PRIMER)  # AGTCAG

# ---------------------------------------------------------------------------
# Amplicons (the sequence BETWEEN fwd and rev primers in the read)
# ---------------------------------------------------------------------------
# 60-nt amplicons (A/T/G/C only)
AMP_PASSES_ALL  = "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"
AMP_FAILS_PROP  = "TTAATTAATTAATTAATTAATTAATTAATTAATTAATTAATTAATTAATTAATTAATTAA"
AMP_FAILS_COUNT = "GCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGC"
AMP_FAILS_LENGTH = "AAATTTCCCGGG"   # 12 nt — well below --l 50 threshold

assert len(AMP_PASSES_ALL)  == 60, "AMP_PASSES_ALL must be 60 nt"
assert len(AMP_FAILS_PROP)  == 60, "AMP_FAILS_PROP must be 60 nt"
assert len(AMP_FAILS_COUNT) == 60, "AMP_FAILS_COUNT must be 60 nt"
assert len(AMP_FAILS_LENGTH) < 50, "AMP_FAILS_LENGTH must be < 50 nt"


# ---------------------------------------------------------------------------
# FASTQ record builder
# ---------------------------------------------------------------------------
def make_read(fwd_tag_name: str, rev_tag_name: str, amplicon: str,
              read_id: str, use_primer_b: bool = False) -> str:
    """
    Build a forward-orientation FASTQ record for DAMe sort.

    Read structure (forward orientation):
        [fwd_tag_fwd_seq] [fwd_primer] [amplicon] [rc(rev_primer)] [rc(rev_tag_fwd_seq)]

    DAMe sort checks:
        tag1_str  == TAGS[tag1]['fwd']         (everything before fwd_primer)
        tag2_str  == rc(TAGS[tag2]['fwd'])      (everything after rc(rev_primer))
    """
    fwd_seq = TAGS[fwd_tag_name]
    rev_seq_rc = rc(TAGS[rev_tag_name])   # rc of tag2's fwd seq goes at 3' end

    fwd_primer = FWD_PRIMER_B if use_primer_b else FWD_PRIMER_A

    seq = fwd_seq + fwd_primer + amplicon + RC_REV_PRIMER + rev_seq_rc
    qual = "I" * len(seq)
    return f"@{read_id}\n{seq}\n+\n{qual}\n"


# ---------------------------------------------------------------------------
# Generate reads for one pool
# ---------------------------------------------------------------------------
def generate_pool(pool_name: str, pairs_and_counts: list) -> list:
    """
    Generate FASTQ records for a pool.

    pairs_and_counts: list of (fwd_tag, rev_tag, amplicon, n_reads) tuples.
    Returns a list of FASTQ record strings.
    """
    records = []
    read_num = 1
    for fwd_tag, rev_tag, amplicon, n_reads in pairs_and_counts:
        for i in range(n_reads):
            use_b = (i % 2 == 1)   # alternate primer variants
            read_id = f"{pool_name}_{fwd_tag}_{rev_tag}_{read_num}"
            records.append(
                make_read(fwd_tag, rev_tag, amplicon, read_id, use_b)
            )
            read_num += 1
    return records


# ---------------------------------------------------------------------------
# "Noise" reads: wrong tag pair or no primer match
# ---------------------------------------------------------------------------
def make_noise_reads(pool_name: str, start_num: int, n: int) -> list:
    """
    Generate reads that should NOT be sorted into any valid pair.
    These have a mismatch between tag1/tag2 (wrong combination).
    We use tag1 + fwd_primer + amplicon + rc(rev_primer) + rc(tag6)
    when the pool expects tag1+tag2 -> this produces reads assigned to
    tag1_tag6 which is not a valid sample combination (it gets sorted
    into tag1_tag6.txt but is not referenced in PSinfo.txt).
    """
    records = []
    for i in range(n):
        use_b = (i % 2 == 1)
        # tag1 (fwd) paired with tag6 (rev) — not a valid PSinfo combination
        fwd_seq = TAGS["tag1"]
        rev_seq_rc = rc(TAGS["tag6"])
        fwd_primer = FWD_PRIMER_B if use_b else FWD_PRIMER_A
        seq = fwd_seq + fwd_primer + AMP_PASSES_ALL + RC_REV_PRIMER + rev_seq_rc
        qual = "I" * len(seq)
        read_id = f"{pool_name}_noise_{start_num + i}"
        records.append(f"@{read_id}\n{seq}\n+\n{qual}\n")
    return records


def make_no_primer_reads(pool_name: str, start_num: int, n: int) -> list:
    """
    Generate reads with no recognizable primer — pure noise.
    These will be counted in DAMe sort's error tally.
    """
    records = []
    bases = "ACGT"
    for i in range(n):
        random.seed(SEED + start_num + i)
        seq = "".join(random.choices(bases, k=80))
        qual = "I" * len(seq)
        read_id = f"{pool_name}_noprim_{start_num + i}"
        records.append(f"@{read_id}\n{seq}\n+\n{qual}\n")
    return records


# ---------------------------------------------------------------------------
# Pool 1: tag1+tag2 (Sample1 Rep1) + tag5+tag6 (Sample2 Rep1)
# ---------------------------------------------------------------------------
#
# Filter parameters used in tutorial: --x 2 --y 2 --t 2 --l 50
#
# Expected outcomes:
#   AMP_PASSES_ALL:   count >= 2 in both reps  -> PASSES all filters
#   AMP_FAILS_PROP:   only in Rep1 (Pool1)     -> FAILS --y 2 (not in both replicates)
#   AMP_FAILS_COUNT:  count=1 in both reps     -> FAILS --t 2 (count below threshold)
#   AMP_FAILS_LENGTH: short (12 nt < 50 nt)    -> FAILS --l 50

pool1_pairs = [
    # Sample1 Rep1 (tag1+tag2)
    ("tag1", "tag2", AMP_PASSES_ALL,   50),
    ("tag1", "tag2", AMP_FAILS_PROP,   30),   # only in Pool1, not Pool2
    ("tag1", "tag2", AMP_FAILS_COUNT,   1),
    ("tag1", "tag2", AMP_FAILS_LENGTH, 50),
    # Sample2 Rep1 (tag5+tag6)
    ("tag5", "tag6", AMP_PASSES_ALL,   30),
    ("tag5", "tag6", AMP_FAILS_COUNT,   1),
    ("tag5", "tag6", AMP_FAILS_LENGTH, 30),
]

pool1_records = generate_pool("pool1", pool1_pairs)

# Add noise reads (wrong tag combination — sorted to tag1_tag6.txt, not in PSinfo)
pool1_records += make_noise_reads("pool1", 1, 100)
# Add no-primer reads (DAMe sort errors)
pool1_records += make_no_primer_reads("pool1", 200, 100)

random.shuffle(pool1_records)

# ---------------------------------------------------------------------------
# Pool 2: tag3+tag4 (Sample1 Rep2) + tag7+tag8 (Sample2 Rep2)
# ---------------------------------------------------------------------------
pool2_pairs = [
    # Sample1 Rep2 (tag3+tag4)
    ("tag3", "tag4", AMP_PASSES_ALL,   45),
    # AMP_FAILS_PROP intentionally absent from Pool2 -> fails --y 2
    ("tag3", "tag4", AMP_FAILS_COUNT,   1),
    ("tag3", "tag4", AMP_FAILS_LENGTH, 45),
    # Sample2 Rep2 (tag7+tag8)
    ("tag7", "tag8", AMP_PASSES_ALL,   20),
    ("tag7", "tag8", AMP_FAILS_COUNT,   1),
    ("tag7", "tag8", AMP_FAILS_LENGTH, 20),
]

pool2_records = generate_pool("pool2", pool2_pairs)
pool2_records += make_noise_reads("pool2", 1, 80)
pool2_records += make_no_primer_reads("pool2", 200, 80)

random.shuffle(pool2_records)

# ---------------------------------------------------------------------------
# Write FASTQ files
# ---------------------------------------------------------------------------
with open("Pool1.fastq", "w") as fh:
    fh.writelines(pool1_records)

with open("Pool2.fastq", "w") as fh:
    fh.writelines(pool2_records)

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
print("=" * 60)
print("DAMe Tutorial Data Generator")
print("=" * 60)
print(f"\nPool1.fastq: {len(pool1_records)} reads")
print(f"Pool2.fastq: {len(pool2_records)} reads")

print("\nDesigned amplicon breakdown:")
print(f"  {'Amplicon':<20} {'Pool1 Sample1':>14} {'Pool2 Sample1':>14} {'Pool1 Sample2':>14} {'Pool2 Sample2':>14}")
print(f"  {'-'*20} {'-'*14} {'-'*14} {'-'*14} {'-'*14}")
print(f"  {'AMP_PASSES_ALL':<20} {'50':>14} {'45':>14} {'30':>14} {'20':>14}")
print(f"  {'AMP_FAILS_PROP':<20} {'30':>14} {'absent':>14} {'absent':>14} {'absent':>14}")
print(f"  {'AMP_FAILS_COUNT':<20} {'1':>14} {'1':>14} {'1':>14} {'1':>14}")
print(f"  {'AMP_FAILS_LENGTH':<20} {'50':>14} {'45':>14} {'30':>14} {'20':>14}")

print("\nNoise reads per pool:")
print("  Pool1: 100 wrong-pair reads + 100 no-primer reads")
print("  Pool2:  80 wrong-pair reads +  80 no-primer reads")

print("\nExpected filter outcomes with --x 2 --y 2 --t 2 --l 50:")
print("  AMP_PASSES_ALL  -> PASSES (count >= 2 in both replicates, length >= 50)")
print("  AMP_FAILS_PROP  -> FAILS  (absent from one replicate, y_count < 2)")
print("  AMP_FAILS_COUNT -> FAILS  (count = 1 < t=2 in both replicates)")
print("  AMP_FAILS_LENGTH-> FAILS  (length 12 < l=50)")
print()
print("Next steps:")
print("  mkdir -p pool1 pool2")
print("  dame sort --fq Pool1.fastq --primers Primers.txt --tags Tags.txt")
print("  mv tag*.txt SummaryCounts.txt pool1/")
print("  dame sort --fq Pool2.fastq --primers Primers.txt --tags Tags.txt")
print("  mv tag*.txt SummaryCounts.txt pool2/")
print("  dame filter --ps-info PSinfo.txt --x 2 --y 2 --t 2 --l 50")
print("  dame rsi Comparisons_2PCRs.txt")
