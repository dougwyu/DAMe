# DAMe Tutorial: DNA Metabarcoding from Sort to RSI

This tutorial walks through the full DAMe workflow using a synthetic but realistic dataset.
You will sort raw amplicon reads, filter them by PCR presence and abundance, compute
Renkonen Similarity Index (RSI) values, and optionally decollapse sequences.

---

## Section 1: Introduction

**DAMe** (DNA Metabarcoding) is a toolkit for processing amplicon sequencing data that
has been multiplexed with combinatorial tags. Each PCR reaction is uniquely identified
by a forward and a reverse tag flanking the primer-amplicon-primer region.

### What DAMe Does

```
dame sort   -> demultiplex reads by tag combination, collapse to unique sequences
dame filter -> apply presence (y), count (t), and length (l) thresholds across replicates
dame rsi    -> compute Renkonen Similarity Index between PCR replicates
dame decollapse -> expand collapsed sequences back to one record per original read
```

### Prerequisites

- Python 3.6+ (for `generate_tutorial_data.py`)
- **Either** the `dame` Rust binary (recommended for speed) **or** `dame-py` (Python port)

`dame` is ~5–8× faster than `dame-py` and accepts both plain and
gzip-compressed FASTQ (`.fastq.gz`) with no extra flags.

To build the Rust binary from source:

```bash
cd /path/to/DAMe/rust
cargo build --release
# Binary will be at target/release/dame
```

To install `dame-py`:

```bash
cd /path/to/DAMe/python
pip install -e .
```

---

## Section 2: Tutorial Dataset

The synthetic dataset models a common metabarcoding experiment:

- **2 samples** (Sample1, Sample2)
- **2 PCR replicates each** (Rep1 in Pool1, Rep2 in Pool2)
- **2 sequencing pools**

### Pool Assignments

| Sample   | Replicate | Pool | Fwd Tag | Rev Tag |
|----------|-----------|------|---------|---------|
| Sample1  | Rep1      | 1    | tag1    | tag2    |
| Sample1  | Rep2      | 2    | tag3    | tag4    |
| Sample2  | Rep1      | 1    | tag5    | tag6    |
| Sample2  | Rep2      | 2    | tag7    | tag8    |

### Designed Amplicons and Expected Filter Outcomes

The dataset contains four amplicons per sample, each designed to exercise a specific
filter criterion when running `dame filter --x 2 --y 2 --t 2 --l 50`:

| Amplicon        | Length | Sample1 Rep1 | Sample1 Rep2 | Sample2 Rep1 | Sample2 Rep2 | Expected Outcome          |
|-----------------|--------|--------------|--------------|--------------|--------------|---------------------------|
| AMP_PASSES_ALL  | 60 nt  | 50 reads     | 45 reads     | 30 reads     | 20 reads     | PASSES all filters        |
| AMP_FAILS_PROP  | 60 nt  | 30 reads     | absent       | absent       | absent       | FAILS `--y 2` (only 1 rep)|
| AMP_FAILS_COUNT | 60 nt  | 1 read       | 1 read       | 1 read       | 1 read       | FAILS `--t 2` (count < 2) |
| AMP_FAILS_LENGTH| 12 nt  | 50 reads     | 45 reads     | 30 reads     | 20 reads     | FAILS `--l 50` (too short)|

The pools also contain noise reads (wrong tag combinations, no-primer reads) to
demonstrate DAMe's error handling.

---

## Section 3: Generate the Data

```bash
cd /path/to/DAMe/tutorial
python generate_tutorial_data.py
```

Expected output:

```
============================================================
DAMe Tutorial Data Generator
============================================================

Pool1.fastq: 392 reads
Pool2.fastq: 292 reads
...
```

This creates `Pool1.fastq` and `Pool2.fastq` (uncompressed, plain FASTQ).
The `dame` Rust binary also accepts gzip-compressed input transparently —
you can pass `Pool1.fastq.gz` directly with no extra flags.

### Read Structure

Each valid read is structured as:

```
[fwd_tag_seq][fwd_primer][amplicon][rc(rev_primer)][rc(rev_tag_seq)]
```

For example, a Sample1 Rep1 read (tag1 + tag2) with primer CO1 (GCATGC / CTGACT):

```
AACCGGT  GCATGC  <amplicon>  AGTCAG  TGGCCAA
^tag1    ^fwdP   ^barcode    ^rc(R)  ^rc(tag2)
```

DAMe sort also handles the reverse-complement orientation automatically.

The primer `GCRTGC` contains an IUPAC ambiguity code (R = A or G). The generator
alternates between the two resolved forms (`GCATGC` and `GCGTGC`) to verify that
DAMe's regex-based primer matching handles ambiguity correctly.

---

## Section 4: Step 1 — Sort

The sort step demultiplexes reads by tag combination, strips the tags and primers,
and collapses identical amplicon sequences with their counts.

```bash
mkdir -p pool1 pool2

# Sort Pool 1
dame sort \
    --fq Pool1.fastq \
    --primers Primers.txt \
    --tags Tags.txt
mv tag*.txt SummaryCounts.txt pool1/

# Sort Pool 2
dame sort \
    --fq Pool2.fastq \
    --primers Primers.txt \
    --tags Tags.txt
mv tag*.txt SummaryCounts.txt pool2/
```

### Sort Output Files

**`tag1_tag2.txt`** — one file per tag combination found in the data.

Format: `PrimerName TAB Tag1 TAB Tag2 TAB Count TAB Sequence`

```
CO1	tag1	tag2	50	ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
CO1	tag1	tag2	50	AAATTTCCCGGG
CO1	tag1	tag2	30	TTAATTAATTAATTAATTAATTAATTAATTAATTAATTAATTAATTAATTAATTAATTAA
CO1	tag1	tag2	1	GCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGC
```

**`SummaryCounts.txt`** — overview of how many unique sequences and total reads per tag pair.

```
#tagName1	tagName2	NumUniqSeqs	SumTotalFreq
tag1	tag6	1	100
tag1	tag2	4	131
tag5	tag6	3	61
```

Note that `tag1_tag6.txt` appears in Pool1's summary: these are the 100 "wrong pair"
noise reads (tag1 fwd, tag6 rev), which get sorted to their own file but are
not referenced in `PSinfo.txt`.

DAMe also prints the count of reads that could not be assigned to any combination:

```
Number of erroneous sequences (with errors in the sequence of primer or tags, or no barcode amplified): 100
```

This corresponds to the 100 completely random no-primer reads in Pool1.

---

## Section 5: Step 2 — Filter

The filter step reads the PSinfo file to know which tag-combination files belong to
which samples and replicates, then applies presence, count, and length thresholds.

```bash
dame filter \
    --ps-info PSinfo.txt \
    --x 2 \
    --y 2 \
    --t 2 \
    --l 50
```

### Filter Parameters

| Flag        | Meaning                                                             | Tutorial value |
|-------------|---------------------------------------------------------------------|----------------|
| `--ps-info` | PSinfo file mapping samples to tag combinations and pools           | PSinfo.txt     |
| `--x`       | Number of PCR replicates per sample                                 | 2              |
| `--y`       | Minimum number of replicates a sequence must appear in              | 2              |
| `--t`       | Minimum read count in any given replicate                           | 2              |
| `--l`       | Minimum amplicon length (nucleotides)                               | 50             |

### PSinfo.txt Format

```
SampleName TAB FwdTagName TAB RevTagName TAB PoolNum
```

```
Sample1	tag1	tag2	1
Sample1	tag3	tag4	2
Sample2	tag5	tag6	1
Sample2	tag7	tag8	2
```

The filter step automatically generates `PS1_files.txt` and `PS2_files.txt`, each
listing the sorted-output files for replicate 1 and replicate 2 respectively:

```
# PS1_files.txt
pool1/tag1_tag2.txt
pool1/tag5_tag6.txt

# PS2_files.txt
pool2/tag3_tag4.txt
pool2/tag7_tag8.txt
```

### Filter Output Files

The filter step produces seven output files:

| File                                              | Contents                                                          |
|---------------------------------------------------|-------------------------------------------------------------------|
| `Comparisons_2PCRs.txt`                           | All sequences seen in any replicate, with per-replicate counts    |
| `Comparisons_2PCRs.fasta`                         | Same, in FASTA format                                             |
| `Comparisons_2outOf2PCRs.txt`                     | Sequences passing `--y` (present in ≥ y replicates)              |
| `FilteredReads_atLeast2.fasta`                    | Same, FASTA                                                       |
| `Comparisons_2outOf2PCRs.countsThreshold2.txt`    | Sequences passing both `--y` and `--t`                           |
| `FilteredReads_atLeast2.threshold.fasta`          | Same, FASTA                                                       |
| `FilteredReads.fna`                               | Final filtered reads: passes `--y`, `--t`, AND `--l`             |

### Comparisons_2PCRs.txt Format

```
SampleName TAB Rep1Tags TAB Rep1Count TAB Rep2Tags TAB Rep2Count TAB Sequence
```

The actual output for the tutorial dataset:

```
Sample1	tag1-tag2	50	tag3-tag4	45	AAATTTCCCGGG
Sample1	tag1-tag2	50	tag3-tag4	45	ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
Sample1	tag1-tag2	1	tag3-tag4	1	GCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGC
Sample1	tag1-tag2	30	tag3-tag4	0	TTAATTAATTAATTAATTAATTAATTAATTAATTAATTAATTAATTAATTAATTAATTAA
Sample2	tag5-tag6	30	tag7-tag8	20	AAATTTCCCGGG
Sample2	tag5-tag6	30	tag7-tag8	20	ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
Sample2	tag5-tag6	1	tag7-tag8	1	GCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGC
```

### Which Amplicons Pass Each Stage

**After `--y 2`** (sequences present in both replicates):
- AMP_FAILS_PROP drops out (count=30 in Rep1, count=0 in Rep2 for Sample1)
- AMP_FAILS_LENGTH remains (both reps have count >= 1)
- AMP_FAILS_COUNT remains (present in both reps with count=1)

**After `--t 2`** (sequences with count >= 2 in at least `y` replicates):
- AMP_FAILS_COUNT drops out (count=1 in both reps, below threshold)

**After `--l 50`** (sequences with length >= 50):
- AMP_FAILS_LENGTH drops out (12 nt < 50 nt)

### FilteredReads.fna (Final Result)

The final filtered output, containing only amplicons that pass all three criteria:

```
>Sample1	tag1-tag2.tag3-tag4_2		50_45
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
>Sample2	tag5-tag6.tag7-tag8_2		30_20
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
```

Only `AMP_PASSES_ALL` survives for both samples.

The FASTA header format is:

```
>SampleName TAB Rep1Tags.Rep2Tags_SeqID TAB TAB Rep1Count_Rep2Count
```

---

## Section 6: Step 3 — RSI

The Renkonen Similarity Index measures how compositionally similar two PCR replicates
are. A value of 0 means the replicates are identical in composition; a value of 1
means no overlap at all.

```bash
dame rsi Comparisons_2PCRs.txt
```

Output (`RSI_output.txt`):

```
Sample	RSI
Sample1	0.2290076335877862
Sample2	0.007996801279488208
```

### Interpreting RSI Values

**Sample2 RSI ≈ 0.008** — very low, indicating highly similar replicates. Sample2 has
a clean composition: one dominant amplicon (AMP_PASSES_ALL) and one rare one
(AMP_FAILS_COUNT) in proportionally similar counts across both replicates.

**Sample1 RSI ≈ 0.229** — higher, because Sample1 has AMP_FAILS_PROP present in Rep1
but absent from Rep2, creating a stronger imbalance between the replicates.

For pairwise RSI between specific replicates (useful with > 2 replicates):

```bash
dame rsi --explicit Comparisons_2PCRs.txt
```

Output format:

```
Sample	ReplicateA	ReplicateB	RSI
```

---

## Section 7: Step 4 — Decollapse (Optional)

The sort step collapses identical reads into a single entry with a count. Decollapse
reverses this: it expands each unique sequence back to one FASTA record per original
read count.

```bash
dame decollapse --input pool1/tag1_tag2.txt --out-fas tag1_tag2_decollapsed.fasta
```

Example output (first few records):

```
>tag1.tag2.50_1
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
>tag1.tag2.50_2
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
...
```

The header format is `>Tag1.Tag2.OriginalCount_RecordIndex`.

This is useful when downstream tools (e.g., OTU clustering) expect one FASTA record
per read rather than a collapsed representation.

---

## Section 8: dame-py Equivalents

All steps above work identically with `dame-py`. The flag names differ slightly:

### Sort

```bash
# dame (Rust)
dame sort --fq Pool1.fastq --primers Primers.txt --tags Tags.txt

# dame-py (Python)
dame-py sort -fq Pool1.fastq -p Primers.txt -t Tags.txt
```

### Filter

```bash
# dame (Rust)
dame filter --ps-info PSinfo.txt --x 2 --y 2 --t 2 --l 50

# dame-py (Python)
dame-py filter -psInfo PSinfo.txt -x 2 -y 2 -t 2 -l 50
```

### RSI

```bash
# dame (Rust)
dame rsi Comparisons_2PCRs.txt
dame rsi --explicit Comparisons_2PCRs.txt

# dame-py (Python)
dame-py rsi Comparisons_2PCRs.txt
dame-py rsi -e Comparisons_2PCRs.txt
```

### Decollapse

```bash
# dame (Rust)
dame decollapse --input pool1/tag1_tag2.txt --out-fas decollapsed.fasta

# dame-py (Python)
dame-py decollapse -input pool1/tag1_tag2.txt -outFas decollapsed.fasta
```

---

## Section 9: Understanding Sort Output in Detail

Each `tagA_tagB.txt` file has five tab-separated columns:

| Column | Description                                          |
|--------|------------------------------------------------------|
| 1      | Primer name (from Primers.txt, e.g., `CO1`)          |
| 2      | Forward tag name (e.g., `tag1`)                      |
| 3      | Reverse tag name (e.g., `tag2`)                      |
| 4      | Read count (number of times this exact sequence appeared) |
| 5      | Amplicon sequence (between primers, primers stripped) |

The sequence is what lies **between** the forward and reverse primers. Tags and primers
are stripped unless `--keep-primers-seq` is passed to `dame sort`.

### What Sort Does Not Output

Sort does NOT output files for reads that fail to match any tag+primer combination.
Those are counted and printed as "erroneous sequences" to stdout. This includes:

- Reads where a primer is found but the flanking tag is unrecognized
- Reads where no primer is found at all
- Reads with partial primer matches

---

## Section 10: Troubleshooting

### "No output files created after sort"

Check that:
1. Tags.txt uses the format `TagSeq TAB TagName` (sequence first, name second)
2. Primers.txt uses the format `Name TAB FwdSeq TAB RevSeq`
3. The FASTQ file contains reads with tags flanking the primers

Run with a small test file and verify the erroneous-sequence count is not 100%.

### "Filter produces empty FilteredReads.fna"

The `FilteredReads.fna` file is the most restrictive output (passes `--y`, `--t`, and
`--l`). Try the less-strict files first:

- `Comparisons_2PCRs.txt` — everything
- `Comparisons_2outOf2PCRs.txt` — passes `--y` only
- `Comparisons_2outOf2PCRs.countsThreshold2.txt` — passes `--y` and `--t`

If those are empty, verify that the `pool{N}/tagX_tagY.txt` paths referenced in
`PS1_files.txt` / `PS2_files.txt` actually exist.

### "sort output files are missing some tag combinations"

Only tag combinations that actually appear in the FASTQ data get output files. If
`tag3_tag4.txt` is missing, no reads in that FASTQ were assigned to that pair.
Check the SummaryCounts.txt to see what was found.

### "RSI returns 'no replicates in the file'"

RSI requires at least 2 replicates (i.e., `--x 2` or higher during filter). If you
used `--x 1`, the Comparisons file only has one count column and RSI cannot be
computed.

### Chimera Checking (Advanced)

Between sort and filter, you can run chimera detection:

```bash
dame chimera --input tag1_tag2.txt --abskew 1.9
# produces tag1_tag2.noChim.txt
```

Then re-run filter with `--chimera-checked`:

```bash
dame filter --ps-info PSinfo.txt --x 2 --y 2 --t 2 --l 50 --chimera-checked
```

When `--chimera-checked` is set, filter looks for `tagX_tagY_pool.noChim.txt` files
instead of `pool{N}/tagX_tagY.txt`.

---

## Quick Reference Card

```bash
# 1. Generate tutorial data (one-time setup)
python generate_tutorial_data.py

# 2. Sort
mkdir -p pool1 pool2
dame sort --fq Pool1.fastq --primers Primers.txt --tags Tags.txt
mv tag*.txt SummaryCounts.txt pool1/
dame sort --fq Pool2.fastq --primers Primers.txt --tags Tags.txt
mv tag*.txt SummaryCounts.txt pool2/

# 3. Filter
dame filter --ps-info PSinfo.txt --x 2 --y 2 --t 2 --l 50

# 4. RSI
dame rsi Comparisons_2PCRs.txt

# 5. (Optional) Decollapse
dame decollapse --input pool1/tag1_tag2.txt --out-fas tag1_tag2_decollapsed.fasta
```

### File Format Summary

| File         | Format                                        |
|--------------|-----------------------------------------------|
| Primers.txt  | `Name TAB FwdSeq TAB RevSeq`                  |
| Tags.txt     | `TagSeq TAB TagName`                          |
| PSinfo.txt   | `SampleName TAB FwdTag TAB RevTag TAB PoolNum`|
| tagA_tagB.txt| `PrimerName TAB Tag1 TAB Tag2 TAB Count TAB Seq` |
