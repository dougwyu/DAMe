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
dame sort       -> demultiplex reads by tag combination, collapse to unique sequences
dame filter     -> apply presence (y), count (t), and length (l) thresholds across replicates
dame convert    -> reformat FilteredReads.fna for USEARCH or sumaclust clustering
dame rsi        -> compute Renkonen Similarity Index between PCR replicates
dame decollapse -> expand collapsed sequences back to one record per original read
```

### Prerequisites

- Python 3.6+ (for `generate_tutorial_data.py`)
- **Either** the `dame` Rust binary (recommended for speed) **or** `dame-py` (Python port)

`dame` is ~5–8× faster than `dame-py`. Both accept plain and gzip-compressed
FASTQ (`.fastq.gz`) input with no extra flags.

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
Both `dame` and `dame-py` also accept gzip-compressed input transparently —
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
mv tag*_*.txt SummaryCounts.txt pool1/

# Sort Pool 2
dame sort \
    --fq Pool2.fastq \
    --primers Primers.txt \
    --tags Tags.txt
mv tag*_*.txt SummaryCounts.txt pool2/
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
It does not realign similar sequences or apply the primer and tag mismatch
allowances: amplicon sequences emitted by `sort` are compared exactly.

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

## Section 6: Step 2b — Convert for Clustering

Convert the filtered reads for downstream clustering with USEARCH or sumaclust:

```bash
# Sumaclust format (default):
dame-py convert -i FilteredReads.fna
# → FilteredReads.forsumaclust.fna

# USEARCH format (adds ;size= tag):
dame-py convert -i FilteredReads.fna -u

# USEARCH with fixed-length N-padding (useful for some USEARCH versions):
dame-py convert -i FilteredReads.fna -u --max-length 313

# Per-sample fastas (creates SampleFastas/ directory):
dame-py convert -i FilteredReads.fna -s
```

| Flag | Description |
|------|-------------|
| `-i` / `--in-fasta` | Input `FilteredReads.fna` (required) |
| `-u` / `--usearch` | USEARCH output (`>Sample;size=N`); default is sumaclust (`>Sample:N count=N`) |
| `--min-length N` | Drop sequences shorter than N |
| `--max-length N` | Drop sequences longer than N; pad to N in USEARCH mode |
| `-s` / `--sample-fastas` | Write per-sample fastas to `SampleFastas/` |
| `--per-pcr` | One FASTA record per sequence per PCR, plus `PCRinfo.txt`; see [Per-PCR OTU tables](#per-pcr-otu-tables-for-occupancy-and-detection-models) |
| `--ps-info FILE` | With `--per-pcr`: `PSinfo.txt`, so samples with no reads still get `PCRinfo.txt` rows; see the same section |

`dame-py` also accepts the original v1.0 spellings: `--inFasta`, `-lmin`, `-lmax`, `--sampleFastas`.

## Per-PCR OTU tables for occupancy and detection models

Occupancy and detection models such as occPlus and occJSDM need one row per PCR reaction, and
they estimate false positives and false negatives themselves. The usual DAMe route gets in their
way twice:

1. **Filtered cells become zeros.** `filter --y 2` keeps a sequence in a sample only if it is in
   at least 2 of that sample's PCRs. A real sequence seen once, in one PCR, is dropped from that
   sample even though it passes elsewhere, so its OTU column survives but the cell is 0.
2. **PCRs are summed.** `convert` adds the replicate counts into one record per sample.

The fix is to keep DAMe's decision about *which* sequences are real, and rebuild the *counts* by
mapping every read, PCR by PCR, onto those sequences. [Ji et al. (2025)](https://doi.org/10.1111/ele.70302)
built their occPlus table this way (see Reference below).

### Example

`tests/fixtures/perpcr/` holds a small synthetic data set (4 samples x 3 PCRs, a 120-bp marker).
After `filter --x 3 --y 2 --t 1 --l 100`, the usual route gives:

```
sample  OTU_A  OTU_B  OTU_C
S1        146      0      0
S2          0     55      0
S3          0      0     14
```

In this table, B is missing from S1 although it was seen once, in S1 PCR2; and C and A, each seen once in S2, are missing from S2.

The per-PCR route below gives one row per PCR, keeps those detections, and lists every PCR,
including S3 PCR2 (no reads) and sample S4 (no reads in any PCR):

```
pcr_id   sample pcr tag_pair pool reads_pre_mapping OTU_seq1 OTU_seq2 OTU_seq3
S1_PCR1  S1     1   t1-t2    1    57                55       0        0
S1_PCR2  S1     2   t3-t4    1    47                46       1        0
S1_PCR3  S1     3   t5-t6    2    45                45       0        0
S2_PCR1  S2     1   t7-t8    1    30                0        30       0
S2_PCR2  S2     2   t9-t10   1    26                1        25       0
S2_PCR3  S2     3   t11-t12  2    4                 0        3        1
S3_PCR1  S3     1   t13-t14  1    8                 0        0        8
S3_PCR2  S3     2   t15-t16  1    0                 0        0        0
S3_PCR3  S3     3   t17-t18  2    8                 0        0        6
S4_PCR1  S4     1   t19-t20  1    0                 0        0        0
S4_PCR2  S4     2   t21-t22  1    0                 0        0        0
S4_PCR3  S4     3   t23-t24  2    0                 0        0        0
```

(OTU_seq1 is A plus its 1-bp variant A2, OTU_seq2 is B, OTU_seq3 is C.)

### Recipe

```bash
# 1. Per-PCR reads from the UNFILTERED comparisons file, plus PCRinfo.txt
dame convert -i Comparisons_3PCRs.fasta --per-pcr --ps-info PSinfo.txt

# 2. Reference: the DAMe-passed sequences, dereplicated
#    (no --max-length here: with -u it pads with N, and padded sequences cannot match exactly;
#     use vsearch --minseqlength/--maxseqlength for a length range)
dame convert -i FilteredReads.fna -u
vsearch --derep_fulllength FilteredReads.forusearch.fna --sizein --sizeout --relabel seq --output passed.fna

# 3. Map: exact, contained matches of at least N bases (N a little below the amplicon length)
vsearch --usearch_global FilteredReads.perpcr.fna --db passed.fna \
    --sizein --id 1.0 --mincols 110 --query_cov 1.0 --otutabout table.tsv

# 4. Cluster the passed sequences into OTUs
vsearch --cluster_size passed.fna --sizein --id 0.97 --uc clusters.uc

# 5. Sum sequences within OTUs, one row per PCR, all PCRs listed
Rscript perpcr_to_occupancy.R table.tsv clusters.uc PCRinfo.txt survey.tsv
```

`table.tsv` has one row per passed sequence, so every step can be inspected and other
clusterings (or LULU-style curation) can be tried without remapping.
`perpcr_to_occupancy.R` is in this directory. Its output's first columns (`pcr_id` to
`reads_pre_mapping`) are the per-PCR covariates; the `OTU_*` columns are the read counts.

### Why these mapping options

Checked with vsearch 2.31 and 2.32:

- `--id 1.0` on its own is not an exact-match test: vsearch ignores end gaps when computing
  identity, so any fragment of a reference scores 100% (a 30-bp fragment was counted as the full
  sequence). `--mincols N` sets the shortest match accepted.
- `--query_cov 1.0` makes every base of a read align, so a read must lie inside its reference
  and cannot be longer than it. Without it, a read with extra bases at one end is accepted.
  `--maxqt 1.0` is not a substitute: a read shifted off one end of its reference passes it.
- `--minseqlength` does not filter reads in a search; it removes reference sequences. To drop
  short reads, use `--min-length` on `convert --per-pcr`.
- For markers whose length varies a lot, `--target_cov F` (fraction of the reference aligned)
  is the length-relative alternative to `--mincols`.
- A read can match two passed sequences equally if they differ only beyond its ends; vsearch
  then picks one. This matters only when the two are in different OTUs.

With these options a read counts only if it matches a DAMe-passed sequence exactly. Sequences
that failed in every sample (errors, chimeras) match nothing. Mapping at a lower identity
(`--id 0.97 --mincols N --query_cov 1.0`) lets error copies add their reads to their parent
sequence, at some risk of also absorbing chimeras or reads of rare relatives.

### Notes

- **Use `Comparisons_<X>PCRs.fasta`, not `FilteredReads.fna`.** `FilteredReads.fna` has already
  been filtered by `--y`/`--t`, so converting it per PCR brings the zeros back. `convert` warns
  when the input file name starts with `FilteredReads`.
- **`reads_pre_mapping`** counts every read written for that PCR, including reads that later
  match nothing. It measures sequencing depth; it is not the row total of the mapped table.
- **PCRs with no reads** have no column in `table.tsv`; `PCRinfo.txt` lists them with
  `reads_pre_mapping = 0` and the R script adds them as all-zero rows. DAMe cannot tell a PCR
  that was run and gave no reads from one that was never run; that is for you to decide.
- **Without `--ps-info`**, samples with no reads in any PCR do not appear at all, empty PCRs have
  `tag_pair = empty`, and there is no `pool` column.
- Tag names containing `.`, `-` or `_` break the header format, as for `convert` and `rsi`.

### Reference

Ji, Y., Diana, A., Li, X., Matechou, E., Griffin, J. E., Liu, S., Luo, M., Wu, C., Bai, R., Yao, C., Yin, T., Dong, F., Wu, F., Wang, K., Yu, Z., Chen, X., Jiang, X., Che, J., Yu, D. W., & Popescu, V. D. (2025). High Quality, Granular, Timely, Trustworthy and Efficient Vertebrate Species Distribution Data Across a 30,000 km² Protected Area Complex. *Ecology Letters*, 28(12), e70302. <https://doi.org/10.1111/ele.70302>

---

## Section 7: Step 3 — RSI

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

## Section 8: Step 4 — Decollapse (Optional)

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

## Section 9: dame-py Equivalents

All steps above work identically with `dame-py`. The flag names differ slightly:

### Sort

```bash
# dame (Rust)
dame sort --fq Pool1.fastq --primers Primers.txt --tags Tags.txt
dame sort --fq Pool1.fastq --primers Primers.txt --tags Tags.txt --primer-mismatches 1
dame sort --fq Pool1.fastq --primers Primers.txt --tags Tags.txt --tag-mismatches 1

# dame-py (Python) — accepts both --long and single-dash flags
dame-py sort -fq Pool1.fastq -p Primers.txt -t Tags.txt
dame-py sort -fq Pool1.fastq -p Primers.txt -t Tags.txt -m 1
dame-py sort -fq Pool1.fastq -p Primers.txt -t Tags.txt -mt 1
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

### Convert

```bash
# dame (Rust)
dame convert -i FilteredReads.fna
dame convert -i FilteredReads.fna -u
dame convert -i FilteredReads.fna -u --max-length 313
dame convert -i FilteredReads.fna -s
dame convert -i Comparisons_3PCRs.fasta --per-pcr --ps-info PSinfo.txt

# dame-py (Python) — also accepts v1.0 spellings: --inFasta, -lmin, -lmax, --sampleFastas
dame-py convert -i FilteredReads.fna
dame-py convert -i FilteredReads.fna -u
dame-py convert -i FilteredReads.fna -u --max-length 313
dame-py convert -i FilteredReads.fna -s
dame-py convert -i Comparisons_3PCRs.fasta --per-pcr --ps-info PSinfo.txt
```

For per-PCR output, `dame` takes `--per-pcr` and `--ps-info`; `dame-py` also accepts the single-dash aliases `-perPCR` and `-psInfo`. See the per-PCR section above.

### Decollapse

```bash
# dame (Rust)
dame decollapse --input pool1/tag1_tag2.txt --out-fas decollapsed.fasta

# dame-py (Python)
dame-py decollapse -input pool1/tag1_tag2.txt -outFas decollapsed.fasta
```

---

## Section 10: Understanding Sort Output in Detail

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

### Sort Options

| Flag (Rust) | Flag (Python) | Meaning | Default |
|---|---|---|---|
| `--keep-primers-seq` | `--keepPrimersSeq` | Retain primer sequences in output instead of stripping them | off |
| `--primer-mismatches N` | `-m N` | Allow up to N substitutions per primer match (IUPAC-aware); primer indels are not tolerated. | 0 |
| `--tag-mismatches N` | `-mt N` | Allow up to N substitutions per tag (tag1 and tag2 independently, IUPAC-aware); tag indels are not tolerated. Triggers the anchored matcher; ambiguous reads are discarded. | 0 |

These options use fixed-length, end-anchored comparisons rather than general
sequence alignment. A variable-length amplicon—and an insertion or deletion
within the amplicon—is supported because the second tag and primer are located
relative to the end of the read. An insertion or deletion within a primer or
tag changes an expected fixed-length region, so the matcher normally rejects
the read instead of treating the indel as one mismatch.

### What Sort Does Not Output

Sort does NOT output files for reads that fail to match any tag+primer combination.
Those are counted and printed as "erroneous sequences" to stdout. This includes:

- Reads where a primer is found but the flanking tag is unrecognized
- Reads where no primer is found at all
- Reads with partial primer matches
- Reads whose primer or tag insertion/deletion prevents a valid fixed-window match

---

## Section 11: Troubleshooting

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
mv tag*_*.txt SummaryCounts.txt pool1/
dame sort --fq Pool2.fastq --primers Primers.txt --tags Tags.txt
mv tag*_*.txt SummaryCounts.txt pool2/

# 3. Filter
dame filter --ps-info PSinfo.txt --x 2 --y 2 --t 2 --l 50

# 3b. Convert for clustering (sumaclust default; use -u for USEARCH)
dame convert -i FilteredReads.fna

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
| FilteredReads.perpcr.fna | `>Seq;size=N;sample=<sample>_PCR<k>;` then sequence (convert `--per-pcr`) |
| PCRinfo.txt  | One row per PCR: `pcr_id`, sample, pcr, tag_pair, pool, `reads_pre_mapping` (convert `--per-pcr`) |
