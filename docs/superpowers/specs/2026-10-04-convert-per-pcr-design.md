# Spec: `convert --per-pcr` (v3.2.0)

**Date:** 2026-10-04
**Branch:** `convert-per-pcr`
**Status:** Draft, awaiting review

---

## Background

DAMe's job is to decide which sequences are real, so that only real OTUs become columns of the
OTU table. `dame filter` does this by keeping a sequence in a sample only if it has at least
`--t` reads in at least `--y` of the sample's `--x` PCR replicates (`MakeComparisonFile` in
`modules_filter.py` / `make_comparison_file` in `filter.rs`).

That test is applied to each (sample, unique sequence) pair, not to each sequence across the
whole study. A pair that fails is simply not written to `FilteredReads.fna`. The usual route from
there is `dame convert`, which sums the per-PCR counts into one `size=` per sample, then
clustering, with OTU table counts taken from cluster membership. Two things go wrong for
occupancy and detection models (occPlus, occJSDM and similar), which need one row per PCR
reaction and model false positives and false negatives themselves:

1. **Filtered cells become zeros.** A real sequence with 1 read in 1 of 3 PCRs in sample A fails
   `--y 2` in sample A. Its OTU column survives because the sequence passes elsewhere, but the
   cell for sample A is 0, although the sequence was detected there. The model then sees a data
   set already thresholded at "y of x PCRs".
2. **PCR replicates are summed.** `convert` collapses the replicate counts into one record per
   sample, so the PCR structure the models need is lost.

Both are avoidable because `dame filter` also writes `Comparisons_<X>PCRs.fasta`, which holds
every sequence in every sample, unfiltered, with per-PCR counts, in the same record format as
`FilteredReads.fna`.

The fix is to separate deciding which sequences are real from counting them. The filtered reads
still define the set of DAMe-passed sequences, so that decision is unchanged. The counts are then
rebuilt by **mapping** every unfiltered read, per PCR, onto those sequences with the user's own
tool. Ji et al. (2025) built the occPlus OTU table this way. DAMe's part is to emit the per-PCR
reads in a tool-agnostic format.

The recommended recipe (tutorial) is map, then cluster:

1. Dereplicate the DAMe-passed sequences into a reference FASTA (`convert -u` on
   `FilteredReads.fna`, then `vsearch --derep_fulllength`). Do not pass `--max-length` here:
   with `-u` it N-pads the reference, and padded sequences cannot match at 100%. Apply the same
   `--min-length`/`--max-length` limits to both files so they cover the same length range.
2. Map the per-PCR FASTA onto that reference at 100% identity (`vsearch --usearch_global --id 1.0
   --otutabout`). This gives a per-PCR table with one column per passed sequence.
3. Cluster the passed sequences into OTUs with any method, and sum the columns within each OTU.

At 100% identity, a read is counted exactly when its sequence was accepted by DAMe in at least one
sample: the motivating 1-read case is counted, and sequences that failed in every sample (errors,
chimeras) match nothing and are dropped. The per-PCR FASTA contains those failed sequences too,
so a user can instead map at a lower identity to let error variants add their reads to their
parent sequence, at some risk of absorbing chimeras or rare relatives. Clustering first and then
mapping onto OTU representatives also works. All of these use the same per-PCR FASTA.

### Scope

- In: a `--per-pcr` mode for `convert`, a companion `PCRinfo.txt`, optional PSinfo
  cross-checking, tests, docs, version 3.2.0.
- Out: clustering or mapping inside DAMe. Users' tools vary by project, and keeping mapping as an
  explicit external step is more transparent and leaves the identity threshold to the user. Also
  out: `--sample-fastas` in per-PCR mode.

---

## Interface

```
dame convert -i Comparisons_3PCRs.fasta --per-pcr [--ps-info PSinfo.txt] [--min-length N] [--max-length N]
```

Identical in `dame` (Rust) and `dame-py` (Python). `dame-py` also accepts single-dash aliases
consistent with the existing flags (`-perPCR`, `-psInfo`).

| Flag | Meaning |
|------|---------|
| `--per-pcr` | Write each PCR's reads as separate records instead of summing counts per sample. A zero count has no reads, so it gets no record; its zero appears in the mapped OTU table, and a PCR with no reads at all is listed in `PCRinfo.txt` |
| `--ps-info FILE` | Optional. Adds `pool` and real tag pairs for empty PCRs to `PCRinfo.txt`, lists samples with no reads, and cross-checks tag pairs |
| `--min-length`, `--max-length` | Drop per-PCR records outside the length range. No N-padding in per-PCR mode |
| `-u` | No effect in per-PCR mode, which always writes the `size=`/`sample=` label; a note is printed to stderr |
| `-s` / `--sample-fastas` | Rejected in combination with `--per-pcr` |

`--ps-info` without `--per-pcr` is rejected.

Without `--per-pcr`, `convert` behaves exactly as in v3.1.0, byte for byte.

---

## Input

Records as written by `dame filter`, two lines each:

```
>SampleA	t1-t2.t3-t4.t5-t6_42	0_1_0
ACGT...
```

Header tokens (whitespace-separated, as `convert` already parses them):

1. `>` + sample name (PSinfo column 1)
2. tag pairs for PCR 1..X joined by `.`, with `_<record id>` appended to the last; a PCR with no
   reads at all appears as `empty-empty`
3. per-PCR counts joined by `_`, in the same order

X is the number of tag pairs in token 2. PCR numbers are positions in that list, 1 to X, which
is the order `filter` assigns from PSinfo.

The intended input is `Comparisons_<X>PCRs.fasta`. `FilteredReads.fna` is accepted but gives
back the zeros this mode exists to avoid (see Warnings).

---

## Outputs

### Per-PCR FASTA

`FilteredReads.perpcr.fna`, a new name, so a per-sample run's output is never overwritten.

There is one label format. The file is input for mapping, and sumaclust does not map, so a
sumaclust variant would have no consumer.

For each input record that passes the length filters, and for each PCR k with count c > 0, write
one record with the PCR ID `<sample>_PCR<k>`:

```
>SampleA_PCR2.17;size=1;sample=SampleA_PCR2;
```

- The trailing number is a running record counter across the whole output, starting at 1, in
  output order.
- The explicit `sample=` annotation makes vsearch and usearch name the table column
  `SampleA_PCR2` whatever characters the sample name contains. Users of other mappers can split
  the label on `;`.
- Sequences are never N-padded. Padding was for old usearch clustering; in mapping, trailing Ns
  are unmatched positions that lower identity and would reject real reads at `--id 1.0`.
- Output order: input record order, and within a record, PCR 1 to X.

### `PCRinfo.txt`

Tab-separated with a header row, one row per PCR, sorted by sample (in first-seen order: PSinfo
order with `--ps-info`, else input order) then PCR number.

| Column | Meaning |
|--------|---------|
| `pcr_id` | `<sample>_PCR<k>`, matching the FASTA `sample=` and the mapped table's column |
| `sample` | PSinfo sample name |
| `pcr` | k, 1 to X |
| `tag_pair` | `Ftag-Rtag`; `empty` for a PCR with no reads when `--ps-info` is not given |
| `pool` | Only with `--ps-info`: PSinfo column 4 |
| `reads` | Total count written to the FASTA for this PCR, after length filters; 0 for empty PCRs |

Empty PCRs are listed with `reads = 0` so users can add the all-zero rows that mapping cannot
create. DAMe cannot tell "run but yielded no reads" from "never run"; the docs say so and leave
it to the user.

With `--ps-info`, a sample present in PSinfo but absent from the input (no reads in any PCR, so
`filter` wrote nothing for it) gets X rows with `reads = 0`. Without PSinfo such a sample is
invisible; the docs say so.

---

## PSinfo handling

A new read-only helper in the filter module of each implementation turns PSinfo into rows of
(sample, PCR number, tag pair, pool). It assigns PCR numbers with exactly the line-number rule
`makePSnumFiles` / `make_ps_num_files` uses (`NR % X`, where blank and short lines still consume
a slot), so `convert` and `filter` cannot disagree about which PCR is which. X comes from the
input headers. `filter` itself is not changed in this release.

The tag pair from PSinfo is `<col2>-<col3>`, compared with the header's tag pair for every PCR
that is not `empty-empty`.

---

## Warnings and errors

**Warning** (stderr, run continues): with `--per-pcr`, if the input file's basename starts with
`FilteredReads`:

```
Warning: --per-pcr input looks like FilteredReads output, which has already been
filtered by --y/--t, so sequences that failed in a sample will appear as zeros.
For per-PCR OTU tables, use Comparisons_<X>PCRs.fasta instead.
```

**Note** (stderr, run continues): with `--per-pcr -u`:

```
Note: -u has no effect with --per-pcr; per-PCR output always uses the
;size=N;sample=<pcr_id>; label and is never padded.
```

**Errors** (non-zero exit; message text identical in both implementations). Outputs are written
to temporary files in the working directory and renamed only on success, so an error never leaves
a partial FASTA or `PCRinfo.txt` under the final names:

- In a header, the number of counts differs from the number of tag pairs.
- The number of PCRs (X) differs between records.
- One sample appears with different tag-pair lists in different records. This includes a PCR
  that is `empty-empty` in some of a sample's records and a real tag pair in others, which
  cannot happen in genuine `filter` output because emptiness is decided per replicate file.
- With `--ps-info`: a sample in the input is not in PSinfo; a sample has a number of PSinfo rows
  other than X; or a PCR's header tag pair differs from PSinfo at the same position. The message
  names the sample, the PCR number and both tag pairs.
- `--sample-fastas` with `--per-pcr`; `--ps-info` without `--per-pcr`.

Records with fewer than three header tokens or a missing sequence line are skipped, as
`convert` does now. Count parts that are not integers are an error in per-PCR mode (they would
silently shift PCR positions), unlike the per-sample mode, which keeps skipping them.

**Known limitation:** the header format uses `.`, `-` and `_` as separators, so tag names
containing them break parsing. This is already true of `convert` and `rsi`; documented, not
fixed here.

---

## Testing

**Unit tests**, Python (`python/tests/test_convert.py`, `test_filter.py`) and Rust (`convert.rs`,
`filter.rs` test modules), against one small fixture: a three-PCR `Comparisons` FASTA plus a
PSinfo, covering:

- a sequence with 1 read in 1 PCR becomes one record with `size=1` under the right PCR ID;
- zero counts produce no record; PCR numbers follow header position;
- an `empty-empty` PCR, and a PSinfo sample with no records (`reads = 0` rows, real tag pair and
  pool with `--ps-info`, `empty` without);
- length filters on per-PCR records, with no padding even when `--max-length` is given; record
  counter numbering;
- `PCRinfo.txt` with and without `--ps-info`, including sort order;
- the PSinfo helper reproduces `makePSnumFiles` slot assignment, including a blank line;
- the `FilteredReads` warning; the `-u` note; rejected flag combinations;
- existing convert tests unchanged (default mode is byte-identical).

**Malformed inputs** in `tests/fixtures/malformed/`, run by `tests/integration/run_malformed.sh`:
counts/tag-pair length mismatch, inconsistent X, inconsistent tag pairs within a sample,
non-integer count, PSinfo tag mismatch, sample missing from PSinfo. Both implementations must
exit non-zero with the same message.

**Parity:** `tests/integration/run_convert.sh` runs both implementations in per-PCR mode (with
and without `--ps-info`, and with length filters) and compares the FASTA and `PCRinfo.txt` byte for
byte.

**End to end:** `tests/integration/run_pipeline.sh` on the tutorial data runs sort, filter with
`--y 2`, then `convert --per-pcr` on `Comparisons_2PCRs.fasta`, and checks that a (sample,
sequence) pair absent from `FilteredReads.fna` is present in the per-PCR output. If `vsearch` is
on PATH it also runs the tutorial recipe (dereplicate the passed sequences, `--usearch_global --id
1.0 --otutabout`) and checks that the table's column names are a subset of `PCRinfo.txt`'s
`pcr_id` and that the motivating pair has a non-zero cell; skipped otherwise, as the chimera tests
treat `usearch`. If the tutorial data has no pair that fails `--y 2`, the fixture generator
(`tutorial/generate_tutorial_data.py`) gains one.

---

## Documentation

- **README.md:** the new flags and output files in the convert section and the pipeline summary;
  changelog entry 16 for v3.2.0.
- **tutorial/README.md:** new section "Per-PCR OTU tables for occupancy and detection models":
  the zeroing problem with a worked example; the map-then-cluster recipe from Background (`convert
  --per-pcr --ps-info` on `Comparisons`; dereplicate the passed sequences; map with vsearch at
  `--id 1.0`; cluster and sum columns; join to `PCRinfo.txt`; add zero rows for empty PCRs); why
  100% is the default and what a lower identity changes; cluster-then-map as an alternative; the
  `FilteredReads` caveat; and a note that Ji et al. (2025) built
  their occPlus table by mapping, with occJSDM as another consumer.
- **Versions:** `python/pyproject.toml` and `rust/Cargo.toml` to 3.2.0.
- **occJSDM (separate PR, later):** a one-line pointer to the tutorial section from Lesson 0's
  input-data section.
