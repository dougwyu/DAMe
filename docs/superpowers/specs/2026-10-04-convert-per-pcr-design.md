# Spec: `convert --per-pcr` (v3.2.0)

**Date:** 2026-10-04
**Branch:** `convert-per-pcr`
**Status:** Approved, proceed to implementation plan

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

The recommended recipe (tutorial) maps every read onto the passed sequences themselves, then
combines sequences into OTUs:

1. Dereplicate the DAMe-passed sequences into a reference FASTA (`convert -u` on
   `FilteredReads.fna`, then `vsearch --derep_fulllength`). Do not pass `--max-length` to
   `convert` here: with `-u` it N-pads the reference, and padded sequences cannot match exactly.
   Apply the length range at dereplication instead (`--minseqlength`/`--maxseqlength`), matching
   the `--min-length`/`--max-length` given to `convert --per-pcr`.
2. Map the per-PCR FASTA onto the reference: `vsearch --usearch_global perpcr.fna --db
   passed.fna --id 1.0 --mincols N --query_cov 1.0 --otutabout table.tsv`, where N is a little
   below the amplicon length (e.g. 300 for a 313-bp marker). Reference labels must be unique and
   carry no `otu=` annotation, so `table.tsv` has one row per passed sequence and one column per
   PCR.
3. Cluster the reference into OTUs with any method that reports membership (e.g. `vsearch
   --cluster_size --id 0.97 --uc clusters.uc`), and write the membership as a two-column table
   (sequence, OTU).
4. In R or Python, join `table.tsv` to the membership and sum the rows within each OTU. Keeping
   the sequence-level table makes every step inspectable, and lets users try other clusterings or
   LULU-style curation without remapping.

(vsearch can also do step 4 itself: if each reference label carries its OTU as `;otu=OTU1;`,
`--otutabout` sums rows by `otu=`. The tutorial mentions this as a shortcut only, because it hides
the sequence-level table. Renaming every member to the bare OTU name instead gives duplicate FASTA
labels, which break other tools such as `makeblastdb -parse_seqids` and `samtools faidx`.)

Why these mapping options (each checked with vsearch 2.31 and 2.32 against a 313-bp reference):

- `--id 1.0` alone is not an exact-match test. vsearch's default identity (`--iddef 2`) ignores
  terminal gaps, so any fragment of a reference scores 100%; a 30-bp fragment was counted as the
  full sequence. `--mincols N` sets the shortest acceptable match: with N = 300, reads of 313 and
  300 bp (from either end or the middle) were counted, 299 bp and a 300-bp read with one internal
  mismatch were not.
- `--query_cov 1.0` requires every base of the read to be aligned, so with `--id 1.0` a read must
  lie entirely within its reference and cannot be longer than it. vsearch has no `--maxcols`;
  this is the equivalent. Without it, overhanging bases count as terminal gaps and are ignored: a
  313-bp match plus 7 extra bases was accepted. `--maxqt 1.0` (query no longer than target) is not
  enough: a 313-bp read shifted 10 bp off one end of a 313-bp reference passed `--maxqt 1.0` but
  was rejected by `--query_cov 1.0`.
- `--minseqlength` is not a read-length filter in searches: it discards reference sequences, not
  queries. A minimum read length can also be set upstream with `convert --per-pcr --min-length`.
- `--target_cov F` (fraction of the reference aligned) is the length-relative alternative to
  `--mincols`, for markers whose length varies widely; `--target_cov 0.958` matched `--mincols
  300` exactly on the 313-bp test.
- A short read can match two passed sequences exactly when they differ only beyond its ends;
  vsearch then picks one. This is irrelevant when both are in the same OTU and rare otherwise.

With these options, a read is counted when it is an exact, contained match of at least N bases to
a DAMe-passed sequence: the motivating 1-read case is counted, and sequences that failed in every
sample (errors, chimeras) match nothing and are dropped. The per-PCR FASTA contains those failed
sequences too, so a user can instead map at a lower identity (`--id 0.97 --mincols N --query_cov
1.0`) to let error variants add their reads to their parent sequence, at some risk of absorbing
chimeras or rare relatives.

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

Header tokens (whitespace-separated, as `convert` already parses them; real `filter` output has
a doubled tab before the counts, which whitespace splitting absorbs):

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
| `reads_pre_mapping` | Total count written to the FASTA for this PCR, after length filters; 0 for empty PCRs. Includes reads that later match no reference sequence (errors, chimeras, short fragments), so it is a sequencing-depth measure, not the row total of the mapped table |

Empty PCRs are listed with `reads_pre_mapping = 0` so users can add the all-zero rows that mapping cannot
create. DAMe cannot tell "run but yielded no reads" from "never run"; the docs say so and leave
it to the user.

With `--ps-info`, a sample present in PSinfo but absent from the input (no reads in any PCR, so
`filter` wrote nothing for it) gets X rows with `reads_pre_mapping = 0`. Without PSinfo such a sample is
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

**Unit tests**, Python (`python/tests/test_convert.py`, `test_filter.py`, `test_perpcr.py`) and Rust
(`rust/tests/perpcr_test.rs`, `convert_perpcr_cli_test.rs`, `filter_test.rs`; none are in-module), against the `Comparisons_3PCRs.fasta` that `dame filter` produces
from the `tests/fixtures/perpcr/` data set (committed alongside it) plus its `PSinfo.txt`, and
small inline inputs where a case needs one, covering:

- a sequence with 1 read in 1 PCR becomes one record with `size=1` under the right PCR ID;
- zero counts produce no record; PCR numbers follow header position;
- an `empty-empty` PCR, and a PSinfo sample with no records (`reads_pre_mapping = 0` rows, real tag pair and
  pool with `--ps-info`, `empty` without);
- length filters on per-PCR records, with no padding even when `--max-length` is given; record
  counter numbering;
- `PCRinfo.txt` with and without `--ps-info`, including sort order;
- the PSinfo helper reproduces `makePSnumFiles` slot assignment, including a blank line;
- the `FilteredReads` warning; the `-u` note; rejected flag combinations;
- existing convert tests unchanged (default mode is byte-identical).

**Malformed inputs** in `tests/fixtures/malformed/`, run by `tests/integration/run_malformed.sh`,
cover three per-PCR error cases: counts/tag-pair length mismatch, non-integer count, and PSinfo
tag mismatch. Both implementations must exit 1 with the same message and write no output. The
full set of eight error messages (those three plus inconsistent X, inconsistent tag pairs within a
sample, sample missing from PSinfo, wrong number of PSinfo rows for a sample, and empty input) is pinned by
identical case lists in `python/tests/test_perpcr.py` and `rust/tests/perpcr_test.rs`.

**Parity:** `tests/integration/run_perpcr.sh` runs both implementations in per-PCR mode (with and
without `--ps-info`, and with length filters) and compares the FASTA and `PCRinfo.txt` byte for
byte, against the committed expected files where they exist and between the two implementations
otherwise.

**End to end:** a new fixture, `tests/fixtures/perpcr/`, holds a small synthetic data set in
`dame sort` output form (`pool1/`, `pool2/` tag-pair files) plus `PSinfo.txt`, and the script that
generated it (seeded, so it can be regenerated). It was first built and run through `dame filter`,
a prototype of this spec, vsearch 2.31 and the R join on 2026-10-04. Design: 4 samples x 3 PCRs,
a 120-bp marker, filtered with `--x 3 --y 2 --t 1 --l 100`. Sequences:

- A, B, C: three real sequences; A2: a 1-bp variant of A that also passes (same OTU at 97%);
- the motivating case: B with 1 read in 1 PCR of S1 (fails `--y 2` there, passes in S2);
- C with 1 read in S2 and A with 1 read in S2 (fail there, pass elsewhere);
- E: a 2-substitution error copy of B, in S1 only (fails everywhere);
- Bt: B trimmed to 112 bp, inside B (fails in S2, counted by mapping with `--mincols 110`);
- Bs: a 90-bp fragment of B (rejected by `--mincols 110`);
- S3 PCR2 has no reads at all (`empty-empty`), and S4 has no reads in any PCR (absent from
  `Comparisons`, present only through `--ps-info`).

`tests/integration/run_perpcr.sh` runs `filter`, then `convert --per-pcr --ps-info` with both
implementations, and compares `FilteredReads.perpcr.fna` and `PCRinfo.txt` with committed
expected files (15 records; 12 `PCRinfo.txt` rows, four of them `reads_pre_mapping = 0`). If
`vsearch` is on PATH it also runs the recipe (dereplicate the passed sequences, `--usearch_global
--id 1.0 --mincols 110 --query_cov 1.0 --otutabout`, `--cluster_size --id 0.97 --uc`) and checks
the sequence-level table against a committed expected table: rows `seq1`..`seq4` (A, B, C, A2,
numbered by `--derep_fulllength --relabel seq` in order of abundance); B in `S1_PCR2` = 1; B in
`S2_PCR3` = 3 (Bt); no column for `S3_PCR2` or S4. If `Rscript` with dplyr and tidyr is also
available, it runs `tutorial/perpcr_to_occupancy.R` and checks the final 12-row table, including the
summed column `OTU_seq1` (55 in `S1_PCR1`) and the zero rows. Each optional step is skipped with a message
when its tool is missing, as the chimera tests treat `usearch`.

---

## Documentation

- **README.md:** the new flags and output files in the convert section and the pipeline summary;
  changelog entry 16 for v3.2.0.
- **tutorial/README.md:** new section "Per-PCR OTU tables for occupancy and detection models":
  the zeroing problem, worked through on the `tests/fixtures/perpcr/` data set (the usual
  per-sample route beside the per-PCR result, so readers see the recovered detections); the recipe from Background (`convert --per-pcr
  --ps-info` on `Comparisons`; dereplicate the passed sequences; map with
  `--usearch_global --id 1.0 --mincols N --query_cov 1.0` to a sequence-level table; cluster and
  sum rows within OTUs in R; transpose to one row per PCR and join to `PCRinfo.txt`, adding
  all-zero rows for PCRs with no column), as the script `tutorial/perpcr_to_occupancy.R` (shared with the end-to-end test) ending in one row per PCR, ready to split into occJSDM-style `info` and `OTU`. The script strips `;size=` annotations, since vsearch drops them from table row names but keeps them in `.uc` files, and stops if any sequence has no OTU;
  the mapping-option findings above, including why `--id 1.0` needs `--mincols` and
  `--query_cov 1.0` and that `--minseqlength` does not filter reads; the `otu=` shortcut; what a lower
  identity changes; the `FilteredReads` caveat; and a note that Ji et al. (2025) built
  their occPlus table by mapping, with occJSDM as another consumer.
- **Versions:** `python/pyproject.toml` and `rust/Cargo.toml` to 3.2.0.
- **occJSDM (separate PR, later):** a one-line pointer to the tutorial section from Lesson 0's
  input-data section.
