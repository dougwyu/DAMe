"""convert --per-pcr: one FASTA record per PCR replicate, plus PCRinfo.txt.

Input is Comparisons_<X>PCRs.fasta (or FilteredReads.fna) from dame filter.
Output reads keep their PCR as `sample=<sample>_PCR<k>` so a mapper such as
vsearch --usearch_global --otutabout builds a per-PCR table. The Rust
implementation (rust/src/perpcr.rs) must produce byte-identical output and
identical messages.
"""
import os
import sys

from dame.modules_filter import readPSinfoRows

FASTA_OUT = "FilteredReads.perpcr.fna"
INFO_OUT = "PCRinfo.txt"

FILTERED_READS_WARNING = (
    "Warning: --per-pcr input looks like FilteredReads output, which has already been\n"
    "filtered by --y/--t, so sequences that failed in a sample will appear as zeros.\n"
    "For per-PCR OTU tables, use Comparisons_<X>PCRs.fasta instead."
)

USEARCH_NOTE = (
    "Note: -u has no effect with --per-pcr; per-PCR output always uses the\n"
    ";size=N;sample=<pcr_id>; label and is never padded."
)


class PerPcrError(Exception):
    """An input problem that stops convert --per-pcr. The message is shown to the user."""


def flag_error(per_pcr, ps_info, sample_fastas):
    """Return the message for an invalid flag combination, or None."""
    if per_pcr and sample_fastas:
        return "--sample-fastas cannot be used with --per-pcr"
    if ps_info is not None and not per_pcr:
        return "--ps-info requires --per-pcr"
    return None


def _records(path):
    """Yield (ordinal, sample, tag_pairs, count_parts, seq) for each record.

    ordinal counts header lines from 1. A header with fewer than three tokens,
    or one not followed by a sequence line, is skipped, as in per-sample convert.
    """
    ordinal = 0
    pending = None
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                ordinal += 1
                pending = None
                toks = line.split()
                if len(toks) >= 3:
                    tags = toks[1].rsplit("_", 1)[0].split(".")
                    pending = (ordinal, toks[0][1:], tags, toks[2].split("_"))
            elif pending is not None:
                yield pending + (line.rstrip(),)
                pending = None


def _load_ps_info(ps_info, x):
    """Return {sample: [(tag_pair, pool) for PCR 1..x]} in PSinfo first-seen order."""
    by_sample = {}
    for sample, pcr, tag_pair, pool in readPSinfoRows(ps_info, x):
        rows = by_sample.setdefault(sample, {})
        if pcr in rows:
            raise PerPcrError(
                "sample %s does not have exactly one PSinfo row for each PCR 1..%d" % (sample, x))
        rows[pcr] = (tag_pair, pool)
    out = {}
    for sample, rows in by_sample.items():
        if sorted(rows) != list(range(1, x + 1)):
            raise PerPcrError(
                "sample %s does not have exactly one PSinfo row for each PCR 1..%d" % (sample, x))
        out[sample] = [rows[k] for k in range(1, x + 1)]
    return out


def _write(in_fasta, ps_info, min_length, max_length, fasta_tmp, info_tmp):
    x = None
    ps_rows = None
    samples = {}  # sample -> {"tags": [...], "reads": [...]}, in first-seen order
    counter = 0
    with open(fasta_tmp, "w") as out:
        for ordinal, sample, tags, parts, seq in _records(in_fasta):
            counts = []
            for part in parts:
                try:
                    counts.append(int(part))
                except ValueError:
                    raise PerPcrError("record %d for sample %s has a non-integer count '%s'"
                                      % (ordinal, sample, part))
            if len(counts) != len(tags):
                raise PerPcrError("record %d for sample %s has %d counts but %d tag pairs"
                                  % (ordinal, sample, len(counts), len(tags)))
            if x is None:
                x = len(tags)
                if ps_info is not None:
                    ps_rows = _load_ps_info(ps_info, x)
            elif len(tags) != x:
                raise PerPcrError("record %d for sample %s has %d PCRs but earlier records have %d"
                                  % (ordinal, sample, len(tags), x))

            state = samples.get(sample)
            if state is None:
                if ps_rows is not None:
                    if sample not in ps_rows:
                        raise PerPcrError("sample %s is in the input but not in %s"
                                          % (sample, ps_info))
                    for k, (tag_pair, _pool) in enumerate(ps_rows[sample], start=1):
                        if tags[k - 1] != "empty-empty" and tags[k - 1] != tag_pair:
                            raise PerPcrError(
                                "sample %s PCR %d: tag pair %s in the input but %s in PSinfo"
                                % (sample, k, tags[k - 1], tag_pair))
                state = samples[sample] = {"tags": tags, "reads": [0] * x}
            elif tags != state["tags"]:
                raise PerPcrError(
                    "sample %s has tag pairs %s in record %d but %s in an earlier record"
                    % (sample, ".".join(tags), ordinal, ".".join(state["tags"])))

            if len(seq) < min_length or (max_length is not None and len(seq) > max_length):
                continue
            for k, c in enumerate(counts, start=1):
                if c > 0:
                    counter += 1
                    out.write(">%s_PCR%d.%d;size=%d;sample=%s_PCR%d;\n%s\n"
                              % (sample, k, counter, c, sample, k, seq))
                    state["reads"][k - 1] += c

    if x is None:
        raise PerPcrError("no records in %s" % in_fasta)

    with open(info_tmp, "w") as f:
        if ps_rows is None:
            f.write("pcr_id\tsample\tpcr\ttag_pair\treads_pre_mapping\n")
            for sample, state in samples.items():
                for k in range(1, x + 1):
                    tag_pair = state["tags"][k - 1]
                    if tag_pair == "empty-empty":
                        tag_pair = "empty"
                    f.write("%s_PCR%d\t%s\t%d\t%s\t%d\n"
                            % (sample, k, sample, k, tag_pair, state["reads"][k - 1]))
        else:
            f.write("pcr_id\tsample\tpcr\ttag_pair\tpool\treads_pre_mapping\n")
            for sample, rows in ps_rows.items():
                reads = samples[sample]["reads"] if sample in samples else [0] * x
                for k, (tag_pair, pool) in enumerate(rows, start=1):
                    f.write("%s_PCR%d\t%s\t%d\t%s\t%s\t%d\n"
                            % (sample, k, sample, k, tag_pair, pool, reads[k - 1]))


def per_pcr_convert(in_fasta, ps_info=None, min_length=0, max_length=None,
                    usearch=False, out_dir="."):
    """Write FilteredReads.perpcr.fna and PCRinfo.txt into out_dir.

    Outputs are written under .tmp names and renamed only when the whole input
    has been read without error. Returns the two output paths.
    """
    if os.path.basename(in_fasta).startswith("FilteredReads"):
        print(FILTERED_READS_WARNING, file=sys.stderr)
    if usearch:
        print(USEARCH_NOTE, file=sys.stderr)

    fasta_out = os.path.join(out_dir, FASTA_OUT)
    info_out = os.path.join(out_dir, INFO_OUT)
    fasta_tmp = fasta_out + ".tmp"
    info_tmp = info_out + ".tmp"
    try:
        _write(in_fasta, ps_info, min_length, max_length, fasta_tmp, info_tmp)
    except BaseException:
        for path in (fasta_tmp, info_tmp):
            if os.path.exists(path):
                os.remove(path)
        raise
    os.replace(fasta_tmp, fasta_out)
    os.replace(info_tmp, info_out)
    return fasta_out, info_out
