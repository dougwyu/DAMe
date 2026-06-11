import argparse
import os
import sys

from dame.utils import smart_open
from dame.modules_sort import (
    readTags, readPrimers, GetPiecesInfo, GetPiecesInfoMismatch, FillHAP,
    PrintSortedCollapsedCountedSeqs, PrintSummaryFile,
    PrintSplitSummaryFile, read_valid_pairs,
)

AMBIG = {
    'A': "A", 'B': "[CGT]", 'C': "C", 'D': "[AGT]", 'G': "G",
    'H': "[ACT]", 'K': "[GT]", 'M': "[AC]", 'N': "[ACGT]", 'R': "[AG]",
    'S': "[CG]", 'T': "T", 'V': "[ACG]", 'W': "[AT]", 'Y': "[CT]",
}


def register_subcommand(subparsers):
    p = subparsers.add_parser(
        "sort",
        description="Sort amplicon sequences tagged on each end by tag combination",
    )
    p.add_argument("-fq", required=True, help="Input fastq with amplicon sequences")
    p.add_argument("-p", required=True,
                   help="Input text file with primer name and sequences [Format: Name\\tForwardSeq\\tReverseSeq]")
    p.add_argument("-t", required=True,
                   help="Input text file with tag names and sequences [Format: TagSeq\\tTagName]")
    p.add_argument("--keepPrimersSeq", action="store_true",
                   help="Keep primer sequences instead of trimming them [default not set]")
    p.add_argument("-psInfo", dest="psinfo", default=None,
                   help="PCRsetsInfo file: tab-separated (SampleName, FwdTagName, RevTagName, "
                        "PoolNumber), same format as filter -psInfo.  Pool number is inferred "
                        "from the current working directory name.  When supplied, writes "
                        "SplitSummary_<pool>.txt categorising reads as valid pair / same-tag "
                        "pair / different-pool pair / one tag only / no tags found.")
    p.add_argument("-m", dest="m", type=int, default=0,
                   help="Max mismatches allowed per primer (F and R each) [default 0]")
    p.add_argument("-mt", dest="mt", type=int, default=0,
                   help="Max mismatches allowed per tag (tag1 and tag2 each) [default 0]")
    p.set_defaults(func=run)


def run(args):
    TAGS = {}
    PRIMERS = {}
    HAP = {}
    HAP_err = {}   # section 3: primer found, exactly one tag matched
    no_tags_seqs = {}  # section 4: no primer or no tags at all
    CountErrors = 0

    TAGS = readTags(args.t, TAGS)
    PRIMERS = readPrimers(args.p, PRIMERS, AMBIG)

    with smart_open(args.fq) as f:
        line = f.readline()  # header line
        while line:
            line = f.readline().rstrip()  # seq line
            if not line:
                break
            if args.m == 0 and args.mt == 0:
                Info = GetPiecesInfo(line, PRIMERS, TAGS, args.keepPrimersSeq)
            else:
                Info = GetPiecesInfoMismatch(line, PRIMERS, TAGS, args.keepPrimersSeq,
                                             args.m, args.mt)
            if Info[0] is None:
                # Primer pair found but not both tags; Info = [None, t1, t2, pname, between]
                t1, t2, pname, between = Info[1], Info[2], Info[3], Info[4]
                if t1 or t2:
                    HAP_err = FillHAP(HAP_err, t1 or 'none', t2 or 'none', pname, between)
                else:
                    no_tags_seqs[between] = no_tags_seqs.get(between, 0) + 1
                CountErrors += 1
            elif len(Info) == 1:
                # No primer found at all
                no_tags_seqs[line] = no_tags_seqs.get(line, 0) + 1
                CountErrors += 1
            else:
                HAP = FillHAP(HAP, Info[0], Info[1], Info[2], Info[3])
            f.readline()  # "+" line
            f.readline()  # qual line
            line = f.readline()  # next header

    # Derive prefix and pool name from input file path
    fq_abs = os.path.abspath(args.fq)
    basename = os.path.basename(fq_abs)
    for ext in ('.fastq.gz', '.fq.gz', '.fastq', '.fq'):
        if basename.endswith(ext):
            basename = basename[:-len(ext)]
            break
    prefix = basename.split('_')[0]
    # Use the current working directory name as the pool identifier so that
    # running "dame-py sort -fq ../file.fq.gz ..." from inside pool1/ correctly
    # yields pool="pool1" even when the fq file lives in the parent directory.
    pool = os.path.basename(os.getcwd())

    # Resolve valid tag pairs from PCRsetsInfo if provided
    valid_pairs = None
    if args.psinfo:
        pool_digits = ''.join(c for c in pool if c.isdigit())
        if pool_digits:
            valid_pairs = read_valid_pairs(args.psinfo, int(pool_digits))

    PrintSortedCollapsedCountedSeqs(HAP)
    PrintSummaryFile(HAP)
    PrintSplitSummaryFile(HAP, HAP_err, no_tags_seqs, prefix, pool, valid_pairs)
    print(f"Number of erroneous sequences (with errors in the sequence of primer or tags, "
          f"or no barcode amplified): {CountErrors}")


def main():
    parser = argparse.ArgumentParser(
        description="Sort amplicon sequences tagged on each end by tag combination"
    )
    parser.add_argument("-fq", required=True)
    parser.add_argument("-p", required=True)
    parser.add_argument("-t", required=True)
    parser.add_argument("--keepPrimersSeq", action="store_true")
    parser.add_argument("-psInfo", dest="psinfo", default=None)
    parser.add_argument("-m", dest="m", type=int, default=0)
    parser.add_argument("-mt", dest="mt", type=int, default=0)
    run(parser.parse_args())


if __name__ == "__main__":
    main()
