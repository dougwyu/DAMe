import sys

from dame.modules_sort import (
    readTags, readPrimers, GetPiecesInfo, GetPiecesInfoMismatch, FillHAP,
    PrintSortedCollapsedCountedSeqs, PrintSummaryFile, min_equal_length_tag_distance,
    compile_primer_regexes,
)
from dame.utils import smart_open

def register_subcommand(subparsers):
    p = subparsers.add_parser(
        "sort",
        description="Sort amplicon sequences tagged on each end by tag combination",
    )
    p.add_argument("-fq", "--fq", dest="fq", required=True,
                   help="Input fastq with amplicon sequences")
    p.add_argument("-p", "--primers", dest="p", required=True,
                   help="Input text file with primer name and sequences [Format: Name\\tForwardSeq\\tReverseSeq]")
    p.add_argument("-t", "--tags", dest="t", required=True,
                   help="Input text file with tag names and sequences [Format: TagSeq\\tTagName]")
    p.add_argument("--keepPrimersSeq", "--keep-primers-seq", dest="keepPrimersSeq", action="store_true",
                   help="Keep primer sequences instead of trimming them [default not set]")
    p.add_argument("-m", "--primer-mismatches", dest="primer_mismatches", type=int, default=0,
                   help="Max substitutions allowed per primer match [default 0]")
    p.add_argument("-mt", "--tag-mismatches", dest="tag_mismatches", type=int, default=0,
                   help="Max substitutions allowed per tag (tag1 and tag2 independently) [default 0]")
    p.set_defaults(func=run)


def run(args):
    TAGS = {}
    PRIMERS = {}
    HAP = {}
    CountErrors = 0

    TAGS = readTags(args.t, TAGS)
    PRIMERS = readPrimers(args.p, PRIMERS)

    if args.tag_mismatches > 0:
        min_d = min_equal_length_tag_distance(TAGS)
        if min_d is not None and 2 * args.tag_mismatches + 1 > min_d:
            safe = (min_d - 1) // 2
            print(f"WARNING: --tag-mismatches {args.tag_mismatches} may misassign reads: "
                  f"minimum tag Hamming distance is {min_d}, so the safe maximum is {safe}.",
                  file=sys.stderr)
    use_anchored = args.primer_mismatches > 0 or args.tag_mismatches > 0
    # Exact path uses the fast compiled-regex matcher; build patterns once.
    PRIMERS_RX = None if use_anchored else compile_primer_regexes(PRIMERS)

    with smart_open(args.fq) as f:
        line = f.readline()  # header line
        while line:
            line = f.readline().rstrip()  # seq line
            if not line:
                break
            if use_anchored:
                Info = GetPiecesInfoMismatch(line, PRIMERS, TAGS, args.keepPrimersSeq,
                                             args.primer_mismatches, args.tag_mismatches)
            else:
                Info = GetPiecesInfo(line, PRIMERS_RX, TAGS, args.keepPrimersSeq)
            if len(Info) == 1:
                f.readline()  # "+" line
                f.readline()  # qual line
                line = f.readline()  # next header
                CountErrors += 1
            else:
                HAP = FillHAP(HAP, Info[0], Info[1], Info[2], Info[3])
                f.readline()  # "+" line
                f.readline()  # qual line
                line = f.readline()  # next header

    PrintSortedCollapsedCountedSeqs(HAP)
    PrintSummaryFile(HAP)
    print(f"Number of erroneous sequences (with errors in the sequence of primer or tags, "
          f"or no barcode amplified): {CountErrors}")
