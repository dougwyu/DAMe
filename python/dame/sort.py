from dame.modules_sort import (
    readTags, readPrimers, GetPiecesInfo, FillHAP,
    PrintSortedCollapsedCountedSeqs, PrintSummaryFile,
)

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
    p.add_argument("-m", "--primer-mismatches", dest="primer_mismatches", type=int, default=0,
                   help="Max substitutions allowed per primer match [default 0]")
    p.set_defaults(func=run)


def run(args):
    TAGS = {}
    PRIMERS = {}
    HAP = {}
    CountErrors = 0

    TAGS = readTags(args.t, TAGS)
    PRIMERS = readPrimers(args.p, PRIMERS)

    with open(args.fq) as f:
        line = f.readline()  # header line
        while line:
            line = f.readline().rstrip()  # seq line
            if not line:
                break
            Info = GetPiecesInfo(line, PRIMERS, TAGS, args.keepPrimersSeq, args.primer_mismatches)
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
