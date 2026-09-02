from dame.modules_filter import (
    makePSnumFiles, ReadPSnumFiles, MakeSampleNameArray,
    ReadHapsForASample, buildReplicateIndexes, allSequences, MakeComparisonFile,
)


def register_subcommand(subparsers):
    p = subparsers.add_parser(
        "filter",
        description="Filter multiplexed sequences by PCR presence, abundance, and length",
    )
    p.add_argument("-psInfo", "--ps-info", dest="psInfo", required=True,
                   help="Text file with tag combination info per PCR reaction per sample")
    p.add_argument("-x", "--x", dest="x", type=int, default=2, help="Number of PCR rxns performed per sample")
    p.add_argument("-y", "--y", dest="y", type=int, default=1, help="Number of PCR rxns sequence must be present in")
    p.add_argument("-p", "--p", dest="p", type=int, default=1,
                   help="Accepted for backward compatibility but ignored; the pool "
                        "number is read from column 4 of PSinfo [default 1]")
    p.add_argument("-t", "--t", dest="t", type=int, default=1, help="Minimum count per unique sequence")
    p.add_argument("-l", "--l", dest="l", type=int, default=100, help="Minimum sequence length")
    p.add_argument("--chimeraChecked", "--chimera-checked", dest="chimeraChecked", action="store_true",
                   help="Use chimera-checked sorted collapsed files [default not set]")
    p.set_defaults(func=run)


def run(args):
    PSinfo = args.psInfo
    X = args.x
    Y = args.y
    P = args.p
    T = args.t
    L = args.l
    chimeraChecked = args.chimeraChecked

    OUT = open("Comparisons_%sPCRs.txt" % X, "w")
    OUTYX = open("Comparisons_%soutOf%sPCRs.txt" % (Y, X), "w")
    OUTthresh = open("Comparisons_%soutOf%sPCRs.countsThreshold%s.txt" % (Y, X, T), "w")
    OUT_fas = open("Comparisons_%sPCRs.fasta" % X, "w")
    OUTYX_fas = open("FilteredReads_atLeast%s.fasta" % Y, "w")
    OUTthresh_fas = open("FilteredReads_atLeast%s.threshold.fasta" % Y, "w")
    OUTthreshLen_fas = open("FilteredReads.fna", "w")

    makePSnumFiles(PSinfo, X, P, chimeraChecked)
    PSinsLines = ReadPSnumFiles(X)
    sampleName = MakeSampleNameArray(PSinfo)

    for i in range(len(PSinsLines["0"])):
        haps = ReadHapsForASample(X, PSinsLines, i)
        indexes = buildReplicateIndexes(X, haps)
        seqsALL = allSequences(indexes)
        MakeComparisonFile(X, seqsALL, haps, indexes,
                           OUT, OUTthresh, OUTYX, OUT_fas, OUTthresh_fas,
                           OUTYX_fas, OUTthreshLen_fas, Y, T, L, sampleName, i)

    for fh in [OUT, OUTYX, OUTthresh, OUT_fas, OUTYX_fas, OUTthresh_fas, OUTthreshLen_fas]:
        fh.close()
