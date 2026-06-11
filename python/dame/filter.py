import argparse
import os
from dame.modules_filter import (
    makePSnumFiles, ReadPSnumFiles, MakeSampleNameArray,
    ReadHapsForASample, getSeqsSetsAndFRcounts, MakeComparisonFile,
)


def _auto_outdir(psinfo_path, Y, T):
    """Build a default output directory name from filter parameters.

    Pattern: Filter_min{Y}PCRs_min{T}copies_{marker}
    where marker is derived from the psInfo filename, e.g.
    PCRsetsInfo_MIFISH.txt -> MIFISH
    """
    basename = os.path.basename(psinfo_path)
    # Strip known extensions
    for ext in ('.txt.gz', '.txt', '.gz'):
        if basename.endswith(ext):
            basename = basename[:-len(ext)]
            break
    # Strip leading "PCRsetsInfo_" prefix if present
    prefix = 'PCRsetsInfo_'
    if basename.startswith(prefix):
        marker = basename[len(prefix):]
    else:
        marker = basename
    return "Filter_min%sPCRs_min%scopies_%s" % (Y, T, marker)


def register_subcommand(subparsers):
    p = subparsers.add_parser(
        "filter",
        description="Filter multiplexed sequences by PCR presence, abundance, and length",
    )
    p.add_argument("-psInfo", required=True,
                   help="PSinfo file: tab-separated columns (SampleName, FwdTagName, "
                        "RevTagName, PoolNumber), one PCR reaction per line (plain or gzip)")
    p.add_argument("-x", type=int, default=2,
                   help="Number of PCR replicates per sample [default 2]")
    p.add_argument("-y", type=int, default=1,
                   help="Minimum number of replicates a sequence must appear in to pass [default 1]")
    p.add_argument("-p", type=int, default=1,
                   help="Number of sequencing pools [default 1]")
    p.add_argument("-t", type=int, default=1,
                   help="Minimum read count per sequence per replicate; counts below this "
                        "threshold increment the below-threshold counter used with -y [default 1]")
    p.add_argument("-l", type=int, default=100,
                   help="Minimum sequence length in nucleotides "
                        "(applied to FilteredReads.fna only) [default 100]")
    p.add_argument("--chimeraChecked", action="store_true",
                   help="Use chimera-checked input files ({tag1}_{tag2}_{pool}.noChim.txt) "
                        "instead of the default sort output (pool{n}/{tag1}_{tag2}.txt)")
    p.add_argument("-o", dest="outdir", default=None,
                   help="Output directory; auto-named Filter_min{Y}PCRs_min{T}copies_{marker} "
                        "if omitted, where marker is derived from the -psInfo filename")
    p.set_defaults(func=run)


def run(args):
    PSinfo = args.psInfo
    X = args.x
    Y = args.y
    P = args.p
    T = args.t
    L = args.l
    chimeraChecked = args.chimeraChecked

    outdir = args.outdir if args.outdir else _auto_outdir(PSinfo, Y, T)
    os.makedirs(outdir, exist_ok=True)

    def out(filename):
        return open(os.path.join(outdir, filename), "w")

    OUT          = out("Comparisons_%sPCRs.txt" % X)
    OUTYX        = out("Comparisons_%soutOf%sPCRs.txt" % (Y, X))
    OUTthresh    = out("Comparisons_%soutOf%sPCRs.countsThreshold%s.txt" % (Y, X, T))
    OUT_fas      = out("Comparisons_%sPCRs.fasta" % X)
    OUTYX_fas    = out("FilteredReads_atLeast%s.fasta" % Y)
    OUTthresh_fas    = out("FilteredReads_atLeast%s.threshold.fasta" % Y)
    OUTthreshLen_fas = out("FilteredReads.fna")

    makePSnumFiles(PSinfo, X, P, chimeraChecked, outdir)
    PSinsLines = ReadPSnumFiles(X, outdir)
    sampleName = MakeSampleNameArray(PSinfo)

    for i in range(len(PSinsLines["0"])):
        haps = ReadHapsForASample(X, PSinsLines, i)
        seqsALL, F, R, seqs = getSeqsSetsAndFRcounts(X, haps)
        MakeComparisonFile(X, seqsALL, haps, F, R, seqs,
                           OUT, OUTthresh, OUTYX, OUT_fas, OUTthresh_fas,
                           OUTYX_fas, OUTthreshLen_fas, Y, T, L, sampleName, i)

    for fh in [OUT, OUTYX, OUTthresh, OUT_fas, OUTYX_fas, OUTthresh_fas, OUTthreshLen_fas]:
        fh.close()

    print("Output written to: %s" % outdir)


def main():
    parser = argparse.ArgumentParser(
        description="Filter multiplexed sequences by PCR presence, abundance, and length"
    )
    parser.add_argument("-psInfo", required=True)
    parser.add_argument("-x", type=int, default=2)
    parser.add_argument("-y", type=int, default=1)
    parser.add_argument("-p", type=int, default=1)
    parser.add_argument("-t", type=int, default=1)
    parser.add_argument("-l", type=int, default=100)
    parser.add_argument("--chimeraChecked", action="store_true")
    parser.add_argument("-o", dest="outdir", default=None)
    run(parser.parse_args())


if __name__ == "__main__":
    main()
