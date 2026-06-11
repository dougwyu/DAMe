from dame.modules_chimera_check import (
    makeTagFiles, makeTagFilesWithPools, MakeSizeOutFastas,
    SortFasta, MakeFasSeqOneLine, MakeNoChimHaps,
)


def register_subcommand(subparsers):
    p = subparsers.add_parser(
        "chimera",
        description="Create necessary files to operate on sequences per PCR reaction",
    )
    p.add_argument("-psInfo", "--ps-info", dest="psInfo", required=True,
                   help="Text file with tag combination info per PCR reaction per sample")
    p.add_argument("-x", "--x", dest="x", required=True, type=int, help="Number of PCR rxns performed per sample")
    p.add_argument("-p", "--p", dest="p", type=int, default=1, help="Number of pools [default 1]")
    p.set_defaults(func=run)


def run(args):
    PSinfo = args.psInfo
    X = args.x
    P = args.p

    if P == 1:
        makeTagFiles(PSinfo, X)
    else:
        makeTagFilesWithPools(PSinfo, X)

    MakeSizeOutFastas(P, X)
    SortFasta(P)
    MakeFasSeqOneLine(P)
    MakeNoChimHaps(P)
