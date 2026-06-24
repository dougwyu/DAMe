# python/dame/convert.py
import os


def _parse_fasta(path):
    """Yield (sample, size, sequence) tuples from a FilteredReads.fna file."""
    with open(path) as fh:
        header = None
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                toks = line.split()
                sample = toks[0][1:]
                size = sum(int(x) for x in toks[2].split("_"))
                header = (sample, size)
            elif header is not None:
                yield header[0], header[1], line
                header = None


def convert(in_fasta, min_length=0, max_length=None, usearch=False, sample_fastas=False):
    """
    Convert FilteredReads.fna to USEARCH or sumaclust format.

    Returns the path of the main output file created.
    """
    out_name = "FilteredReads.forusearch.fna" if usearch else "FilteredReads.forsumaclust.fna"

    if sample_fastas:
        os.makedirs("SampleFastas", exist_ok=True)

    sample_handles = {}
    counter = 1

    with open(out_name, "w") as out:
        for sample, size, seq in _parse_fasta(in_fasta):
            if len(seq) < min_length:
                continue
            if max_length is not None and len(seq) > max_length:
                continue

            if usearch:
                hdr = f">{sample};size={size}"
                out_seq = seq.ljust(max_length, "N") if max_length is not None else seq
            else:
                hdr = f">{sample}:{counter} count={size}"
                out_seq = seq
                counter += 1

            out.write(hdr + "\n" + out_seq + "\n")

            if sample_fastas:
                if sample not in sample_handles:
                    sample_handles[sample] = open(
                        f"SampleFastas/{sample}.fixed.fasta", "w"
                    )
                sample_handles[sample].write(hdr + "\n" + out_seq + "\n")

    for fh in sample_handles.values():
        fh.close()

    return out_name


def register_subcommand(subparsers):
    p = subparsers.add_parser(
        "convert",
        description="Convert FilteredReads.fna to USEARCH or sumaclust input format",
    )
    p.add_argument(
        "-i", "--in-fasta", "--inFasta",
        dest="in_fasta", required=True, metavar="FILE",
        help="Input FilteredReads.fna file",
    )
    p.add_argument(
        "--min-length", "-lmin", "--minLength",
        dest="min_length", type=int, default=0, metavar="N",
        help="Drop sequences shorter than N [default 0]",
    )
    p.add_argument(
        "--max-length", "-lmax", "--maxLength",
        dest="max_length", type=int, default=None, metavar="N",
        help="Drop sequences longer than N; pad to N in USEARCH mode",
    )
    p.add_argument(
        "-u", "--usearch",
        dest="usearch", action="store_true",
        help="Write USEARCH output format (default: sumaclust)",
    )
    p.add_argument(
        "-s", "--sample-fastas", "--sampleFastas",
        dest="sample_fastas", action="store_true",
        help="Write per-sample fastas to SampleFastas/",
    )
    p.set_defaults(func=run)


def run(args):
    convert(
        in_fasta=args.in_fasta,
        min_length=args.min_length,
        max_length=args.max_length,
        usearch=args.usearch,
        sample_fastas=args.sample_fastas,
    )
