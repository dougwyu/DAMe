# python/dame/convert.py
import os
import sys

from dame.perpcr import PerPcrError, flag_error, per_pcr_convert


def _sum_counts(field):
    """Sum the underscore-separated per-replicate counts in a header field.

    Unparseable parts are skipped rather than raising, mirroring the Rust
    implementation's filter_map(|x| x.parse().ok()).
    """
    total = 0
    for part in field.split("_"):
        try:
            total += int(part)
        except ValueError:
            continue
    return total


def _parse_fasta(path):
    """Yield (sample, size, sequence) tuples from a FilteredReads.fna file."""
    with open(path) as fh:
        header = None
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                toks = line.split()
                # A header needs at least sample, tag pair and counts. Skip the
                # record instead of raising IndexError, as the Rust
                # implementation does (convert.rs: `if toks.len() < 3`). The
                # sequence line that follows is then ignored, because `header`
                # stays None.
                if len(toks) < 3:
                    header = None
                    continue
                sample = toks[0][1:]
                size = _sum_counts(toks[2])
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
        try:
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
        finally:
            for fh in sample_handles.values():
                fh.close()

    return out_name


def register_subcommand(subparsers):
    p = subparsers.add_parser(
        "convert",
        description="Convert FilteredReads.fna to USEARCH or sumaclust input format, or write per-PCR reads for mapping (--per-pcr)",
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
    p.add_argument(
        "--per-pcr", "-perPCR",
        dest="per_pcr", action="store_true",
        help="Write one record per PCR replicate (for mapping) and PCRinfo.txt, "
             "instead of summing counts per sample. Use Comparisons_<X>PCRs.fasta as input",
    )
    p.add_argument(
        "--ps-info", "-psInfo",
        dest="ps_info", default=None, metavar="FILE",
        help="With --per-pcr: PSinfo file, adding pool and real tag pairs to PCRinfo.txt "
             "and cross-checking tag pairs",
    )
    p.set_defaults(func=run)


def run(args):
    message = flag_error(args.per_pcr, args.ps_info, args.sample_fastas)
    if message is not None:
        sys.exit("Error: " + message)
    if args.per_pcr:
        try:
            per_pcr_convert(
                in_fasta=args.in_fasta,
                ps_info=args.ps_info,
                min_length=args.min_length,
                max_length=args.max_length,
                usearch=args.usearch,
            )
        except PerPcrError as e:
            sys.exit("Error: %s" % e)
        return
    convert(
        in_fasta=args.in_fasta,
        min_length=args.min_length,
        max_length=args.max_length,
        usearch=args.usearch,
        sample_fastas=args.sample_fastas,
    )
