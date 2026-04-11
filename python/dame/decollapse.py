import argparse


def register_subcommand(subparsers):
    p = subparsers.add_parser(
        "decollapse",
        description="Expand unique sequences to individual reads by frequency",
    )
    p.add_argument("-input", required=True,
                   help="Text file with tag combination and freq of each unique seq")
    p.add_argument("-outFas", default="Decollapsed.fasta",
                   help='Output fasta file [default "Decollapsed.fasta"]')
    p.set_defaults(func=run)


def run(args):
    seq_id = 0
    with open(args.input) as IN, open(args.outFas, "w") as OUT:
        for line in IN:
            line = line.rstrip().split()
            count = 1
            while count <= int(line[3]):
                seq_id += 1
                OUT.write(">" + line[1] + "." + line[2] + "." + line[3]
                          + "_" + str(seq_id) + "\n" + line[4] + "\n")
                count += 1


def main():
    parser = argparse.ArgumentParser(
        description="Expand unique sequences to individual reads by frequency"
    )
    parser.add_argument("-input", required=True)
    parser.add_argument("-outFas", default="Decollapsed.fasta")
    run(parser.parse_args())


if __name__ == "__main__":
    main()
