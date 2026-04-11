import sys
import argparse
import numpy as np


def compare(matrix, j, a, b):
    print(f"Comparing replicates {a!r} and {b!r} from sample {j}")
    total = matrix.sum(axis=0)
    if total[0] == 0:
        total[0] = 1
        print(f"This sample gave zero in replicate: {a!r}")
    if total[1] == 0:
        total[1] = 1
        print(f"This sample gave zero in replicate: {b!r}")
    col_a = np.array([row[0] / float(total[0]) for row in matrix])
    col_b = np.array([row[1] / float(total[1]) for row in matrix])
    percent = np.column_stack((col_a, col_b))
    return 1 - np.sum([row.min() for row in percent])


def register_subcommand(subparsers):
    p = subparsers.add_parser(
        "rsi",
        description="Compute Renkonen Similarity Index between PCR replicates",
    )
    p.add_argument("input", help="Input comparison file")
    p.add_argument("-e", "--explicit", action="store_true",
                   help="Output explicit RSI for every pairwise comparison")
    p.add_argument("-o", "--output", dest="outfile", metavar="FILE",
                   help="Write output to FILE [default RSI_output.txt]")
    p.set_defaults(func=run)


def run(args):
    data = []
    with open(args.input) as f:
        for line in f:
            parts = line.split()
            data.append(parts)
    data = np.array(data)
    names = set(row[0] for row in data)
    no_rep = (len(data[0]) - 2) // 2
    if no_rep < 2:
        print("There are no replicates in the file.")
        sys.exit(0)

    rkn = []
    for i in names:
        subset = np.array([row for row in data if row[0] == i])
        if args.explicit:
            for A in range(1, no_rep):
                for B in range(A + 1, no_rep + 1):
                    sample = subset[:, [A * 2, B * 2]].astype(int)
                    output = compare(sample, i, A, B)
                    rkn.append([i, A, B, output])
        else:
            rep = 0
            output = 0
            for A in range(1, no_rep):
                for B in range(A + 1, no_rep + 1):
                    sample = subset[:, [A * 2, B * 2]].astype(int)
                    output += compare(sample, i, A, B)
                    rep += 1
            rkn.append([i, output / rep])

    outfile = args.outfile if args.outfile else "RSI_output.txt"
    with open(outfile, "w") as out:
        if args.explicit:
            out.write("Sample\tReplicateA\tReplicateB\tRSI\n")
            for row in rkn:
                out.write("%s\t%s\t%s\t%s\n" % (row[0], row[1], row[2], row[3]))
        else:
            out.write("Sample\tRSI\n")
            for row in rkn:
                out.write("%s\t%s\n" % (row[0], row[1]))


def main():
    parser = argparse.ArgumentParser(
        description="Compute Renkonen Similarity Index between PCR replicates"
    )
    parser.add_argument("input")
    parser.add_argument("-e", "--explicit", action="store_true")
    parser.add_argument("-o", "--output", dest="outfile", metavar="FILE")
    run(parser.parse_args())


if __name__ == "__main__":
    main()
