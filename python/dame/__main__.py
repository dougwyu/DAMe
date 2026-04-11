import argparse

from dame import sort, chimera_check, decollapse, rsi
import dame.filter as filter_mod


def main():
    parser = argparse.ArgumentParser(
        prog="dame-py",
        description="DAMe: DNA Metabarcoding pipeline toolkit",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    sort.register_subcommand(subparsers)
    chimera_check.register_subcommand(subparsers)
    filter_mod.register_subcommand(subparsers)
    decollapse.register_subcommand(subparsers)
    rsi.register_subcommand(subparsers)

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
