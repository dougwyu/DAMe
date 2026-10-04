# Per-PCR fixture

Synthetic data set for `convert --per-pcr`, in `dame sort` output form. See the docstring in `make_data.py` for what each sequence and sample exercises.

Regenerate with `python3 make_data.py`, then `dame filter --ps-info PSinfo.txt --x 3 --y 2 --t 1 --l 100`, keeping only `Comparisons_3PCRs.fasta` and `FilteredReads.fna` from the filter output.

`expected/` holds the outputs the tests compare against: the per-PCR FASTA and `PCRinfo.txt` (with and without `--ps-info`), the vsearch 2.31 sequence-level table, and the survey table written by `tutorial/perpcr_to_occupancy.R`.
