"""Generate the synthetic per-PCR fixture in `dame sort` output form.

Run from this directory: python3 make_data.py
Then: dame filter --ps-info PSinfo.txt --x 3 --y 2 --t 1 --l 100

Sequences (120-bp marker):
  A, B, C  real sequences
  A2       1-bp variant of A that also passes the filter (same OTU at 97%)
  E        2-substitution error copy of B (fails everywhere)
  Bt       B trimmed to 112 bp, inside B (counted by mapping with --mincols 110)
  Bs       90-bp fragment of B (rejected by --mincols 110)
S3 PCR2 gets no reads at all; S4 gets no reads in any PCR.
"""
import os
import random

random.seed(11)


def rs(n):
    return "".join(random.choice("ACGT") for _ in range(n))


def sub(s, i):
    return s[:i] + ("A" if s[i] != "A" else "C") + s[i + 1:]


A = rs(120)
A2 = sub(A, 60)
B = rs(120)
C = rs(120)
E = sub(sub(B, 30), 90)
Bt = B[4:116]
Bs = B[:90]
seqs = dict(A=A, A2=A2, B=B, C=C, E=E, Bt=Bt, Bs=Bs)

# PSinfo: sample, forward tag, reverse tag, pool. 4 samples x 3 PCRs.
ps = []
tag = 1
for s in ["S1", "S2", "S3", "S4"]:
    for k in range(3):
        ps.append((s, f"t{tag}", f"t{tag + 1}", 1 if k < 2 else 2))
        tag += 2
with open("PSinfo.txt", "w") as f:
    f.write("".join("\t".join(map(str, r)) + "\n" for r in ps))

# Reads per sample, per sequence, per PCR
counts = {
    "S1": {"A": (50, 40, 45), "A2": (5, 6, 0), "B": (0, 1, 0), "E": (2, 0, 0)},
    "S2": {"B": (30, 25, 0), "Bt": (0, 0, 3), "C": (0, 0, 1), "A": (0, 1, 0)},
    "S3": {"C": (8, 0, 6), "Bs": (0, 0, 2)},
    "S4": {},
}
for s, ft, rt, pool in ps:
    k = [r for r in ps if r[0] == s].index((s, ft, rt, pool))
    rows = [(n, c[k]) for n, c in counts[s].items() if c[k] > 0]
    if not rows:
        continue  # sort writes no file for a tag pair with no reads
    os.makedirs(f"pool{pool}", exist_ok=True)
    with open(f"pool{pool}/{ft}_{rt}.txt", "w") as f:
        for n, c in rows:
            f.write(f"COI\t{ft}\t{rt}\t{c}\t{seqs[n]}\n")
