import os


def makePSnumFiles(PSinfo, X, P, chimeraChecked):
    PSouts = [open("PS%s_files.txt" % (i + 1), "w") for i in range(X)]
    with open(PSinfo) as f:
        PS = f.readlines()
    for NR, psinfo in enumerate(PS):
        NR = NR + 1
        # Skip blank and short lines rather than raising IndexError, matching
        # the Rust implementation (filter.rs). NR still counts skipped lines,
        # exactly as Rust's enumerate index does, so a blank line consumes a
        # replicate slot in both and the PS file assignment stays in step.
        psinfo = psinfo.strip()
        if not psinfo:
            continue
        psinfo = psinfo.split()
        if len(psinfo) < 4:
            continue
        residue = NR % X
        idx = residue - 1 if residue != 0 else X - 1
        if not chimeraChecked:
            PSouts[idx].write("pool%s/%s_%s.txt\n" % (psinfo[3], psinfo[1], psinfo[2]))
        else:
            PSouts[idx].write("%s_%s_%s.noChim.txt\n" % (psinfo[1], psinfo[2], psinfo[3]))
    for out in PSouts:
        out.close()


def ReadPSnumFiles(X):
    PSinsLines = {}
    for i in range(X):
        with open("PS%s_files.txt" % (i + 1)) as f:
            PSinsLines[str(i)] = f.readlines()
    return PSinsLines


def MakeSampleNameArray(PSinfo):
    sampleName = []
    with open(PSinfo) as f:
        for line in f:
            # Skip blank lines rather than raising IndexError, matching
            # make_sample_name_array in filter.rs.
            line = line.strip()
            if not line:
                continue
            name = line.split()[0]
            if name not in sampleName:
                sampleName.append(name)
    return sampleName


def ReadHapsForASample(X, PSinsLines, i):
    haps = {}
    for j in range(X):
        haps[str(j)] = []
        path = PSinsLines[str(j)][i].rstrip()
        if path != "empty" and os.path.exists(path):
            with open(path) as f:
                for line in f:
                    haps[str(j)].append(line.split())
    return haps


def buildReplicateIndexes(X: int, haps: dict) -> dict:
    indexes = {}
    for j in range(X):
        rows = haps[str(j)]
        if not rows:
            continue
        counts_by_sequence = {}
        for row in rows:
            counts_by_sequence.setdefault(row[4], row[3])
        indexes[str(j)] = {
            "forward_tag": rows[0][1],
            "reverse_tag": rows[0][2],
            "counts_by_sequence": counts_by_sequence,
        }
    return indexes


def allSequences(indexes: dict) -> set[str]:
    return {
        sequence
        for replicate in indexes.values()
        for sequence in replicate["counts_by_sequence"]
    }


def MakeComparisonFile(X, seqsALL, haps, indexes,
                       OUT, OUTthresh, OUTYX, OUT_fas, OUTthresh_fas,
                       OUTYX_fas, OUTthreshLen_fas, Y, T, L, sampleName, i):
    idnum = 1
    # Sort for deterministic output (and to match the Rust implementation's
    # ordering): a Python set iterates in hash-seed-dependent order, which makes
    # line order and the per-sequence _N ids non-reproducible run to run.
    for seq in sorted(seqsALL):
        line = sampleName[i] + "\t"
        lineFasIDs = ">" + sampleName[i] + "\t"
        lineFasCounts = "\t"
        y = 0
        t = 0
        for j in range(X):
            if len(haps[str(j)]) != 0:
                replicate = indexes[str(j)]
                count = replicate["counts_by_sequence"].get(seq)
                if count is None:
                    count = 0
                else:
                    y += 1
                    if int(count) < T:
                        t += 1
                F = replicate["forward_tag"]
                R = replicate["reverse_tag"]
                line = line + F + "-" + R + "\t" + str(count) + "\t"
                if j < (X - 1):
                    lineFasIDs = lineFasIDs + F + "-" + R + "."
                    lineFasCounts = lineFasCounts + str(count) + "_"
                else:
                    lineFasIDs = lineFasIDs + F + "-" + R + "_" + str(idnum) + "\t"
                    lineFasCounts = lineFasCounts + str(count) + "\n" + seq
            if len(haps[str(j)]) == 0:
                line = line + "empty\t0\t"
                if j < (X - 1):
                    lineFasIDs = lineFasIDs + "empty-empty."
                    lineFasCounts = lineFasCounts + "0_"
                else:
                    lineFasIDs = lineFasIDs + "empty-empty_" + str(idnum) + "\t"
                    lineFasCounts = lineFasCounts + "0\n" + seq
        line = line + seq + "\n"
        lineFas = lineFasIDs + lineFasCounts + "\n"
        OUT.write(line)
        OUT_fas.write(lineFas)
        if y >= Y:
            OUTYX.write(line)
            OUTYX_fas.write(lineFas)
        if (y - t) >= Y:
            OUTthresh.write(line)
            OUTthresh_fas.write(lineFas)
            if len(seq) >= L:
                OUTthreshLen_fas.write(lineFas)
        idnum += 1
