import os

from dame.utils import smart_open


def makePSnumFiles(PSinfo, X, P, chimeraChecked, outdir="."):
    from collections import OrderedDict
    PSouts = [open(os.path.join(outdir, "PS%s_files.txt" % (i + 1)), "w") for i in range(X)]
    sample_lines = OrderedDict()
    with smart_open(PSinfo) as f:
        for line in f:
            line = line.rstrip()
            if not line:
                continue
            psinfo = line.split()
            if len(psinfo) < 4:
                continue
            sample = psinfo[0]
            if sample not in sample_lines:
                sample_lines[sample] = []
            if not chimeraChecked:
                entry = "pool%s/%s_%s.txt\n" % (psinfo[3], psinfo[1], psinfo[2])
            else:
                entry = "%s_%s_%s.noChim.txt\n" % (psinfo[1], psinfo[2], psinfo[3])
            sample_lines[sample].append(entry)
    for entries in sample_lines.values():
        for k in range(X):
            if k < len(entries):
                PSouts[k].write(entries[k])
            else:
                PSouts[k].write("empty\n")
    for out in PSouts:
        out.close()


def ReadPSnumFiles(X, outdir="."):
    PSinsLines = {}
    for i in range(X):
        with open(os.path.join(outdir, "PS%s_files.txt" % (i + 1))) as f:
            PSinsLines[str(i)] = f.readlines()
    return PSinsLines


def MakeSampleNameArray(PSinfo):
    sampleName = []
    with smart_open(PSinfo) as f:
        for line in f:
            line = line.rstrip()
            if not line:
                continue
            parts = line.split()
            if not parts:
                continue
            name = parts[0]
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
                    row = line.split()
                    if len(row) >= 5:
                        haps[str(j)].append(row)
    return haps


def getSeqsSetsAndFRcounts(X, haps):
    F = {}
    R = {}
    seqs = {}   # dict[replicate_str -> dict[seq -> count_str]]
    seqsALL = []
    for j in range(X):
        if len(haps[str(j)]) != 0:
            seqs[str(j)] = {}
            F[str(j)] = haps[str(j)][0][1]
            R[str(j)] = haps[str(j)][0][2]
            for k in range(len(haps[str(j)])):
                seq = haps[str(j)][k][4]
                seqs[str(j)][seq] = haps[str(j)][k][3]
                seqsALL.append(seq)
    seqsALL = set(seqsALL)
    return (seqsALL, F, R, seqs)


def MakeComparisonFile(X, seqsALL, haps, F, R, seqs,
                       OUT, OUTthresh, OUTYX, OUT_fas, OUTthresh_fas,
                       OUTYX_fas, OUTthreshLen_fas, Y, T, L, sampleName, i):
    idnum = 1
    for seq in sorted(seqsALL):
        line = sampleName[i] + "\t"
        lineFasIDs = ">" + sampleName[i] + "\t"
        lineFasCounts = "\t"
        y = 0
        t = 0
        for j in range(X):
            if len(haps[str(j)]) != 0:
                count_str = seqs[str(j)].get(seq, "0")   # O(1) dict lookup
                if count_str != "0":
                    y += 1
                    try:
                        count_val = int(count_str)
                    except ValueError:
                        count_val = 0
                    if count_val < T:
                        t += 1
                line = line + F[str(j)] + "-" + R[str(j)] + "\t" + count_str + "\t"
                if j < (X - 1):
                    lineFasIDs = lineFasIDs + F[str(j)] + "-" + R[str(j)] + "."
                    lineFasCounts = lineFasCounts + count_str + "_"
                else:
                    lineFasIDs = lineFasIDs + F[str(j)] + "-" + R[str(j)] + "_" + str(idnum) + "\t"
                    lineFasCounts = lineFasCounts + count_str + "\n" + seq
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
