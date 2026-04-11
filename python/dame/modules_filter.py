import os


def makePSnumFiles(PSinfo, X, P, chimeraChecked):
    PSouts = [open("PS%s_files.txt" % (i + 1), "w") for i in range(X)]
    with open(PSinfo) as f:
        PS = f.readlines()
    for NR, psinfo in enumerate(PS):
        NR = NR + 1
        psinfo = psinfo.rstrip().split()
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


def getSeqsSetsAndFRcounts(X, haps):
    F = {}
    R = {}
    counts = {}
    seqs = {}
    seqsALL = []
    for j in range(X):
        if len(haps[str(j)]) != 0:
            seqs[str(j)] = []
            F[str(j)] = haps[str(j)][0][1]
            R[str(j)] = haps[str(j)][0][2]
            counts[str(j)] = []
            for k in range(len(haps[str(j)])):
                counts[str(j)].append(haps[str(j)][k][3])
                seqs[str(j)].append(haps[str(j)][k][4])
                seqsALL.append(haps[str(j)][k][4])
    seqsALL = set(seqsALL)
    return (seqsALL, F, R, counts, seqs)


def MakeComparisonFile(X, seqsALL, haps, F, R, counts, seqs,
                       OUT, OUTthresh, OUTYX, OUT_fas, OUTthresh_fas,
                       OUTYX_fas, OUTthreshLen_fas, Y, T, L, sampleName, i):
    idnum = 1
    for seq in seqsALL:
        line = sampleName[i] + "\t"
        lineFasIDs = ">" + sampleName[i] + "\t"
        lineFasCounts = "\t"
        y = 0
        t = 0
        for j in range(X):
            if len(haps[str(j)]) != 0:
                pos = [pos for pos, s in enumerate(seqs[str(j)]) if seq == s]
                if len(pos) == 0:
                    count = 0
                else:
                    y += 1
                    count = counts[str(j)][pos[0]]
                    if int(count) < T:
                        t += 1
                line = line + F[str(j)] + "-" + R[str(j)] + "\t" + str(count) + "\t"
                if j < (X - 1):
                    lineFasIDs = lineFasIDs + F[str(j)] + "-" + R[str(j)] + "."
                    lineFasCounts = lineFasCounts + str(count) + "_"
                else:
                    lineFasIDs = lineFasIDs + F[str(j)] + "-" + R[str(j)] + "_" + str(idnum) + "\t"
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
