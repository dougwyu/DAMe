import re
import os
import sys
import subprocess


def makeTagFiles(PSinfo, X):
    PSouts = [open("PS%s.tags.txt" % (i + 1), "w") for i in range(X)]
    with open(PSinfo) as f:
        PS = f.readlines()
    for NR, psinfo in enumerate(PS):
        NR = NR + 1
        psinfo = psinfo.rstrip().split()
        residue = NR % X
        idx = residue - 1 if residue != 0 else X - 1
        PSouts[idx].write("%s\t%s\n" % (psinfo[1], psinfo[2]))
    for out in PSouts:
        out.close()


def makeTagFilesWithPools(PSinfo, X):
    PSouts = [open("PS%s.tags.txt" % (i + 1), "w") for i in range(X)]
    with open(PSinfo) as f:
        PS = f.readlines()
    for NR, psinfo in enumerate(PS):
        NR = NR + 1
        psinfo = psinfo.rstrip().split()
        residue = NR % X
        idx = residue - 1 if residue != 0 else X - 1
        PSouts[idx].write("%s\t%s\t%s\n" % (psinfo[1], psinfo[2], psinfo[3]))
    for out in PSouts:
        out.close()


def MakeSizeOutFastas(P, X):
    OUTS = [open("Pool%s.fasta" % (pool + 1), "w") for pool in range(P)]
    for num in range(X):
        with open("PS%s.tags.txt" % (num + 1)) as f:
            line = f.readline()
            while line:
                line = line.rstrip().split()
                hap = "_".join([line[0], line[1]]) + ".txt"
                if P > 1:
                    hap = "./pool" + str(line[2]) + "/" + hap
                if not os.path.exists(hap):
                    line = f.readline()
                    continue
                with open(hap) as hap_f:
                    idNum = 1
                    for seq in hap_f:
                        seq = seq.rstrip().split()
                        a = (">" + "_".join([seq[0], line[0], line[1], str(idNum)])
                             + ";size=" + str(seq[3]) + "\n" + seq[4])
                        pool_idx = int(line[2]) - 1 if P > 1 else 0
                        OUTS[pool_idx].write("%s\n" % a)
                        idNum += 1
                line = f.readline()
    for out in OUTS:
        out.close()


def SortFasta(P):
    for pool in range(P):
        input_file = "Pool%s.fasta" % (pool + 1)
        output = "Pool%s.sort.fasta" % (pool + 1)
        cmd = "usearch --sortsize " + input_file + " --output " + output
        p_core = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = p_core.communicate()
        with open("sort%s.out" % (pool + 1), "wb") as fh:
            fh.write(stdout)
        with open("sort%s.err" % (pool + 1), "wb") as fh:
            fh.write(stderr)
        input_file = output
        output1 = "Pool%s.Chim.fasta" % (pool + 1)
        output2 = "Pool%s.noChim.fasta" % (pool + 1)
        cmd = "usearch -uchime " + input_file + " -chimeras " + output1 + " -nonchimeras " + output2
        p_core = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = p_core.communicate()
        with open("chimeraCheck%s.out" % (pool + 1), "wb") as fh:
            fh.write(stdout)
        with open("chimeraCheck%s.err" % (pool + 1), "wb") as fh:
            fh.write(stderr)


def MakeFasSeqOneLine(P):
    for pool in range(P):
        with open("Pool%s.noChim.fasta" % (pool + 1)) as fasta:
            with open("Pool%s.noChim.oneLiner.fasta" % (pool + 1), "w") as fastaOne:
                seq = ""
                for line in fasta:
                    line = line.rstrip()
                    if line.startswith(">"):
                        if seq:
                            fastaOne.write("%s\n" % seq)
                        fastaOne.write("%s\n" % line)
                        seq = ""
                    else:
                        seq += line
                if seq:
                    fastaOne.write("%s\n" % seq)


def MakeNoChimHaps(P):
    HAP = {}
    for pool in range(P):
        with open("Pool%s.noChim.oneLiner.fasta" % (pool + 1)) as fasta:
            primerName = tagName1 = tagName2 = freq = tagHapKey = ""
            for line in fasta:
                line = line.rstrip()
                if line.startswith(">"):
                    line = line[1:]
                    primerName = line.split("_")[0]
                    tagName1 = line.split("_")[1]
                    tagName2 = line.split("_")[2]
                    freq = line.split("=")[1]
                    tagHapKey = "_".join([tagName1, tagName2, str(pool + 1)])
                    if tagHapKey not in HAP:
                        HAP[tagHapKey] = [[primerName], [tagName1], [tagName2], [freq], []]
                    else:
                        HAP[tagHapKey][0].append(primerName)
                        HAP[tagHapKey][1].append(tagName1)
                        HAP[tagHapKey][2].append(tagName2)
                        HAP[tagHapKey][3].append(freq)
                else:
                    HAP[tagHapKey][4].append(line + "\n")
    for TagComb in HAP:
        with open("%s.noChim.txt" % TagComb, "w") as out:
            for i in range(len(HAP[TagComb][0])):
                a = "\t".join([HAP[TagComb][0][i], HAP[TagComb][1][i],
                               HAP[TagComb][2][i], HAP[TagComb][3][i],
                               HAP[TagComb][4][i]])
                out.write(a)
