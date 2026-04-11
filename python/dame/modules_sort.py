import re
import os
import sys


def RC(seq):
    seq = seq[::-1]
    transtab = str.maketrans('ACGTMRWSYKVHDB', 'TGCAKYWSRMBDHV')
    return seq.translate(transtab)


def readTags(tags, TAGS):
    with open(tags) as f:
        for line in f:
            line = line.rstrip().split()
            if not line:
                continue
            if line[1] not in TAGS:
                TAGS[line[1]] = []
            TAGS[line[1]].append(line[0])
            TAGS[line[1]].append(RC(line[0]))
    return TAGS


def readPrimers(primers, PRIMERS, AMBIG):
    with open(primers) as f:
        for line in f:
            line = line.rstrip().split()
            if not line:
                continue
            if line[0] not in PRIMERS:
                PRIMERS[line[0]] = [[], []]
            Frc = RC(line[1])
            Rrc = RC(line[2])
            F = line[1]
            R = line[2]
            for key in AMBIG:
                Frc = re.sub(key, AMBIG[key], Frc)
                Rrc = re.sub(key, AMBIG[key], Rrc)
                F = re.sub(key, AMBIG[key], F)
                R = re.sub(key, AMBIG[key], R)
            PRIMERS[line[0]][0].append(F)
            PRIMERS[line[0]][0].append(R)
            PRIMERS[line[0]][1].append(Frc)
            PRIMERS[line[0]][1].append(Rrc)
    return PRIMERS


def GetPiecesInfo(line, PRIMERS, TAGS, keepPrimersSeq):
    for key in PRIMERS:
        primIniPos = [(m.start(0), m.end(0)) for m in re.finditer(PRIMERS[key][0][0], line)]
        if len(primIniPos) > 0:
            if keepPrimersSeq:
                primIniPosPrim = primIniPos[0][0]
                primIniPosTags = primIniPos[0][0]
            else:
                primIniPosPrim = primIniPos[0][1]
                primIniPosTags = primIniPos[0][0]
            primFinPos = [(m.start(0), m.end(0)) for m in re.finditer(PRIMERS[key][1][1], line)]
            if len(primFinPos) > 0:
                if keepPrimersSeq:
                    primFinPosPrim = primFinPos[0][1]
                    primFinPosTags = primFinPos[0][1]
                else:
                    primFinPosPrim = primFinPos[0][0]
                    primFinPosTags = primFinPos[0][1]
                PrimerName = key
                between = line[primIniPosPrim:primFinPosPrim]
                if len(between) == 0:
                    return [1]
                tag1 = line[:primIniPosTags]
                tag2 = line[primFinPosTags:]
                tagName1 = [t for t in TAGS if TAGS[t][0] == tag1]
                tagName2 = [t for t in TAGS if TAGS[t][1] == tag2]
                if len(tagName1) > 0 and len(tagName2) > 0:
                    return [tagName1[0], tagName2[0], PrimerName, between]
                return [1]
            return [1]
        else:
            primIniPos = [(m.start(0), m.end(0)) for m in re.finditer(PRIMERS[key][0][1], line)]
            if len(primIniPos) > 0:
                if keepPrimersSeq:
                    primIniPosPrim = primIniPos[0][0]
                    primIniPosTags = primIniPos[0][0]
                else:
                    primIniPosPrim = primIniPos[0][1]
                    primIniPosTags = primIniPos[0][0]
                primFinPos = [(m.start(0), m.end(0)) for m in re.finditer(PRIMERS[key][1][0], line)]
                if len(primFinPos) > 0:
                    if keepPrimersSeq:
                        primFinPosPrim = primFinPos[0][1]
                        primFinPosTags = primFinPos[0][1]
                    else:
                        primFinPosPrim = primFinPos[0][0]
                        primFinPosTags = primFinPos[0][1]
                    PrimerName = key
                    between = line[primIniPosPrim:primFinPosPrim]
                    if len(between) == 0:
                        return [1]
                    between = RC(between)
                    tag1 = line[:primIniPosTags]
                    tag2 = line[primFinPosTags:]
                    tagName2 = [t for t in TAGS if TAGS[t][0] == tag1]
                    tagName1 = [t for t in TAGS if TAGS[t][1] == tag2]
                    if len(tagName1) > 0 and len(tagName2) > 0:
                        return [tagName1[0], tagName2[0], PrimerName, between]
                    return [1]
                return [1]
    return [1]


def FillHAP(HAP, tagName1, tagName2, PrimerName, between):
    tagHapKey = "_".join([tagName1, tagName2])
    if tagHapKey not in HAP:
        HAP[tagHapKey] = [tagName1, tagName2, {}]
    if between not in HAP[tagHapKey][2]:
        HAP[tagHapKey][2][between] = [1, PrimerName]
    else:
        HAP[tagHapKey][2][between][0] += 1
    return HAP


def PrintSortedCollapsedCountedSeqs(HAP):
    for TagComb in HAP:
        with open("%s.txt" % TagComb, "w") as out:
            tagName1 = HAP[TagComb][0]
            tagName2 = HAP[TagComb][1]
            for Seq in HAP[TagComb][2]:
                a = "\t".join([HAP[TagComb][2][Seq][1], tagName1, tagName2,
                               str(HAP[TagComb][2][Seq][0]), Seq])
                out.write("%s\n" % a)


def PrintSummaryFile(HAP):
    with open("SummaryCounts.txt", "w") as out:
        out.write("%s\n" % "\t".join(("#tagName1", "tagName2", "NumUniqSeqs", "SumTotalFreq")))
        for TagComb in HAP:
            tagName1 = HAP[TagComb][0]
            tagName2 = HAP[TagComb][1]
            NumUniqSeqs = len(HAP[TagComb][2])
            SumTotalFreq = sum(HAP[TagComb][2][Seq][0] for Seq in HAP[TagComb][2])
            out.write("%s\n" % "\t".join((tagName1, tagName2, str(NumUniqSeqs), str(SumTotalFreq))))
