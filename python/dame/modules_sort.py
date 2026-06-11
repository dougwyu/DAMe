import re
import os
import sys

from dame.utils import smart_open


def RC(seq):
    seq = seq[::-1]
    transtab = str.maketrans('ACGTMRWSYKVHDB', 'TGCAKYWSRMBDHV')
    return seq.translate(transtab)


IUPAC_SETS = {
    'A': {'A'}, 'C': {'C'}, 'G': {'G'}, 'T': {'T'},
    'M': {'A', 'C'}, 'R': {'A', 'G'}, 'W': {'A', 'T'},
    'S': {'C', 'G'}, 'Y': {'C', 'T'}, 'K': {'G', 'T'},
    'V': {'A', 'C', 'G'}, 'H': {'A', 'C', 'T'},
    'D': {'A', 'G', 'T'}, 'B': {'C', 'G', 'T'},
    'N': {'A', 'C', 'G', 'T'},
}


def iupac_match(p, r):
    """Return True if read base r satisfies IUPAC primer base p."""
    return r.upper() in IUPAC_SETS.get(p.upper(), set())


def hamming_iupac(primer, region):
    """Count positions where region does not satisfy the IUPAC primer pattern.

    Expects len(primer) == len(region); comparison is limited to the shorter sequence.
    """
    return sum(0 if iupac_match(p, r) else 1 for p, r in zip(primer, region))


def readTags(tags, TAGS):
    with smart_open(tags) as f:
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
    with smart_open(primers) as f:
        for line in f:
            line = line.rstrip().split()
            if not line:
                continue
            if line[0] not in PRIMERS:
                PRIMERS[line[0]] = [[], [], [], []]
            Frc = RC(line[1])
            Rrc = RC(line[2])
            F = line[1]
            R = line[2]
            # Store raw IUPAC sequences before regex substitution (indices 2, 3)
            PRIMERS[line[0]][2].append(F)
            PRIMERS[line[0]][2].append(R)
            PRIMERS[line[0]][3].append(Frc)
            PRIMERS[line[0]][3].append(Rrc)
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
                t1 = tagName1[0] if tagName1 else ''
                t2 = tagName2[0] if tagName2 else ''
                return [None, t1, t2, PrimerName, between]
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
                    t1 = tagName1[0] if tagName1 else ''
                    t2 = tagName2[0] if tagName2 else ''
                    return [None, t1, t2, PrimerName, between]
                return [1]
    return [1]


def GetPiecesInfoMismatch(line, PRIMERS, TAGS, keepPrimersSeq, max_primer_mm, max_tag_mm):
    """Tag-anchored primer search allowing up to max_primer_mm / max_tag_mm mismatches.

    Checks each known tag as a prefix/suffix (Hamming <= max_tag_mm), then verifies
    the adjacent primer at the anchored position (IUPAC-aware Hamming <= max_primer_mm).
    Returns the combination with the fewest total mismatches; ties are discarded.
    """
    read_len = len(line)
    best_mm = None
    best = None
    ambiguous = False

    # Precompute tag candidates once per read — O(48) scan instead of O(48²) per primer.
    # Both orientations use TAGS[t][0] for t1 (prefix) and TAGS[t][1] for t2 (suffix).
    t1_cands = []  # (tag_name, tag_mm, tag_len)
    t2_cands = []
    for t in TAGS:
        t1_seq = TAGS[t][0]
        t1_len = len(t1_seq)
        if t1_len <= read_len:
            mm = hamming_iupac(t1_seq, line[:t1_len])
            if mm <= max_tag_mm:
                t1_cands.append((t, mm, t1_len))
        t2_seq = TAGS[t][1]
        t2_len = len(t2_seq)
        if t2_len <= read_len:
            mm = hamming_iupac(t2_seq, line[read_len - t2_len:])
            if mm <= max_tag_mm:
                t2_cands.append((t, mm, t2_len))

    if not t1_cands or not t2_cands:
        return [1]

    for key in PRIMERS:
        F_raw, R_raw = PRIMERS[key][2][0], PRIMERS[key][2][1]
        Frc_raw, Rrc_raw = PRIMERS[key][3][0], PRIMERS[key][3][1]

        # ── forward orientation: [tag1_fwd][F][amplicon][RC(R)][tag2_rc] ──
        for (t1, t1_mm, t1_len) in t1_cands:
            f_end = t1_len + len(F_raw)
            if f_end > read_len:
                continue
            f_mm = hamming_iupac(F_raw, line[t1_len:f_end])
            if f_mm > max_primer_mm:
                continue
            for (t2, t2_mm, t2_len) in t2_cands:
                rrc_end = read_len - t2_len
                rrc_start = rrc_end - len(Rrc_raw)
                if rrc_start < f_end or rrc_start < 0:
                    continue
                rrc_mm = hamming_iupac(Rrc_raw, line[rrc_start:rrc_end])
                if rrc_mm > max_primer_mm:
                    continue
                between = line[t1_len:rrc_end] if keepPrimersSeq else line[f_end:rrc_start]
                if not between:
                    continue
                total = t1_mm + f_mm + rrc_mm + t2_mm
                if best_mm is None or total < best_mm:
                    best_mm, best, ambiguous = total, [t1, t2, key, between], False
                elif total == best_mm:
                    ambiguous = True

        # ── reverse orientation: [tag1_fwd][R][RC(amplicon)][RC(F)][tag2_rc] ──
        for (t1, t1_mm, t1_len) in t1_cands:
            r_end = t1_len + len(R_raw)
            if r_end > read_len:
                continue
            r_mm = hamming_iupac(R_raw, line[t1_len:r_end])
            if r_mm > max_primer_mm:
                continue
            for (t2, t2_mm, t2_len) in t2_cands:
                frc_end = read_len - t2_len
                frc_start = frc_end - len(Frc_raw)
                if frc_start < r_end or frc_start < 0:
                    continue
                frc_mm = hamming_iupac(Frc_raw, line[frc_start:frc_end])
                if frc_mm > max_primer_mm:
                    continue
                between = RC(line[t1_len:frc_end] if keepPrimersSeq else line[r_end:frc_start])
                if not between:
                    continue
                total = t1_mm + r_mm + frc_mm + t2_mm
                if best_mm is None or total < best_mm:
                    best_mm, best, ambiguous = total, [t1, t2, key, between], False
                elif total == best_mm:
                    ambiguous = True

    if best is None or ambiguous:
        return [1]
    return best


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


def read_valid_pairs(psinfo_file, pool_num):
    """Return set of (tag1_lower, tag2_lower) valid pairs for the given pool number.

    PCRsetsInfo format: SampleName  Tag1  Tag2  PoolNumber
    """
    valid = set()
    with smart_open(psinfo_file) as f:
        for line in f:
            parts = line.rstrip().split()
            if len(parts) < 4:
                continue
            try:
                if int(parts[3]) == pool_num:
                    valid.add((parts[1].lower(), parts[2].lower()))
            except ValueError:
                continue
    return valid


def PrintSplitSummaryFile(HAP, HAP_err, no_tags_seqs, prefix, pool, valid_pairs=None):
    """Write splitSummaryByPSInfo_{prefix}_{pool}.txt with a 4-section breakdown."""

    def _stats(d):
        total = sum(sum(v[2][s][0] for s in v[2]) for v in d.values())
        uniq = sum(len(v[2]) for v in d.values())
        return total, uniq

    if valid_pairs is not None:
        same = {k: v for k, v in HAP.items()
                if (v[0].lower(), v[1].lower()) in valid_pairs}
        diff = {k: v for k, v in HAP.items()
                if (v[0].lower(), v[1].lower()) not in valid_pairs}
    else:
        same = {k: v for k, v in HAP.items() if v[0].lower() == v[1].lower()}
        diff = {k: v for k, v in HAP.items() if v[0].lower() != v[1].lower()}

    same_total, same_uniq = _stats(same)
    diff_total, diff_uniq = _stats(diff)
    onetag_total, onetag_uniq = _stats(HAP_err)
    notag_total = sum(no_tags_seqs.values()) if no_tags_seqs else 0
    notag_uniq = len(no_tags_seqs)

    grand_total = same_total + diff_total + onetag_total + notag_total
    grand_uniq = same_uniq + diff_uniq + onetag_uniq + notag_uniq

    def pct(num, denom):
        if denom == 0 or num == 0:
            return "0.0"
        return str(round(num / denom * 100, 2))

    def detail_rows(d, out):
        for k in d:
            t1, t2 = d[k][0], d[k][1]
            n_uniq = len(d[k][2])
            n_total = sum(d[k][2][s][0] for s in d[k][2])
            out.write("%s\t%s\t%s\t%s\n" % (t1, t2, n_uniq, n_total))

    filename = "splitSummaryByPSInfo_%s_%s.txt" % (prefix, pool)
    with open(filename, "w") as out:
        out.write("%-52s\tTotal seqs\tTotal unique seqs\t%% Total seqs\t%%Total unique seqs\n" % "")
        out.write("%-52s\t%s\t%s\t%s\t%s\n" % (
            "Tag combinations where the tag pair was used",
            same_total, same_uniq, pct(same_total, grand_total), pct(same_uniq, grand_uniq)))
        out.write("Tag combinations where both tags used\n")
        out.write("%-52s\t%s\t%s\t%s\t%s\n" % (
            "but not in this combination",
            diff_total, diff_uniq, pct(diff_total, grand_total), pct(diff_uniq, grand_uniq)))
        out.write("%-52s\t%s\t%s\t%s\t%s\n" % (
            "Tag combinations where only one of the tags was used",
            onetag_total, onetag_uniq, pct(onetag_total, grand_total), pct(onetag_uniq, grand_uniq)))
        out.write("%-52s\t%s\t%s\t%s\t%s\n" % (
            "Tag combinations where neither tag was used",
            notag_total, notag_uniq, pct(notag_total, grand_total), pct(notag_uniq, grand_uniq)))
        out.write("%-52s\t%s\t%s\t100.00\t100.00\n" % ("Total", grand_total, grand_uniq))
        out.write("\n")

        out.write("Tag combinations where the tag pair was used.\n")
        out.write("---------------------------------------------\n")
        out.write("Tag1Name\tTag2Name\tNumUniqSeqs\tSumTotalFreq\n")
        detail_rows(same, out)
        out.write("\n")

        out.write("Tag combinations where both tags used - but not in this combination.\n")
        out.write("--------------------------------------------------------------------\n")
        out.write("Tag1Name\tTag2Name\tNumUniqSeqs\tSumTotalFreq\n")
        detail_rows(diff, out)
        out.write("\n")

        out.write("Tag combinations where only one of the tags was used.\n")
        out.write("-----------------------------------------------------\n")
        out.write("Tag1Name\tTag2Name\tNumUniqSeqs\tSumTotalFreq\n")
        detail_rows(HAP_err, out)
        out.write("\n")

        out.write("Tag combinations where neither tag was used.\n")
        out.write("--------------------------------------------\n")
        out.write("Tag1Name\tTag2Name\tNumUniqSeqs\tSumTotalFreq\n")


def PrintSummaryFile(HAP):
    with open("SummaryCounts.txt", "w") as out:
        out.write("%s\n" % "\t".join(("#tagName1", "tagName2", "NumUniqSeqs", "SumTotalFreq")))
        for TagComb in HAP:
            tagName1 = HAP[TagComb][0]
            tagName2 = HAP[TagComb][1]
            NumUniqSeqs = len(HAP[TagComb][2])
            SumTotalFreq = sum(HAP[TagComb][2][Seq][0] for Seq in HAP[TagComb][2])
            out.write("%s\n" % "\t".join((tagName1, tagName2, str(NumUniqSeqs), str(SumTotalFreq))))
