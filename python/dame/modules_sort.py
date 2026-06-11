
def RC(seq):
    seq = seq[::-1]
    transtab = str.maketrans('ACGTMRWSYKVHDB', 'TGCAKYWSRMBDHV')
    return seq.translate(transtab)


_IUPAC = {
    'A': frozenset('A'), 'C': frozenset('C'), 'G': frozenset('G'),
    'T': frozenset('T'), 'R': frozenset('AG'), 'Y': frozenset('CT'),
    'S': frozenset('CG'), 'W': frozenset('AT'), 'K': frozenset('GT'),
    'M': frozenset('AC'), 'B': frozenset('CGT'), 'D': frozenset('AGT'),
    'H': frozenset('ACT'), 'V': frozenset('ACG'), 'N': frozenset('ACGT'),
}


def iupac_matches(primer_base, read_base):
    """True if read_base (A/C/G/T) satisfies the primer's IUPAC code."""
    allowed = _IUPAC.get(primer_base)
    return allowed is not None and read_base in allowed


def find_primer(primer, seq, max_mismatches=0):
    """Leftmost window in seq matching primer with <= max_mismatches
    substitutions (IUPAC-aware). Returns (start, end) or None."""
    plen = len(primer)
    slen = len(seq)
    if plen > slen:
        return None
    for i in range(slen - plen + 1):
        mismatches = 0
        ok = True
        for p, s in zip(primer, seq[i:i + plen]):
            if not iupac_matches(p, s):
                mismatches += 1
                if mismatches > max_mismatches:
                    ok = False
                    break
        if ok:
            return (i, i + plen)
    return None


def hamming_iupac(pattern, region):
    """Count positions where region fails pattern's IUPAC constraint.
    Returns a sentinel greater than len(pattern) when lengths differ."""
    if len(pattern) != len(region):
        return len(pattern) + 1
    return sum(0 if iupac_matches(p, r) else 1 for p, r in zip(pattern, region))


def min_equal_length_tag_distance(TAGS):
    """Minimum pairwise Hamming distance among equal-length forward tag
    sequences (TAGS[name][0]). Returns None if no two tags share a length."""
    seqs = [TAGS[t][0] for t in TAGS]
    best = None
    for i in range(len(seqs)):
        for j in range(i + 1, len(seqs)):
            a, b = seqs[i], seqs[j]
            if len(a) != len(b):
                continue
            d = sum(1 for x, y in zip(a, b) if x != y)
            best = d if best is None else min(best, d)
    return best


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


def readPrimers(primers, PRIMERS):
    with open(primers) as f:
        for line in f:
            line = line.rstrip().split()
            if not line:
                continue
            if line[0] not in PRIMERS:
                PRIMERS[line[0]] = [[], []]
            F = line[1]
            R = line[2]
            Frc = RC(F)
            Rrc = RC(R)
            PRIMERS[line[0]][0].append(F)
            PRIMERS[line[0]][0].append(R)
            PRIMERS[line[0]][1].append(Frc)
            PRIMERS[line[0]][1].append(Rrc)
    return PRIMERS


def GetPiecesInfo(line, PRIMERS, TAGS, keepPrimersSeq, maxMismatches=0):
    # Normalise soft-masked / lowercase bases so they are not silently dropped.
    line = line.upper()
    for key in PRIMERS:
        # Forward orientation: F at start, RC(R) at end
        primIni = find_primer(PRIMERS[key][0][0], line, maxMismatches)
        if primIni is not None:
            if keepPrimersSeq:
                primIniPosPrim = primIni[0]
                primIniPosTags = primIni[0]
            else:
                primIniPosPrim = primIni[1]
                primIniPosTags = primIni[0]
            primFin = find_primer(PRIMERS[key][1][1], line, maxMismatches)
            if primFin is not None:
                if keepPrimersSeq:
                    primFinPosPrim = primFin[1]
                    primFinPosTags = primFin[1]
                else:
                    primFinPosPrim = primFin[0]
                    primFinPosTags = primFin[1]
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
            # Reverse orientation: R at start, RC(F) at end
            primIni = find_primer(PRIMERS[key][0][1], line, maxMismatches)
            if primIni is not None:
                if keepPrimersSeq:
                    primIniPosPrim = primIni[0]
                    primIniPosTags = primIni[0]
                else:
                    primIniPosPrim = primIni[1]
                    primIniPosTags = primIni[0]
                primFin = find_primer(PRIMERS[key][1][0], line, maxMismatches)
                if primFin is not None:
                    if keepPrimersSeq:
                        primFinPosPrim = primFin[1]
                        primFinPosTags = primFin[1]
                    else:
                        primFinPosPrim = primFin[0]
                        primFinPosTags = primFin[1]
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


def GetPiecesInfoMismatch(line, PRIMERS, TAGS, keepPrimersSeq, max_primer_mm, max_tag_mm):
    """Anchored matcher: tolerate <= max_primer_mm per primer and <= max_tag_mm
    per tag. Returns [tag1, tag2, primer, between] for the unique best-scoring
    assembly, or [1] on no match / ambiguous tie. Mirrors the Rust
    get_pieces_info_anchored."""
    line = line.upper()
    slen = len(line)
    best = None
    best_score = None
    tied = False

    for key in PRIMERS:
        F, R = PRIMERS[key][0][0], PRIMERS[key][0][1]
        Frc, Rrc = PRIMERS[key][1][0], PRIMERS[key][1][1]
        for orientation in (0, 1):
            if orientation == 0:
                start_primer, end_primer = F, Rrc
            else:
                start_primer, end_primer = R, Frc

            for tag1_name in TAGS:
                tag1_seq = TAGS[tag1_name][0]
                t1l = len(tag1_seq)
                if t1l + len(start_primer) > slen:
                    continue
                tag1_mm = hamming_iupac(tag1_seq, line[0:t1l])
                if tag1_mm > max_tag_mm:
                    continue
                p_start = t1l
                p_end = t1l + len(start_primer)
                start_mm = hamming_iupac(start_primer, line[p_start:p_end])
                if start_mm > max_primer_mm:
                    continue

                for tag2_name in TAGS:
                    tag2_seq = TAGS[tag2_name][1]
                    t2l = len(tag2_seq)
                    if t2l + len(end_primer) > slen:
                        continue
                    tag2_start = slen - t2l
                    tag2_mm = hamming_iupac(tag2_seq, line[tag2_start:])
                    if tag2_mm > max_tag_mm:
                        continue
                    ep_start = tag2_start - len(end_primer)
                    ep_end = tag2_start
                    if ep_start < p_end:
                        continue
                    end_mm = hamming_iupac(end_primer, line[ep_start:ep_end])
                    if end_mm > max_primer_mm:
                        continue

                    if keepPrimersSeq:
                        b_start, b_end = p_start, ep_end
                    else:
                        b_start, b_end = p_end, ep_start
                    if b_start >= b_end:
                        continue
                    between = line[b_start:b_end]
                    if orientation == 1:
                        between = RC(between)

                    if orientation == 0:
                        tn1, tn2 = tag1_name, tag2_name
                    else:
                        tn1, tn2 = tag2_name, tag1_name

                    score = tag1_mm + start_mm + end_mm + tag2_mm
                    if best_score is None or score < best_score:
                        best_score = score
                        best = [tn1, tn2, key, between]
                        tied = False
                    elif score == best_score:
                        tied = True

    if tied or best is None:
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


def PrintSummaryFile(HAP):
    with open("SummaryCounts.txt", "w") as out:
        out.write("%s\n" % "\t".join(("#tagName1", "tagName2", "NumUniqSeqs", "SumTotalFreq")))
        for TagComb in HAP:
            tagName1 = HAP[TagComb][0]
            tagName2 = HAP[TagComb][1]
            NumUniqSeqs = len(HAP[TagComb][2])
            SumTotalFreq = sum(HAP[TagComb][2][Seq][0] for Seq in HAP[TagComb][2])
            out.write("%s\n" % "\t".join((tagName1, tagName2, str(NumUniqSeqs), str(SumTotalFreq))))
