#! /usr/bin/env python
################################################ 
# Convert the FilteredReads.fna file to format #
# accepted by Usearch or sumaclust.            #
# Split by samples and add the size tag. Also  #
# filter by length.                            #
################################################

import argparse
import numpy as np

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Take fasta file and fix it for input to USEARCH/Sumaclust")
    parser.add_argument("-i", "--inFasta", dest="fasta", type=str, metavar="InputFasta", required=True, help='Input fasta file.')
#    parser.add_argument("-n", "--numPCRs", dest="pcrs", type=int, metavar="NumPCRs", required=False, help='Number of PCRs', default=1)
    parser.add_argument("-lmin", "--minLength", dest="lmin", type=int, metavar="minLength", required=False, help='Minimum length of read', default=0)
    parser.add_argument("-lmax", "--maxLength", dest="lmax", type=int, metavar="maxLength", required=False, help='Maximum length of read', default=1e6)
#    parser.add_argument("-o", "--outfile", dest="output", type=str, metavar="OutputFile", required=False, help='Combined output file name', default='')
#    parser.add_argument("-d", "--dereplicate", action="store_true")
#    parser.add_argument("-c", "--combined", action="store_true")
    parser.add_argument("-u", "--usearch", dest="usearch", help="Generate Usearch input (default: sumaclust input)", action="store_true")
    parser.add_argument("-s", "--sampleFastas", dest="samp", help="Make sample level fastas in SampleFastas directory", action="store_true")
    args = parser.parse_args()
    args.lmin -= 1
    args.lmax += 1

## Read the fasta file one line at a time and output if ok.
## If not, skip.
fasta = open(args.fasta)
curreadinfo = ""
sample = ""
curcnt = 1
if not args.usearch:
    outfile = open("FilteredReads.forsumaclust.fna", "w")
else:
    outfile = open("FilteredReads.forusearch.fna", "w")
if args.samp:
    outfiles = {}

for line in fasta:
    toks = line.strip().split()
    if line[0] == '>': ## read info line
        toks[2] = str(sum([int(x) for x in toks[2].split("_")]))
        if args.usearch:
            curreadinfo = toks[0]+";size="+toks[2]+"\n"
        else:
            curreadinfo = toks[0]+":"+str(curcnt)+" count="+toks[2]+"\n"
        sample = toks[0][1:]
        curcnt += 1
    elif len(toks[0]) > args.lmin and len(toks[0]) < args.lmax:
        if args.samp:
            if sample not in outfiles:
                outfiles[sample] = open("SampleFastas/"+sample+".fixed.fasta", "w")
            outfiles[sample].write(curreadinfo)
            outfiles[sample].write(line)
        if args.usearch:
            outfile.write(curreadinfo)
            line = line.strip()
            if len(line) < args.lmax:
                line = line+"N"*(args.lmax - len(line))
            outfile.write(line+"\n")
        else:
            outfile.write(curreadinfo)
            outfile.write(line)
        

outfile.close()
if args.samp:
    for f in outfiles:
        outfiles[f].close()
fasta.close()
        
