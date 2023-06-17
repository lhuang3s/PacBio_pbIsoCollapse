#!/usr/bin/python3
#ReOrderIsoCollapseOutput.py
#This script takes collapsed isoforms in sorted bed12 format and match the abundance of them using the Prefix.Original2Final.matching file.
#Output: A new abundance file with the same isoform order as the bed12.
#Version: Yu Sun, 2022-10-23

import sys
import re

def ReOrderIsoCollapseOutput():
    fi=open(sys.argv[1],'r')
    fmatching=open(sys.argv[2],'r')
    fo=open(sys.argv[3],'w')
    
    AbundanceDict={}
    for line in fmatching:
        CurrLine=line.strip().split()
        CurrIsoform=CurrLine[1]
        CurrAbd=CurrLine[2:]
        if CurrIsoform not in AbundanceDict.keys():
            AbundanceDict[CurrIsoform]=CurrAbd

    for isoform in fi:
        CurrIso=isoform.strip().split()
        Name=CurrIso[3]
        fo.write(Name+"\t"+"\t".join(AbundanceDict[Name])+"\n")

    fi.close()
    fmatching.close()
    fo.close()

if len(sys.argv) != 4:
    print("This script takes collapsed isoforms in sorted bed12 format and match the abundance of them using the Prefix.Original2Final.matching file.")
    print("Output: A new abundance file with the same isoform order as the bed12.")
    print("Usage: [ReOrderIsoCollapseOutput.py] [Data.sorted.bed12] [Data.Original2Final.matching] [Data.abundance.txt]")
else:
    ReOrderIsoCollapseOutput()
