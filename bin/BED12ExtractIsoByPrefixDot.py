#!/usr/bin/python3
#BED12ExtractIsoByPrefixDot.py
#This script extracts bed12 isoform by the name matching a prefix plus a dot
#Version: 2021-08-07
                                                                          
import sys
import re

def BED12ExtractIsoByPrefixDot():
    fi=open(sys.argv[1],'r')
    fpre=sys.argv[2]
    fo=open(sys.argv[3],'w')
    
    for line in fi:
        Name=line.strip().split()[3]
        if re.match(fpre+"\.", Name):
            fo.write(line)

    fi.close()
    fo.close()

if len(sys.argv) != 4:
    print("This script extracts bed12 isoform by the name matching a prefix plus a dot")
    print("Usage: [BED12ExtractIsoByPrefixDot.py] [Input.bed|at least 4 columns] [PatternPrefix] [OutputName]")
else:
    BED12ExtractIsoByPrefixDot()
