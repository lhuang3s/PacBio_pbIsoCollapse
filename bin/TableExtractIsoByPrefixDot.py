#!/usr/bin/python3
#TableExtractIsoByPrefixDot.py
#This script takes a multi-column table and a prefix as input, extract items based on the prefix, and output a new subset of the table
#For pbIsoCollapse.sh
#Version: 2021-08-09
                                                                          
import sys
import re

def TableExtractIsoByPrefixDot():
    fi=open(sys.argv[1],'r')
    fpre=sys.argv[2]
    fo=open(sys.argv[3],'w')

    for line in fi:
        Name=line.strip().split()[0]
        if re.match(fpre+"\.", Name):
            fo.write(line)

    fi.close()
    fo.close()

if len(sys.argv) != 4:
    print("This script takes a multi-column table and a prefix as input, extract items based on the prefix, and output a new subset of the table")
    print("For pbIsoCollapse.sh")
    print("Usage: [TableExtractIsoByPrefixDot.py] [Table.2col] [Prefix] [OutputFile]")
else:
    TableExtractIsoByPrefixDot()
