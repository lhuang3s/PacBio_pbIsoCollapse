#!/usr/bin/python3
#BED14Correct3endByLastIntron.py
#This script takes a bed12 + intron col2-col3 in col13-col14 as input, correct the 3ends if they have same last intron splicing junction.
#Then it will output bed12 corrected isoforms, keeping the names of all original isoforms.
#This script will take all input isoforms for calculation, so please make sure they are on the same strand, and can be compared.
#Version: 2022-10-23

import sys
import re

def Counter():
    fi=open(sys.argv[1],'r')
    strand=sys.argv[2]
    fo=open(sys.argv[3],'w')
    
    LastExonSite2Longest3end={}

    #Determine the longest last exon and save them into a dict
    for line in fi:
        Data=line.strip().split()
        if strand=="+":
            IntronSite=int(Data[13])
            Tx3end=int(Data[2])
        else:
            IntronSite=int(Data[12])
            Tx3end=int(Data[1])
        if len(LastExonSite2Longest3end) == 0:
            LastExonSite2Longest3end[IntronSite]=Tx3end
        else:
            if IntronSite not in LastExonSite2Longest3end.keys():
                LastExonSite2Longest3end[IntronSite]=Tx3end
            else:
                if strand=="+":
                    if Tx3end > LastExonSite2Longest3end[IntronSite]:
                        LastExonSite2Longest3end[IntronSite]=Tx3end
                else:
                    if Tx3end < LastExonSite2Longest3end[IntronSite]:
                        LastExonSite2Longest3end[IntronSite]=Tx3end

    #print(LastExonSite2Longest3end)

    fi.seek(0)
    #Compare each 3'end to the longest 3end sharing the same last exon junction, and correct it if needed.
    for line in fi:
        Data=line.strip().split()
        if strand=="+":
            IntronSite=int(Data[13])
        else:
            IntronSite=int(Data[12])
        Longest3end=LastExonSite2Longest3end[IntronSite]
        if strand=="+":
            CurrentTx3end=int(Data[2])
            MoveRight=Longest3end-CurrentTx3end
            ExonBlocks=re.sub(',$','',Data[10]).split(",")
            ExonStarts=re.sub(',$','',Data[11]).split(",")
            ExonBlocks[len(ExonBlocks)-1]=int(ExonBlocks[len(ExonBlocks)-1])+MoveRight
            fo.write(Data[0]+"\t"+Data[1]+"\t"+str(Longest3end)+"\t"+Data[3]+"\t"+Data[4]+"\t"+Data[5]+"\t"+Data[6]+"\t"+Data[7]+"\t"+Data[8]+"\t"+Data[9]+"\t"+",".join(str(x) for x in ExonBlocks)+",\t"+",".join(str(x) for x in ExonStarts)+",\n")
        else:
            CurrentTx3end=int(Data[1])
            MoveRight=CurrentTx3end-Longest3end
            ExonBlocks=re.sub(',$','',Data[10]).split(",")
            ExonStarts=re.sub(',$','',Data[11]).split(",")
            ExonBlocks[0]=int(ExonBlocks[0])+MoveRight
            for i in range(len(ExonStarts)):
                if i!=0:
                    ExonStarts[i]=int(ExonStarts[i])+MoveRight
            fo.write(Data[0]+"\t"+str(Longest3end)+"\t"+Data[2]+"\t"+Data[3]+"\t"+Data[4]+"\t"+Data[5]+"\t"+Data[6]+"\t"+Data[7]+"\t"+Data[8]+"\t"+Data[9]+"\t"+",".join(str(x) for x in ExonBlocks)+",\t"+",".join(str(x) for x in ExonStarts)+",\n")

    fi.close()
    fo.close()

if len(sys.argv) != 4:
    print("This script takes a bed12 + intron col2-col3 in col13-col14 as input, correct the 3ends if they have same last intron splicing junction.")
    print("Then it will output bed12 corrected isoforms, keeping the names of all original isoforms.")
    print("This script will take all input isoforms for calculation, so please make sure they are on the same strand, and can be compared.")
    print("Usage: [BED14Correct3endByLastIntron.py] [Data.bed14] [Strand] [Output.bed12]")
else:
    Counter()
