#!/usr/bin/python3
#BED14Correct5endByFirstIntron.py
#This script takes a bed12 + intron col2-col3 in col13-col14 as input, correct the 5ends if they have same first intron splicing junction.
#Then it will output bed12 corrected isoforms, keeping the names of all original isoforms.
#This script will take all input isoforms for calculation, so please make sure they are on the same strand, and can be compared.
#Version: 2022-10-23

import sys
import re

def Counter():
    fi=open(sys.argv[1],'r')
    strand=sys.argv[2]
    fo=open(sys.argv[3],'w')
    
    FirstExonSite2Longest5end={}

    for line in fi:
        Data=line.strip().split()
        if strand=="+":
            IntronSite=int(Data[12])
            Tx5end=int(Data[1])
        else:
            IntronSite=int(Data[13])
            Tx5end=int(Data[2])
        if len(FirstExonSite2Longest5end) == 0:
            FirstExonSite2Longest5end[IntronSite]=Tx5end
        else:
            if IntronSite not in FirstExonSite2Longest5end.keys():
                FirstExonSite2Longest5end[IntronSite]=Tx5end
            else:
                if strand=="+":
                    if Tx5end < FirstExonSite2Longest5end[IntronSite]:
                        FirstExonSite2Longest5end[IntronSite]=Tx5end
                else:
                    if Tx5end > FirstExonSite2Longest5end[IntronSite]:
                        FirstExonSite2Longest5end[IntronSite]=Tx5end

#    print(FirstExonSite2Longest5end)

    fi.seek(0)
    for line in fi:
        Data=line.strip().split()
        if strand=="+":
            IntronSite=int(Data[12])
        else:
            IntronSite=int(Data[13])
        Longest5end=FirstExonSite2Longest5end[IntronSite]
        if strand=="+":
            CurrentTx5end=int(Data[1])
            MoveLeft=CurrentTx5end-Longest5end
            ExonBlocks=re.sub(',$','',Data[10]).split(",")
            ExonStarts=re.sub(',$','',Data[11]).split(",")
            ExonBlocks[0]=int(ExonBlocks[0])+MoveLeft
            for i in range(len(ExonStarts)):
                if i!=0:
                    ExonStarts[i]=int(ExonStarts[i])+MoveLeft
            if MoveLeft==0:
                fo.write(Data[0]+"\t"+str(Longest5end)+"\t"+Data[2]+"\t"+Data[3]+"\t"+Data[4]+"\t"+Data[5]+"\t"+Data[6]+"\t"+Data[7]+"\t"+Data[8]+"\t"+Data[9]+"\t"+",".join(str(x) for x in ExonBlocks)+",\t"+",".join(str(x) for x in ExonStarts)+",\n")
            else:
                fo.write(Data[0]+"\t"+str(Longest5end)+"\t"+Data[2]+"\t"+Data[3]+"\t"+Data[4]+"\t"+Data[5]+"\t"+Data[6]+"\t"+Data[7]+"\t"+Data[8]+"\t"+Data[9]+"\t"+",".join(str(x) for x in ExonBlocks)+",\t"+",".join(str(x) for x in ExonStarts)+",\n")
        else:
            CurrentTx5end=int(Data[2])
            MoveLeft=Longest5end-CurrentTx5end
            ExonBlocks=re.sub(',$','',Data[10]).split(",")
            ExonStarts=re.sub(',$','',Data[11]).split(",")
            ExonBlocks[len(ExonBlocks)-1]=int(ExonBlocks[len(ExonBlocks)-1])+MoveLeft
            if MoveLeft==0:
                fo.write(Data[0]+"\t"+Data[1]+"\t"+str(Longest5end)+"\t"+Data[3]+"\t"+Data[4]+"\t"+Data[5]+"\t"+Data[6]+"\t"+Data[7]+"\t"+Data[8]+"\t"+Data[9]+"\t"+",".join(str(x) for x in ExonBlocks)+",\t"+",".join(str(x) for x in ExonStarts)+",\n")
            else:
                fo.write(Data[0]+"\t"+Data[1]+"\t"+str(Longest5end)+"\t"+Data[3]+"\t"+Data[4]+"\t"+Data[5]+"\t"+Data[6]+"\t"+Data[7]+"\t"+Data[8]+"\t"+Data[9]+"\t"+",".join(str(x) for x in ExonBlocks)+",\t"+",".join(str(x) for x in ExonStarts)+",\n")

    fi.close()
    fo.close()

if len(sys.argv) != 4:
    print("This script takes a bed12 + intron col2-col3 in col13-col14 as input, correct the 5ends if they have same first intron splicing junction.")
    print("Then it will output bed12 corrected isoforms, keeping the names of all original isoforms.")
    print("This script will take all input isoforms for calculation, so please make sure they are on the same strand, and can be compared.")
    print("Usage: [BED14Correct5endByFirstIntron.py] [Data.bed14] [Strand] [Output.bed12]")
else:
    Counter()
