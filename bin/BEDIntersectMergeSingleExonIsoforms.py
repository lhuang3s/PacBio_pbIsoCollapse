#!/usr/bin/python3
#BEDIntersectMergeSingleExonIsoforms.py
#This script processes bedtools intersect results for PacBio isoform collapsing, for single-exon isoforms
#Version: Yu Sun, 2022-10-23, Updated to support multiple abundance columns
                                                                          
import sys

def Convert2Int(List):
    for i in range(len(List)):
        List[i] = int(List[i])
    return(List)

def SumTwoLists(ListA, ListB):
    ListC=[]
    for i in range(len(ListA)):
        ListC.append(int(ListA[i]) + int(ListB[i]))
    return(ListC)

def Counter():
    fi=open(sys.argv[1],'r')
    fabd=open(sys.argv[2],'r')
    fomatching=open(sys.argv[3]+"_temp",'w')
    fomatchingfinal=open(sys.argv[3],'w')
    fobed12=open(sys.argv[4],'w')

    Locations={}

    for line in fi:
        LineElement=line.strip().split()
        Strand=LineElement[8]
        Corr=LineElement[0]+";"+LineElement[1]+";"+LineElement[2]+";"+Strand
        if len(Locations) == 0:
            Locations[Corr]=LineElement[6]
            fomatching.write(LineElement[6]+"\t"+LineElement[6]+"\n")
        else:
            if Corr not in Locations.keys():
                Locations[Corr]=LineElement[6]
                fomatching.write(LineElement[6]+"\t"+LineElement[6]+"\n")
            else:
                fomatching.write(Locations[Corr]+"\t"+LineElement[6]+"\n")

    Locations_rev={v: k for k, v in Locations.items()}
    fomatching.close()

    AbundanceDic={}
    for abd in fabd:
        tx=abd.strip().split()[0]
        abundance=abd.strip().split()[1:]
        AbundanceDic[tx]=abundance

    Old2NewName={}
    FinalAbdSum={}   #final matching to the updated name, with reads summed
    fiupdate=open(sys.argv[3]+"_temp",'r')
    for tx in fiupdate:
        txname=tx.strip().split()[0]
        subname=tx.strip().split()[1]
        if subname not in Old2NewName.keys():
            Old2NewName[subname]=txname
        if txname not in FinalAbdSum.keys():
            FinalAbdSum[txname]=Convert2Int(AbundanceDic[subname])
        else:
            FinalAbdSum[txname]=SumTwoLists(FinalAbdSum[txname], AbundanceDic[subname])

    fabd.seek(0)    #Go back to the beginning, and loop the file again
    for tx in fabd:
        txname=tx.strip().split()[0]
        NewName=Old2NewName[txname]
        NewAbd=FinalAbdSum[NewName]
        fomatchingfinal.write(txname+"\t"+NewName+"\t"+"\t".join([str(x) for x in NewAbd])+"\n")

    for name in Locations.keys():
        Corrf=name.split(";")
        fobed12.write(Corrf[0]+"\t"+Corrf[1]+"\t"+Corrf[2]+"\t"+Locations[name]+"\t0\t"+Corrf[3]+"\t"+Corrf[1]+"\t"+Corrf[2]+"\t0\t1\t"+str(int(Corrf[2])-int(Corrf[1]))+",\t0,\n")

    fobed12.close()
    fiupdate.close()

if len(sys.argv) != 5:
    print("This script processes bedtools intersect results for PacBio isoform collapsing, for single-exon isoforms")
    print("Usage: [BEDIntersectMergeSingleExonIsoforms.py] [Data.intersect.bed] [Data.abundance] [OutputMatching] [OutputBed12]")
else:
    Counter()
