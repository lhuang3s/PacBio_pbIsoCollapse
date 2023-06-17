#!/usr/bin/python3
#BED12Collapse3endIsoforms.py
#This script takes the 5end corrected isoforms using BED14Correct5endByFirstIntron.py, and collapse the isoforms that have exactly the same, 
#    or a subset of exon structure, and merge them into the longest one
#Single-exon genes need to be removed! This script only supports multi-exon transcript processing
#Version: 2022-10-23, Updated to support multiple abundance columns

import sys
import re

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
    folist=open(sys.argv[3],'w')
    fobed12=open(sys.argv[4],'w')

    LongestIsoList=[]

    #Loop through the bed12 isoform annotations with corrected 5end, and find the longest isoforms
    for line in fi:
        Data=line.strip().split()
        if len(LongestIsoList) == 0:
            LongestIsoList.append(Data)
        else:
            if Data[5] == "+":
                found=0
                for i in range(len(LongestIsoList)):
                    if Data[1] == LongestIsoList[i][1]:
                        LongestIsoBED12=LongestIsoList[i][11]
                        DataBED12=Data[11]
                        if (LongestIsoBED12 in DataBED12) or (DataBED12 in LongestIsoBED12):
                            LongestIsoBED11=re.sub(',$','',LongestIsoList[i][10]).split(",")
                            DataBED11=re.sub(',$','',Data[10]).split(",")
                            if len(LongestIsoBED11) > len(DataBED11):
                                DataBED11_firstNminus1=",".join([str(x) for x in DataBED11[:-1]])+","
                                if (DataBED11_firstNminus1 in LongestIsoList[i][10]) and (int(DataBED11[len(DataBED11)-1]) <= int(LongestIsoBED11[len(DataBED11)-1])):
#                                    print("This is a subset of the current longest one")
                                    found=1
                            elif len(LongestIsoBED11) < len(DataBED11):
                                LongestIsoBED11_firstNminus1=",".join([str(x) for x in LongestIsoBED11[:-1]])+","
                                if (LongestIsoBED11_firstNminus1 in Data[10]) and (int(LongestIsoBED11[len(LongestIsoBED11)-1]) <= int(DataBED11[len(LongestIsoBED11)-1])):
#                                    print("Update the current longest anno")
                                    LongestIsoList[i]=Data
                                    found=1
                            else:
                                DataBED11_firstNminus1=",".join([str(x) for x in DataBED11[:-1]])+","
                                LongestIsoBED11_firstNminus1=",".join([str(x) for x in LongestIsoBED11[:-1]])+","
                                if DataBED11_firstNminus1==LongestIsoBED11_firstNminus1:
                                    if int(DataBED11[-1]) <= int(LongestIsoBED11[-1]):
#                                        print("This is a subset of the current longest one")
                                        found=1
                                    else:
#                                        print("Update the current longest anno")
                                        LongestIsoList[i]=Data
                                        found=1
#                print(found)
                if found==0:
                    LongestIsoList.append(Data)  #This is a new pattern
            else:
                found=0
                for i in range(len(LongestIsoList)):
                    if Data[2] == LongestIsoList[i][2]:
                        Data3end=int(Data[1])
                        LongestIso3end=int(LongestIsoList[i][1])
                        LongestIsoBED12=LongestIsoList[i][11]
                        LongestIsoBED11=LongestIsoList[i][10]
                        DataBED12=Data[11]
                        DataBED11=Data[10]
                        LongestIsoBED12list=re.sub(',$','',LongestIsoBED12).split(",")
                        LongestIsoBED11list=re.sub(',$','',LongestIsoBED11).split(",")
                        DataBED12list=re.sub(',$','',DataBED12).split(",")
                        DataBED11list=re.sub(',$','',DataBED11).split(",")
                        MergedLongestBED12list=[]
                        MergedDataBED12list=[]
#                        print(LongestIsoList[i][2])
                        for j in range(len(LongestIsoBED12list)):
                            MergedLongestBED12list.append(int(LongestIsoBED12list[j])+int(LongestIsoBED11list[j])+LongestIso3end)
                        for j in range(len(DataBED12list)):
                            MergedDataBED12list.append(int(DataBED12list[j])+int(DataBED11list[j])+Data3end)
#                        print(MergedLongestBED12list)
#                        print(MergedDataBED12list)
                        MergedLongestBED12=",".join([str(x) for x in MergedLongestBED12list])
                        MergedDataBED12=",".join([str(x) for x in MergedDataBED12list])
#                        print(MergedLongestBED12,MergedDataBED12)
                        if (MergedLongestBED12 in MergedDataBED12) or (MergedDataBED12 in MergedLongestBED12):
                            LongestIsoBED11listrev=LongestIsoBED11list[::-1]
                            LongestIsoBED11listrev_cat=",".join([str(x) for x in LongestIsoBED11listrev])+","
                            DataBED11listrev=DataBED11list[::-1]
                            DataBED11listrev_cat=",".join([str(x) for x in DataBED11listrev])+","
                            if len(LongestIsoBED11listrev) > len(DataBED11listrev):
                                DataBED11listrev_firstNminus1=",".join([str(x) for x in DataBED11listrev[:-1]])+","
                                if (DataBED11listrev_firstNminus1 in LongestIsoBED11listrev_cat) and (int(DataBED11listrev[-1]) <= int(LongestIsoBED11listrev[len(DataBED11listrev)-1])):
#                                    print("This is a subset of the current longest one")
                                    found=1
                            elif len(LongestIsoBED11listrev) < len(DataBED11listrev):
                                LongestIsoBED11listrev_firstNminus1=",".join([str(x) for x in LongestIsoBED11listrev[:-1]])+","
                                if (LongestIsoBED11listrev_firstNminus1 in DataBED11listrev_cat) and (int(LongestIsoBED11listrev[-1]) <= int(DataBED11listrev[len(LongestIsoBED11listrev)-1])):
#                                    print("Update the current longest anno")
                                    LongestIsoList[i]=Data
                                    found=1
                            else:
                                DataBED11listrev_firstNminus1=",".join([str(x) for x in DataBED11listrev[:-1]])+","
                                LongestIsoBED11listrev_firstNminus1=",".join([str(x) for x in LongestIsoBED11listrev[:-1]])+","
                                if DataBED11listrev_firstNminus1 == LongestIsoBED11listrev_firstNminus1:
                                    if int(DataBED11listrev[-1]) <= int(LongestIsoBED11listrev[-1]):
#                                        print("This is a subset of the current longest one")
                                        found=1
                                    else:
#                                        print("Update the current longest anno")
                                        LongestIsoList[i]=Data
                                        found=1
                if found==0:
                    LongestIsoList.append(Data)

#    print(len(LongestIsoList))
#    print(LongestIsoList)
    Matching={} #Old isoforms -> Collapsed longest isoform names
    fi.seek(0)  #Go back and loop the file again. Now we know the longest isoforms, we just need to map all isoforms back to them.
    for line in fi:
        Data=line.strip().split()
#        print(Data)
        for i in range(len(LongestIsoList)):
            if Data[5] == "+":
                if Data[1] == LongestIsoList[i][1]:
                    LongestIsoBED12=LongestIsoList[i][11]
                    DataBED12=Data[11]
                    if DataBED12 in LongestIsoBED12:
                        LongestIsoBED11=re.sub(',$','',LongestIsoList[i][10]).split(",")
                        DataBED11=re.sub(',$','',Data[10]).split(",")
                        if len(LongestIsoBED11) > len(DataBED11):
                            DataBED11_firstNminus1=",".join([str(x) for x in DataBED11[:-1]])+","
                            if (DataBED11_firstNminus1 in LongestIsoList[i][10]) and (int(DataBED11[len(DataBED11)-1]) <= int(LongestIsoBED11[len(DataBED11)-1])):
                                Matching[Data[3]]=LongestIsoList[i][3]
                        elif len(LongestIsoBED11) == len(DataBED11):
                            DataBED11_firstNminus1=",".join([str(x) for x in DataBED11[:-1]])+","
                            LongestIsoBED11_firstNminus1=",".join([str(x) for x in LongestIsoBED11[:-1]])+","
                            if DataBED11_firstNminus1==LongestIsoBED11_firstNminus1:
                                if int(DataBED11[-1]) <= int(LongestIsoBED11[-1]):
                                    Matching[Data[3]]=LongestIsoList[i][3]
            else:
                if Data[2] == LongestIsoList[i][2]:
                    Data3end=int(Data[1])
                    LongestIso3end=int(LongestIsoList[i][1])
                    LongestIsoBED12=LongestIsoList[i][11]
                    LongestIsoBED11=LongestIsoList[i][10]
                    DataBED12=Data[11]
                    DataBED11=Data[10]
                    LongestIsoBED12list=re.sub(',$','',LongestIsoBED12).split(",")
                    LongestIsoBED11list=re.sub(',$','',LongestIsoBED11).split(",")
                    DataBED12list=re.sub(',$','',DataBED12).split(",")
                    DataBED11list=re.sub(',$','',DataBED11).split(",")
                    MergedLongestBED12list=[]
                    MergedDataBED12list=[]
#                   print(LongestIsoList[i][2])                                                                                                                                                                                                         
                    for j in range(len(LongestIsoBED12list)):
                        MergedLongestBED12list.append(int(LongestIsoBED12list[j])+int(LongestIsoBED11list[j])+LongestIso3end)
                    for j in range(len(DataBED12list)):
                        MergedDataBED12list.append(int(DataBED12list[j])+int(DataBED11list[j])+Data3end)
                    MergedLongestBED12=",".join([str(x) for x in MergedLongestBED12list])
                    MergedDataBED12=",".join([str(x) for x in MergedDataBED12list])
                    if (MergedDataBED12 in MergedLongestBED12):
                        LongestIsoBED11listrev=LongestIsoBED11list[::-1]
                        LongestIsoBED11listrev_cat=",".join([str(x) for x in LongestIsoBED11listrev])+","
                        DataBED11listrev=DataBED11list[::-1]
                        DataBED11listrev_cat=",".join([str(x) for x in DataBED11listrev])+","
                        if len(LongestIsoBED11listrev) > len(DataBED11listrev):
                            DataBED11listrev_firstNminus1=",".join([str(x) for x in DataBED11listrev[:-1]])+","
                            if (DataBED11listrev_firstNminus1 in LongestIsoBED11listrev_cat) and (int(DataBED11listrev[-1]) <= int(LongestIsoBED11listrev[len(DataBED11listrev)-1])):
                                Matching[Data[3]]=LongestIsoList[i][3]
                        elif len(LongestIsoBED11listrev) == len(DataBED11listrev):
                            DataBED11listrev_firstNminus1=",".join([str(x) for x in DataBED11listrev[:-1]])+","
                            LongestIsoBED11listrev_firstNminus1=",".join([str(x) for x in LongestIsoBED11listrev[:-1]])+","
                            if DataBED11listrev_firstNminus1==LongestIsoBED11listrev_firstNminus1:
                                if int(DataBED11listrev[-1]) <= int(LongestIsoBED11listrev[-1]):
                                    Matching[Data[3]]=LongestIsoList[i][3]

    Adundance={}  #Corrected 5ends -> Abundance
    for abd in fabd:
        data=abd.strip().split()
        Adundance[data[0]]=data[1:]

    FinalAbd_ma={}  #Final names -> final abundance, summing up all sub-isoforms
    for name in Matching.keys():
        CollapsedIso=Matching[name]
        if CollapsedIso not in FinalAbd_ma.keys():
            FinalAbd_ma[CollapsedIso]=Convert2Int(Adundance[name])
        else:
            FinalAbd_ma[CollapsedIso]=SumTwoLists(FinalAbd_ma[CollapsedIso], Adundance[name])

    fabd.seek(0)
    for line in fabd:
        Data=line.strip().split()
        DataMatching=Matching[Data[0]]
        folist.write(Data[0]+"\t"+DataMatching+"\t"+"\t".join([str(x) for x in FinalAbd_ma[DataMatching]])+"\n")

    for finalname in FinalAbd_ma.keys():
        for longiso in LongestIsoList:
            if longiso[3] == finalname:
                fobed12.write("\t".join(longiso)+"\n")

    fi.close()
    fabd.close()
    folist.close()
    fobed12.close()

if len(sys.argv) != 5:
    print("This script takes the 5end corrected isoforms using BED14Correct5endByFirstIntron.py, and collapse the isoforms that have exactly the same")
    print("    or a subset of exon structure, and merge them into the longest one")
    print("Single-exon genes need to be removed! This script only supports multi-exon transcript processing")
    print("Usage: [BED12Collapse3endIsoformsGeneral.py] [Input.bed12] [Abundance|two or more columns] [Output.list] [Output.bed12]")
else:
    Counter()
