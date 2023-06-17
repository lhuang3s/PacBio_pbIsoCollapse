#!/usr/bin/env bash
#pbIsoCollapse.sh
#This pipeline collapses PacBio Iso-seq defined isoforms by collapsing redundant or degraded isoforms, for exon-usage level studies.
#This script is for any cDNA_Cupcake refined isoform bed12 annotations
#Version: Libin Huang, Yu Sun, 2023

if [ ! -n "$4" ]
then
  echo "    This pipeline collapses PacBio Iso-seq defined isoforms by collapsing redundant or degraded isoforms, for exon-usage level studies."
  echo "    This script is for any cDNA_Cupcake refined isoform bed12 annotations"
  echo "    Input isoform bed12: col4: PB.1.1, or with any suffix such as PB.1.1_Both_2_6_8"
  echo "    Input list file: only the base part of the isoforms such as: PB.id"
  echo "    Input abundance file: two columns, isoform name same as the bed12 file, and a number"
  echo "    Usage: `basename $0` [Annotation.bed12] [List|PB.id] [Abundance] [OutputPrefix]"
  echo "    Output: sorted bed12 annotation of collapsed isoforms OutputPrefix.bed12"
  echo "       and a: OriginalName|CorrectedName|SumAbundance matching file, OutputPrefix.Original2Final.matching"
  echo "       and a: CorrectedName|Abundance file, OutputPrefix.abundance.txt"
else
  Data=$1
  List=$2
  Abd=$3
  OutputPrefix=$4
  
  rm -rf $OutputPrefix.bed12 $OutputPrefix.bed12.ttt $OutputPrefix.Original2Final.matching
  
  for gene in `cat $List`;
  do
      echo ">>Processing gene locus: "$gene
      BED12ExtractIsoByPrefixDot.py $Data $gene ${Data}.${gene}
      TableExtractIsoByPrefixDot.py $Abd $gene ${Abd}.${gene}

      #Check single-exon isoforms
      awk '$10==1' ${Data}.${gene} > ${Data}.${gene}.single
      SingeDataLine=`wc -l ${Data}.${gene}.single|awk '{print $1}'`
      if [[ "$SingeDataLine" -gt "0" ]];then
	  		echo "  Processing single exon isoforms"
	  		if [ $SingeDataLine == "1" ];then
	    		  echo "    1 isoform, no update"
	    		  cat ${Data}.${gene}.single > ${Data}.${gene}.single.corrected
	    		  IsoName=`awk '{print $4}' ${Data}.${gene}.single.corrected`
	    		  awk -v t=$IsoName '{OFS="\t";if ($1==t) print $1,$0}' ${Abd}.${gene} > ${Data}.${gene}.single.matching
	  	  elif [[ "$SingeDataLine" -gt "1" ]];then
	  		    awk '{OFS="\t";print $1,$2,$3,$4,$5,$6,$3-$2}' ${Data}.${gene}.single |sort -k7,7nr > ${Data}.${gene}.single.cut
	  		    bedtools merge -i ${Data}.${gene}.single > ${Data}.${gene}.single.merged
	 		      bedtools intersect -a ${Data}.${gene}.single.merged -b ${Data}.${gene}.single.cut -wa -wb > ${Data}.${gene}.single.intersect     #we ignore strands here
	 		      for name in `awk '{print $4}' ${Data}.${gene}.single`;do awk -v t=$name '{if ($1==t) print $0}' ${Abd}.${gene};done > ${Abd}.${gene}.single
	 		      BEDIntersectMergeSingleExonIsoforms.py ${Data}.${gene}.single.intersect ${Abd}.${gene}.single ${Data}.${gene}.single.matching ${Data}.${gene}.single.corrected
	 		      FinalDataLineSingle=`wc -l ${Data}.${gene}.single.corrected|awk '{print $1}'`
	 		      echo "    Collapsed "$SingeDataLine" isoforms to "$FinalDataLineSingle
	 		      #rm -rf ${Data}.${gene}.single.matching_temp ${Data}.${gene}.single.cut ${Data}.${gene}.single.intersect
			  fi
	  	  cat ${Data}.${gene}.single.matching >> $OutputPrefix.matching
	 		  cat ${Data}.${gene}.single.corrected >> $OutputPrefix.bed12
	 		  rm -rf ${Data}.${gene}.single.*
      fi

      #Check multi-exon isoforms:
      awk '$10>1' ${Data}.${gene} > ${Data}.${gene}.multi
      MultiDataLine=`wc -l ${Data}.${gene}.multi|awk '{print $1}'`
      if [[ "$MultiDataLine" -gt "0" ]];then
	  		echo "  Processing multi-exon isoforms"
	  		AllStrands=`awk '{print $6}' ${Data}.${gene}.multi|sort|uniq|wc -l`
	  		if [[ "$AllStrands" -gt "1" ]];then
	  				echo "Warning! Please double check the isoforms on opposite strands!"
	  		fi	
	  		Strand=`head -1 ${Data}.${gene}.multi|awk '{print $6}'`

		  	if [ $MultiDataLine == "1" ];then
		      echo "    1 isoform, no update"
	  	    cat ${Data}.${gene}.multi > ${Data}.${gene}.multi.corrected
	    	  IsoNameM=`awk '{print $4}' ${Data}.${gene}.multi.corrected`
	      	awk -v t=$IsoNameM '{OFS="\t";if ($1==t) print $1,$0}' ${Abd}.${gene} > ${Data}.${gene}.multi.matching
	    	elif [[ "$MultiDataLine" -gt "1" ]];then
	      	echo "    Correcting 5'ends"
		      BED12Extractor.sh -a intron -i ${Data}.${gene}.multi -o ${Data}.${gene}.multi.intron > /dev/null
		      awk '{OFS="\t";print $2,$3}' ${Data}.${gene}.multi.intron > ${Data}.${gene}.multi.intron.col23
	  	    paste ${Data}.${gene}.multi ${Data}.${gene}.multi.intron.col23 > ${Data}.${gene}.multi.withintronpos
	    	  BED14Correct5endByFirstIntron.py ${Data}.${gene}.multi.withintronpos $Strand ${Data}.${gene}.multi.5endcorrected
	      	for name in `awk '{print $4}' ${Data}.${gene}.multi`;do awk -v t=$name '{if ($1==t) print $0}' ${Abd}.${gene};done > ${Abd}.${gene}.multi
		      #paste ${Data}.${gene}.multi.5endcorrected ${Abd}.${gene} |awk '{OFS="\t";print $4,$14}' > ${Data}.${gene}.multi.5endcorrected.abd
		      echo "    Collapsing isoforms"
	  	    BED12Collapse3endIsoforms.py ${Data}.${gene}.multi.5endcorrected ${Abd}.${gene}.multi ${Data}.${gene}.multi.matching ${Data}.${gene}.multi.corrected.t
	    	  echo "    Correcting 3'ends"
	      	BED12Extractor.sh -a intron -i ${Data}.${gene}.multi.corrected.t -o ${Data}.${gene}.multi.correct3end.intron > /dev/null
	      	awk '{OFS="\t";print $2,$3}' ${Data}.${gene}.multi.correct3end.intron > ${Data}.${gene}.multi.correct3end.intron.col23
	      	paste ${Data}.${gene}.multi.corrected.t ${Data}.${gene}.multi.correct3end.intron.col23 > ${Data}.${gene}.multi.correct3end.withintronpos
	      	BED14Correct3endByLastIntron.py ${Data}.${gene}.multi.correct3end.withintronpos $Strand ${Data}.${gene}.multi.corrected

	      	FinalDataLineMulti=`wc -l ${Data}.${gene}.multi.corrected|awk '{print $1}'`
	      	echo "    Collapsed "$MultiDataLine" isoforms to "$FinalDataLineMulti
	    	fi
	  		cat ${Data}.${gene}.multi.matching >> $OutputPrefix.Original2Final.matching
	  		cat ${Data}.${gene}.multi.corrected >> $OutputPrefix.bed12.ttt
	  		rm -rf ${Data}.${gene}.multi.*
      fi
      rm -rf ${Abd}.${gene}.* ${Data}.${gene}.*
      rm -rf ${Abd}.${gene} ${Data}.${gene}
  echo ""
  done

  #Final processing of bed12 and abundance file
  sort -k1,1 -k2,2n $OutputPrefix.bed12.ttt > $OutputPrefix.bed12
  rm -rf $OutputPrefix.bed12.ttt
  ReOrderIsoCollapseOutput.py $OutputPrefix.bed12 $OutputPrefix.Original2Final.matching $OutputPrefix.abundance.txt
  echo "Output generated: "$OutputPrefix.bed12
  echo "Output generated: "$OutputPrefix.Original2Final.matching
  echo "Output generated: "$OutputPrefix.abundance.txt
  echo "Done pbIsoCollapse pipeline."
fi
