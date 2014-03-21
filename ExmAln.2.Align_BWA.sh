#!/bin/bash
#$ -cwd
#This script aligns a single/pair of fastq file(s) against the provided reference using BWA-MEM
#Inputs:
#	FqFil - the primary input is a table containing the fastq file name and the RG read header for the output SAM file.
#      - Columns should be as follows:
#            For single end - <fastq files> <readgroup headers>
#            For paired end - <Read 1 fastq files> <readgroup headers> <Read 2 fastq files>
#	FqDir - The directory containing the fastq files(s)
#	Settings - Global settings file
#	LogFil - Log file

ChaIn="no"

while getopts i:f:s:l:c: opt; do
  case "$opt" in
      i) FqFil="$OPTARG";;
	  f) FqDir="$OPTARG";;
      s) Settings="$OPTARG";;
      l) LogFil="$OPTARG";;
	  c) ChaIn="$OPTARG";;
  esac
done

#load settings file
. $Settings

#set local variables
SamNum=$SGE_TASK_ID
TmpLog=$LogFil.$SamNum
NCOL=$(head -n1 $FqFil | wc -w | cut -d" " -f1)
fastq1=$(tail -n+$SamNum $FqFil | head -n 1 | cut -f1)
rgheader=$(tail -n+$SamNum $FqFil | head -n 1 | cut -f2)
rgID=${rgheader##*ID:}
rgID=${rgID%%\\tSM*}
AlignDir=$rgID.align
mkdir -p $AlignDir
if [ $NCOL -eq 3 ]; then
	fastq2=$(tail -n+$SamNum $FqFil | head -n 1 | cut -f3)
fi

#start log
cat $LogFil > $TmpLog
LogFil=${LogFil/.log/}.$rgID.log
cat $TmpLog > $LogFil
uname -a >> $LogFil
echo "Start Align with BWA - $0:`date`" >> $LogFil
echo " Job name: "$JOB_NAME >> $LogFil
echo " Job ID: "$JOB_ID >> $LogFil
echo " Array task ID: "$SGE_TASK_ID >> $LogFil
echo " fastq file 1: $fastq1" >> $LogFil
echo " fastq file 2: $fastq2" >> $LogFil
echo " Sample Name: "$rgID >> $LogFil
echo " Sample File Size: "$(wc $FqDir/$fastq1)" lines" >> $LogFil
echo " RG header: "$rgheader >> $LogFil
echo "----------------------------------------------------------------" >> $LogFil

#Run Jobs
#Align using BWA mem algorithm
echo "- Align with BWA mem `date`...">> $LogFil
if [ $NCOL -eq 3 ]; then
	echo "    $BWA mem -M -R \"$rgheader\" -t 4 $REF $FqDir/$fastq1 $FqDir/$fastq2 > $rgID.sam" >> $LogFil
	$BWA mem -M -t 4 -R "$rgheader" $REF $FqDir/$fastq1 $FqDir/$fastq2 > $AlignDir/$rgID.sam
else
	echo "    $BWA mem -M -R \"$rgheader\" -t 2 $REF $FqDir/$fastq1 > $rgID.sam" >> $LogFil
	$BWA mem -M -t 2 -R "$rgheader" $REF $FqDir/$fastq1 > $AlignDir/$rgID.sam
fi

if [[ $? == 1 ]]; then
	echo "----------------------------------------------------------------" >> $LogFil
    echo "Align with BWA mem failed `date` " >> $LogFil
    qstat -j $JOB_ID | grep -E "usage *$SGE_TASK_ID" >> $LogFil
	exit 1
fi

#move into the Sample alignment directory
mv $LogFil $AlignDir/$LogFil
cd $AlignDir
mkdir stdostde
echo "----------------------------------------------------------------" >> $LogFil

#Call next jobs if chain
if [[ $ChaIn = "chain" ]]; then
	echo "- Call Convert SAM to BAM and deduplication `date`:" >> $LogFil
	JobNm=${JOB_NAME#*.}
	cmd="qsub -l $ConvS2BAlloc -N Sm2Bm.$JobNm -o stdostde/ -e stdostde/ $EXOMSCR/ExmAln.3.ConvSamtoBam.sh -i $rgID -s $Settings -l $LogFil -c chain"
	echo "    "$cmd  >> $LogFil
	$cmd
	echo "----------------------------------------------------------------" >> $LogFil
fi

#End log
echo "End Align with BWA $0:`date`" >> $LogFil
qstat -j $JOB_ID | grep -E "usage *$SGE_TASK_ID" >> $LogFil
echo "===========================================================================================" >> $LogFil
echo "" >> $LogFil

#remove temporary files
rm ../$TmpLog