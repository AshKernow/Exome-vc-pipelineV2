#!/bin/bash
#$ -cwd

ChaIn="no"

while getopts i:s:l:c: opt; do
  case "$opt" in
      i) BamFil="$OPTARG";;
      s) Settings="$OPTARG";;
      l) LogFil="$OPTARG";;
	  c) ChaIn="$OPTARG";;
  esac
done

#Load settings file
. $Settings

#Set local variables
TmpDir=$BamFil.Stampytempdir
mkdir -p $TmpDir

#Start Log
uname -a >> $LogFil
echo "Align BAM with Stampy - $0:`date`" >> $LogFil
echo " Job name: "$JOB_NAME >> $LogFil
echo " Job ID: "$JOB_ID >> $LogFil
echo "----------------------------------------------------------------" >> $LogFil

#Run Jobs
#Align with Stampy
echo "- Align BAM with Stampy `date`..."  >> $LogFil
cmd="$STAMPY -g $STIDX -h $STHSH --bamkeepgoodreads --bwamark -t 4 -M $BamFil.bam -v 3 -o $BamFil.stampy.sam"
echo "    "$cmd >> $LogFil
$cmd
if [[ $? == 1 ]]; then
	echo "----------------------------------------------------------------" >> $LogFil
    echo "Align BAM with Stampy failed `date`" >> $LogFil
	qstat -j $JOB_ID | grep -E "usage" >> $LogFil
    exit 1
fi
SamFil=$BamFil.stampy

#Call Next Job if chain
if [[ $ChaIn = "chain" ]]; then
	echo "- Call Convert SAM to BAM and dedup `date`:" >> $LogFil
	JobNm=${JOB_NAME#*.}
	cmd="qsub -l $ConvS2BAlloc -N Sm2Bm2.$JobNm -o stdostde/ -e stdostde/ $EXOMSCR/ExmAln.TEST.ConvSamtoBamandDedup.sh -i $SamFil -s $Settings -l $LogFil -c chain"
	echo "    "$cmd  >> $LogFil
	$cmd
	echo "----------------------------------------------------------------" >> $LogFil
fi

#End Log
echo "End Align BAM with Stampy $0:`date`" >> $LogFil
qstat -j $JOB_ID | grep -E "usage" >> $LogFil
echo "===========================================================================================" >> $LogFil
echo "" >> $LogFil

#Remove temp files
#rm -r $TmpDir $BamFil.bam $BamFil.bai 
