#!/bin/bash
#$ -cwd

ChaIn="no"

while getopts i:s:l:c: opt; do
  case "$opt" in
      i) SamFil="$OPTARG";;
      s) Settings="$OPTARG";;
      l) LogFil="$OPTARG";;
	  c) ChaIn="$OPTARG";;
  esac
done

#Load settings file
. $Settings

#Set local variables
TmpDir=$SamFil.Conversiontempdir
mkdir -p $TmpDir

#Start Log
uname -a >> $LogFil
echo "Start Convert SAM to BAM and Dedup - $0:`date`" >> $LogFil
echo " Job name: "$JOB_NAME >> $LogFil
echo " Job ID: "$JOB_ID >> $LogFil
echo "----------------------------------------------------------------" >> $LogFil

#Run Jobs
#convert to bam and sort
echo "- Convert SAM to BAM and reorder using PICARD `date`..."  >> $LogFil
cmd="$JAVA7BIN -Xmx4G -Djava.io.tmpdir=$TmpDir -jar $PICARD/SortSam.jar INPUT=$SamFil.sam OUTPUT=$SamFil.bam SORT_ORDER=coordinate CREATE_INDEX=TRUE"
echo "    "$cmd >> $LogFil
$cmd
if [[ $? == 1 ]]; then
	echo "----------------------------------------------------------------" >> $LogFil
    echo "Convert SAM to BAM and reorder using PICARD failed `date`" >> $LogFil
	qstat -j $JOB_ID | grep -E "usage" >> $LogFil
    exit 1
fi
#Mark the duplicates
echo "- Mark PCR Duplicates using PICARD `date`..." >> $LogFil
BamFil=$SamFil.dedup
cmd="$JAVA7BIN -Xmx4G -Djava.io.tmpdir=$TmpDir -jar $PICARD/MarkDuplicates.jar INPUT=$SamFil.bam OUTPUT=$BamFil.bam METRICS_FILE=$SamFil.Dup.metrics.txt CREATE_INDEX=TRUE"
echo "    "$cmd >> $LogFil
$cmd
if [[ $? == 1 ]]; then
	echo "----------------------------------------------------------------" >> $LogFil
    echo "Mark PCR Duplicates using PICARD failed `date`" >> $LogFil
	qstat -j $JOB_ID | grep -E "usage" >> $LogFil
    exit 1
fi

# raw flagstat
echo "- Output flag stats using `date`..." >> $LogFil
echo "    $SAMTOOLS flagstat $BamFil.bam > $BamFil.flagstat" >> $LogFil
$SAMTOOLS flagstat $BamFil.bam > $BamFil.flagstat
echo "- Output idx stats using `date`..." >> $LogFil
echo "    $SAMTOOLS idxstats $BamFil.bam > $BamFil.idxstats" >> $LogFil
$SAMTOOLS idxstats $BamFil.bam > $BamFil.idxstats
echo "----------------------------------------------------------------" >> $LogFil

#Call Next Job if chain
if [[ $ChaIn = "chain" ]]; then
	echo "- Call GC Metrics `date`:" >> $LogFil
	JobNm=${JOB_NAME#*.}
	BamLst=$BamFil"_bam_files.list"
	echo $BamFil.bam > $BamLst
	cmd="qsub -l $GCstatAlloc -N GCstat.$JobNm -o stdostde/ -e stdostde/ $EXOMSCR/ExmAln.4a.GC_metrics.sh -i $BamFil -s $Settings -l $LogFil"
	echo "    "$cmd  >> $LogFil
	$cmd
	echo "- Call GATK realign, 1 job for each Chr `date`:" >> $LogFil
	cmd="qsub -pe smp $NumCores -t 1-24 -l $realnAlloc -N realn.$JobNm -o stdostde/ -e stdostde/ $EXOMSCR/ExmAln.4.LocalRealignment.sh -i $BamFil -b $BamLst -s $Settings -l $LogFil -c chain"
	echo "    "$cmd >> $LogFil
	$cmd
	echo "----------------------------------------------------------------" >> $LogFil
fi

#End Log
echo "End Convert SAM to BAM and dedup $0:`date`" >> $LogFil
qstat -j $JOB_ID | grep -E "usage" >> $LogFil
echo "===========================================================================================" >> $LogFil
echo "" >> $LogFil

#Remove temp files
rm -r $TmpDir $SamFil.sam $SamFil.bam 
