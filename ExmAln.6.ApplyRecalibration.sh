#!/bin/bash
#$ -cwd -l mem=12G,time=4:: -N AppRcl


#This script takes a bam file and uses a previously generated base quality score recalibration (BQSR) table to recalibrate them using GATK. If the bam file has previously been split into chromosomes (default 24, i.e. 1-22, X, Y) a list can be provided. The filename of list MUST end ".list"
#	InpFil - (required) - Path to Bam file or a list of BamFiles to be recalibrated
#	RclTab - (required) - Previously generated BQSR table
#	RefFiles - (required) - shell file to export variables with locations of reference files, jar files, and resource directories; see list below
#	TgtBed - (required) - Exome capture kit targets bed file (must end .bed for GATK compatability)
#	LogFil - (optional) - File for logging progress
#	Flag - A - AllowMisencoded - see GATK manual, causes GATK to ignore abnormally high quality scores that would otherwise indicate that the quality score encoding was incorrect
#	Flag - P - PipeLine - call the next step in the pipeline at the end of the job
#	Flag - B - BadET - prevent GATK from phoning home
#	Help - H - (flag) - get usage information

#list of required vairables in reference file:
# $REF - reference genome in fasta format - must have been indexed using 'bwa index ref.fa'
# $EXOMPPLN - directory containing exome analysis pipeline scripts
# $GATK - GATK jar file 
# $ETKEY - GATK key file for switching off the phone home feature, only needed if using the B flag

#list of required tools:
# java <http://www.oracle.com/technetwork/java/javase/overview/index.html>
# GATK <https://www.broadinstitute.org/gatk/> <https://www.broadinstitute.org/gatk/download>

## This file also require exome.lib.sh - which contains various functions used throughout my Exome analysis scripts; this file should be in the same directory as this script

###############################################################

#set default arguments
usage="
ExmAln.6.ApplyRecalibration.sh -i <InputFile> -x <GATK BQSR table> -r <reference_file> -t <targetfile> -l <logfile> -PABH

	 -i (required) - Path to Bam file or \".list\" file containing a multiple paths
	 -x (required) - Previously generated BQSR table
	 -r (required) - shell file to export variables with locations of reference files and resource directories
	 -t (required) - Exome capture kit targets bed file (must end .bed for GATK compatability)
	 -l (optional) - Log file
	 -P (flag) - Call next step of exome analysis pipeline after completion of script
	 -A (flag) - AllowMisencoded - see GATK manual
	 -B (flag) - Prevent GATK from phoning home
	 -H (flag) - echo this message and exit
"

AllowMisencoded="false"
PipeLine="false"
BadET="false"

#get arguments
while getopts i:x:r:t:l:PABH opt; do
	case "$opt" in
		i) InpFil="$OPTARG";;
		x) RclTab="$OPTARG";;
		r) RefFil="$OPTARG";; 
		t) TgtBed="$OPTARG";; 
		l) LogFil="$OPTARG";;
		P) PipeLine="true";;
		A) AllowMisencoded="true";;
		B) BadET="true";;
		H) echo "$usage"; exit;;
	esac
done

#load settings file
source $RefFil

#Load script library
source $EXOMPPLN/exome.lib.sh #library functions begin "func"

#Set Local Variables
ArrNum=$SGE_TASK_ID
funcGetTargetFile
InpNam=`basename $InpFil | sed s/.bam// | sed s/.list//`
RclLst=Recalibrated.$InpNam.list #File listing paths to recalibrated bams
funcFilfromList #if the input is a list get the appropriate input file for this job of the array --> $InpFil
BamFil=`readlink -f $InpFil` #resolve absolute path to bam
BamNam=`basename $BamFil | sed s/.bam//` #a name to use for the various files
if [[ -z "$LogFil" ]];then LogFil=$BamNam.appBQSR.log; fi # a name for the log file
if [[ "$SGE_TASK_ID" == "undefined" ]]; then
	RclFil=$BamNam.recalibrated.bam #file to output recalibrated bam to
else
	RclDir=$InpNam.recalibrated ; mkdir -p $RclDir # a directory to collect individual recalibrated files in 
	RclFil=$RclDir/$BamNam.recalibrated.bam #file to output recalibrated bam to
	Chr=$(echo $SGE_TASK_ID | sed s/24/Y/ | sed s/23/X/)
	if [[ "$BUILD" = "hg19" ]]; then Chr=chr$Chr; fi
fi
GatkLog=$BamNam.AppBQSR.gatklog #a log for GATK to output to, this is then trimmed and added to the script log
TmpLog=$BamNam.AppBQSR.temp.log #temporary log file
TmpDir=$BamNam.AppBQSR.tempdir; mkdir -p $TmpDir #temporary directory
TmpTar=TmpTarFil.$Chr.bed #temporary target file

#Start Log
ProcessName="Recalibrate Base Quality Scores with GATK" # Description of the script - used in log
funcWriteStartLog

#Make chromosome specific exome target file if necessary
if [[ "$SGE_TASK_ID" == "undefined" ]]; then
	TmpTar=$TgtBed
else
	grep -E "^$Chr[[:blank:]]" $TgtBed > $TmpTar
fi

#Apply Recalibration
StepName="Apply recalibration using GATK PrintReads" # Description of this step - used in log
StepCmd="java -Xmx7G -Djava.io.tmpdir=$TmpDir -jar $GATKJAR
 -T PrintReads 
 -R $REF 
 -I $BamFil 
 -BQSR $RclTab
 -L $TmpTar
 -ip 50
 -o $RclFil 
 --filter_mismatching_base_and_quals
 -log $GatkLog" #command to be run
funcGatkAddArguments
funcRunStep

#Call next step
#if the original input was a list then only first job of the array calls the script to merge all the bams and puts it on hold until entire array has finished
#if there was only one bam supplied the depth of coverage and quality distribution scripts are called
if [[ "$SGE_TASK_ID" -eq 1 ]]; then
	#generate realigned file list
	find `pwd` | grep -E bam$ | grep $RclDir | sort -V > $RclLst
	NextJob="Merge Recalibrated Bams"
	QsubCmd="qsub -hold_jid $JOB_ID -o stdostde/ -e stdostde/ $EXOMPPLN/ExmAln.7.MergeBams.sh -i $RclLst -r $RefFil -t $TgtBed -l $LogFil -P"
	if [[ "$AllowMisencoded" == "true" ]]; then QsubCmd=$QsubCmd" -A"; fi
	if [[ "$BadET" == "true" ]]; then QsubCmd=$QsubCmd" -B"; fi
	funcPipeLine
elif [[ "$SGE_TASK_ID" == "undefined" ]]; then
	NextJob="Get Depth of Coverage Statistics"
	QsubCmd="qsub -o stdostde/ -e stdostde/ $EXOMPPLN/ExmAln.8a.DepthofCoverage.sh -i $RclFil -r $RefFil -t $TgtBed -l $LogFil"
	if [[ "$AllowMisencoded" == "true" ]]; then QsubCmd=$QsubCmd" -A"; fi
	if [[ "$BadET" == "true" ]]; then QsubCmd=$QsubCmd" -B"; fi
	funcPipeLine
	NextJob="Get basic bam metrics"
	QsubCmd="qsub -o stdostde/ -e stdostde/ $EXOMPPLN/ExmAln.3a.Bam_metrics.sh -i $RclFil -r $RefFil -l $LogFil -Q"
	funcPipeLine
fi

#End Log
funcWriteEndLog

#Clean up
if [[ -e $RclFil ]]; then rm $BamFil ${BamFil/bam/bai}; fi
if [[ "$TmpTar" != "$TgtBed" ]]; then rm $TmpTar; fi
