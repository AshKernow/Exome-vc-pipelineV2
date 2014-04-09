#!/bin/bash
#$ -cwd -l mem=12G,time=2:: -N MrgBam


#This script a list of bam files and merges them into a single file. The filename of list MUST end ".list"
#	InpFil - (required) - A list of BamLstes to be merged
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
ExmAln.7.MergeBams.sh -i <InputFile> -r <reference_file> -t <targetfile> -l <logfile> -GIQH

	 -i (required) - Path to Bam file or \".list\" file containing a multiple paths
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
while getopts i:r:t:l:PABH opt; do
	case "$opt" in
		i) InpFil="$OPTARG";;
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
BamLst=`readlink -f $InpFil` #resolve absolute path to bam
BamNam=`basename $BamFil | sed s/.bam//` #a name to use for the various files
if [[ -z "$LogFil" ]];then LogFil=$BamNam.MrgBam.log; fi # a name for the log file
MrgFil=$BamNam.merged.bam #file to output
GatkLog=$BamNam.MrgBam.gatklog #a log for GATK to output to, this is then trimmed and added to the script log
TmpLog=$BamNam.MrgBam.temp.log #temporary log file
TmpDir=$BamNam.MrgBams.tempdir; mkdir -p $TmpDir #temporary directory

#Start Log
ProcessName="Merge Bams with GATK" # Description of the script - used in log
funcWriteStartLog

#Merge Bams
StepName="Merge Bams using GATK PrintReads" # Description of this step - used in log
StepCmd="java -Xmx7G -Djava.io.tmpdir=$TmpDir -jar $GATKJAR
 -T PrintReads
 -R $REF
 -I $BamLst 
 -L $TgtBed
 -ip 50
 -o $MrgFil
 -log $GatkLog" #command to be run
funcGatkAddArguments
funcRunStep

#Call next step
NextJob="Get Depth of Coverage Statistics"
QsubCmd="qsub -o stdostde/ -e stdostde/ $EXOMPPLN/ExmAln.8a.DepthofCoverage.sh -i $MrgFil -r $RefFil -t $TgtBed -l $LogFil"
if [[ "$AllowMisencoded" == "true" ]]; then QsubCmd=$QsubCmd" -A"; fi
if [[ "$BadET" == "true" ]]; then QsubCmd=$QsubCmd" -B"; fi
funcPipeLine
NextJob="Get basic bam metrics"
QsubCmd="qsub -o stdostde/ -e stdostde/ $EXOMPPLN/ExmAln.3a.Bam_metrics.sh -i $MrgFil -r $RefFil -l $LogFil -Q"
funcPipeLine

#End Log
funcWriteEndLog

#Clean up
if [[ -s $MrgFil ]]; then rm $(cat $BamLst); fi
