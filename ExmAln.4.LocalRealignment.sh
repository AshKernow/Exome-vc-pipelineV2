#!/bin/bash
#$ -cwd -t 1-24 -pe smp 6 -l mem=2G,time=2:: -N LocRln


#This script takes a bam file andperfomrs local indel realignment using GATK, the script runs as an array across 24 Chromosomes
#	InpFil - (required) - Path to Bam file to be realigned
#	RefFiles - (required) - shell file to export variables with locations of reference files, jar files, and resource directories; see list below
#	LogFil - (optional) - File for logging progress
#	Flag - A - AllowMisencoded - see GATK manual (https://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_CommandLineGATK.html#--allow_potentially_misencoded_quality_scores), causes GATK to ignore abnormally high quality scores that would otherwise indicate that the quality score encoding was incorrect
#	Flag - P - PipeLine - call the next step in the pipeline at the end of the job
#	Flag - B - BadET - prevent GATK from phoning home
#	Help - H - (flag) - get usage information

#list of required vairables in reference file:
# $REF - reference genome in fasta format - must have been indexed using 'bwa index ref.fa'
# $INDEL - Gold standard INDEL reference from GATK
# $INDEL1KG - INDEL reference from 1000 genomes
# $EXOMPPLN - directory containing exome analysis pipeline scripts
# $GATK - GATK jar file 
# $ETKEY - GATK key file for switching off the phone home feature, only needed if using the B flag

#list of required tools:
# java
# GATK

## This file also require exome.lib.sh - which contains various functions used throughout my Exome analysis scripts; this file should be in the same directory as this script

###############################################################

#set default arguments
usage="
ExmAln.8a.DepthofCoverage.sh -i <InputFile> -r <reference_file> -t <targetfile> -l <logfile> -GIQH

	 -i (required) - Path to Bam file to be aligned or \".list\" file containing a multiple paths
	 -r (required) - shell file to export variables with locations of reference files and resource directories
	 -l (optional) - Log file
	 -P (flag) - Call next step of exome analysis pipeline after completion of script
	 -A (flag) - AllowMisencoded - see GATK manual
	 -B (flag) - Prevent GATK from phoning home
	 -H (flag) - echo this message and exit
"

AllowMisencoded="false"
PipeLine="false"
BadEt="false"

#get arguments
while getopts i:r:l:PABH opt; do
	case "$opt" in
		i) InpFil="$OPTARG";;
		r) RefFil="$OPTARG";; 
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
Chr=$SGE_TASK_ID
Chr=${Chr/23/X}
Chr=${Chr/24/Y}
if [[ "$BUILD" = "hg19" ]]; then
	Chr=chr$Chr
fi
BamFil=`readlink -f $InpFil` #resolve absolute path to bam
BamNam=`basename ${BamFil/.bam/}` #a name to use for the various files
TmpDir=$BamNam.$Chr.LocRealignjavdir #temp directory for java machine
mkdir -p $TmpDir
if [[ -z $LogFil ]];then
	LogFil=$BamNam.LocReal.log # a name for the log file
fi
TmpLog=$LogFil.LocReal.$Chr.log #temporary log file 
RalDir=LocRealign.$BamNam #directory to collect individual chromosome realignments
mkdir -p $RalDir
RalLst=LocRealign.$BamNam.list #File listing paths to individual chromosome realignments
StatFil=RealignStatus_$Chr.$JOB_ID.LocReal.stat #Status file to check if all chromosome are complete
TgtFil=$RalDir/$BamNam.$Chr.target_intervals.list #temporary file for target intervals for the CHR
realignedFile=$RalDir/realigned.$BamNam.Chr_$Chr.bam # the output - a realigned bam file for the CHR
GatkLog=$BamNam.$Chr.LocRealign.gatklog #a log for GATK to output to, this is then trimmed and added to the script log

#Start Log
ProcessName="Start Local Realignment around InDels on Chromosome $Chr" # Description of the script - used in log
funcWriteStartLog

#Generate target file
StepName="Create target interval file using GATK RealignerTargetCreator" # Description of this step - used in log
StepCmd="java -Xmx7G -Djava.io.tmpdir=$TmpDir -jar $GATKJAR
 -T RealignerTargetCreator
 -R $REF
 -I $BamFil
 -L $Chr
 -known $INDEL
 -known $INDEL1KG
 -o $TgtFil
 -nt 6
 -log $GatkLog" #command to be run
funcGatkAddArguments
funcRunStep

#Realign InDels
StepName="Realign InDels file using GATK IndelRealigner" # Description of this step - used in log
StepCmd="java -Xmx7G -Djava.io.tmpdir=$TmpDir -jar $GATKJAR
 -T IndelRealigner
 -R $REF
 -I $BamFil
 -targetIntervals $TgtFil
 -L $Chr
 -known $INDEL
 -known $INDEL1KG
 -o $realignedFile
 -log $GatkLog" #command to be run
funcGatkAddArguments
funcRunStep

#generate realigned file list
find `pwd` | grep -E bam$ | grep $RalDir | sort -V > $RalLst

#End Log
funcWriteEndLog
#Clean up
rm -r $TmpDir $TmpLog $TgtFil
