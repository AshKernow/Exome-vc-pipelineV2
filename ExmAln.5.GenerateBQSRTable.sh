#!/bin/bash
#$ -cwd -pe smp 6 -l mem=2G,time=6:: -N GenBQSR


#This script takes a bam file or a list of bam files (filename must end ".list") and generates the base quality score recalibration table using GATK
#	InpFil - (required) - Path to Bam file or a list of BamFiles to be recalibrated
#	RefFiles - (required) - shell file to export variables with locations of reference files, jar files, and resource directories; see list below
#	TgtBed - (required) - Exome capture kit targets bed file (must end .bed for GATK compatability)
#	LogFil - (optional) - File for logging progress
#	Flag - A - AllowMisencoded - see GATK manual, causes GATK to ignore abnormally high quality scores that would otherwise indicate that the quality score encoding was incorrect
#	Flag - P - PipeLine - call the next step in the pipeline at the end of the job
#	Flag - B - BadET - prevent GATK from phoning home

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
AllowMisencoded="false"
FixMisencoded="false"
PipeLine="false"
BadEt="false"

#get arguments
while getopts i:r:t:l:PAFB opt; do
	case "$opt" in
		i) InpFil="$OPTARG";;
		r) RefFil="$OPTARG";; 
		t) TgtBed="$OPTARG";; 
		l) LogFil="$OPTARG";;
		P) PipeLine="true";;
		A) AllowMisencoded="true";;
		F) FixMisencoded="true";;
		B) BadET="true";;
	esac
done

#load settings file
source $RefFil

#Load script library
source $EXOMPPLN/exome.lib.sh #library functions begin "func"

#Set Local Variables
BamFil=`readlink -f $InpFil` #resolve absolute path to bam
BamNam=`basename ${BamFil/.bam/}` #a name to use for the various files
BamNam=${BamNam/.list/} 
TmpDir=$BamNam.GenBQSRjavdir #temp directory for java machine
mkdir -p $TmpDir
if [[ -z $LogFil ]];then
	LogFil=$BamNam.GenBQSR.log # a name for the log file
fi
TmpLog=$LogFil.GenBQSR.log #temporary log file 
RclTable=$BamFil.recal.table # output - base quality score recalibration table
GatkLog=$BamNam.gatklog #a log for GATK to output to, this is then trimmed and added to the script log

#Start Log
ProcessName="Start Generate Base Quality Score Recalibration Table with GATK" # Description of the script - used in log
funcWriteStartLog

#Generate target file
StepName="Create recalibration data file using GATK BaseRecalibrator" # Description of this step - used in log
StepCmd="java -Xmx7G -Djava.io.tmpdir=$TmpDir -jar $GATKJAR
 -T BaseRecalibrator 
 -R $REF 
 -L $TARGET 
 -I $BamFil 
 -knownSites $DBSNP 
 -knownSites $INDEL 
 -knownSites $INDEL1KG 
 -o $RclTable 
 -nct 6
 --filter_mismatching_base_and_quals
 -log $GatkLog" #command to be run
funcGatkAddArguments
funcRunStep

#End Log
funcWriteEndLog
#Clean up
rm -r $TmpDir $TmpLog
