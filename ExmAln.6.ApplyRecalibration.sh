#!/bin/bash
#$ -cwd -t 1-24 -pe smp 6 -l mem=2G,time=6:: -N LocRln


#This script takes a bam file and uses a previously generated base quality score recalibration (BQSR) table to recalibrate them using GATK. If the bam file has previously been split into chromosomes (default 24, i.e. 1-22, X, Y) a list can be provided, ensuring each the files are listed in chromosome order (1-22 then X then Y). The filename of list MUST end ".list"
#	InpFil - (required) - Path to Bam file or a list of BamFiles to be recalibrated
#	RclTab - (required) - Previously generated BQSR table
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
while getopts i:x:r:t:l:PAB opt; do
	case "$opt" in
		i) InpFil="$OPTARG";;
		x) RclTab="$OPTARG";;
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
Chr=$SGE_TASK_ID
Chr=${Chr/23/X}
Chr=${Chr/24/Y}
if [[ "$BUILD" = "hg19" ]]; then
	Chr=chr$Chr
fi
funcFilfromList
BamFil=`readlink -f $InpFil` #resolve absolute path to bam
BamNam=`basename ${BamFil/.bam/}` #a name to use for the various files
BamNm=${BamNam/.list/} 
TmpDir=$BamNam.$Chr.appBQSRjavdir #temp directory for java machine
mkdir -p $TmpDir
if [[ -z $LogFil ]];then
	LogFil=$BamNam.recal.log # a name for the log file
fi
TmpLog=$LogFil.recal.$Chr.log #temporary log file 
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
 -I $RalLst 
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
