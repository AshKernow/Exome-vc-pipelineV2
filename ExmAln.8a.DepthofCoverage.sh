#!/bin/bash
#$ -cwd -l mem=8G,time=4:: -N DepOfCov


#This script takes a bam file and generates depth of coverage statistics using GATK
#	InpFil - (required) - Path to Bam file to be aligned or a file containing a list of bam files one per line (file names must end ".list")
#			if it is a list then call the job as an array job with -t 1:n where n is the number of bams
#	RefFiles - (required) - shell file to export variables with locations of reference files, jar files, and resource directories; see list below
#	TgtBed - (required) - Exome capture kit targets bed file (must end .bed for GATK compatability)
#	LogFil - (optional) - File for logging progress
#	Flag - A - AllowMisencoded - see GATK manual, causes GATK to ignore abnormally high quality scores that would otherwise indicate that the quality score encoding was incorrect

#list of required vairables in reference file:
# $TARGET - exome capture intervals bed file or other target file (must end ".bed")
# $EXOMPPLN - directory containing exome analysis pipeline scripts
# $GATK - GATK jar file 

#list of required tools:
# java
# GATK

## This file also require exome.lib.sh - which contains various functions used throughout my Exome analysis scripts; this file should be in the same directory as this script

###############################################################

#set default arguments
AllowMisencoded="false"

#get arguments
while getopts i:r:t:l:A opt; do
	case "$opt" in
		i) InpFil="$OPTARG";;
		r) RefFil="$OPTARG";; 
		t) TgtBed="$OPTARG";; 
		l) LogFil="$OPTARG";;
		A) AllowMisencoded="true";;
	esac
done

#load settings file
source $RefFil

#Load script library
source $EXOMPPLN/exome.lib.sh #library functions begin "func"

#Set Local Variables
funcFilfromList
BamFil=`readlink -f $InpFil` #resolve absolute path to bam
TmpDir=$BamFil.DoC.tempdir #temp directory for java machine
mkdir -p $TmpDir
if [[ -z $LogFil ]];then
	LogFil=$BamFil.DoC.log # a name for the log file
fi
TmpLog=$BamFil.DoC.temp.log #temporary log file 
OutFil=$BamFil.DoC #prefix used in names of output files
GatkLog=$BamFil.gatklog #a log for GATK to output to, this is then trimmed and added to the script log

#Start Log
ProcessName="Depth of Coverage with GATK" # Description of the script - used in log
funcWriteStartLog

#Calculate depth of coverage statistics
StepName="Calculate depth of coverage statistics using GATK DepthOfCoverage" # Description of this step - used in log
StepCmd="java -Xmx5G -Djava.io.tmpdir=$TmpDir -jar $GATKJAR -T DepthOfCoverage -R $REF -I $BamFil -L $TgtBed -o $OutFil -ct 1  -ct 5 -ct 10 -ct 15 -ct 20 -omitIntervals -log $GatkLog" #command to be run
funcGatkAddArguments
funcRunStep

#End Log
funcWriteEndLog
#Clean up
rm -r $TmpDir $TmpLog
