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
	 -t (required) - Exome capture kit targets or other genomic intervals bed file (must end .bed for GATK compatability)
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
BamFil=`readlink -f $InpFil` #resolve absolute path to bam
BamNam=`basename ${BamFil/.bam/}` #a name to use for the various files
BamNam=${BamNam/.list/} 
TmpDir=$BamNam.GenBQSRjavdir #temp directory for java machine
mkdir -p $TmpDir
if [[ -z $LogFil ]];then LogFil=$BamNam.GenBQSR.log; fi # a name for the log file
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
 -L $TgtBed 
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

#Call next step
NextJob="Apply Base Quality Score Recalibration"
QsubCmd="qsub $EXOMPPLN/ExmAln.6.ApplyRecalibration.sh -i $BamFil -x $RclTable -r $RefFil -t $TgtBed -l $LogFil -P"
if [[ $AllowMisencoded == "true" ]]; then QsubCmd=$QsubCmd" -A"; fi
if [[ $BadET == "true" ]]; then QsubCmd=$QsubCmd" -B"; fi
funcPipeLine

#End Log
funcWriteEndLog
#Clean up
rm -r $TmpDir $TmpLog
