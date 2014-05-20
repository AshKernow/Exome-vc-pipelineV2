#!/bin/bash
#$ -cwd -l mem=10G,time=6:: -N GenBQSR

#This script takes a bam file or a list of bam files  and generates the base quality score recalibration table using GATK
#	InpFil - (required) - Path to Bam file or a list of BamFiles to be recalibrated, if a list filename must end ".list"
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
# $DBSNP - dbSNP vcf from GATK
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
ExmAln.5.GenerateBQSRTable.sh -i <InputFile> -r <reference_file> -t <targetfile> -l <logfile> -PABH

	 -i (required) - Path to Bam file or \".list\" file containing a multiple paths
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
RefFil=`readlink -f $RefFil`
source $RefFil

#Load script library
source $EXOMPPLN/exome.lib.sh #library functions begin "func"

#Set Local Variables
ArrNum=$SGE_TASK_ID
funcGetTargetFile
BamFil=`readlink -f $InpFil` #resolve absolute path to bam
BamNam=`basename $BamFil | sed s/.bam// | sed s/.list//` #a name to use for the various files
if [[ -z "$LogFil" ]];then LogFil=$BamNam.GenBQSR.log; fi # a name for the log file
RclTable=$BamNam.recal.table # output - base quality score recalibration table
GatkLog=$BamNam.gatklog #a log for GATK to output to, this is then trimmed and added to the script log
TmpLog=$BamNam.GenBQSR.temp.log #temporary log file 
TmpDir=$BamNam.GenBQSR.tempdir; mkdir -p $TmpDir #temporary directory

#Start Log
ProcessName="Generate Base Quality Score Recalibration Table with GATK" # Description of the script - used in log
funcWriteStartLog

#Generate target file
StepName="Create recalibration data file using GATK BaseRecalibrator" # Description of this step - used in log
StepCmd="java -Xmx7G -Djava.io.tmpdir=$TmpDir -jar $GATKJAR
 -T BaseRecalibrator 
 -R $REF 
 -I $BamFil 
 -L $TgtBed 
 -ip 50
 -knownSites $DBSNP 
 -knownSites $INDEL 
 -knownSites $INDEL1KG 
 -o $RclTable 
 --filter_mismatching_base_and_quals
 -log $GatkLog" #command to be run
funcGatkAddArguments
funcRunStep


#Call next step; if original file was a list of bams then need to call as an array job
ChecList=${InpFil##*.}
if [[ "$ChecList" == "list" ]];then
   echo $ChecList
   nJobs=$(cat $BamFil | wc -l)
   NextJob="Apply Base Quality Score Recalibration"
   QsubCmd="qsub -t 1-$nJobs -o stdostde/ -e stdostde/ $EXOMPPLN/ExmAln.6.ApplyRecalibration.sh -i $BamFil -x $RclTable -r $RefFil -t $TgtBed -l $LogFil -P -K"
   if [[ "$AllowMisencoded" == "true" ]]; then QsubCmd=$QsubCmd" -A"; fi
   if [[ "$BadET" == "true" ]]; then QsubCmd=$QsubCmd" -B"; fi
   funcPipeLine
else
    NextJob="Apply Base Quality Score Recalibration"
    QsubCmd="qsub -o stdostde/ -e stdostde/ $EXOMPPLN/ExmAln.6.ApplyRecalibration.sh -i $BamFil -x $RclTable -r $RefFil -t $TgtBed -l $LogFil -P -K"
    if [[ "$AllowMisencoded" == "true" ]]; then QsubCmd=$QsubCmd" -A"; fi
    if [[ "$BadET" == "true" ]]; then QsubCmd=$QsubCmd" -B"; fi
    funcPipeLine
fi

#End Log
funcWriteEndLog
