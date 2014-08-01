#!/bin/bash
#$ -cwd -l mem=12G,time=9:: -N AppRcl

#This script takes a bam file and uses a previously generated base quality score recalibration (BQSR) table to recalibrate them using GATK.
#    InpFil - (required) - Path to Bam file
#    RclTab - (required) - Previously generated BQSR table
#    RefFil - (required) - shell file containing variables with locations of reference files, jar files, and resource directories; see list below for required variables
#    TgtBed - (required) - Exome capture kit targets bed file (must end .bed for GATK compatability); may be specified using a code corresponding to a variable in the RefFil giving the path to the target file
#    LogFil - (optional) - File for logging progress
#    Flag - P - PipeLine - call the next step in the pipeline at the end of the job
#    Flag - K - KillFile - this will cause the script to delete the original bam file once the recalibration has successfully completed
#    Flag - A - AllowMisencoded - see GATK manual, causes GATK to ignore abnormally high quality scores that would otherwise indicate that the quality score encoding was incorrect
#    Flag - B - BadET - prevent GATK from phoning home
#    Help - H - (flag) - get usage information

#list of required vairables in reference file:
# $REF - reference genome in fasta format - must have been indexed using 'bwa index ref.fa'
# $EXOMPPLN - directory containing exome analysis pipeline scripts
# $GATK - GATK jar file 
# $ETKEY - GATK key file for switching off the phone home feature, only needed if using the B flag

#list of required tools:
# java <http://www.oracle.com/technetwork/java/javase/overview/index.html>
# GATK <https://www.broadinstitute.org/gatk/> <https://www.broadinstitute.org/gatk/download>

## This file also requires exome.lib.sh - which contains various functions used throughout the Exome analysis scripts; this file should be in the same directory as this script

###############################################################

#set default arguments
usage="
ExmAln.6.ApplyRecalibration.sh -i <InputFile> -x <GATK BQSR table> -r <reference_file> -t <targetfile> -l <logfile> -PABH

     -i (required) - Path to Bam file
     -x (required) - Previously generated BQSR table
     -r (required) - shell file containing variables with locations of reference files and resource directories
     -t (required) - Exome capture kit targets bed file (must end .bed for GATK compatability)
     -l (optional) - Log file
     -P (flag) - Call next step of exome analysis pipeline after completion of script
     -K (flag) - Causes the script to delete the original bam file once the recalibration has successfully completed
     -A (flag) - AllowMisencoded - see GATK manual
     -B (flag) - Prevent GATK from phoning home
     -H (flag) - echo this message and exit
"

AllowMisencoded="false"
PipeLine="false"
BadET="false"
KillFile="false"

#get arguments
while getopts i:x:r:t:l:PKABH opt; do
    case "$opt" in
        i) InpFil="$OPTARG";;
        x) RclTab="$OPTARG";;
        r) RefFil="$OPTARG";; 
        t) TgtBed="$OPTARG";; 
        l) LogFil="$OPTARG";;
        P) PipeLine="true";;
        K) KillFile="true";;
        A) AllowMisencoded="true";;
        B) BadET="true";;
        H) echo "$usage"; exit;;
    esac
done

#check all required paramaters present
if [[ ! -e "$InpFil" ]] || [[ ! -e "$RclTab" ]] || [[ ! -e "$RefFil" ]] || [[ -z "$TgtBed" ]]; then echo "Missing/Incorrect required arguments"; echo "$usage"; exit; fi

#Call the RefFil to load variables
RefFil=`readlink -f $RefFil`
source $RefFil

#Load script library
source $EXOMPPLN/exome.lib.sh #library functions begin "func" #library functions begin "func"

#Set Local Variables
funcGetTargetFile #If the target file has been specified using a code, get the full path from the exported variable
InpNam=`basename $InpFil | sed s/.bam//`
BamFil=`readlink -f $InpFil` #resolve absolute path to bam
BamNam=`basename $BamFil | sed s/.bam//` #a name to use for the various files
if [[ -z "$LogFil" ]];then LogFil=$BamNam.appBQSR.log; fi # a name for the log file
RclFil=$BamNam.recalibrated.bam #file to output recalibrated bam to
FlgStat=${RclFil/bam/flagstat} # file to output samtools flagstats on the final file to
GatkLog=$BamNam.AppBQSR.gatklog #a log for GATK to output to, this is then trimmed and added to the script log
TmpLog=$BamNam.AppBQSR.temp.log #temporary log file
TmpDir=$BamNam.AppBQSR.tempdir; mkdir -p $TmpDir #temporary directory
TmpTar=TmpTarFil.$Chr.bed #temporary target file

#Start Log
ProcessName="Recalibrate Base Quality Scores with GATK" # Description of the script - used in log
funcWriteStartLog

#Apply Recalibration
StepName="Apply recalibration using GATK PrintReads" # Description of this step - used in log
StepCmd="java -Xmx7G -Djava.io.tmpdir=$TmpDir -jar $GATKJAR
 -T PrintReads 
 -R $REF 
 -I $BamFil 
 -BQSR $RclTab
 -L $TgtBed
 -ip 50
 -o $RclFil 
 --emit_original_quals 
 --filter_mismatching_base_and_quals
 -log $GatkLog" #command to be run
funcGatkAddArguments # Adds additional parameters to the GATK command depending on flags (e.g. -B or -F)
funcRunStep

#Get flagstat
StepName="Output flag stats using Samtools"
StepCmd="samtools flagstat $RclFil > $FlgStat"
funcRunStep

#Call next steps
NextJob="Get Depth of Coverage Statistics"
QsubCmd="qsub -o stdostde/ -e stdostde/ $EXOMPPLN/ExmAln.8a.DepthofCoverage.sh -i $RclFil -r $RefFil -t $TgtBed -l $LogFil"
if [[ "$AllowMisencoded" == "true" ]]; then QsubCmd=$QsubCmd" -A"; fi
if [[ "$BadET" == "true" ]]; then QsubCmd=$QsubCmd" -B"; fi
funcPipeLine

NextJob="Get basic bam metrics"
QsubCmd="qsub -o stdostde/ -e stdostde/ $EXOMPPLN/ExmAln.3a.Bam_metrics.sh -i $RclFil -r $RefFil -l $LogFil -Q"
funcPipeLine

#End Log
funcWriteEndLog

#Clean up
if [[ -e $RclFil ]] && [[ "$KillFile" == "true" ]]; then rm $BamFil ${BamFil/bam/bai}; fi
