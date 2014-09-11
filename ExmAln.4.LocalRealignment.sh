#!/bin/bash
#$ -cwd -l mem=10G,time=12:: -N LocRln

# This script takes a bam file and performs local indel realignment using GATK
#    InpFil - (required) - Path to Bam file to be realigned or a file containing a list of bam files one per line (file name must end ".list")
#            if it is a list then call the job as an array job with -t 1:n where n is the number of bams
#    RefFil - (required) - shell file containing variables with locations of reference files, jar files, and resource directories; see list below for required variables
#    LogFil - (optional) - File for logging progress
#    TgtBed - (optional) - Exome capture kit targets bed file (must end .bed for GATK compatability) ; may be specified using a code corresponding to a variable in the RefFil giving the path to the target file
#    Flag - A - AllowMisencoded - see GATK manual (https://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_CommandLineGATK.html#--allow_potentially_misencoded_quality_scores), causes GATK to ignore abnormally high quality scores that would otherwise indicate that the quality score encoding was incorrect
#    Flag - P - PipeLine - call the next step in the pipeline at the end of the job
#    Flag - K - KillFile - this will cause the script to delete the original bam file once the recalibration has successfully completed
#    Flag - B - BadET - prevent GATK from phoning home
#    Flag - F - Fix mis-encoded base quality scores - see GATK manual. GATK will subtract 31 from all quality scores; used to fix encoding in some datasets (especially older Illumina ones) which starts at Q64 (see https://en.wikipedia.org/wiki/FASTQ_format#Encoding)
#    Help - H - (flag) - get usage information

#list of required vairables in reference file:
# $REF - reference genome in fasta format - must have been indexed using 'bwa index ref.fa'
# $INDEL - Gold standard INDEL reference from GATK
# $INDEL1KG - INDEL reference from 1000 genomes
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
 (-t <X>-<Y> [if providing a list]) ExmAln.4.LocalRealignment.sh -i <InputFile> -r <reference_file> -t <targetfile> -l <logfile> -PABH

     -i (required) - Path to Bam file or \".list\" file containing a multiple paths
     -r (required) - shell file containing variables with locations of reference files and resource directories
     -l (optional) - Log file
     -t (optional) - Exome capture kit targets or other genomic intervals bed file (must end .bed for GATK compatability); this file is required if calling the pipeline but otherwise can be omitted
     -P (flag) - Call next step of exome analysis pipeline after completion of script
     -K (flag) - Causes the script to delete the original bam file once the recalibration has successfully completed
     -A (flag) - AllowMisencoded - see GATK manual
     -B (flag) - Prevent GATK from phoning home
     -F (flag) - Fix mis-encoded base quality scores - see GATK manual
     -H (flag) - echo this message and exit
"

AllowMisencoded="false"
PipeLine="false"
BadET="false"
FixMisencoded="false"
KillFile="false"

#get arguments
while getopts i:r:l:t:PKABFH opt; do
    case "$opt" in
        i) InpFil="$OPTARG";;
        r) RefFil="$OPTARG";; 
        l) LogFil="$OPTARG";;
        t) TgtBed="$OPTARG";; 
        P) PipeLine="true";;
        K) KillFile="true";;
        A) AllowMisencoded="true";;
        B) BadET="true";;
        F) FixMisencoded="true";;
        H) echo "$usage"; exit;;
    esac
done

#check all required paramaters present
if [[ ! -e "$InpFil" ]] || [[ ! -e "$RefFil" ]]; then echo "Missing/Incorrect required arguments"; echo "$usage"; exit; fi

#Call the RefFil to load variables
RefFil=`readlink -f $RefFil`
source $RefFil

#Load script library
source $EXOMPPLN/exome.lib.sh #library functions begin "func" #library functions begin "func"

#Set Local Variables
ArrNum=$SGE_TASK_ID
funcGetTargetFile #If the target file has been specified using a code, get the full path from the exported variable
funcFilfromList #if the input is a list get the appropriate input file for this job of the array --> $InpFil
BamFil=`readlink -f $InpFil` #resolve absolute path to bam
BamNam=`basename $BamFil | sed s/.bam//` #a name to use for the various files
if [[ -z "$LogFil" ]];then LogFil=$BamNam.LocReal.log; fi # a name for the log file
RalFil=$BamNam.realigned$ChrNam.bam # the output - a realigned bam file for the CHR
FlgStat=$BamNam.realigned$ChrNam.flagstat # file to output samtools flagstats on the realigned file to
TgtFil=$BamNam.target_intervals$ChrNam.list #target intervals file created by GATK RealignerTargetCreator
GatkLog=$BamNam.LocReal$ChrNam.gatklog #a log for GATK to output to, this is then trimmed and added to the script log
TmpLog=$LogFil.LocReal$ChrNam.temp.log #temporary log file 
TmpDir=$BamNam.LocReal$ChrNam.tempdir; mkdir -p $TmpDir #temporary directory

#Start Log
ProcessName="Local Realignment around InDels" # Description of the script - used in log
funcWriteStartLog

#Generate target file
StepName="Create target interval file using GATK RealignerTargetCreator" # Description of this step - used in log
StepCmd="java -Xmx7G -Djava.io.tmpdir=$TmpDir -jar $GATKJAR
 -T RealignerTargetCreator
 -R $REF
 -I $BamFil
 -known $INDEL
 -known $INDEL1KG
 -L $TgtBed
 -o $TgtFil
 --filter_mismatching_base_and_quals
 -log $GatkLog" #command to be run
funcGatkAddArguments # Adds additional parameters to the GATK command depending on flags (e.g. -B or -F)
funcRunStep

#Realign InDels
StepName="Realign InDels file using GATK IndelRealigner" # Description of this step - used in log
StepCmd="java -Xmx7G -Djava.io.tmpdir=$TmpDir -jar $GATKJAR
 -T IndelRealigner
 -R $REF
 -I $BamFil
 -targetIntervals $TgtFil
 -known $INDEL
 -known $INDEL1KG
 -o $RalFil
 --filter_mismatching_base_and_quals
 -log $GatkLog" #command to be run
funcGatkAddArguments # Adds additional parameters to the GATK command depending on flags (e.g. -B or -F)
funcRunStep

#Get flagstat
StepName="Output flag stats using Samtools"
StepCmd="samtools flagstat $RalFil > $FlgStat"
funcRunStep

#Call next step in pipeline if requested 
NextJob="Run Base Quality Score Recalibration"
QsubCmd="qsub -o stdostde/ -e stdostde/ $EXOMPPLN/ExmAln.5.BaseQualityScoreRecalibration.sh -i $RalFil -r $RefFil -t $TgtBed -l $LogFil -P"
if [[ "$AllowMisencoded" == "true" ]]; then QsubCmd=$QsubCmd" -A"; fi
if [[ "$KillFile" == "true" ]]; then QsubCmd=$QsubCmd" -K"; fi
if [[ "$BadET" == "true" ]]; then QsubCmd=$QsubCmd" -B"; fi
funcPipeLine

#End Log
funcWriteEndLog

#Clean up
rm $TgtFil
if [[ -e $RalFil ]] && [[ "$KillFile" == "true" ]]; then rm $BamFil ${BamFil/bam/bai}; fi
