#!/bin/bash
#$ -cwd -l mem=4G,time=2:: -N MergeVCF

#This script concatenates multiple vcfs into a single vcf, for example vcfs that have been split by chromosome. 
#    VcfDir - (required) - A driectory containging vcf files to be concatenated - they should all contain the same samples 
#    RefFil - (required) - shell file containing variables with locations of reference files, jar files, and resource directories; see list below for required variables
#    LogFil - (optional) - File for logging progress
#    Flag - P - PipeLine - call the next step in the pipeline at the end of the job
#    Flag - B - BadET - prevent GATK from phoning home
#    Help - H - (flag) - get usage information

#list of required vairables in reference file:
# $EXOMPPLN - directory containing exome analysis pipeline scripts

#list of required tools:
# java <http://www.oracle.com/technetwork/java/javase/overview/index.html>

## This file also requires exome.lib.sh - which contains various functions used throughout the Exome analysis scripts; this file should be in the same directory as this script

###############################################################

#set default arguments
usage="-t 1-NumberofJobs
Z:\Exome_Seq\scripts\Exome_pipeline_scripts_GATKv3\ExmVC.2ug.MergeVCF.sh -i <InputFile> -r <reference_file> -t <targetfile> -l <logfile> -PABH

     -i (required) - Directory containing vcf files to be merged - all vcfs in the directory will be merged
     -r (required) - shell file containing variables with locations of reference files and resource directories
     -l (optional) - Log file
     -P (flag) - Call next step of exome analysis pipeline after completion of script
     -B (flag) - Prevent GATK from phoning home - only if calling pipeline
     -H (flag) - echo this message and exit
"

AllowMisencoded="false"
PipeLine="false"
BadET="false"

PipeLine="false"
while getopts i:r:l:PBH opt; do
    case "$opt" in
        i) VcfDir="$OPTARG";;
        r) RefFil="$OPTARG";; 
        l) LogFil="$OPTARG";;
        P) PipeLine="true";;
        B) BadET="true";;
        H) echo "$usage"; exit;;
  esac
done

#check all required paramaters present
if [[ ! -e "$VcfDir" ]] || [[ ! -e "$RefFil" ]]; then echo "Missing/Incorrect required arguments"; echo "$usage"; exit; fi

#Call the RefFil to load variables
RefFil=`readlink -f $RefFil`
source $RefFil

#Load script library
source $EXOMPPLN/exome.lib.sh #library functions begin "func" #library functions begin "func"

##Set local parameters
VcfDir=${VcfDir%/} # remove trailing slash
VcfNam=${VcfDir%%.*}
VcfFil=$VcfNam.rawvariants.vcf #Outputfile
if [[ -z "$LogFil" ]];then LogFil=$VcfNam.MergeVCF.log; fi # a name for the log file
TmpLog=$VcfNam.MergeVCFtemp.log #temporary log file 

#Start Log File
ProcessName="Merge & Sort individual chromosome VCFs with vcftools" # Description of the script - used in log
funcWriteStartLog

##Merge and sort variant files 
StepName="Merge & sort with vcftools" # Description of this step - used in log
StepCmd="vcf-concat -p $VcfDir/*vcf | vcf-sort -c > $VcfFil"
funcRunStep

#Call next job
NextJob="Annotate with Annovar"
QsubCmd="qsub -o stdostde/ -e stdostde/ $EXOMPPLN/ExmVC.3.AnnotatewithANNOVAR.sh -i $VcfFil -r $RefFil -l $LogFil -P"
funcPipeLine

#End Log
funcWriteEndLog
rm -r $VcfDir
