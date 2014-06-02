#!/bin/bash
#$ -cwd -l mem=4G,time=2:: -N MergeVCF

#This script concatenates multiple vcfs into a single vcf, for example vcfs that have been split by chromosome. 
#	InpFil - (required) - A driectory containging vcf files to be concatenated - they should all contain the same samples 
#	RefFiles - (required) - shell file to export variables with locations of reference files, jar files, and resource directories; see list below
#	LogFil - (optional) - File for logging progress
#	Flag - A - AllowMisencoded - see GATK manual, causes GATK to ignore abnormally high quality scores that would otherwise indicate that the quality score encoding was incorrect
#	Flag - P - PipeLine - call the next step in the pipeline at the end of the job
#	Flag - B - BadET - prevent GATK from phoning home
#	Help - H - (flag) - get usage information

#list of required vairables in reference file:
# $EXOMPPLN - directory containing exome analysis pipeline scripts

#list of required tools:
# java <http://www.oracle.com/technetwork/java/javase/overview/index.html>

## This file also require exome.lib.sh - which contains various functions used throughout my Exome analysis scripts; this file should be in the same directory as this script

###############################################################

#set default arguments
usage="-t 1-NumberofJobs
Z:\Exome_Seq\scripts\Exome_pipeline_scripts_GATKv3\ExmVC.2ug.MergeVCF.sh -i <InputFile> -r <reference_file> -t <targetfile> -l <logfile> -PABH

	 -i (required) - Directory containing vcf files to be merged - all vcfs in the directory will be merged
	 -r (required) - shell file to export variables with locations of reference files and resource directories
	 -l (optional) - Log file
	 -P (flag) - Call next step of exome analysis pipeline after completion of script
	 -A (flag) - AllowMisencoded - see GATK manual - only if calling pipeline
	 -B (flag) - Prevent GATK from phoning home - only if calling pipeline
	 -H (flag) - echo this message and exit
"

AllowMisencoded="false"
PipeLine="false"
BadET="false"

PipeLine="false"
while getopts i:r:l:PABH opt; do
	case "$opt" in
		i) VcfDir="$OPTARG";;
		r) RefFil="$OPTARG";; 
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

#Get VCF stats with python script
StepNam="Get VCF stats"
StepCmd="python $EXOMPPLN/VCF_summary_Stats.py -v $VcfFil -o ${VcfFil/vcf/stats.tsv}"
funcRunStep

#Call next steps
NextJob="Recalibrate Variant Quality"
QsubCmd="qsub -o stdostde/ -e stdostde/ $EXOMPPLN/ExmVC.3.RecalibrateVariantQuality.sh -i $VcfFil -r $RefFil -l $LogFil -P"
if [[ "$AllowMisencoded" == "true" ]]; then QsubCmd=$QsubCmd" -A"; fi
if [[ "$BadET" == "true" ]]; then QsubCmd=$QsubCmd" -B"; fi 
funcPipeLine

#End Log
funcWriteEndLog
rm -r $VcfDir
