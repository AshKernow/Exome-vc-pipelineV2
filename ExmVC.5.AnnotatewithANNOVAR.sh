#!/bin/bash
#$ -cwd -l mem=12G,time=2:: -N AnnVC

#This script takes a bam file or a list of bam files (filename must end ".list") and runs variant calling using the HaplotypeCaller in gVCF mode
#	InpFil - (required) - Path to Bam file to be aligned. Alternatively a file with a list of bams can be provided and the script run as an array job. List file name must end ".list"
#	RefFiles - (required) - shell file to export variables with locations of reference files, jar files, and resource directories; see list below
#	LogFil - (optional) - File for logging progress
#	Flag - P - PipeLine - call the next step in the pipeline at the end of the job
#	Help - H - (flag) - get usage information

#list of required vairables in reference file:
# $EXOMPPLN - directory containing exome analysis pipeline scripts

#list of required tools:
# annovar <

## This file also require exome.lib.sh - which contains various functions used throughout my Exome analysis scripts; this file should be in the same directory as this script

###############################################################

#set default arguments
usage="
ExmVC.1.HaplotypeCaller_GVCFmode.sh -i <InputFile> -r <reference_file> -t <targetfile> -l <logfile> -PABH

	 -i (required) - Path to Bam file for variant calling or \".list\" file containing a multiple paths
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

PipeLine="false"
while getopts i:r:l:t:PABH opt; do
	case "$opt" in
		i) InpFil="$OPTARG";;
		r) RefFil="$OPTARG";; 
		l) LogFil="$OPTARG";;
		P) PipeLine="true";;
		H) echo "$usage"; exit;;
  esac
done

#load settings file
source $RefFil

#Load script library
source $EXOMPPLN/exome.lib.sh #library functions begin "func"

#Set local Variables
##Set local parameters
AnnFil=${InpFil/.vcf/}
TmpLog=$AnnFil.AnnVC.temp.log #temporary log file
TmpDir=$AnnFil.AnnVC; mkdir -p $TmpDir


#Start Log File
ProcessName="Annotate VCF with ANNOVAR" # Description of the script - used in log
funcWriteStartLog

##Convert VCF to ANNOVAR input file using ANNOVAR - use old vcf method
StepNam="Convert VCF to ANNOVAR input file using ANNOVAR"
StepCmd="convert2annovar.pl $InpFil -format vcf4old -outfile $TmpDir/$AnnFil" #command to be run
funcGatkAddArguments
funcRunStep

#End Log
funcWriteEndLog