#!/bin/bash
#$ -cwd  -l mem=12G,time=4:: -N HCgVCF

#This script takes a bam file or a list of bam files (filename must end ".list") and runs variant calling using the HaplotypeCaller in gVCF mode
#	InpFil - (required) - Path to Bam file to be aligned. Alternatively a file with a list of bams can be provided and the script run as an array job. List file name must end ".list"
#	RefFiles - (required) - shell file to export variables with locations of reference files, jar files, and resource directories; see list below
#	TgtBed - (optional) - Exome capture kit targets bed file (must end .bed for GATK compatability) - only required if calling pipeline
#	LogFil - (optional) - File for logging progress
#	Flag - A - AllowMisencoded - see GATK manual, causes GATK to ignore abnormally high quality scores that would otherwise indicate that the quality score encoding was incorrect
#	Flag - P - PipeLine - call the next step in the pipeline at the end of the job
#	Flag - B - BadET - prevent GATK from phoning home
#	Help - H - (flag) - get usage information

#list of required vairables in reference file:
# $REF - reference genome in fasta format - must have been indexed using 'bwa index ref.fa'
# $DBSNP - dbSNP vcf from GATK
# $HAPMAP - hapmap vcf from GATKf
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
		t) TgtBed="$OPTARG";; 
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


#Set local Variables
funcGetTargetFile
if [[ "$ArrNum" != "undefined" ]]; then 
	InpNam=`basename ${InpNam/.list/}`
	VcfLst=HCgVCF.$InpNam.list #File listing paths to gVCF
fi
ArrNum=$SGE_TASK_ID
funcFilfromList #if the input is a list get the appropriate input file for this job of the array --> $InpFil
BamFil=`readlink -f $InpFil` #resolve absolute path to bam
BamNam=`basename $BamFil | sed s/.bam//`
BamNam=${BamNam/.bam/} # a name for the output files
if [[ -z $LogFil ]]; then LogFil=$BamNam.HCgVCF.log; fi # a name for the log file
VcfFil=$BamNam.g.vcf #Output File
GatkLog=$BamNam.HCgVCF.gatklog #a log for GATK to output to, this is then trimmed and added to the script log
TmpLog=$BamNam.HCgVCF.temp.log #temporary log file
TmpDir=$BamNam.HCgVCF.tempdir; mkdir -p $TmpDir #temporary directory
infofields="-A AlleleBalance -A BaseQualityRankSumTest -A Coverage -A HaplotypeScore -A HomopolymerRun -A MappingQualityRankSumTest -A MappingQualityZero -A QualByDepth -A RMSMappingQuality -A SpanningDeletions -A FisherStrand -A InbreedingCoeff" #Annotation fields to output into vcf files

#Start Log File
ProcessName="Genomic VCF generatation with GATK HaplotypeCaller" # Description of the script - used in log
funcWriteStartLog

##Run genomic VCF generation
StepNam="gVCF generation with GATK HaplotypeCaller"
StepCmd="java -Xmx7G -Djava.io.tmpdir=$TmpDir -jar $GATKJAR
 -T HaplotypeCaller
 -R $REF
 -L $TgtBed
 -I $BamFil
 --genotyping_mode DISCOVERY
 -stand_emit_conf 10
 -stand_call_conf 30
 --emitRefConfidence GVCF
 --variant_index_type LINEAR
 --variant_index_parameter 128000
 -o $VcfFil
 -D $DBSNP
 --comp:HapMapV3 $HpMpV3 
 -pairHMM VECTOR_LOGLESS_CACHING
 -rf BadCigar
 $infofields
 -log $GatkLog" #command to be run
funcGatkAddArguments
funcRunStep

#Call next step
#NextJob="Get basic bam metrics"
#QsubCmd="qsub -hold_jid $JOB_ID -o stdostde/ -e stdostde/ $EXOMPPLN/ExmAln.3a.Bam_metrics.sh -i $RclLst -r $RefFil -l $LogFil -Q"
#funcPipeLine

#End Log
funcWriteEndLog
