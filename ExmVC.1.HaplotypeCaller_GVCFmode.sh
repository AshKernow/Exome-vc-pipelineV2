#!/bin/bash
#$ -cwd  -l mem=12G,time=8:: -N HCgVCF

#This script takes a bam file or a list of bam files (filename must end ".list") and runs variant calling using the HaplotypeCaller in gVCF mode
#    InpFil - (required) - Path to Bam file to be aligned. Alternatively a file with a list of bams can be provided and the script run as an array job. List file name must end ".list"
#    RefFil - (required) - shell file containing variables with locations of reference files, jar files, and resource directories; see list below for required variables
#    TgtBed - (optional) - Exome capture kit targets bed file (must end .bed for GATK compatability) ; may be specified using a code corresponding to a variable in the RefFil giving the path to the target file- only required if calling pipeline
#    VcfNam - (optional) - A name for the analysis - to be used for naming output files. Will be derived from input filename if not provided; only used if calling pipeline
#    LogFil - (optional) - File for logging progress
#    Flag - P - PipeLine - call the next step in the pipeline at the end of the job
#    Flag - B - BadET - prevent GATK from phoning home
#    Help - H - (flag) - get usage information

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

## This file also requires exome.lib.sh - which contains various functions used throughout the Exome analysis scripts; this file should be in the same directory as this script

###############################################################

#set default arguments
usage="
ExmVC.1.HaplotypeCaller_GVCFmode.sh -i <InputFile> -r <reference_file> -t <targetfile> -l <logfile> -PABH

     -i (required) - Path to Bam file for variant calling or \".list\" file containing a multiple paths
     -r (required) - shell file containing variables with locations of reference files and resource directories
     -t (required) - Exome capture kit targets or other genomic intervals bed file (must end .bed for GATK compatability)
     -n (optional) - Analysis/output VCF name - will be derived from input filename if not provided; only used if calling pipeline
     -l (optional) - Log file
     -P (flag) - Call next step of exome analysis pipeline after completion of script
     -B (flag) - Prevent GATK from phoning home
     -H (flag) - echo this message and exit
"

AllowMisencoded="false"
PipeLine="false"
BadET="false"

PipeLine="false"
while getopts i:r:n:l:t:PBH opt; do
    case "$opt" in
        i) InpFil="$OPTARG";;
        r) RefFil="$OPTARG";; 
        n) VcfNam="$OPTARG";;
        l) LogFil="$OPTARG";;
        t) TgtBed="$OPTARG";; 
        P) PipeLine="true";;
        B) BadET="true";;
        H) echo "$usage"; exit;;
  esac
done

#check all required paramaters present
if [[ ! -e "$InpFil" ]] || [[ ! -e "$RefFil" ]] || [[ -z "$TgtBed" ]]; then echo "Missing/Incorrect required arguments"; echo "$usage"; exit; fi

#Call the RefFil to load variables
RefFil=`readlink -f $RefFil`
source $RefFil

#Load script library
source $EXOMPPLN/exome.lib.sh #library functions begin "func" #library functions begin "func"


#Set local Variables
funcGetTargetFile #If the target file has been specified using a code, get the full path from the exported variable
if [[ "$ArrNum" != "undefined" ]]; then 
    if [[ -z "$VcfNam" ]];then VcfNam=`basename $InpFil`; VcfNam=${VcfNam%%.*}; fi  # a name for the output files
    VcfLst=HCgVCF.$VcfNam.list #File listing paths to gVCF
    if [[ -z $LogFil ]]; then LogFil=$VcfNam.HCgVCF.log; fi # a name for the log file
    if [[ "$PipeLine" == "true" ]]; then awk '{ gsub(/.*\//, ""); gsub(/.bam$/, ".g.vcf"); print }' $InpFil > $VcfLst; fi #make a list of the output gVCF files for passing to next step of pipeline
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
funcGatkAddArguments # Adds additional parameters to the GATK command depending on flags (e.g. -B or -F)
funcRunStep

#Call next step
if [[ "$ArrNum" -eq 1 ]]; then
    mkdir -p stdostde
    NextJob="Genotype gVCFs"
    QsubCmd="qsub -hold_jid $JOB_ID -o stdostde/ -e stdostde/ $EXOMPPLN/ExmVC.2.GenotypeGVCFs.sh -i $VcfLst -r $RefFil -t $TgtBed -l $LogFil -P"
    if [[ "$BadET" == "true" ]]; then QsubCmd=$QsubCmd" -B"; fi
    funcPipeLine
fi
#End Log
funcWriteEndLog
