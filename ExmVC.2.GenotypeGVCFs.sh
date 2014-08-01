#!/bin/bash
#$ -cwd -l mem=10G,time=4:: -N GgVCFs

#This script takes a list of gVCF files generated by the HaplotypeCaller (filename must end ".list") and performs the multi-sample joint aggregation step and merges the records together.
#    InpFil - (required) - List of gVCF files. List file name must end ".list"
#    RefFil - (required) - shell file containing variables with locations of reference files, jar files, and resource directories; see list below for required variables
#    TgtBed - (optional) - Exome capture kit targets bed file (must end .bed for GATK compatability) ; may be specified using a code corresponding to a variable in the RefFil giving the path to the target file- only required if calling pipeline
#    VcfNam - (optional) - A name for the analysis - to be used for naming output files. Will be derived from input filename if not provided
#    LogFil - (optional) - File for logging progress
#    Flag - P - PipeLine - call the next step in the pipeline at the end of the job
#    Flag - B - BadET - prevent GATK from phoning home
#    Help - H - (flag) - get usage information

#list of required vairables in reference file:
# $REF - reference genome in fasta format - must have been indexed using 'bwa index ref.fa'
# $DBSNP - dbSNP vcf from GATK
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
ExmVC.2.GenotypeGVCFs.sh -i <InputFile> -r <reference_file> -t <targetfile> -l <logfile> -PABH

     -i (required) - List of gVCF files. List file name must end \".list\"
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
while getopts i:r:t:n:l:PBH opt; do
    case "$opt" in
        i) InpFil="$OPTARG";;
        r) RefFil="$OPTARG";;
        t) TgtBed="$OPTARG";;
        n) VcfNam="$OPTARG";;
        l) LogFil="$OPTARG";;
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
if [[ -z "$VcfNam" ]];then VcfNam=`basename $InpFil`; VcfNam=${VcfNam/.list/}; fi # a name for the output files
if [[ -z $LogFil ]]; then LogFil=$VcfNam.GgVCF.log; fi # a name for the log file
VcfFil=$VcfNam.vcf #Output File
VcfAnnFil=$VcfNam.ann.vcf
GatkLog=$VcfNam.GgVCF.gatklog #a log for GATK to output to, this is then trimmed and added to the script log
TmpLog=$VcfNam.GgVCF.temp.log #temporary log file
TmpDir=$VcfNam.GgVCF.tempdir; mkdir -p $TmpDir #temporary directory
infofields="-A AlleleBalance -A BaseQualityRankSumTest -A Coverage -A HomopolymerRun -A MappingQualityRankSumTest -A MappingQualityZero -A QualByDepth -A RMSMappingQuality -A SpanningDeletions -A FisherStrand -A InbreedingCoeff" #Annotation fields to output into vcf files

#Start Log File
ProcessName="Joint calling of gVCFs with GATK GenotypeGVCFs" # Description of the script - used in log
funcWriteStartLog

##Run Joint Variant Calling
StepNam="Joint call gVCFs with GATK" >> $TmpLog
StepCmd="java -Xmx7G -Djava.io.tmpdir=$TmpDir -jar $GATKJAR
 -T GenotypeGVCFs
 -R $REF
 -L $TgtBed
 -V $InpFil
 -o $VcfFil
 -D $DBSNP
  $infofields
 -log $GatkLog" #command to be run
funcGatkAddArguments # Adds additional parameters to the GATK command depending on flags (e.g. -B or -F)
funcRunStep

##Annotate VCF with GATK
StepNam="Joint call gVCFs" >> $TmpLog
StepCmd="java -Xmx7G -Djava.io.tmpdir=$TmpDir -jar $GATKJAR
 -T VariantAnnotator 
 -R $REF
 -L $VcfFil
 -V $VcfFil
 -o $VcfAnnFil
 -D $DBSNP
  $infofields
 -log $GatkLog" #command to be run
funcGatkAddArguments # Adds additional parameters to the GATK command depending on flags (e.g. -B or -F)
funcRunStep
mv -f $VcfAnnFil $VcfFil

#Call next job
NextJob="Annotate with Annovar"
QsubCmd="qsub -o stdostde/ -e stdostde/ $EXOMPPLN/ExmVC.3.AnnotateVCF.sh -i $VcfFil -r $RefFil -l $LogFil -P"
funcPipeLine

#End Log
funcWriteEndLog
