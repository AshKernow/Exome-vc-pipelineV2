#!/bin/bash
#$ -cwd -l mem=16G,time=6:: -N MergeGVCF

# This script can be used to merge multiple gVCF files in. It requires the a list of gVCFs to be merged.
# Usage notes:
#    InpFil - (required) - List of gVCF files to be merged
#    RefFil - (required) - shell file containing variables with locations of reference files, jar files, and resource directories; see list below for required variables
#    OutFil - (optional) - Output file name - derived from InpFil otherwise
#    LogFil - (optional) - File for logging progress
#    Help - H - (flag) - get usage information

#list of required tools:
# java <http://www.oracle.com/technetwork/java/javase/overview/index.html>
# GATK <https://www.broadinstitute.org/gatk/> <https://www.broadinstitute.org/gatk/download>

#list of required variables in reference file:
# $REF - reference genome in fasta format - must have been indexed using 'bwa index ref.fa'
# $DBSNP - dbSNP vcf from GATK
# $EXOMPPLN - directory containing exome analysis pipeline scripts
# $GATK - GATK jar file 
# $ETKEY - GATK key file for switching off the phone home feature, only needed if using the B flag

## This file also requires exome.lib.sh - which contains various functions used throughout the Exome analysis scripts; this file should be in the same directory as this script

###############################################################

#set default arguments
usage="
ExmAdHoc.15.MergeGVCFs.sh -i <InputFile> -r <reference_file> -o <OutputFilename> -l <logfile> -H

     -i (required) - Aligned bam file
     -r (required) - shell file containing variables with locations of reference files and resource directories
     -o (optional) - Output file name
     -l (optional) - Log file
     -H (flag) - echo this message and exit
"

#get arguments
while getopts i:r:o:l:H opt; do
    case "$opt" in
        i) InpFil="$OPTARG";; 
        r) RefFil="$OPTARG";;
        o) OutFil="true";;
        l) LogFil="$OPTARG";;
        H) echo "$usage"; exit;;
    esac
done

#Load script library
source $EXOMPPLN/exome.lib.sh #library functions begin "func" #library functions begin "func"

#Call the RefFil to load variables
RefFil=`readlink -f $RefFil`
source $RefFil

#check all required paramaters present
if [[ ! -e "$InpFil" ]]; then 
    echo "Missing/Incorrect required arguments"
    echo "$usage"
    echo "Provided arguments:"
    echo "    InpFil: "$InpFil
    exit
fi

#set local variables
if [[ -z "$OutFil" ]]; then OutFil=`basename $InpFil | sed s/.list$/`; fi # a name for the output file
if [[ -z "$LogFil" ]]; then LogFil=${OutFil%bam}log; fi # a name for the log file
echo "    InpFil: "$InpFil"    OutFil: "$OutFil
TmpLog=$InpFil.mrgGVCF.temp.log #temporary log file

#Start Log File
ProcessName="Merging gVCFs with GATK CombineGVCFs " # Description of the script - used in log
funcWriteStartLog

##Run Joint Variant Calling
StepName="Joint call gVCFs with GATK"
GVCFargument=`sed s/^/ -V /g $InpFil`
StepCmd="java -Xmx20G -Djava.io.tmpdir=$TmpDir -jar $GATKJAR
 -T CombineGVCFs 
 -R $REF
 $GVCFargument
 -o $VcfFil
 -log $GatkLog" #command to be run
#funcRunStep

#End Log
funcWriteEndLog
