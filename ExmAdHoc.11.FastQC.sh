#!/bin/bash
#$ -cwd -l mem=8G,time=4:: -N FastqQC

#This script takes a fastq file and generates QC information using fastqc
#    InpFil - (required) - Path to fastq file to be QC'ed or a file containing a list of fastq files one per line (file names must end ".list")
#    Help - H - (flag) - get usage information

#list of required tools:
# java <http://www.oracle.com/technetwork/java/javase/overview/index.html>
# fastqc <http://finwww.bioinformatics.babraham.ac.uk/projects/fastqc/>

## This file also requires exome.lib.sh - which contains various functions used throughout the Exome analysis scripts; this file should be in the same directory as this script

###############################################################

#set default arguments
usage="
 (-t <X>-<Y> [if providing a list]) ExmAdHoc.11.FastQC.sh -i <InputFile> -H

     -i (required) - Path to Bam file or \".list\" file containing a multiple paths
     -H (flag) - echo this message and exit
"



#get arguments
while getopts i:H opt; do
    case "$opt" in
        i) InpFil="$OPTARG";;
        H) echo "$usage"; exit;;
    esac
done


#Load script library
source $EXOMPPLN/exome.lib.sh #library functions begin "func" #library functions begin "func"

#Set Local Variables
ArrNum=$SGE_TASK_ID
funcFilfromList #if the input is a list get the appropriate input file for this job of the array --> $InpFil
FQFil=`readlink -f $InpFil`

#Run fastqc
echo $InpFil
echo $FQFil
fastqc --outdir ./ $FQFil