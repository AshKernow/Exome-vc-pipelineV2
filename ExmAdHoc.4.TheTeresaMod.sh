#!/bin/bash
#$ -l mem=4G,time=2:: -cwd -S /bin/bash -N TerModBam

# This script is to fix a problem with paired-end Bam files that were aligned with old verions of bwa in which "/3" or "/4" was added to the mate in each read pair rather than /1 and /2, as in:
#D8GSQ5P1:4:1102:10621:68803#0    73    1    10009    0    101M    =    10009    0    ACCCTAACCCTAACCCTAACCCTA....
#D8GSQ5P1:4:1102:10621:68803#0/3    133    1    10009    0    *    =    10009    0    GTTAGGGTTAGGGTTAGGGCTGGG....
# This causes problems with the bam2fq-->bwa pipeline - the bam2fq does not recoginise the pairs and doesn't interleave them properly and hence bwa then sees them as all singletons.
#    InpBam - (required) - The bam file to be modified
#    OutBam - (optional) - A name for the output file. If this is not provided it will be derived from the input name.
#    Help - H - (flag) - get usage information

#list of required tools:
# samtools <http://samtools.sourceforge.net/> <http://sourceforge.net/projects/samtools/files/>

## This file also requires exome.lib.sh - which contains various functions used throughout the Exome analysis scripts; this file should be in the same directory as this script

###############################################################

#set default arguments
usage="
ExmAdHoc.4.TheTeresaMod.sh -i <InputName> -o <OutputName>

     -i (required) - Path to directory containing fastq.gz files
     -o (optional) - Output filename - if not provided \"fixed\" will be added into the original filename
     -H (flag) - echo this message and exit
"

#get arguments
while getopts i:o:H opt; do
    case "$opt" in
        i) InpBam="$OPTARG";;
        o) OutBam="$OPTARG";;
        H) echo "$usage"; exit;;
    esac
done

#check all required paramaters present
if [[ ! -e "$InpBam" ]]; then echo "Missing/Incorrect required arguments"; echo "$usage"; exit; fi

#hpc workaround
if [[ /bin/hostname==*.hcp ]]; then source $HOME/.bash_profile; fi

if [[ -z $OutBam ]]; then
    OutBam=`basename $InpBam | sed 's/bam$/fixed.bam/'`
fi

# bam --> sam | use awk to make modification | sam --> bam
samtools view -h $InpBam | awk 'BEGIN { OFS = "\t" } { gsub(/#0\/3$/,"#0",$1); gsub(/#0\/4$/,"#0",$1); print $0 }' | samtools view -bS - > $OutBam

echo "Completed bam modification"
