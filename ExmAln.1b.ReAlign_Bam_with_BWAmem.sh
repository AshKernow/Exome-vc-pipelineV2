#!/bin/bash
#$ -cwd -pe smp 6 -l mem=2G,time=1:: -N BamBWABam
# -cwd -l mem=1G,time=:10: -N BamBWABam
#This script takes a bam file and reverts it to sam format and then realigns with BWA mem
#	BamFil - (required) - Path to Bam file to be aligned
#	RefFiles - (required) - shell file to export variables with locations of reference files and resource directories; see list below
#	BamNam - (optional) - A base name for output files and directories
#	LogFil - (optional) - File for logging progress
#	Chain - (flag) - will start the GATK realign and recalibrate pipeline using the files generated by this script

#list of required reference files:
# $REF - reference genome in fasta format - must have been indexed using 'bwa index ref.fa'
# $EXOMPPLN - directory containing exome analysis pipeline scripts, 
#list of required tools:
# samtools
# bwa mem
# picard

## This file also require exome.lib.sh - which contains various functions used throughout my Exome analysis scripts; this file should be in the same directory as this script

###############################################################

#get arguments
Chain="false"
BamNam="none"
LogFil="none"
while getopts i:r:n:l:C opt; do
	case "$opt" in
		i) BamFil="$OPTARG";;
		r) RefFil="$OPTARG";; 
		n) BamNam="$OPTARG";; 
		l) LogFil="$OPTARG";;
		C) Chain="true";;
	esac
done

#load RefFil file
source $RefFil #load the references required
#Load script library
source $EXOMPPLN/exome.lib.sh

#set local variables
BamFil=`readlink -f $BamFil` #resolve absolute path to bam
if [[ $BamNam == "none" ]];then
	BamNam=`basename $BamFil` 
	BamNam=${BamNam/.bam/} # a name for the bam
fi
if [[ $LogFil == "none" ]];then
	LogFil=$BamNam.BbB.log # a name for the log file
fi
AlnDir=$BamNam.align # directory in which processing will be done
AlnFil=$BamNam.bwamem.bam #filename for output file
DdpFil=$BamNam.bwamem.mkdup.bam #filename for file with PCR duplicates marked
mkdir -p $AlnDir # create working directory
cd $AlnDir # move into working directory
TmpLog=$LogFil.BbB.temp.log #use TmpLog for LogFil - for consistency with other exome analysis scripts

#start log
ProcessName="Align with BWA"
WriteStartLog
echo " Bam file: $BamFil" >> $TmpLog
echo " Base name for outputs: $BamNam" >> $TmpLog
echo " Build of reference files: "$BUILD >> $TmpLog
echo "----------------------------------------------------------------" >> $TmpLog

#get ReadGroupHeader from input BAM
RgHeader=$(samtools view -H $BamFil | grep ^@RG | awk '{ gsub("\t","\\t") } { print }')
echo "ReadGroup header: $RgHeader" >> $TmpLog
if [[ $RgHeader == "" ]]||[[ $(echo "$RgHeader" | wc -l) -gt 1 ]]; then #check that we have a  RG header and if not write a warning to the log file
	echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" >> $TmpLog
	echo "     Problem with ReadGroup header" >> $TmpLog
	echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" >> $TmpLog
fi

###Align using BWA mem algorithm
StepName="Align with BWA mem"
StepCmd="samtools bamshuf -Ou $BamFil | samtools bam2fq - > $AlnFil.fq"
#| bwa mem -M -R \"$RgHeader\" -t 6 -p $REF - > $AlnFil.sam"
#| samtools view -bS - > $AlnFil"
funcRunStep

#Final CleanUp
cat $TmpLog >> $LogFil
rm $TmpLog
