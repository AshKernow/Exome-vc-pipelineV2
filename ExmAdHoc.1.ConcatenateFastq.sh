#!/bin/bash
#$ -l mem=4G,time=4:: -N CatFastq -S /bin/bash -cwd 

# This script concatenates fastq.gz files into a single file
# The main input it is directory containing fastq files
# All fastq's in the directory will be combined into 1 (single end) or 2 (paired end) files, so the direcotry should contain fastqs for only 1 sample
# Need to also specifiy "SE" for single end or "PE" for paired end
# If paired end, the files for paired ends should be marked "_R1_" and "_R2_"
# Generally file names should be in the format <SAMPLENAMELANEETCETC>_R<PAIRNUMBER>_<Section_NUMBER>.fastq.gz
#	FqDir - (required) - A directory containing fastq files
#	Type - (required) - "PE" for paired end, "SE" for single end
#	OutNam - (optional) - A name for the output file. If this is not provided it will be derived from the directory name.
#	Help - H - (flag) - get usage information

#list of required reference files:
# $EXOMPPLN - directory containing exome analysis pipeline scripts

#list of required tools:
# samtools <http://samtools.sourceforge.net/> <http://sourceforge.net/projects/samtools/files/>
# java <http://www.oracle.com/technetwork/java/javase/overview/index.html>
# picard <http://picard.sourceforge.net/> <http://sourceforge.net/projects/picard/files/>

## This file also require exome.lib.sh - which contains various functions used throughout my Exome analysis scripts; this file should be in the same directory as this script

###############################################################

#set default arguments
usage="
ExmAdHoc.1.ConcatenateFastq.sh -i <InputDirectory> -t <Type>

	 -i (required) - Path to directory containing fastq.gz files
	 -t (required) - \"PE\" for paired-end \"SE\" for single-end
	 -o (optional) - Output filename - if not provided the directory name will be used
	 -H (flag) - echo this message and exit
"

#get arguments
while getopts i:t:o:H opt; do
	case "$opt" in
		i) FqDir="$OPTARG";;
		t) Type="$OPTARG";; 
		o) OutNam="$OPTARG";;
		H) echo "$usage"; exit;;
	esac
done

#check directory exists
if [[ ! -d $FqDir ]]; then
	echo  "Need provide a directory"
	echo $usage
	exit
fi

#check for PE/SE specification
if [[ "$Type" != "PE" ]] && [[ "$Type" != "SE" ]]; then
	echo  "Need to specify paired-end or single-end"
	echo $usage
	exit
fi

if [[ -z "$OutNam" ]]; then OutNam=$FqDir; fi #Output file name if not provided

if [[ "$Type" == "SE" ]]; then
	echo "Single End"
	#If single end concatenate all fastq
	FqFils=$(find $FqDir | grep "fastq.gz" | uniq | sort)
	echo "----Fastq List---"
	echo "$FqFils"
	zcat $FqFils | gzip > $OutNam".fastq.gz"
	echo "Done"
elif [[ "$Type" == "PE" ]]; then
	echo "Paired End"
	#if paired end do R1 first then R2
	FqFils=$(find $FqDir | grep "_R1_" | grep "fastq.gz" | uniq | sort)
	echo "----Fastq List R1 ---"
	echo "$FqFils"
	zcat $FqFils | gzip > $OutNam"_R1.fastq.gz"
	echo "Done R1"
	FqFils=$(find $FqDir | grep "_R2_" | grep "fastq.gz" | uniq | sort)
	echo "----Fastq List R2 ---"
	echo "$FqFils"
	zcat $FqFils | gzip > $OutNam"_R2.fastq.gz"
	echo "Done R2"
fi
