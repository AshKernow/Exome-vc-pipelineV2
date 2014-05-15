#!/bin/bash
#$ -l mem=4G,time=4:: -N CatFastq -S /bin/bash -cwd 

# This script concatenates fastq.gz files into a single file
# The main input it is directory containing fastq files
# All fastq's in the directory will be combined into 1 (single end) or 2 (paired end) files, so the direcotry should contain fastqs for only 1 sample
# Need to also specifiy "SE" for single end or "PE" for paired end
# If paired end, the files for paired ends should be marked "_R1_" and "_R2_"
# Generally file names should be in the format <SAMPLENAMELANEETCETC>_R<PAIRNUMBER>_<Section_NUMBER>.fastq.gz

usage="
ExmAdHoc.1.ConcatenateFastq.sh -i <InputDirectory> -t <Type>

	 -i (required) - Path to directory containing fastq.gz files
	 -t (required) - \"PE\" for paired-end \"SE\" for single-end
	 -o (optional) - Output filename - if not specified the directory name will be used
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
	echo  "Need specify paired-end or single-end"
	echo $usage
	exit
fi

if [[ -z "$OutNam" ]]; then OutNam=$FqDir; fi #Output file name

if [[ "$Type" == "SE" ]]; then
	echo "Single End"
	#If single end concatenate all fastq
	FqFils=$(find $FqDir | grep "fastq.gz" | uniq | sort)
	echo "----Fastq List---"
	echo "$FqFils"
	zcat $FqFils | gzip > $OutNam".fastq.gz"
elif [[ "$Type" == "PE" ]]; then
	echo "Paired End"
	#if paired end do R1 first then R2
	FqFils=$(find $FqDir | grep "_R1_" | grep "fastq.gz" | uniq | sort)
	echo "----Fastq List R1 ---"
	echo "$FqFils"
	zcat $FqFils | gzip > $OutNam"_R1.fastq.gz"
	FqFils=$(find $FqDir | grep "_R2_" | grep "fastq.gz" | uniq | sort)
	echo "----Fastq List R2 ---"
	echo "$FqFils"
	zcat $FqFils | gzip > $OutNam"_R2.fastq.gz"
fi
