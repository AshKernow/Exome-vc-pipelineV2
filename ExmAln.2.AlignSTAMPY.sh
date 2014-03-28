#!/bin/bash
#$ -cwd -l mem=8G,time=1:: -N Stampy

#This script takes a bam file and reverts it to sam format and then realigns with BWA mem
#	InpFil - (required) - Path to Bam file to be aligned. Alternatively a file with a list of bams can be provided and the script run in an array job. List file name must end ".list"
#	RefFiles - (required) - shell file to export variables with locations of reference files and resource directories; see list below
#	LogFil - (optional) - File for logging progress
#	PipeLine - P -(flag) - will start the GATK realign and recalibrate pipeline using the files generated by this script
#	Help - H - (flag) - get usage information

#list of required reference files:
# $STHSH - hash file for Stampy
# $STIDX - genome index file for Stampy
# above to files are created by running Stampy on the reference genome, which should be in the same directory. See Stampy documenation for further details

#list of required tools:
# picard <http://picard.sourceforge.net/> <http://sourceforge.net/projects/picard/files/>
# stampy <http://www.well.ox.ac.uk/project-stampy>
# HTSlib <https://github.com/samtools/htslib> <https://github.com/samtools/htslib/archive/master.zip>

## This file also require exome.lib.sh - which contains various functions used throughout my Exome analysis scripts; this file should be in the same directory as this script

###############################################################

usage="
ExmAln.2.AlignSTAMPY.sh -i <InputFile> -r <reference_file> -l <logfile> -PH

	 -i (required) - Path to Bam file to be aligned or \".list\" file containing a multiple paths
	 -r (required) - shell file to export variables with locations of reference files and resource directories
	 -l (optional) - Log file
	 -P (flag) - Initiate next step of exome analysis pipeline after completion of script
	 -H (flag) - echo this message and exit
"

#get arguments
PipeLine="false"

while getopts i:r:l:PH opt; do
  case "$opt" in
      i) InpFil="$OPTARG";;
      r) RefFil="$OPTARG";;
      l) LogFil="$OPTARG";;
	  P) PipeLine="$OPTARG";;
	  H) echo "$usage"; exit;;
  esac
done

#load RefFil file
source $RefFil #load the references required
#Load script library
source $EXOMPPLN/exome.lib.sh

#check all required paramaters present
if [[ ! -e $InpFil ]] || [[ ! -e $RefFil ]]; then echo "Missing required arguments"; echo "$usage"; exit; fi

#Set local variablesfuncFilfromList
funcFilfromList
BamFil=`readlink -f $InpFil` #resolve absolute path to bam
BamNam=`basename $BamFil` 
BamNam=${BamNam/.bam/} # a name for the output files
BamNam=${BamNam/.list/} 
if [[ -z $LogFil ]]; then LogFil=$BamNam.Stampy.log; fi # a name for the log file
TmpLog=$LogFil.temp.log #temporary log file
AlnFil=$BamNam.Stampy.bam #filename for stampy aligned file
SrtFil=$BamNam.Stampy.sorted.bam #output file for sorted bam

#Start Log
ProcessName="Align BAM with Stampy"
funcWriteStartLog

#Run Jobs
#Align with Stampy
StepName="Align BAM with Stampy"
StepCmd="stampy.py -g $STIDX -h $STHSH
 --bamkeepgoodreads 
 --bwamark 
 -t 4 
 -M $BamFil 
 -v 3 | 
 htscmd samview -bS - > $AlnFil"
funcRunStep

#Sort the bam file by coordinate
StepName="Sort Bam using PICARD"
StepCmd="java -Xmx4G -Djava.io.tmpdir=$TmpDir -jar $PICARD/SortSam.jar
 INPUT=$AlnFil
 OUTPUT=$SrtFil
 SORT_ORDER=coordinate
 CREATE_INDEX=TRUE
 VERBOSITY=DEBUG"
funcRunStep
#rm $AlnFil

#Index with HTSlib
#StepName="Index BAM with HTSlib"
#StepCmd="htscmd bamidx $SrtFil"
#funcRunStep

#Final CleanUp
funcWriteEndLog
rm $TmpLog
