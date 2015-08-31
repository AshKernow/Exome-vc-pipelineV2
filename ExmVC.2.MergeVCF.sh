#!/bin/bash
#$ -cwd -l mem=4G,time=2:: -N MergeVCF

#This script concatenates multiple vcfs into a single vcf, for example vcfs that have been split by chromosome. 
#    VcfDir - (required) - A driectory containging vcf files to be concatenated - they should all contain the same samples 
#    RefFil - (required) - shell file containing variables with locations of reference files, jar files, and resource directories; see list below for required variables
#    LogFil - (optional) - File for logging progress
#    Flag - C - ChecPrg - Only used by pipeline. The script will check that the number of files in the "Progress directory" is the same as the number in the VcfDir before commencing the merge. If the parameter is not provided the script will proceed regardless.
#    Flag - P - PipeLine - call the next step in the pipeline at the end of the job
#    Flag - B - BadET - prevent GATK from phoning home
#    Help - H - (flag) - get usage information

#list of required vairables in reference file:
# $EXOMPPLN - directory containing exome analysis pipeline scripts

#list of required tools:
# vcftools - http://vcftools.sourceforge.net/index.html

## This file also requires exome.lib.sh - which contains various functions used throughout the Exome analysis scripts; this file should be in the same directory as this script

###############################################################

#set default arguments
usage="-t 1-NumberofJobs
Z:\Exome_Seq\scripts\Exome_pipeline_scripts_GATKv3\ExmVC.2ug.MergeVCF.sh -i <InputFile> -r <reference_file> -t <targetfile> -p <progressdirectory> -l <logfile> -PABH

     -i (required) - Directory containing vcf files to be merged - all vcfs in the directory will be merged
     -r (required) - shell file containing variables with locations of reference files and resource directories
     -l (optional) - Log file
     -C (flag)  - Check the directory of \"progess files\" - only used by pipeline to check that the all genotyping jobs completed successfully
     -P (flag) - Call next step of exome analysis pipeline after completion of script
     -X (flag) - Do not run Variant Quality Score Recalibration - only if calling pipeline
     -B (flag) - Prevent GATK from phoning home - only if calling pipeline
     -H (flag) - echo this message and exit
"

ChecPrg="false"
PipeLine="false"
NoRecal="false"
BadET="false"

while getopts i:r:l:CPXBH opt; do
    case "$opt" in
        i) VcfDir="$OPTARG";;
        r) RefFil="$OPTARG";; 
        l) LogFil="$OPTARG";;
        C) ChecPrg="true";;
        P) PipeLine="true";;
        X) NoRecal="true";;
        B) BadET="true";;
        H) echo "$usage"; exit;;
  esac
done

#check all required paramaters present
if [[ ! -e "$VcfDir" ]] || [[ ! -e "$RefFil" ]]; then echo "Missing/Incorrect required arguments"; echo "$usage"; exit; fi

#Call the RefFil to load variables
RefFil=`readlink -f $RefFil`
source $RefFil

#Load script library
source $EXOMPPLN/exome.lib.sh #library functions begin "func" #library functions begin "func"

##Set local parameters
VcfDir=${VcfDir%/} # remove trailing slash
PrgDir=${VcfDir/splitfiles/progfiles}
VcfNam=${VcfDir%%.*}
VcfFil=$VcfNam.rawvariants.vcf #Outputfile
if [[ -z "$LogFil" ]];then LogFil=$VcfNam.MergeVCF.log; fi # a name for the log file
TmpLog=$VcfNam.MergeVCFtemp.log #temporary log file 

#Start Log File
ProcessName="Merge & Sort individual chromosome VCFs with vcftools" # Description of the script - used in log
funcWriteStartLog

#Check progress directory if provided
CountVCF=$(ls $VcfDir/ | grep vcf$ | wc -l)
echo $ChecPrg
if [[ "$ChecPrg" == "true" ]]; then
    CountPrg=$(ls $PrgDir/ | grep genotypingcomplete$ | wc -l)
    if [[ $CountPrg -ne $CountVCF ]]; then 
        echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" >> $TmpLog
        echo "     Check for completion of genotyping: Failed `date`" >> $TmpLog
        echo "     There are $CountVCF vcfs but only $CountPrg progress files. " >> $TmpLog
        qstat -j $JOB_ID | grep -E "usage *:" >> $TmpLog #get cluster usage stats
        echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" >> $TmpLog
        echo "=================================================================" >> $TmpLog
        cat $TmpLog >> $LogFil
        rm $TmpLog
        exit 1
    else
        echo "     Check for completion of genotyping: Successful `date`" >> $TmpLog
        echo "     There are $CountVCF vcfs and $CountPrg progress files. " >> $TmpLog
    fi
    rm -rf $PrgDir
fi


##Merge and sort variant files 
StepName="Merge & sort with vcftools" # Description of this step - used in log
echo "Merging ... "$CountVCF" ... vcfs" >> $TmpLog
StepCmd="vcf-concat -p $VcfDir/*vcf | vcf-sort -c > $VcfFil"
funcRunStep

#Create a list of sample names used in the vcf
StepName="Output sample list"
StepCmd="for i in \$(grep -m 1 ^#CHROM $VcfFil | cut -f 10-); do echo \$i >> $VcfNam.vcfheaderline.txt; done"
funcRunStep

#gzip and index
StepName="Gzip the vcf and index" # Description of this step - used in log
StepCmd="bgzip $VcfFil; tabix -f -p vcf $VcfFil.gz"
funcRunStep
rm $VcfFil $VcfFil.idx
VcfFil=$VcfFil.gz

#Call next job
NextJob="Annotate with Annovar"
QsubCmd="qsub -o stdostde/ -e stdostde/ $EXOMPPLN/ExmVC.3.AnnotateVCF.sh -i $VcfFil -r $RefFil -l $LogFil -P"
if [[ "$BadET" == "true" ]]; then QsubCmd=$QsubCmd" -B"; fi 
if [[ "$NoRecal" == "true" ]]; then QsubCmd=$QsubCmd" -X"; fi
funcPipeLine

#End Log
funcWriteEndLog
rm -r $VcfDir
