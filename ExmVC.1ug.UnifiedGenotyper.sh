#!/bin/bash
#$ -cwd -l mem=12G,time=12:: -N UG.VC

#This script takes a bam file or a list of bam files (filename must end ".list") and runs variant calling using the UnifiedGenotyper
# The job should be called as an array. The the variants calling will then be split into X jobs, where X is the number of jobs in the array. This should be larger for bigger jobs (more samples).
#    InpFil - (required) - A list of bams for variant calling. List file name must end ".list". 
#    RefFil - (required) - shell file containing variables with locations of reference files, jar files, and resource directories; see list below for required variables
#    TgtBed - (optional) - Exome capture kit targets bed file (must end .bed for GATK compatability) ; may be specified using a code corresponding to a variable in the RefFil giving the path to the target file
#    VcfNam - (optional) - A name for the analysis - to be used for naming output files. Will be derived from input filename if not provided
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
usage="-t 1-<NumberofJobs> ExmVC.1ug.UnifiedGenotyper.sh -i <InputFile> -r <reference_file> -t <targetfile> -l <logfile> -PABH

    <NumberofJobs> - number of jobs to split the variant calling across
     -i (required) - Path to list of Bam files for variant calling
     -r (required) - shell file containing variables with locations of reference files and resource directories
     -t (required) - Exome capture kit targets or other genomic intervals bed file (must end .bed for GATK compatability)
     -n (optional) - Analysis/output VCF name - will be derived from input filename if not provided
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
# The target file needs to be divided evenly between all the jobs. i.e. if the target file is 1000 lines long and there are 40 jobs, each job should have 25 lines of the target file
# bash arithmetic division actually gives the quotient, so if there are 1010 lines and 40 jobs the division would still give 25 lines per a job and the last 10 lines would be lost
# to compensate for this we will find the remainder (RemTar) and then add an extra line to the first $RemTar jobs
ArrNum=$SGE_TASK_ID
NumJobs=$SGE_TASK_LAST
echo $NumJobs
TarLen=$(cat $TgtBed | wc -l) 
RemTar=$(( TarLen % NumJobs )) # get remainder of target file length and number of jobs
QuoTar=$(( TarLen / NumJobs )) # get quotient of target file length and number of jobs
SttLn=1
DivLen=0
echo $RemTar
echo $SttLn
for ((i=1; i <= $ArrNum; i++)); do
    SttLn=$(( SttLn + DivLen ))
    if [[ $i -le $RemTar ]]; then
        DivLen=$(( QuoTar + 1 ))
        else
        DivLen=$QuoTar
    fi
done

###

BamFil=`readlink -f $InpFil` #resolve absolute path to bam
if [[ -z "$VcfNam" ]];then VcfNam=`basename $BamFil | sed s/_recalibrated// | sed s/_bams// | sed s/.bam// | sed s/.list//` ; fi  #analysis name to use for the various files
if [[ -z "$LogFil" ]];then LogFil=$VcfNam.CallVC.log; fi # a name for the log file
VcfDir=$VcfNam.splitfiles; mkdir -p $VcfDir # Directory to output slices to
VcfFil=$VcfDir/$VcfNam.$ArrNum.raw_variants.vcf #Output File
GatkLog=$VcfNam.$ArrNum.gatklog #a log for GATK to output to, this is then trimmed and added to the script log
TmpLog=$VcfNam.CallVC.temp.$ArrNum.log #temporary log file 
TmpDir=$VcfNam.CallVC.$ArrNum.tempdir; mkdir -p $TmpDir #temporary directory
TgtFil=$TmpDir/Range.$VcfNam.$ArrNum.bed #exome capture range
tail -n+$SttLn $TgtBed | head -n $DivLen > $TgtFil #get exome capture range

infofields="-A AlleleBalance -A BaseQualityRankSumTest -A Coverage -A HaplotypeScore -A HomopolymerRun -A MappingQualityRankSumTest -A MappingQualityZero -A QualByDepth -A RMSMappingQuality -A SpanningDeletions" #Annotation fields to output into vcf files

#Start Log File
ProcessName="Variant Calling on Chromosome $CHR with GATK UnifiedGenotyper" # Description of the script - used in log
funcWriteStartLog
echo "Target file line range: $SttLn - $(( $SttLn + $DivLen ))" >> $TmpLog

##Run Joint Variant Calling
StepName="Variant Calling with GATK UnifiedGenotyper" # Description of this step - used in log
StepCmd="java -Xmx7G -Djava.io.tmpdir=$TmpDir -jar $GATKJAR
 -T UnifiedGenotyper
 -R $REF
 -L $TgtFil
 -I $BamFil
 -stand_emit_conf 10
 -stand_call_conf 30
 -o $VcfFil
 -glm BOTH
 -D $DBSNP
 --comp:HapMapV3 $HAPMAP 
 $infofields
 -rf BadCigar
 -log $GatkLog" #command to be run
funcGatkAddArguments # Adds additional parameters to the GATK command depending on flags (e.g. -B or -F)
funcRunStep

#Need to wait for all UnifiedGenotyper jobs to finish and then remerge all the vcfs
if [[ "$ArrNum" -eq 1 ]]; then
    NextJob="Remerge vcf files"
    QsubCmd="qsub -hold_jid $JOB_ID -o stdostde/ -e stdostde/ $EXOMPPLN/ExmVC.2.MergeVCF.sh -i $VcfDir -r $RefFil -l $LogFil -P"
    if [[ "$AllowMisencoded" == "true" ]]; then QsubCmd=$QsubCmd" -A"; fi
    if [[ "$BadET" == "true" ]]; then QsubCmd=$QsubCmd" -B"; fi
    funcPipeLine
fi

#End Log
funcWriteEndLog

#Clean up
