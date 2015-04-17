#!/bin/bash
#$ -cwd -pe smp 6 -l mem=2G,time=12:: -N BamBWABam

# This script takes a bam file and reverts it to sam format and then realigns with BWA mem
# The process is:
#  InputBam --bamshuf--> Randomly Shuffled Bam --bam2fq--> Fastq --bwa mem--> Sam --> Bam
# See http://gatkforums.broadinstitute.org/discussion/2908/howto-revert-a-bam-file-to-fastq-format and associated discussion for rationale, particularly with regard to shuffling the bam, not clear if it is essential, but it probably can't hurt (though it is a little time consuming)
# bam2fq generates and interleaved fastq, if the Bam contains singletons these disrupt the order of the the interleaved read pairs and bwa mem throws an error, therefore they are removed with -s flag and dumped into a separate fastq.gz file
# IMPORTANT READ GROUP CONSIDERATION - if the bam files contains multiple read groups these must be split up and processed separately if you wish to maintain this information as the bam --> bam2fq --> bwa mem --> sam is read group illiterate and will reassign all reads in the bam to the same read group. The script will by default do this for you by calling a new set of jobs, one for each readgroup. Alternatively you can force this with the -R flag.
#   You can use "samtools view -h BAMFILE -r ReadGroupID | samtools -bS > BAMFILE_RGID" to split the file
# Note regarding Quality Scores:
#   The bam2fq call will use the OQ field (old as in unrecalibrated quality scores) if it is present.
#   It is a good idea to check the QS and see what coding they are in. If they are in pre-Illumina 1.8 coding then use the -F flag if calling the pipeline so that GATK adjusts them to the new coding scale.
# Usage notes:
#    InpFil - (required) - Path to Bam file to be aligned. Alternatively a file with a list of bams can be provided and the script run in an array job. List file name must end ".list"
#    RefFil - (required) - shell file containing variables with locations of reference files and resource directories; see list below for required variables
#    LogFil - (optional) - File for logging progress
#    TgtBed - (optional) - Exome capture kit targets bed file (must end .bed for GATK compatability) ; may be specified using a code corresponding to a variable in the RefFil giving the path to the target file- only required if calling pipeline
#    PipeLine - P -(flag) - will start the GATK realign and recalibrate pipeline using the files generated by this script
#    Flag - F - Fix mis-encoded base quality scores - Rewrite the qulaity scores from Illumina 1.5+ to Illumina 1.8+, will subtract 31 from all quality scores; used to fix encoding in some datasets (especially older Illumina ones) which starts at Q64 (see https://en.wikipedia.org/wiki/FASTQ_format#Encoding)
#    Flag - S - Split readgroups and process separately. The script must be fed a single Bam file as input, and called using a array of X jobs, where X is the number or read groups in the bam file.
#    Help - H - (flag) - get usage information

#list of variables required in reference file:
# $REF - reference genome in fasta format - must have been indexed using 'bwa index ref.fa'
# $EXOMPPLN - directory containing exome analysis pipeline scripts, 

#list of required tools:
# samtools <http://samtools.sourceforge.net/> <http://sourceforge.net/projects/samtools/files/>
# bwa mem <http://bio-bwa.sourceforge.net/> <http://sourceforge.net/projects/bio-bwa/files/>
# java <http://www.oracle.com/technetwork/java/javase/overview/index.html>
# picard <http://picard.sourceforge.net/> <http://sourceforge.net/projects/picard/files/>
# seqtk <https://github.com/lh3/seqtk>

## This file also requires exome.lib.sh - which contains various functions used throughout the Exome analysis scripts; this file should be in the same directory as this script

###############################################################

#set default arguments
usage="
ExmAln.1b.ReAlign_Bam_with_BWAmem.sh -i <InputFile> -r <reference_file> -t <target intervals file> -l <logfile> -PRFH

     -i (required) - Path to Bam file or \".list\" file containing a multiple paths
     -r (required) - shell file containing variables with locations of reference files and resource directories
     -l (optional) - Log file
     -t (optional) - Exome capture kit targets or other genomic intervals bed file (must end .bed for GATK compatability); this file is required if calling the pipeline but otherwise can be omitted
     -P (flag) - Initiate exome analysis pipeline after completion of script
     -F (flag) - Fix mis-encoded base quality scores - Rewrite the qulaity scores from Illumina 1.5+ to Illumina 1.8+
     -S (flag) - split by read group prior to processing
     -H (flag) - echo this message and exit
"

PipeLine="false"
FixMisencoded="false"
SplitReads="false"

#get arguments
while getopts i:r:l:t:PFSH opt; do
    case "$opt" in
        i) InpFil="$OPTARG";;
        r) RefFil="$OPTARG";; 
        l) LogFil="$OPTARG";;
        t) TgtBed="$OPTARG";; 
        P) PipeLine="true";;
        F) FixMisencoded="true";;
        S) SplitReads="true";;
        H) echo "$usage"; exit;;
    esac
done

#check all required paramaters present
if [[ ! -e "$InpFil" ]] || [[ ! -e "$RefFil" ]]; then echo "Missing/Incorrect required arguments"; echo "$usage"; exit; fi

#Call the RefFil to load variables
RefFil=`readlink -f $RefFil`
source $RefFil 

#Load script library
source $EXOMPPLN/exome.lib.sh #library functions begin "func"

#set local variables
ArrNum=$SGE_TASK_ID
funcFilfromList #if the input is a list get the appropriate input file for this job of the array --> $InpFil
BamFil=`readlink -f $InpFil` #resolve absolute path to bam
BamNam=`basename $BamFil` 
BamNam=`basename $BamFil | sed s/.bam$//` # a name for the output files
if [[ -z "$LogFil" ]]; then LogFil=$BamNam.BbB.log; fi # a name for the log file
AlnDir=wd.$BamNam.align # directory in which processing will be done
RgIdLst=$BamNam.byRGID.list # a file containing a list of RG ID specific bamfiles
AlnFil=$BamNam.bwamem.bam #filename for bwa-mem aligned file
SngFil=$BamNam.singletons.gz #output file for the singletons to be dumped to
SrtFil=$BamNam.bwamem.sorted.bam #output file for sorted bam
DdpFil=$BamNam.bwamem.mkdup.bam #output file with PCR duplicates marked
FlgStat=$BamNam.bwamem.flagstat #output file for bam flag stats
IdxStat=$BamNam.idxstats #output file for bam index stats
AlnDir=`readlink -f $AlnDir` #resolve absolute path to working directory
mkdir -p $AlnDir # create working directory
cd $AlnDir # move into working directory
TmpLog=$BamNam.BwB$ArrNum.temp.log #temporary log file
TmpDir=$BamNam.BwB$ArrNum.tempdir; mkdir -p $TmpDir #temporary directory
echo $SplitReads

#start log
ProcessName="Align with BWA"
if [[ $SplitReads == "true" ]]; then  ProcessName="Split by readgroups and Align with BWA"; fi
funcWriteStartLog
echo " Build of reference files: "$BUILD >> $TmpLog
echo "----------------------------------------------------------------" >> $TmpLog

#get ReadGroupHeader from input BAM and check readgroup header
RgHeader=$(samtools view -H $BamFil | grep -m 1 ^@RG | awk '{ gsub("\t","\\t") } { print }')
if [[ $RgHeader == "" ]]; then #check that we have a  RG header and if not write a warning to the log file and exit
    echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" >> $TmpLog
    echo "     Problem with ReadGroup header" >> $TmpLog
    echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" >> $TmpLog
    cat $TmpLog >> $LogFil
    rm $TmpLog
    exit 1
fi
#check for multiple readgroups and split if necessary
if [[ $SplitReads == "false" ]]; then
    NumOfRG=$(samtools view -H $BamFil | grep ^@RG | wc -l)
    if [[ $NumOfRG -gt 1 ]]; then 
        echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" >> $TmpLog
        echo "There are "$NumOfRG" Readgroups in the bam:" >> $TmpLog
        samtools view -H $BamFil | grep ^@RG >> $TmpLog
        echo >> $TmpLog
        echo "They will be processed separately and remerged....">> $TmpLog
        cd ../ #move back out of working directory
        mkdir -p stdostde
        QsubCmd="qsub -t 1-$NumOfRG -o stdostde/ -e stdostde/ $EXOMPPLN/ExmAln.1b.ReAlign_Bam_with_BWAmem.sh -i $BamFil -r $RefFil -l $LogFil -S"
    if [[ $TgtBed ]]; then QsubCmd=$QsubCmd" -t $TgtBed"; fi 
    if [[ "$FixMisencoded" == "true" ]]; then QsubCmd=$QsubCmd" -F"; fi
        if [[ "$PipeLine" == "true" ]]; then QsubCmd=$QsubCmd" -P"; fi
        PipeLine="true"
    LogFil=$AlnDir/$LogFil
    TmpLog=$AlnDir/$TmpLog
    TmpDir=$AlnDir/$TmpDir
    funcPipeLine
    #End Log
    funcWriteEndLog
        exit 0
    fi
else
    echo "Splitting Bam by readgroups..." >> $TmpLog
    RGID=`samtools view -H $BamFil | grep ^@RG | head -n $ArrNum | tail -n 1 | sed "s/.*\<ID://" | sed -e 's/\s\+.*//' `
    echo "ReadGroup: " $RGID >> $TmpLog
    RgHeader=$(samtools view -H $BamFil | grep ^@RG | head -n $ArrNum | tail -n 1 | awk '{ gsub("\t","\\t") } { print }')
    AlnFil=$BamNam.RGID_$RGID.bwamem.bam
    SngFil=$BamNam.RGID_$RGID.singletons.gz #output file for the singletons to be dumped to
    SrtFil=$BamNam.RGID_$RGID.bwamem.sorted.bam #output file for sorted bam
    DdpFil=$BamNam.RGID_$RGID.bwamem.mkdup.bam #output file with PCR duplicates marked
    FlgStat=$BamNam.RGID_$RGID.bwamem.flagstat #output file for bam flag stats
    IdxStat=$BamNam.RGID_$RGID.idxstats #output file for bam index stats
    if [[ $ArrNum -eq 1 ]]; then
        echo "Job 1, creating mergelist and sending merge command with hold."  >> $TmpLog
        samtools view -H $BamFil | grep ^@RG | sed "s/.*\<ID:/$BamNam.RGID_/" | sed -e 's/\s\+.*/.bwamem.mkdup.bam/' | awk -v PWD=$PWD '{ print PWD"/"$0 }' > $RgIdLst
        QsubCmd="qsub -hold_jid $JOB_ID -o stdostde/ -e stdostde/ $EXOMPPLN/ExmAdHoc.2.MergeBams.sh -i $RgIdLst -r $RefFil -t $TgtBed -l $LogFil"
        if [[ "$PipeLine" == "true" ]]; then 
            QsubCmd=$QsubCmd" -P"
        else
            QsubCmd=$QsubCmd" -M"
        fi
        PipeLine="true"
        funcPipeLine
        PipeLine="false"
    fi
fi
echo "ReadGroup header: $RgHeader" >> $TmpLog

###Align using BWA mem algorithm
# use HTSlib to shuffle the bam | then tranform it to an interleaved fastq discarding an singletons to a separate file as they mess up the interleaving | transform sam back to bam
#StepCmd="samtools bamshuf -uOn 128 $BamFil tmp | samtools bam2fq -s $SngFil -O - | gzip | bwa mem -M -R \"$RgHeader\" -t 6 -p $REF - | samtools view -bS - > $AlnFil"

StepName="Align with BWA mem"
StepCmd="samtools bamshuf -uOn 128 $BamFil tmp |
 samtools bam2fq -s $SngFil -O - |"
if [[ $FixMisencoded == "true" ]]; then 
    StepCmd=$StepCmd" seqtk seq -Q64 -V - |"
fi
StepCmd=$StepCmd" gzip | bwa mem -M -R \"$RgHeader\" -t 6 -p $REF - |
 samtools view -bS - > $AlnFil"
if [[ $SplitReads == "true" ]]; then 
    StepCmd=${StepCmd#*|}
    mkdir $RGID.split
    SplitFil=$BamNam.$RGID.bam
    StepCmd="samtools view -r $RGID -b $BamFil > $SplitFil; 
    samtools bamshuf -uOn 128 $SplitFil $RGID.split/tmp |"$StepCmd
fi
funcRunStep

#Sort the bam file by coordinate
StepName="Sort Bam using PICARD"
StepCmd="java -Xmx4G -Djava.io.tmpdir=$TmpDir -jar $PICARD/SortSam.jar
 INPUT=$AlnFil
 OUTPUT=$SrtFil
 SORT_ORDER=coordinate
 CREATE_INDEX=TRUE"
funcRunStep
rm $AlnFil #removed the "Aligned bam"

#Mark the duplicates
StepName="Mark PCR Duplicates using PICARD"
StepCmd="java -Xmx4G -Djava.io.tmpdir=$TmpDir -jar $PICARD/MarkDuplicates.jar
 INPUT=$SrtFil
 OUTPUT=$DdpFil
 METRICS_FILE=$DdpFil.dup.metrics.txt
 CREATE_INDEX=TRUE"
funcRunStep
rm $SrtFil ${SrtFil/bam/bai} #removed the "Sorted bam"


if [[ $SplitReads == "false" ]]; then
    #Call next steps of pipeline if requested
    NextJob="Run Genotype VCF"
    QsubCmd="qsub -o stdostde/ -e stdostde/ $EXOMPPLN/ExmAln.2.HaplotypeCaller_GVCFmode.sh -i $DdpFil -r $RefFil -t $TgtBed -l $LogFil -B"
    funcPipeLine
    NextJob="Get basic bam metrics"
    QsubCmd="qsub -o stdostde/ -e stdostde/ $EXOMPPLN/ExmAln.3a.Bam_metrics.sh -i $DdpFil -r $RefFil -l $LogFil"
    funcPipeLine

    #Get flagstat
    StepName="Output flag stats using Samtools"
    StepCmd="samtools flagstat $DdpFil > $FlgStat"
    funcRunStep

    #get index stats
    StepName="Output idx stats using Samtools"
    StepCmd="samtools idxstats $DdpFil > $IdxStat"
    funcRunStep
fi

#End Log
funcWriteEndLog
