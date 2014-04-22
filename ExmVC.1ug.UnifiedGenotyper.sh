#!/bin/bash
#$ -cwd  -l mem=12G,time=4:: -N HCgVCF

#This script takes a bam file or a list of bam files (filename must end ".list") and runs variant calling using the HaplotypeCaller in gVCF mode
#	InpFil - (required) - Path to Bam file to be aligned. Alternatively a file with a list of bams can be provided and the script run as an array job. List file name must end ".list"
#	RefFiles - (required) - shell file to export variables with locations of reference files, jar files, and resource directories; see list below
#	TgtBed - (optional) - Exome capture kit targets bed file (must end .bed for GATK compatability) - only required if calling pipeline
#	LogFil - (optional) - File for logging progress
#	Flag - A - AllowMisencoded - see GATK manual, causes GATK to ignore abnormally high quality scores that would otherwise indicate that the quality score encoding was incorrect
#	Flag - P - PipeLine - call the next step in the pipeline at the end of the job
#	Flag - B - BadET - prevent GATK from phoning home
#	Help - H - (flag) - get usage information

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

## This file also require exome.lib.sh - which contains various functions used throughout my Exome analysis scripts; this file should be in the same directory as this script

###############################################################

#set default arguments
usage="
ExmVC.1.HaplotypeCaller_GVCFmode.sh -i <InputFile> -r <reference_file> -t <targetfile> -l <logfile> -PABH

	 -i (required) - Path to Bam file for variant calling or \".list\" file containing a multiple paths
	 -r (required) - shell file to export variables with locations of reference files and resource directories
	 -t (required) - Exome capture kit targets or other genomic intervals bed file (must end .bed for GATK compatability)
	 -l (optional) - Log file
	 -P (flag) - Call next step of exome analysis pipeline after completion of script
	 -A (flag) - AllowMisencoded - see GATK manual
	 -B (flag) - Prevent GATK from phoning home
	 -H (flag) - echo this message and exit
"

AllowMisencoded="false"
PipeLine="false"
BadET="false"

PipeLine="false"
while getopts i:r:l:t:PABH opt; do
	case "$opt" in
		i) InpFil="$OPTARG";;
		r) RefFil="$OPTARG";; 
		l) LogFil="$OPTARG";;
		t) TgtBed="$OPTARG";; 
		P) PipeLine="true";;
		A) AllowMisencoded="true";;
		B) BadET="true";;
		H) echo "$usage"; exit;;
  esac
done

#load settings file
source $RefFil

#Load script library
source $EXOMPPLN/exome.lib.sh #library functions begin "func"

#Set local Variables
JobNum=$SGE_TASK_ID
JobNm=${JOB_NAME#*.}
TmpLog=$LogFil.CallVC.$JobNum.log
TmpDir=$JobNm.$JobNum.VC 
mkdir -p $TmpDir
# The target file needs to be divided evenly between all the jobs. i.e. if the target file is 1000 lines long and there are 40 jobs, each job should have 25 lines of the target file
# bash arithmetic division actually gives the quotient, so if there are 1010 lines and 40 jobs the division would still give 25 lines per a job and the last 10 lines would be lost
# to compensate for this we will find the remainder (RemTar) and then add an extra line to the first $RemTar jobs
TarLen=$(cat $TARGET | wc -l) 
RemTar=$(( TarLen % NumJobs )) # get remainder of target file length and number of jobs
QuoTar=$(( TarLen / NumJobs )) # get quotient of target file length and number of jobs
FinLn=0
for ((i=1; i <= $JobNum; i++)); do
	SttLn=$(( FinLn + 1 ))
	if [[ $i -le $RemTar ]]; then
		DivLen=$(( QuoTar + 1 ))
		else
		DivLen=$QuoTar
	fi
	FinLn=$(( FinLn + DivLen ))
done
###
Range=$TmpDir/Range$JobNm$JobNum.bed #exome capture range
tail -n+$SttLn $TARGET | head -n $DivLen > $Range #get exome capture range
VcfDir=$JobNm"_VCF_final" #Output Directory
VcfFil=$VcfDir$JobNm.$JobNum.raw_variants.vcf #Output File
mkdir -p $VcfDir

infofields="-A AlleleBalance -A BaseQualityRankSumTest -A Coverage -A HaplotypeScore -A HomopolymerRun -A MappingQualityRankSumTest -A MappingQualityZero -A QualByDepth -A RMSMappingQuality -A SpanningDeletions -A VariantType" #Annotation fields to output into vcf files

#Start Log File
uname -a >> $TmpLog
echo "Start Variant Calling on Chromosome $CHR with GATK UnifiedGenotyper - $0:`date`" >> $TmpLog
echo "" >> $TmpLog
echo "Job name: "$JOB_NAME >> $TmpLog
echo "Job ID: "$JOB_ID >> $TmpLog
echo "Output Directory: "$VcfDir >> $TmpLog
echo "Output File: "$VcfFil >> $TmpLog
echo "Target file line range: $SttLn - $FinLn" >> $TmpLog

##Run Joint Variant Calling
mkdir -p $TmpDir

echo "Variant Calling with GATK UnifiedGenotyper..." >> $TmpLog
cmd="$JAVA7BIN -Xmx7G -Djava.io.tmpdir=$TmpDir -jar $GATKJAR  -T UnifiedGenotyper -R $REF -L $Range -nct $NumCores -I $BamLst -stand_emit_conf 10 -stand_call_conf 30 -o $VcfDir/$VcfFil -glm BOTH --dbsnp $DBSNP --comp:HapMapV3 $HpMpV3 $infofields -rf BadCigar"
echo $cmd >> $TmpLog
$cmd
if [[ $? == 1 ]]; then
	echo "----------------------------------------------------------------" >> $TmpLog
    echo "Variant Calling with GATK UnifiedGenotyper $JOB_NAME $JOB_ID failed `date`" >> $TmpLog
	qstat -j $JOB_ID | grep -E "usage *$SGE_TASK_ID:" >> $TmpLog
	cat $TmpLog >> $LogFil
	#rm $TmpLog $TmpDir
    exit 1
fi
echo "" >> $TmpLog
echo "Variant Calling done." >> $TmpLog
echo "" >> $TmpLog

#Need to wait for all HaplotypeCaller jobs to finish and then remerge all the vcfs
if [ $SGE_TASK_ID -eq 1 ]; then
	echo "Call Merge with vcftools with hold ...:" >> $TmpLog
	echo ""
	cmd="qsub -l $RmgVCFAlloc -N RmgVCF.$JobNm -o stdostde/ -e stdostde/ -hold_jid $JOB_NAME $EXOMSCR/ExmVC.3.MergeVCF.sh -d $VcfDir -s $Settings -l $LogFil"
	echo $cmd  >> $TmpLog
	$cmd
	echo "" >> $TmpLog
fi

echo "End Variant Calling with GATK UnifiedGenotyper on Chromosome $Chr $0:`date`" >> $TmpLog
echo ""
qstat -j $JOB_ID | grep -E "usage *$SGE_TASK_ID:" >> $TmpLog
echo "" >> $TmpLog
echo "===========================================================================================" >> $TmpLog
echo "" >> $TmpLog
cat $TmpLog >> $LogFil
rm -r $TmpLog $TmpDir
