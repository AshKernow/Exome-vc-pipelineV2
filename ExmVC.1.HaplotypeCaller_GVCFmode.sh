#!/bin/bash
#$ -cwd 
# This script is an array job that will run the HaplotypeCaller in GVCF mode on a list of bam files to generate genomic VCFs, from which variants can then be called as per the GATK 3.0 pipeline. 
#It is intended to be the beginning of the variant calling pipeline for GATK3.0, however, it can be used in isolation on a list of bam files just to generate gVCF files. 
#It will run the HC on each bam file (one array job per file), and in pipeline mode it will then call the HaplotypeCaller again to joint call variants on all of the gVCFs produced. 
#To use it to start the pipeline add the -P flag.

#load function librarys
$(qstat -j $JOB_ID | grep "script_file" | awk '{print $2}')/exome.lib.sh

PipeLine="false"
while getopts i:s:j:l:P opt; do
  case "$opt" in
      i) BamLst="$OPTARG";;   # tab delimited, two colums: <Full path to BamFil>[TAB]<Name for Sample/BamFil>
      s) Settings="$OPTARG";; # settings file containing paths to tools and references
	  j) JobNam="$OPTARG";; #name for entire job array
      l) LogFil="$OPTARG";;   # name for log file for entire job
	  P) PipeLine="true";; # call the next step in the pipeline at the end of the job
  esac
done

#load settings file
. $Settings

#Set local Variables
JobNum=$SGE_TASK_ID
NumJobs=$(cat $BamLst | wc -l)
JobNam=${JOB_NAME#*.}
VcfDir=$JobNam"_gVCF" #Output Directory
mkdir -p $VcfDir
BamFil=$(tail -n+$JobNum $BamLst | head -n 1 | cut -f1)
BamNam=$(tail -n+$JobNum $BamLst | head -n 1 | cut -f2)
TmpLog=$LogFil.$BamNam.CallgVC.log
TmpDir=$JobNam.$BamNam.gVC
mkdir -p $TmpDir
VcfFil=$VcfDir/$JobNam.$JobNum.gvcf #Output File


infofields="-A AlleleBalance -A BaseQualityRankSumTest -A Coverage -A HaplotypeScore -A HomopolymerRun -A MappingQualityRankSumTest -A MappingQualityZero -A QualByDepth -A RMSMappingQuality -A SpanningDeletions -A FisherStrand" #Annotation fields to output into vcf files

#Start Log File
uname -a >> $TmpLog
echo "Start genomic VCF generatation with GATK HaplotypeCaller - $0:`date`" >> $TmpLog
echo "" >> $TmpLog
echo "Job name: "$JOB_NAME >> $TmpLog
echo "Job ID: "$JOB_ID >> $TmpLog
echo "Output Directory: "$VcfDir >> $TmpLog
echo "Input File: " $BamFil >> $TmpLog
echo "Output File: "$VcfFil >> $TmpLog


##Run Joint Variant Calling
StepNam"gVCF generation with GATK HaplotypeCaller..." >> $TmpLog
StepCmd="$JAVA7BIN -Xmx7G -Djava.io.tmpdir=$TmpDir -jar $GATKJAR  -T HaplotypeCaller -R $REF -L $TARGET -nct $NumCores -I $BamFil --genotyping_mode DISCOVERY -stand_emit_conf 10 -stand_call_conf 30 --emitRefConfidence GVCF --variant_index_type LINEAR --variant_index_parameter 128000 -o $VcfDir/$VcfFil --dbsnp $DBSNP --comp:HapMapV3 $HpMpV3 $infofields -rf BadCigar"
funcRunStep

#End log
echo "Genomic VCF generatation with GATK HaplotypeCaller - $0:`date`" >> $TmpLog
echo ""
qstat -j $JOB_ID | grep -E "usage *$SGE_TASK_ID:" >> $TmpLog
echo "" >> $TmpLog
echo "===========================================================================================" >> $TmpLog
echo "" >> $TmpLog
cat $TmpLog >> $LogFil

#Clean up
rm -r $TmpLog $TmpDir
