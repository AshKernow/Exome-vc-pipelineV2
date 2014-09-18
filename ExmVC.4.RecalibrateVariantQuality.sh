#!/bin/bash
#$ -cwd -l mem=24G,time=12:: -N VQSRFilt

#This script takes a raw VCF file and performs GATK's variant quality score recalibration
#    InpFil - (required) - Path to VCF file or a list of VCF Files to be recalibrated
#    RefFil - (required) - shell file containing variables with locations of reference files, jar files, and resource directories; see list below for required variables
#    LogFil - (optional) - File for logging progress
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
usage="ExmVC.3.RecalibrateVariantQuality.sh -i <InputFile> -r <reference_file> -t <targetfile> -l <logfile> -PABH

     -i (required) - Path to list of Bam files for variant calling
     -r (required) - shell file containing variables with locations of reference files and resource directories
     -l (optional) - Log file
     -B (flag) - Prevent GATK from phoning home
     -H (flag) - echo this message and exit
"

BadET="false"

PipeLine="false"
while getopts i:r:l:BH opt; do
    case "$opt" in
        i) InpFil="$OPTARG";;
        r) RefFil="$OPTARG";; 
        l) LogFil="$OPTARG";;
        B) BadET="true";;
        H) echo "$usage"; exit;;
  esac
done

#check all required paramaters present
if [[ ! -e "$InpFil" ]] || [[ ! -e "$RefFil" ]]; then echo "Missing/Incorrect required arguments"; echo "$usage"; exit; fi

#Call the RefFil to load variables
RefFil=`readlink -f $RefFil`
source $RefFil

#Load script library
source $EXOMPPLN/exome.lib.sh #library functions begin "func" #library functions begin "func"

#Set local Variables
ArrNum=$SGE_TASK_ID
funcFilfromList #if the input is a list get the appropriate input file for this job of the array --> $InpFil
VcfFil=`readlink -f $InpFil` #resolve absolute path to bam
HapChec=$(head -n20 $VcfFil | grep "HaplotypeCaller" | wc -l) #check which VC tool was used
if [[ $HapChec -eq 1 ]]; then
    InfoFields="-an DP -an QD -an FS -an MQRankSum -an ReadPosRankSum"
else
    InfoFields="-an DP -an QD -an FS -an MQRankSum -an ReadPosRankSum -an HaplotypeScore"
fi
VcfNam=`basename $VcfFil | sed s/.vcf// | sed s/.list//` #a name to use for the various files
if [[ -z "$LogFil" ]];then LogFil=$VcfNam.VQSR.log; fi # a name for the log file
GatkLog=$VcfNam.gatklog #a log for GATK to output to, this is then trimmed and added to the script log
TmpLog=$VcfNam.VQSR.temp.log #temporary log file 
TmpDir=$VcfNam.VQSR.tempdir; mkdir -p $TmpDir #temporary directory

#Start Log File
ProcessName="Variant Quality Score Recalibration GATK" # Description of the script - used in log
funcWriteStartLog

##Build the SNP recalibration model
StepName="Build the SNP recalibration model with GATK VariantRecalibrator" # Description of this step - used in log
StepCmd="java -Xmx9G -Djava.io.tmpdir=$TmpDir -jar $GATKJAR
 -T VariantRecalibrator 
 -R $REF
 -input $VcfFil
 -resource:hapmap,known=false,training=true,truth=true,prior=15.0 $HAPMAP
 -resource:omni,known=false,training=true,truth=true,prior=12.0 $TGVCF
 -resource:1000G,known=false,training=true,truth=false,prior=10.0 $ONEKG
 -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $DBSNP
  $InfoFields
 -mode SNP
 -tranche 100.0
 -tranche 99.9
 -tranche 99.0
 -tranche 90.0
 -recalFile $VcfNam.recalibrate_SNP.recal
 -tranchesFile $VcfNam.recalibrate_SNP.tranches
 -rscriptFile $VcfNam.recalibrate_SNP_plots.R
 -log $GatkLog" #command to be run
funcGatkAddArguments # Adds additional parameters to the GATK command depending on flags (e.g. -B or -F)
funcRunStep

##Apply SNP recalibration
StepName="Apply SNP recalibration with GATK ApplyRecalibration" # Description of this step - used in log
StepCmd="java -Xmx9G -Djava.io.tmpdir=$TmpDir -jar $GATKJAR
 -T ApplyRecalibration
 -R $REF
 -input $VcfFil
 -mode SNP
 --ts_filter_level 99.0
 -recalFile $VcfNam.recalibrate_SNP.recal
 -tranchesFile $VcfNam.recalibrate_SNP.tranches
 -o $VcfNam.recal_snps.vcf
 -log $GatkLog" #command to be run
funcGatkAddArguments # Adds additional parameters to the GATK command depending on flags (e.g. -B or -F)
funcRunStep
VcfFil=$VcfNam.recal_snps.vcf

##Build the InDel recalibration model
StepName="Build the InDel recalibration model with GATK VariantRecalibrator" # Description of this step - used in log
InfoFields="-an DP -an FS -an MQRankSum -an ReadPosRankSum"
StepCmd="java -Xmx9G -Djava.io.tmpdir=$TmpDir -jar $GATKJAR
 -T VariantRecalibrator
 -R $REF
 -input $VcfFil
 -resource:mills,known=true,training=true,truth=true,prior=12.0 $INDEL
 -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $DBSNP 
  $InfoFields
 -mode INDEL
 -tranche 100.0
 -tranche 99.9
 -tranche 99.0
 -tranche 90.0
 --maxGaussians 4
 -recalFile $VcfNam.recalibrate_INDEL.recal
 -tranchesFile $VcfNam.recalibrate_INDEL.tranches
 -rscriptFile $VcfNam.recalibrate_INDEL_plots.R
 -log $GatkLog" #command to be run
funcGatkAddArguments # Adds additional parameters to the GATK command depending on flags (e.g. -B or -F)
funcRunStep

##Apply InDel recalibration
StepName="Apply InDel recalibration with GATK ApplyRecalibration" # Description of this step - used in log
StepCmd="java -Xmx9G -Djava.io.tmpdir=$TmpDir -jar $GATKJAR
 -T ApplyRecalibration
 -R $REF
 -input $VcfFil
 -mode INDEL
 --ts_filter_level 99.0
 -recalFile $VcfNam.recalibrate_INDEL.recal
 -tranchesFile $VcfNam.recalibrate_INDEL.tranches
 -o $VcfNam.recalibrated.vcf
 -log $GatkLog" #command to be run
funcGatkAddArguments # Adds additional parameters to the GATK command depending on flags (e.g. -B or -F)
funcRunStep
rm $VcfFil
VcfFil=$VcfNam.recalibrated.vcf

#Apply Hard Filters to VCF
StepName="Apply Hard Filters with GATK Variant Filtration" # Description of this step - used in log
StepCmd="java -Xmx9G -Djava.io.tmpdir=$TmpDir -jar $GATKJAR
 -T VariantFiltration
 -R $REF
 --variant $VcfFil
 -o $VcfNam.hardfiltered.vcf
 --clusterWindowSize 10
 --filterExpression \"QUAL<30.0||QD<2.0\"
 --filterName \"StandardFilters\"
 --filterExpression \"MQ0>=4&&(MQ0/DP)>0.1\"
 --filterName \"HARD_TO_VALIDATE\"
 --missingValuesInExpressionsShouldEvaluateAsFailing
 -log $GatkLog" #command to be run
funcGatkAddArguments # Adds additional parameters to the GATK command depending on flags (e.g. -B or -F)
funcRunStep
rm $VcfFil
VcfFil=$VcfNam.hardfiltered.vcf

#Filter VCF with vcftools
StepName="Create filtered vcf using vcftools" # Description of this step - used in log
StepCmd="vcftools --vcf $VcfFil --remove-filtered-all --recode --recode-INFO-all --out $VcfNam.filtered"
funcRunStep
mv $VcfNam.filtered.recode.vcf $VcfNam.filtered.vcf
VcfFil=$VcfNam.filtered.vcf

#Get VCF stats with python script
StepNam="Get VCF stats"
StepCmd="python $EXOMPPLN/ExmPY.VCF_summary_Stats.py -v $VcfFil -o ${VcfFil/vcf/stats.tsv}"
funcRunStep

#End Log
funcWriteEndLog

#Clean up
