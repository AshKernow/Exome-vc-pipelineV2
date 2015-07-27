#!/bin/bash
#$ -cwd -l mem=8G,time=2:: -N MakeKin

#This script takes a VCF file and generates files to check the familial relationships and sex of the samples
#    InpFil - (required) - Path to VCF file or a list of VCF Files to be recalibrated
#    OutNam - (optional) - Name for output files
#    LogFil - (optional) - File for logging progress
#    Help - H - (flag) - get usage info#rmation

#list of required vairables in reference file:
# None

#list of required tools:
# samtools
# vcftools
# R

## This file also requires exome.lib.sh - which contains various functions used throughout the Exome analysis scripts; this file should be in the same directory as this script

###############################################################

#set default arguments
usage="ExmVC.3.RecalibrateVariantQuality.sh -i <InputFile> -o <output_file> -l <logfile> -PABH

     -i (required) - Path to VCF file
     -o (optional) - Base for output file names <default is VCF file name>
     -l (optional) - Log file
     -H (flag) - echo this message and exit
"

while getopts i:o:l:H opt; do
    case "$opt" in
        i) InpFil="$OPTARG";;
        o) OutNam="$OPTARG";; 
        l) LogFil="$OPTARG";;
        H) echo "$usage"; exit;;
  esac
done

#check all required paramaters present
if [[ ! -e "$InpFil" ]]; then echo "Missing/Incorrect required arguments"; echo "$usage"; exit; fi

#Load script library
source $EXOMPPLN/exome.lib.sh #library functions begin "func" #library functions begin "func"

#set local variables
VcfFil=`readlink -f $InpFil` #resolve absolute path to vcf
if [[ ! $OutNam ]]; then OutNam=`basename $InpFil | sed s/.vcf// `; fi
if [[ ! $LogFil ]];then LogFil=$OutNam.kinshipfiles.log; fi
TmpLog=Temp.$LogFil
TmpRscript=Temp.rPlotHist.$OutNam.R
SampList=Temp.$OutNam.`date "+%j%y%H%M%N"`.samplist #file to hold list of samples in vcf
MrgList=Temp.$OutNam.splitplink.merge.list #a file to hold the merge list for plink if it is necessary to convert the vcf in steps

#Start Log File
ProcessName="Generate pedigree analysis files" # Description of the script - used in log
funcWriteStartLog

#get relatedness via KING algorithm
StepName="Get relatedness via KING algorithm in vcftools" # Description of this step - used in log
StepCmd="qsss -c \"vcftools --vcf $VcfFil --relatedness2 --out $OutNam\""
funcRunStep

#get missingness for each individual
StepName="Get missingness with vcftools" # Description of this step - used in log
StepCmd="qsss -c \"vcftools --vcf $VcfFil --missing-indv --out $OutNam.vcf\""
funcRunStep

#convert to plink - the vcftools vcf --> plink tool cannot run with > 1000 samples on the cluster due to limitations on temporary files. Therefore, for larger cohorts, it is necessary to split and remerge the plink files.
grep -m 1 "^#CHROM" $VcfFil | cut -f 10- | tr "\t" "\n" > $SampList
LEN=`cat $SampList | wc -l`
echo "There are $LEN samples in the vcf" >> $TmpLog
if [[ $LEN -lt 1000 ]]; then
    StepName="Convert vcf to plink ped/map using vcftools"
    StepCmd="vcftools --vcf $VcfFil --plink --out Temp.$OutNam"
    funcRunStep
    StepName="Convert ped/map to bed/bim/fam using plink"
    StepCmd="plink --file Temp.$OutNam --make-bed --out $OutNam"
    funcRunStep
else
    split -l 900 $SampList $SampList.split.
    echo "Convert the vcf in "`ls | grep $SampList.split | wc -l`" steps then merge:"
    for i in $SampList.split*; do
        NAM=${i##*split.}
        StepName="Convert vcf to plink ped/map using vcftools - step $NAM"
        StepCmd="vcftools --vcf $VcfFil --keep $i --plink --out Temp.$OutNam.splitplink.$NAM"
        funcRunStep
    done
    ls -r | grep splitplink | grep -E "ped$|map$" | paste - - > $MrgList
    echo "Merge list:" >> $TmpLog
    cat $MrgList >> $TmpLog
    StepName="Merge ped/map split files and convert to bed/bim/fam using plink"
    StepCmd="plink --merge-list $MrgList --make-bed --out $OutNam"
    funcRunStep
fi

#split the pseudo-autosomal X into "XY" for sex check
StepName="Split the pseudo-autosomal X into XY for sex check with plink"
StepCmd="plink --bfile $OutNam --split-x hg19 --make-bed --out $OutNam"
funcRunStep

#run sex check
StepName="Run imputation of sex using plink"
StepCmd="plink --bfile $OutNam --impute-sex --make-bed --out $OutNam.sexcheck"
funcRunStep
mv $OutNam.sexcheck.sexcheck $OutNam.sexcheck


StepName="Plot histogram of F statistics"
echo "pdf(\"$OutNam.plink_sexcheck_histogram.pdf\")
 hist(read.table(\"$OutNam.sexcheck\", header=T)[,\"F\"], main=\"$OutNam plink sexcheck - F statistic\", xlab=\"F\")
 dev.off()" > $TmpRscript
StepCmd="Rscript $TmpRscript"
funcRunStep



#Run IBD with rare variants
StepName="Run IBD with rare variants using plink"
StepCmd="plink --bfile $OutNam --maf 0.01 --genome --out $OutNam"
funcRunStep

#get missingness
StepName="Get missingness with plink"
StepCmd="plink --bfile $OutNam --missing --out $OutNam"
funcRunStep

#End Log
funcWriteEndLog

#Clean up
#rm -f Temp.$OutNam.* $OutNam*.sexcheck.* $OutNam*nosex $OutNam.ped $OutNam.map  $OutNam.log $OutNam.lmiss $OutNam*~
