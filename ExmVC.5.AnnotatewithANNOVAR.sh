#!/bin/bash
#$ -cwd -l mem=4G,time=2:: -N AnnVC

#This script takes a bam file or a list of bam files (filename must end ".list") and runs variant calling using the HaplotypeCaller in gVCF mode
#	InpFil - (required) - Path to Bam file to be aligned. Alternatively a file with a list of bams can be provided and the script run as an array job. List file name must end ".list"
#	RefFiles - (required) - shell file to export variables with locations of reference files, jar files, and resource directories; see list below
#	LogFil - (optional) - File for logging progress
#	Flag - P - PipeLine - call the next step in the pipeline at the end of the job
#	Help - H - (flag) - get usage information

#list of required vairables in reference file:
# $EXOMPPLN - directory containing exome analysis pipeline scripts
# $ANNHDB - directory containing databases for annovar

#list of required tools:
# annovar <http://www.openbioinformatics.org/annovar/> <http://www.openbioinformatics.org/annovar/annovar_download_form.php>

## This file also require exome.lib.sh - which contains various functions used throughout my Exome analysis scripts; this file should be in the same directory as this script

###############################################################

#set default arguments
usage="
ExmVC.5.AnnotatewithANNOVAR.sh -i <InputFile> -r <reference_file> -l <logfile> -PH

	 -i (required) - Path to VCF file or \".list\" file containing a multiple paths
	 -r (required) - shell file to export variables with locations of reference files and resource directories
	 -l (optional) - Log file
	 -P (flag) - Call next step of exome analysis pipeline after completion of script
	 -H (flag) - echo this message and exit
"

PipeLine="false"

while getopts i:r:l:PH opt; do
	case "$opt" in
		i) InpFil="$OPTARG";;
		r) RefFil="$OPTARG";; 
		l) LogFil="$OPTARG";;
		P) PipeLine="true";;
		H) echo "$usage"; exit;;
  esac
done

#load settings file
source $RefFil

#Load script library
source $EXOMPPLN/exome.lib.sh #library functions begin "func"

#Set local Variables
##Set local parameters
ArrNum=$SGE_TASK_ID
funcFilfromList #if the input is a list get the appropriate input file for this job of the array --> $InpFil
AnnNam=${InpFil/.vcf/}
if [[ -z $LogFil ]]; then LogFil=$AnnNam.AnnVC.log; fi # a name for the log file
TmpLog=$AnnNam.AnnVC.temp.log #temporary log file
AnnDir=$AnnNam.AnnVC.tempdir; mkdir -p $AnnDir
TmpVar=$AnnDir/$AnnNam.tempvar
AnnFil=$AnnNam.annovar

#Start Log File
ProcessName="Annotate VCF with ANNOVAR" # Description of the script - used in log
funcWriteStartLog

##Convert VCF to ANNOVAR input file using ANNOVAR - use old vcf method
StepNam="Convert VCF to ANNOVAR input file using ANNOVAR"
StepCmd="convert2annovar.pl $InpFil -format vcf4old -includeinfo | cut -f 1-10 > $TmpVar" #command to be run
funcRunStep

##Generate Annotation table
StepNam="Build Annotation table using ANNOVAR"
StepCmd="table_annovar23.pl $TmpVar $ANNHDB --buildver hg19 --remove -protocol refGene,esp6500si_all,esp6500si_aa,esp6500si_ea,1000g2012apr_all,1000g2012apr_eur,1000g2012apr_amr,1000g2012apr_asn,1000g2012apr_afr,ljb23_all,caddgt10 -operation g,f,f,f,f,f,f,f,f,f,f -otherinfo  -nastring \"\"  --outfile $AnnFil"
funcRunStep
AnnFil=$AnnFil.hg19_multianno.txt

##sort, replace spaces, zip and index
#regarding spaces - annovar output has these in the refgene function, e.g. "synonymous SNV", but they are not permitted in vcf format and other tools (e.g. GATK) will throw an error if they encounter them
head -n 1 $AnnFil > $AnnFil.tempheader
tail -n+2 $AnnFil | sort -V | awk '{gsub( / /, ""); print}' | awk '{gsub( /;/, ","); print}' >> $AnnFil.tempheader
mv $AnnFil.tempheader $AnnFil
bgzip $AnnFil
tabix  -s 1 -b 2 -e 3 $AnnFil.gz

#Annotate with vcftools
StepNam="Annotate with vcftools"
StepCmd="cat $InpFil | vcf-annotate -a $AnnFil.gz
 -c -,-,-,-,-,INFO/Func.refGene,INFO/Gene.refGene,INFO/ExonicFunc.refGene,INFO/AAChange.refGene,INFO/esp6500si_all,-,-,INFO/1000g2012apr_all,-,-,-,-,INFO/LJB23_SIFT_score,-,INFO/LJB23_SIFT_pred,-,-,INFO/LJB23_Polyphen2_HVAR_score,INFO/LJB23_Polyphen2_HVAR_pred,-,-,-,INFO/LJB23_MutationTaster_score,-,INFO/LJB23_MutationTaster_pred,-,-,-,-,-,-,-,-,-,-,-,INFO/LJB23_GERP++,INFO/LJB23_PhyloP,INFO/LJB23_SiPhy,INFO/caddgt,CHROM,POS,-,REF,ALT
 -d key=INFO,ID=Func.refGene,Number=1,Type=String,Description='function'
 -d key=INFO,ID=Gene.refGene,Number=1,Type=String,Description='refGene GeneName'
 -d key=INFO,ID=ExonicFunc.refGene,Number=1,Type=String,Description='functional Class'
 -d key=INFO,ID=AAChange.refGene,Number=1,Type=String,Description='Amino Acid change'
 -d key=INFO,ID=esp6500si_all,Number=1,Type=Float,Description='ESP6500 alternative allele frequency'
 -d key=INFO,ID=1000g2012apr_all,Number=1,Type=Float,Description='1000 genome alternative allele frequency'
 -d key=INFO,ID=LJB23_SIFT_score,Number=1,Type=Float,Description='SIFT score'
 -d key=INFO,ID=LJB23_SIFT_pred,Number=1,Type=String,Description='SIFT Prediction'
 -d key=INFO,ID=LJB23_Polyphen2_HVAR_score,Number=1,Type=Float,Description='PolyPhen2 HVAR Score'
 -d key=INFO,ID=LJB23_Polyphen2_HVAR_pred,Number=1,Type=Character,Description='PolyPhen2 HVAR Prediction'
 -d key=INFO,ID=LJB23_MutationTaster_score,Number=1,Type=Float,Description='MutationTaster Score'
 -d key=INFO,ID=LJB23_MutationTaster_pred,Number=1,Type=Character,Description='MutationTaster Prediction'
 -d key=INFO,ID=LJB23_GERP++,Number=1,Type=Float,Description='GERP++ score'
 -d key=INFO,ID=LJB23_PhyloP,Number=1,Type=Float,Description='PhyloP score'
 -d key=INFO,ID=LJB23_SiPhy,Number=1,Type=Float,Description='SiPhy scores'
 -d key=INFO,ID=caddgt,Number=1,Type=Float,Description='Whole-genome CADD score'
 > $InpFil.annotated.vcf"
funcRunStep


#End Log
funcWriteEndLog

#clean up
rm -rf $AnnDir
