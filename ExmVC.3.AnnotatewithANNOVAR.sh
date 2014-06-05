#!/bin/bash
#$ -cwd -l mem=8G,time=4:: -N AnnVC

#This script takes a bam file or a list of bam files (filename must end ".list") and runs variant calling using the HaplotypeCaller in gVCF mode
#    InpFil - (required) - Path to Bam file to be aligned. Alternatively a file with a list of bams can be provided and the script run as an array job. List file name must end ".list"
#    RefFil - (required) - shell file containing variables with locations of reference files, jar files, and resource directories; see list below for required variables
#    LogFil - (optional) - File for logging progress
#    Flag - P - PipeLine - call the next step in the pipeline at the end of the job
#    Flag - B - BadET - prevent GATK from phoning home
#    Help - H - (flag) - get usage information

#list of required vairables in reference file:
# $EXOMPPLN - directory containing exome analysis pipeline scripts
# $ANNHDB - directory containing databases for annovar

#list of required tools:
# annovar <http://www.openbioinformatics.org/annovar/> <http://www.openbioinformatics.org/annovar/annovar_download_form.php>

## This file also requires exome.lib.sh - which contains various functions used throughout the Exome analysis scripts; this file should be in the same directory as this script

###############################################################

#set default arguments
usage="
ExmVC.5.AnnotatewithANNOVAR.sh -i <InputFile> -r <reference_file> -l <logfile> -PH

     -i (required) - Path to VCF file or \".list\" file containing a multiple paths
     -r (required) - shell file containing variables with locations of reference files and resource directories
     -l (optional) - Log file
     -P (flag) - Call next step of exome analysis pipeline after completion of script
     -B (flag) - Prevent GATK from phoning home - only if calling pipeline
     -H (flag) - echo this message and exit
"

PipeLine="false"

while getopts i:r:l:PBH opt; do
    case "$opt" in
        i) InpFil="$OPTARG";;
        r) RefFil="$OPTARG";; 
        l) LogFil="$OPTARG";;
        P) PipeLine="true";;
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
##Set local parameters
ArrNum=$SGE_TASK_ID
funcFilfromList #if the input is a list get the appropriate input file for this job of the array --> $InpFil
AnnNam=${InpFil/.vcf/}
if [[ -z $LogFil ]]; then LogFil=$AnnNam.AnnVC.log; fi # a name for the log file
TmpLog=$AnnNam.AnnVC.temp.log #temporary log file
AnnDir=$AnnNam.AnnVC.tempdir; mkdir -p $AnnDir
TmpVar=$AnnDir/$AnnNam.tempvar
AnnFil=$AnnNam.annovar
VcfFil=${InpFil/vcf/annotated.vcf}

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

##sort, replace spaces and semi-colons, zip and index
# - annovar output has spaces in the RefSeq function code, e.g. "synonymous SNV", but they are not permitted in vcf format and other tools (e.g. GATK) will throw an error if they encounter them
# - annovar separates multiple gene names in the RefSeq gene name field with semi-colons, this causes and error in the vcf
head -n 1 $AnnFil > $AnnFil.tempheader
tail -n+2 $AnnFil | sort -V | awk '{gsub( / /, ""); print}' | awk '{gsub( /;/, ","); print}' >> $AnnFil.tempheader
mv $AnnFil.tempheader $AnnFil
bgzip $AnnFil
tabix  -s 1 -b 2 -e 3 $AnnFil.gz

#Annotate with vcftools
StepNam="Annotate with vcftools"
StepCmd="cat $InpFil | vcf-annotate -a $AnnFil.gz
 -c -,-,-,-,-,INFO/VarFunc,INFO/GeneName,INFO/VarClass,INFO/AAChange,INFO/ESPfreq,-,-,INFO/1KGfreq,-,-,-,-,INFO/SIFTscr,-,INFO/SIFTprd,-,-,INFO/PP2scr,INFO/PP2prd,-,-,-,INFO/MutTscr,-,INFO/MutTprd,-,-,-,-,-,-,-,-,-,-,-,INFO/GERP,INFO/PhyloP,INFO/SiPhy,INFO/CADD,CHROM,POS,-,REF,ALT
 -d key=INFO,ID=VarFunc,Number=1,Type=String,Description='Genomic region/Sequence Function'
 -d key=INFO,ID=GeneName,Number=1,Type=String,Description='refGene GeneName'
 -d key=INFO,ID=VarClass,Number=1,Type=String,Description='Mutational Class'
 -d key=INFO,ID=AAChange,Number=1,Type=String,Description='Amino Acid change'
 -d key=INFO,ID=ESPfreq,Number=1,Type=Float,Description='Exome Sequencing Project 6500 alternative allele frequency'
 -d key=INFO,ID=1KGfreq,Number=1,Type=Float,Description='1000 genome alternative allele frequency'
 -d key=INFO,ID=SIFTscr,Number=1,Type=Float,Description='SIFT score'
 -d key=INFO,ID=SIFTprd,Number=1,Type=String,Description='SIFT prediction'
 -d key=INFO,ID=PP2scr,Number=1,Type=Float,Description='PolyPhen2 HVAR score'
 -d key=INFO,ID=PP2prd,Number=1,Type=Character,Description='PolyPhen2 HVAR prediction'
 -d key=INFO,ID=MutTscr,Number=1,Type=Float,Description='MutationTaster score'
 -d key=INFO,ID=MutTprd,Number=1,Type=Character,Description='MutationTaster prediction'
 -d key=INFO,ID=GERP,Number=1,Type=Float,Description='GERP++ score'
 -d key=INFO,ID=PhyloP,Number=1,Type=Float,Description='PhyloP score'
 -d key=INFO,ID=SiPhy,Number=1,Type=Float,Description='SiPhy scores'
 -d key=INFO,ID=CADD,Number=1,Type=Float,Description='Whole-genome CADD score'
 > $VcfFil"
funcRunStep

#Get VCF stats with python script
StepNam="Get VCF stats"
StepCmd="python $EXOMPPLN/VCF_summary_Stats.py -v $VcfFil -o ${InpFil/vcf/stats.tsv}"
funcRunStep

#Call next steps of pipeline if requested
NextJob="Recalibrate Variant Quality"
QsubCmd="qsub -o stdostde/ -e stdostde/ $EXOMPPLN/ExmVC.4.RecalibrateVariantQuality.sh -i $VcfFil -r $RefFil -l $LogFil -P"
if [[ "$AllowMisencoded" == "true" ]]; then QsubCmd=$QsubCmd" -A"; fi
if [[ "$BadET" == "true" ]]; then QsubCmd=$QsubCmd" -B"; fi 
funcPipeLine

#End Log
funcWriteEndLog

#clean up
rm -rf $AnnDir
