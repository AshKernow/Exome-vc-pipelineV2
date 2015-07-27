#!/bin/bash


dbToDownload=$1

AnnovarDirectory=/ifs/scratch/c2b2/af_lab/ads2202/bin/annovar

latest1KG=1000g2014oct
latestdbSNP=snp138
latestLJB=ljb26_all
latestESPall=esp6500siv2_all
latestESPaa=esp6500siv2_aa
latestESPea=esp6500siv2_ea
latestCOMIC=cosmic70

cd $AnnovarDirectory

case $dbToDownload in
    1)
        nam="Download: refseq hg19 gene reference"
        cmd="perl annotate_variation.pl -downdb -buildver hg19 -webfrom annovar refGene humandb/"
        ;;
    2)
        nam="Download: 1000 genomes reference"
        cmd="perl annotate_variation.pl -downdb -buildver hg19 -webfrom annovar $latest1KG humandb/"
        ;;
    3)
        nam="Download: dbSNP reference"
        cmd="perl annotate_variation.pl -downdb -buildver hg19 -webfrom annovar $latestdbSNP humandb/"
        ;;
    4)
        nam="Download: LJB"
        cmd="perl annotate_variation.pl -downdb -buildver hg19 -webfrom annovar $latestLJB humandb/"
        ;;
    5)
        nam="Download: ESP alternative allele frequency - all"
        cmd="perl annotate_variation.pl -downdb -buildver hg19 -webfrom annovar $latestESPall humandb/"
        ;;
    6)
        nam="Download: ESP alternative allele frequency - African Americans"
        cmd="perl annotate_variation.pl -downdb -buildver hg19 -webfrom annovar $latestESPaa humandb/"
        ;;
    7)
        nam="Download: SP alternative allele frequency - European Americans"
        cmd="perl annotate_variation.pl -downdb -buildver hg19 -webfrom annovar $latestESPea humandb/"
        ;;
    8)
        nam="Download: superdups "
        cmd="perl annotate_variation.pl -downdb -buildver hg19 genomicSuperDups humandb/"
        ;;
    9)
        nam="Download: Cadd top 10% "
        cmd="perl annotate_variation.pl -downdb -buildver hg19 -webfrom annovar caddgt10 humandb/"
        ;;
    10)
        nam="Download: Cadd indel "
        cmd="perl annotate_variation.pl -downdb -buildver hg19 -webfrom annovar caddindel humandb/"
        ;;
    11)
        nam="Download: exac03"
        cmd="perl annotate_variation.pl -downdb -buildver hg19 -webfrom annovar  exac03 humandb/"
        ;;
    12)
        nam="Download: Cadd full "
        cmd="perl annotate_variation.pl -downdb -buildver hg19 -webfrom annovar cadd humandb/"
        ;;
    13)
        nam="Download: Cadd full "
        cmd="perl annotate_variation.pl -downdb -buildver hg19 -webfrom annovar clinvar_20150330 humandb/"
        ;;
    14)
        nam="Download: COSMIC "
        cmd="perl annotate_variation.pl -downdb -buildver hg19 -webfrom annovar $latestCOMIC humandb/"
        ;;
esac

echo "Start $nam - `date`"
echo $cmd
eval $cmd
if [[ $? == 0 ]]; then
    echo "Finish $nam - `date`"
else
    echo "Error during "${nam/load:/loading}
fi