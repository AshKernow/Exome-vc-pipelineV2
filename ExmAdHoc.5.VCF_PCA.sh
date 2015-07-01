#!/bin/bash
#$ -cwd -l mem=10G,time=4:: -N VcfPCA

#This script takes VCF file and runs PCA on the variants using Eigenstrat.
#    InpFil - (required) - A vcf file or plink bed/bim/fam trio of files. If plink, supply the bed file name (e.g. myfile.bed)
#    OutNam - (optional) - A name for the analysis - to be used for naming output files. Will be derived from input filename if not provided
#    Help - H - (flag) - get usage information

#list of required tools:
# 


###############################################################

usage="
ExmAdHoc.5.VCF_PCA.sh -i <InputFile> -l <logfile> -H

     -i (required) - A vcf input file
     -o (optional) - output name - will be derived from input filename if not provided
     -S (optional) - only analyse samples - omit HapMap data. Default is to include the HapMap data for comparison.
     -H (flag) - echo this message and exit
"

SamOnly="false"
while getopts i:o:SH opt; do
    case "$opt" in
        i) InpFil="$OPTARG";;
        o) OutNam="$OPTARG";;
        S) SamOnly="true";;
        H) echo "$usage"; exit;;
  esac
done

#some variables
EXOMFILT=/ifs/scratch/c2b2/af_lab/ads2202/Exome_Seq/scripts/Filtering_scripts
HapMapDir=/ifs/scratch/c2b2/af_lab/ads2202/Exome_Seq/resources/hapmap3_pop/hg19

InpFil=`readlink -f $InpFil`
if [[ -z "$OutNam" ]];then OutNam=`basename $InpFil`; OutNam=${OutNam/.bed/}; OutNam=${OutNam/.vcf/}; fi # a name for the output files

#check for vcf, if vcf convert to plink format
if [[ "${InpFil##*.}" != "bed" ]]; then
    VcfFil=`readlink -f $InpFil`
    BbfNam=$OutNam
    #filter the VCF for common variants for 1KG and GO-ESP and within cohort alternate allele frequency > 0.05
    $EXOMFILT/ExmFilt.1.FilterbyAlleleFrequency.py -v $VcfFil -o $BbfNam --maf 0.05 -G -W
    VcfFil=$BbfNam.filter.aaf.vcf
    # convert vcf --> plink using vcftools; if there are more than 1000 samples in the vcf vcftools cannot generate the necesary number of temporary files (cluster limitiation) so we need to do it in multiple batches and then remerge
    SamLst=$OutNam.samples.list
    grep -m 1 "^#CHROM" $VcfFil | cut -f 10- | tr '\t' '\n' > $SamLst
    SamLen=`cat $SamLst | wc -l`
    if [[ $SamLen -lt 1000 ]]; then
        vcftools --vcf $VcfFil --plink --out $BbfNam.tmp
        plink --file $BbfNam.tmp --make-bed --out $BbfNam
    else
        while [[ $SamLen -gt 0 ]]; do
            TmpSams=$SamLst.tmp
            head -n 1000 $SamLst > $TmpSams
            vcftools --vcf $VcfFil --keep $TmpSams --plink --out $BbfNam.tmp.$SamLen
            plink --file $BbfNam.tmp.$SamLen --make-bed --out $BbfNam.tmp.$SamLen
            tail -n +1001 $SamLst > $SamLst.2
            mv -f $SamLst.2 $SamLst
            SamLen=`cat $SamLst | wc -l`
            rm -f $TmpSams
        done
        ls $BbfNam.tmp.*bed | sed s/.bed//g > $BbfNam.tmp.splitlist
        plink --merge-list $BbfNam.tmp.splitlist --make-bed --out $BbfNam
        rm -f $BbfNam.tmp.splitlist
    fi
    rm -f $BbfNam.tmp* $SamLst
    echo
    echo "------------------------------------------------------------------------"
    echo
    #change -9 in the fam to 2
    awk '{ gsub( /-9$/, "2"); print }' $BbfNam.fam > $BbfNam.fam.temp
    mv -f $BbfNam.fam.temp $BbfNam.fam
else
    BbfNam=`readlink -f $InpFil`
    BbfNam=${BbfNam/.bed/}
    awk '{ gsub( /-9$/, "2"); print }' $BbfNam.fam > $BbfNam.fam.temp
    mv $BbfNam.fam.temp $BbfNam.fam
fi

# Get HapMap data unless excepted
if [[ "$SamOnly" == "false" ]]; then

    SnpList=$OutNam.snplist
    cut -f 2 $BbfNam.bim > $SnpList


    for i in $(find $HapMapDir | grep bed); do
        HapFil=${i/.bed/}
        HapTmp=${HapFil#*.}
        HapTmp=HAPMAPTEMP_${HapTmp%%.*}
        plink --bfile $HapFil --extract $SnpList --make-bed --noweb --out $HapTmp
        echo
        echo "------------------------------------------------------------------------"
        echo
    done

    ls HAPMAPTEMP*bed | sed s/.bed//g > HAPMAPTEMP_merge-list.list

    HapMapDat=$OutNam"_HapMapData"
    plink --merge-list HAPMAPTEMP_merge-list.list --make-bed --noweb --out $HapMapDat
    echo
    echo "------------------------------------------------------------------------"
    echo
    
    rm -f *HAPMAPTEMP*

    ###HapMap data is b36 so update map before merging
    cut -f 2,4 $BbfNam.bim > update_map.tab
    plink --bfile $HapMapDat --update-map update_map.tab --make-bed --noweb --out $HapMapDat
    echo
    echo "------------------------------------------------------------------------"
    echo
    
    #change -9 in fam to 1 
    awk '{ gsub( /-9$/, "1"); print }' $HapMapDat.fam > $HapMapDat.fam.temp
    mv -f $HapMapDat.fam.temp $HapMapDat.fam


    # Merge data:
    EigDat=Combined_data_for_eigenstrat
    plink --bfile $BbfNam --bmerge $HapMapDat.bed $HapMapDat.bim $HapMapDat.fam --geno 0.05 --make-bed --noweb --out $EigDat
    echo
    echo "------------------------------------------------------------------------"
    echo
    
    #check for mismatched snps and multiple position/chr snps and exclude and remerge if necessary
    grep "Warning: Multiple [cp]" Combined_data_for_eigenstrat.log | sed s/.*\ \'//g | sed s/\'.*//g > ExcludeSNPs.list
    
    cat $EigDat-merge.missnp >> ExcludeSNPs.list
    
    if [[ -s ExcludeSNPs.list ]]; then
        echo
        echo "------------------------------------------------------------------------"
        echo
        plink --bfile $HapMapDat --exclude ExcludeSNPs.list --make-bed --noweb --out $HapMapDat
        echo
        echo "------------------------------------------------------------------------"
        echo
        plink --bfile $BbfNam --bmerge $HapMapDat.bed $HapMapDat.bim $HapMapDat.fam --geno 0.05 --make-bed --allow-no-sex --out $EigDat
    fi
else
    EigDat=$BbfNam
fi

# Convert data to Eigenstrat format

cp $EigDat.fam $EigDat.pedind


echo genotypename: $EigDat.bed > par.BBF.EIGENSTRAT
echo snpname: $EigDat.bim >> par.BBF.EIGENSTRAT
echo indivname: $EigDat.pedind >> par.BBF.EIGENSTRAT
echo outputformat: EIGENSTRAT >> par.BBF.EIGENSTRAT
echo genooutfilename: $OutNam.eigenstratgeno >> par.BBF.EIGENSTRAT
echo snpoutfilename: $OutNam.snp >> par.BBF.EIGENSTRAT
echo indoutfilename: $OutNam.ind >> par.BBF.EIGENSTRAT

convertf -p par.BBF.EIGENSTRAT

# run EigenStrat

perl /ifs/home/c2b2/af_lab/ads2202/scratch/bin/eigenstrat/bin/smartpca.perl -i $OutNam.eigenstratgeno -a $OutNam.snp -b $OutNam.ind -k 10 -o $OutNam.plus.HapMap.pca -p $OutNam.plus.HapMap.plot -e $OutNam.plus.HapMap.eval -l $OutNam.plus.HapMap.log -m 5 -t 2 -s 6.0
