#!/usr/bin/env python
#$ -cwd -l mem=2G,time=:15: -N AAFFilt

# The purpose of this script is to filter a VCF file for variants with delterious exonic alternate alleles. Ouputs a filtered VCF file, an plink/Seq variant mask file, and a log file.
#    -v/--vcf     <required>    The script requires an input vcf
#    -o/--out     <required>    The use should specify a base name for the output files.The script outputs the filtered results as a vcf file. The script also outputs a log file. 


from optparse import OptionParser
import os

parser = OptionParser()
# Basic Input Files, required
parser.add_option("-v", "--vcf", dest="VCFfile",help="input VCF file", metavar="VCFfile")
parser.add_option("-o", "--output", dest="OutputFileName",help="user specified name of file for output to be written", metavar="OutputFileName")


parser.add_option("-T", "--trunc", action='store_true', dest="trunc", help="Filter only truncating")
parser.add_option("-z", "--debug", action='store_true', dest="debug", help="Output a boolean decision file for debugging purposes")

(options, args) = parser.parse_args()

VCF=open(options.VCFfile,'r')
BaseName=str(options.OutputFileName)

#open input and output files
VCF=open(options.VCFfile,'r')
BaseName=str(options.OutputFileName)
OutputFilename=BaseName+'.deleterious.vcf'
PseqFilename=BaseName+'.deleterious.pseqmask'
LogOutputFilename=BaseName+'.deleterious.log'
Trunc=options.trunc
DeBug=options.debug
if Trunc:
    OutputFilename=BaseName+'.deleterious_truncating.vcf'
    PseqFilename=BaseName+'.deleterious_truncating.pseqmask'
    LogOutputFilename=BaseName+'.deleterious_truncating.log'
OutVcf=open(OutputFilename,'w')
OutPseq=open(PseqFilename,'w')
Outlog=open(LogOutputFilename,'w')

if DeBug:
    PassOutputFilename=BaseName+'.deletrious.boolean.log'
    OutPass=open(PassOutputFilename,'w')
    
    
    passheader=['Chromosome','Position','PassFilt','PassFunction','PassClass','PassDamaging','VariantFunction','VariantClass''SIFTprediction','PP2prediction','CADDscore']    
    OutPass.write("\t".join(passheader)+"\n")

#start log file
import datetime
TimeNow=str(datetime.datetime.now().strftime("%Y-%m-%d %H:%M"))
Outlog.write("Deleterious Alternate Allele Filtering log: "+TimeNow+"\n")
Outlog.write("Input VCF: "+os.path.abspath(options.VCFfile)+"\n")
    
OrigCount=0
FiltTrunc=0
FiltMissense=0
ChrCount=0

TruncatingList=['frameshiftdeletion','frameshiftinsertion','frameshiftsubstitution','stopgainSNV','stoplossSNV']

MissenseList=['nonsynonymousSNV','nonframeshiftdeletion','nonframeshiftinsertion','nonframeshiftsubstitution']


print "Start filtering..."
for line in VCF:
    OrigCount=OrigCount+1
    # Output vcf header to new vcf
    if '#' in line:
        OutVcf.write(line)
    # Start filtering variants
    if '#' not in line:
        # Variant must first pass 1KG and GO-ESP frequencies, MQ0 threshold, and be exonic
        linelist=line.split("\t")
        INFOstring=linelist[7]
        INFOcolumnList=INFOstring.split(";")
        INFOdict={}
        ChrPresent=linelist[0]
        if ChrCount != ChrPresent:
            print "Chromosome "+str(ChrPresent)
            ChrCount=ChrPresent
            
        for element in INFOcolumnList:
            if '=' in element:
                FieldName,FieldValue=element.split('=',1)
                INFOdict[FieldName]=FieldValue
        
        # Get values for later
        VariantFunction=INFOdict.get('VarFunc','none')
        VariantFunctionList=VariantFunction.split(',')
        VariantClass=INFOdict.get('VarClass','none')
        VariantClassList=VariantClass.split(',')
        SIFTprediction=INFOdict.get('SIFTprd','NA')
        PP2prediction=INFOdict.get('PP2prd','NA')
        CADDscore=float(INFOdict.get('CADDphred',0))
        

        # Check if Function passes
        PassFunction=False
        if VariantFunction=='exonic' or VariantFunction=='splicing' or VariantFunction=="exonic,splicing" or VariantFunction=='none':
            PassFunction=True
        
        # Check if Variant class passes
        PassTrunc=False
        if VariantClass in TruncatingList:
            if SIFTprediction!='T' or PP2prediction!='B':
                PassTrunc=True
        
        # Check if missense and damaging
        PassMissense=False
        if VariantClass in MissenseList:
            if SIFTprediction=='D' and PP2prediction=='D':
                PassMissense=True
            if CADDscore >= 30:
                PassMissense=True
                
        
        Position='chr'+str(linelist[0])+":"+str(linelist[1])
        if PassFunction and PassTrunc:
            OutVcf.write(line)
            MaskLine=['VAR','Truncating',Position]
            OutPseq.write("\t".join(MaskLine)+"\n")
            MaskLine=['VAR','Deleterious',Position]
            OutPseq.write("\t".join(MaskLine)+"\n")
            FiltTrunc=FiltTrunc+1
        
        
        if PassFunction and PassMissense and not Trunc:
            OutVcf.write(line)
            MaskLine=['VAR','DamagingMissense',Position]
            OutPseq.write("\t".join(MaskLine)+"\n")
            MaskLine=['VAR','Deleterious',Position]
            OutPseq.write("\t".join(MaskLine)+"\n")
            FiltMissense=FiltMissense+1
        
        if DeBug:
            passLine=linelist[0:1]+[str(PassFunction),str(PassTrunc),str(PassMissense),VariantFunction,VariantClass,SIFTprediction,PP2prediction,str(CADDscore)]
            OutPass.write("\t".join(passLine)+"\n")
FiltCount= FiltMissense + FiltTrunc
Outlog.write("  Number of variants in original VCF: "+str(OrigCount)+"\n")
Outlog.write("  Number of damaging variants - Truncating: "+str(FiltTrunc)+"\n")
Outlog.write("  Number of damaging variants - Missense: "+str(FiltMissense)+"\n")
Outlog.write("  Number of damaging variants - Total: "+str(FiltCount)+"\n")
TimeNow=str(datetime.datetime.now().strftime("%Y-%m-%d %H:%M"))
Outlog.write("  Filtering Finished: "+TimeNow+"\n")
print "Done"