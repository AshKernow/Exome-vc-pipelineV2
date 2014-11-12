#!/usr/bin/env python
#$ -cwd -l mem=2G,time=:15: -N AAFFilt

# The purpose of this script is to filter a VCF file by minor/alternate allele frequency as provided by 1KG and GO-ESP
#    -v/--vcf     <required>    The script requires an input vcf
#    -o/--out     <required>    The use should specify a base name for the output files.The script outputs the filtered results as a vcf file. The script also outputs a log file. 
#    -m/--maf     <optional>    Minor allele frequency for filtering. Default is 0.01
#    -G/--greater <flag>        Filter for variants with maf greater than or equal to the filter level. Default is less than or equal to.
#    -W/--within  <flag>        Filter for allele frequency within the cohort. Default is just 1KG and ESP frequencies.


from optparse import OptionParser
import os

parser = OptionParser()
# Basic Input Files, required
parser.add_option("-v", "--vcf", dest="VCFfile",help="input VCF file", metavar="VCFfile")
parser.add_option("-o", "--output", dest="OutputFileName",help="user specified name of file for output to be written", metavar="OutputFileName")
parser.set_defaults(MAFFilter="0.01")
parser.add_option("-m", "--maf", dest="MAFFilter",help="user specified allele frequency", metavar="MAFFilter")
parser.set_defaults(GREATER=False)
parser.add_option("-G", "--greater", action='store_true', dest="GREATER", help="Filter for variants with maf greater than or equal to the filter level. Default is less than or equal to.")
parser.set_defaults(WITHIN=False)
parser.add_option("-W", "--within", action='store_true', dest="WITHIN", help="Filter for allele frequency within the cohort. Default is just 1KG and ESP frequencies.")


(options, args) = parser.parse_args()

VCF=open(options.VCFfile,'r')
BaseName=str(options.OutputFileName)
MafCutOff=float(options.MAFFilter)
GreaterThan=options.GREATER
WithinCohort=options.WITHIN

#open input and output files
VcfOutputFilename=BaseName+'.filter.aaf.vcf'
LogOutputFilename=BaseName+'.filter.aaf.log'
Outvcf=open(VcfOutputFilename,'w')
Outlog=open(LogOutputFilename,'w')

#start log file
import datetime
TimeNow=str(datetime.datetime.now().strftime("%Y-%m-%d %H:%M"))
Outlog.write("Allele Frequency Filtering log: "+TimeNow+"\n")
Outlog.write("Input VCF: "+os.path.abspath(options.VCFfile)+"\n")
Outlog.write("  Alternate allele frequency maximum: "+str(MafCutOff)+"\n")
if GreaterThan:
    Outlog.write("  Comparison: Greater than or equal to.\n")
else:
    Outlog.write("  Comparison: Less than or equal to.\n")
if WithinCohort:
    Outlog.write("  Filter applied to: 1000 Genome, GO-ESP, within Cohort frequency.\n")
else:
    Outlog.write("  Filter applied to: 1000 Genome and GO-ESP.\n")
    
OrigCount=0
FiltCount=0
ChrCount=0
print "Start filtering..."
for line in VCF:
    OrigCount=OrigCount+1
    # Output vcf header to new vcf
    if '#' in line:
        Outvcf.write(line)
    # Start filtering variants
    if '#' not in line:
        linelist=line.split("\t")
        ChrPresent=linelist[0]
        if ChrCount != ChrPresent:
            print "Chromosome "+str(ChrPresent)
            ChrCount=ChrPresent
        # Variant must first pass 1KG and GO-ESP frequencies, MQ0 threshold, and be exonic
        INFOstring=linelist[7]
        INFOcolumnList=INFOstring.split(";")
        INFOdict={}
        for element in INFOcolumnList:
            if '=' in element:
                FieldName,FieldValue=element.split('=',1)
                INFOdict[FieldName]=FieldValue
        
        # Get values for later
        KGscore=INFOdict.get('1KGfreq',0)
        if type(KGscore) is str:
            KGscore=KGscore.split(",")
            KGscore=float(KGscore[0])
            
        # Get withing cohort values for later
        AFscore=INFOdict.get('AF',0)
        if type(AFscore) is str:
            AFscore=AFscore.split(",")
            AFscore=float(AFscore[0])
        
        ESPscore=INFOdict.get('ESPfreq',0)
        if type(ESPscore) is str:
            ESPscore=ESPscore.split(",")
            ESPscore=float(ESPscore[0])
        
        
        # Check if KG passes threshold
        PassMAF=True
        #First if less than/equal to (default)
        if ( float(KGscore) < MafCutOff or float(ESPscore) < MafCutOff ) and GreaterThan:
            PassMAF=False
        
        #First if less than/equal to (default)
        if ( float(KGscore) > MafCutOff or float(ESPscore) > MafCutOff ) and not GreaterThan:
            PassMAF=False
        
        #check within cohort if requested
        PassCHT=True
        if WithinCohort:
            if AFscore < MafCutOff and GreaterThan:
                PassCHT=False
            if AFscore > MafCutOff and not GreaterThan:
                PassCHT=False
        
        if PassMAF and PassCHT:
            Outvcf.write(line)
            FiltCount=FiltCount+1

Outlog.write("  Number of variants in original VCF: "+str(OrigCount)+"\n")
Outlog.write("  Number of de novo variants: "+str(FiltCount)+"\n")
TimeNow=str(datetime.datetime.now().strftime("%Y-%m-%d %H:%M"))
Outlog.write("  Filtering Finished: "+TimeNow+"\n")
print "Done"

