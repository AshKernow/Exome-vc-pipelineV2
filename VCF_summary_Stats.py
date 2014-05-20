from optparse import OptionParser
parser = OptionParser()
# Basic Input Files, required
parser.add_option("-v", "--vcf", dest="VCFfile",help="input VCF file", metavar="VCFfile")
parser.add_option("-o", "--output", dest="OutputFileName",help="user specified name of file for output to be written", metavar="OutputFileName")

(options, args) = parser.parse_args()

Output=open(options.OutputFileName,'w')

Nucleotides=['A','C','G','T']
MutUnknown=['unknown']
MutSilent=['synonymousSNV']
MutMissense=['nonsynonymousSNV']
MutNonsense=['stopgainSNV', 'stoploss SNV']
MutNonframeshift=['nonframeshiftinsertion', 'nonframeshiftdeletion']
MutFrameshift=['frameshiftinsertion', 'frameshiftdeletion']
CodingCodes=['splicing', 'exonic', 'exonic;splicing']
CodingTypes=['Coding variants', 'Non-coding Variants', 'All variants']

for WhichCode in CodingTypes:
	print WhichCode
	VCF=open(options.VCFfile,'r')
	# Read VCF file
	for line in VCF:
		# Map column name to number, and then find column numbers of each set of trios
		if '#CHROM' in line:
			Samples=line.split("\t")
			Samples=Samples[8:]
			Samples[0]='Total'
			SamLength=len(Samples)
			SNVcount=[0]*SamLength
			InDelCount=[0]*SamLength
			KnownCount=[0]*SamLength
			NovelCount=[0]*SamLength
			knownTiCount=[0]*SamLength
			knownTvCount=[0]*SamLength
			novelTiCount=[0]*SamLength
			novelTvCount=[0]*SamLength
			silentCount=[0]*SamLength
			missenseCount=[0]*SamLength
			nonsenseCount=[0]*SamLength
			unknownCount=[0]*SamLength
			frameShiftCount=[0]*SamLength
			homozygousCount=[0]*SamLength
			heterozygousCount=[0]*SamLength
			referenceCount=[0]*SamLength
			KGRareCount=[0]*SamLength
			ESPRareCount=[0]*SamLength
			RareCount=[0]*SamLength
			NoAAF=[0]*SamLength
			KnownTiTvRat=[0]*SamLength
			NovelTiTvRat=[0]*SamLength
			TotalTiTvRat=[0]*SamLength
			UnCalled=[0]*SamLength
		# Start Countin
		if '#' not in line:
			# Variant must first pass 1KG and GO-ESP frequencies, MQ0 threshold, and be exonic
			linelist=line.split("\t")
			INFOstring=linelist[7]
			INFOcolumnList=INFOstring.split(";")
			INFOdict={}
			for element in INFOcolumnList:
				if '=' in element:
					FieldName,FieldValue=element.split('=',1)
					INFOdict[FieldName]=FieldValue
			
			# Get values for later
			KGscore=float(INFOdict.get('1KGfreq',2))
			ESPscore=float(INFOdict.get('ESPfreq',2))
			MutationFunct=str(INFOdict.get('VarFunc','none'))
			MutationClass=str(INFOdict.get('VarClass','none'))
			ID=str(linelist[2])
			REF=linelist[3].split(",")
			REF=[ str(i) for i in REF ]
			ALT=linelist[4].split(",")
			ALT=[ str(i) for i in ALT ]
			#Known or novel
			Known=False
			if 'DB' in INFOcolumnList:
				Known=True
			if ID!=".":
				Known=True
			#Snp or Indel
			InDel=True
			if all(i in Nucleotides for i in REF) and all(i in Nucleotides for i in ALT):
				InDel=False
			#Ti or Tv
			Transversion=False
			Transition=False
			if (REF[0]=='A' and ALT[0]=='G') or (REF[0]=='G' and ALT[0]=='A') or (REF[0]=='C' and ALT[0]=='T') or (REF[0]=='T' and ALT[0]=='C'):
				Transition=True
			if (REF[0]=='C' and ALT[0]=='A') or (REF[0]=='A' and ALT[0]=='C') or (REF[0]=='G' and ALT[0]=='C') or (REF[0]=='C' and ALT[0]=='G') or (REF[0]=='T' and ALT[0]=='A') or (REF[0]=='A' and ALT[0]=='T') or (REF[0]=='T' and ALT[0]=='G') or (REF[0]=='G' and ALT[0]=='T'):
				Transversion=True
			codingPass=False
			if WhichCode=="Coding variants" and MutationFunct in CodingCodes:
				codingPass=True
			if WhichCode=="Non-coding Variants" and MutationFunct not in CodingCodes:
				codingPass=True
			if WhichCode=="All variants":
				codingPass=True
			
			if codingPass:
				for i in range(0, SamLength):
					ColNum=i+8
					if i!=0:
						InfoList=linelist[ColNum].strip()
						if "./." not in InfoList:
							InfoListSPL=InfoList.split(":")
							GT=InfoListSPL[0]
							if str(GT)!= "0/0":
								if str(GT)=="0/1":
									heterozygousCount[i]=heterozygousCount[i]+1
								if str(GT)=="1/1":
									homozygousCount[i]=homozygousCount[i]+1
						else:
							UnCalled[i]=UnCalled[i]+1
					else:
						GT="./."
					if str(GT)=="0/1" or str(GT)=="1/1" or i==0:
						if InDel:
								InDelCount[i]=InDelCount[i]+1
						else:
							SNVcount[i]=SNVcount[i]+1
						if Known:
							KnownCount[i]=KnownCount[i]+1
							if Transition:
								knownTiCount[i]=knownTiCount[i]+1
							if Transversion:
								knownTvCount[i]=knownTvCount[i]+1
						else:
							NovelCount[i]=NovelCount[i]+1
							if Transition:
								novelTiCount[i]=novelTiCount[i]+1
							if Transversion:
								novelTvCount[i]=novelTvCount[i]+1
						if MutationClass in MutSilent:
							silentCount[i]=silentCount[i]+1
						if MutationClass in MutMissense:
							missenseCount[i]=missenseCount[i]+1
						if MutationClass in MutNonsense:
							nonsenseCount[i]=nonsenseCount[i]+1
						if MutationClass in MutUnknown:
							unknownCount[i]=unknownCount[i]+1
						if MutationClass in MutFrameshift:
							frameShiftCount[i]=frameShiftCount[i]+1
						if KGscore<=0.01:
							KGRareCount[i]=KGRareCount[i]+1
						if ESPscore<=0.01:
							ESPRareCount[i]=ESPRareCount[i]+1
						if KGscore<=0.01 or ESPscore<=0.01:
							RareCount[i]=RareCount[i]+1
						if KGscore==2 and ESPscore==2:
							NoAAF[i]=NoAAF[i]+1
	for i in range(0, SamLength):
		if int(knownTvCount[i])>0:
			KnownTiTvRat[i]=float(knownTiCount[i])/float(knownTvCount[i])
		if int(novelTvCount[i])>0:
			NovelTiTvRat[i]=float(novelTiCount[i])/float(novelTvCount[i])
		TotalTvCount=float(novelTvCount[i])+float(knownTvCount[i])
		TotalTiCount=float(novelTiCount[i])+float(knownTiCount[i])
		if TotalTvCount>0:
			TotalTiTvRat[i]=float(TotalTiCount)/float(TotalTvCount)
	#write output
	Output.write("\t".join([str(WhichCode)])+"\n")
	Output.write("\t".join(['Samples:']+[ str(i) for i in Samples ]))
	Output.write("\t".join(['SNVs']+[ str(i) for i in SNVcount ])+"\n")
	Output.write("\t".join(['InDels']+[ str(i) for i in InDelCount ])+"\n")
	Output.write("\t".join(['Known']+[ str(i) for i in KnownCount ])+"\n")
	Output.write("\t".join(['Novel']+[ str(i) for i in NovelCount ])+"\n")
	Output.write("\t".join(['known Ti']+[ str(i) for i in knownTiCount ])+"\n")
	Output.write("\t".join(['known TV']+[ str(i) for i in knownTvCount ])+"\n")
	Output.write("\t".join(['known Ti/TV ratio']+[ str(i) for i in KnownTiTvRat ])+"\n")
	Output.write("\t".join(['novel Ti']+[ str(i) for i in novelTiCount ])+"\n")
	Output.write("\t".join(['novel TV']+[ str(i) for i in novelTvCount ])+"\n")
	Output.write("\t".join(['novel Ti/TV ratio']+[ str(i) for i in NovelTiTvRat ])+"\n")
	Output.write("\t".join(['total Ti/TV ratio']+[ str(i) for i in TotalTiTvRat ])+"\n")
	if str(CodingTypes)!='Non-coding Variants':
		Output.write("\t".join(['Silent']+[ str(i) for i in silentCount ])+"\n")
		Output.write("\t".join(['Missense']+[ str(i) for i in missenseCount ])+"\n")
		Output.write("\t".join(['Nonsense']+[ str(i) for i in nonsenseCount ])+"\n")
		Output.write("\t".join(['Unknown']+[ str(i) for i in unknownCount ])+"\n")
		Output.write("\t".join(['Frameshift']+[ str(i) for i in frameShiftCount ])+"\n")
	Output.write("\t".join(['Homozygous']+[ str(i) for i in homozygousCount ])+"\n")
	Output.write("\t".join(['Heterozygous']+[ str(i) for i in heterozygousCount ])+"\n")
	Output.write("\t".join(['Rare (AAF<0.01) - 1KG']+[ str(i) for i in KGRareCount ])+"\n")
	Output.write("\t".join(['Rare (AAF<0.01) - ESP']+[ str(i) for i in ESPRareCount ])+"\n")
	Output.write("\t".join(['Rare (AAF<0.01) - 1KG or ESP']+[ str(i) for i in RareCount ])+"\n")
	Output.write("\t".join(['No AAF']+[ str(i) for i in NoAAF ])+"\n")
	Output.write("\t".join(['not called']+[ str(i) for i in UnCalled ])+"\n")
	Output.write("\t\n")


