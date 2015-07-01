## Resource Directories
export EXOMPPLN="/ifs/scratch/c2b2/af_lab/ads2202/Exome_Seq/scripts/Exome_pipeline_scripts_GATKv3" # Directory containing pipeline shell scripts
export EXOMRES="/ifs/scratch/c2b2/af_lab/ads2202/Exome_Seq/resources" # Directory containing resources/references for pipeline

#jar files  and directories for software
GATKJAR="/ifs/scratch/c2b2/af_lab/ads2202/jar/GenomeAnalysisTK_Current/GenomeAnalysisTK.jar" #Current GATK jar file
PICARD="/ifs/scratch/c2b2/af_lab/ads2202/jar/picard-tools/" #directory containing Picard jar files
SNPEFF="/ifs/scratch/c2b2/af_lab/ads2202/jar/snpEff/snpEff.jar" # Current snpEff jar file

## References
export BUILD="b37" # shorthand for build
export DBSNP="$EXOMRES/references/b37/dbsnp_137.b37.vcf" # dbSNP vcf from GATK
export INDEL="$EXOMRES/references/b37/Mills_and_1000G_gold_standard.indels.b37.vcf" # Gold standard INDEL reference from GATK
export INDEL1KG="$EXOMRES/references/b37/1000G_phase1.indels.b37.vcf" # INDEL reference from 1000 genomes
export REF="$EXOMRES/references/yshen_hg19fasta/hg19.fasta" # UCSC hg19 build from Yufeng
export HAPMAP="$EXOMRES/references/b37/hapmap_3.3.b37.vcf" # hapmap vcf from GATK
export TGVCF="$EXOMRES/references/b37/1000G_omni2.5.b37.vcf" 
export ONEKG="$EXOMRES/references/b37/1000G_phase1.snps.high_confidence.b37.vcf" # 1000 genome SNPs vcf
export STHSH="$EXOMRES/references/b37/stampy_b37" # hash file for Stampy - omit ".sthash" extension for compatibility with Stampy
export STIDX="$EXOMRES/references/b37/stampy_b37" # genome index file for Stampy - omit ".stidx" extension for compatibility with Stampy
export ANNHDB="/ifs/home/c2b2/af_lab/ads2202/scratch/bin/annovar/humandb/" #Location of annovar databases

#GATK no-phone-home key
export ETKEY="$EXOMRES/ads2202_c2b2.columbia.edu.key"

#Capture Kit Target Files
export AgtV2="$EXOMRES/CaptureKitBeds/SureSelect_All_Exon_V2.b37.ordered.bed"
export AgtV4="$EXOMRES/CaptureKitBeds/SureSelect_All_Exon_V4_b37.ordered.bed"
export AgtV5="$EXOMRES/CaptureKitBeds/SureSelect_Human_All_Exon_V5_Covered.ordered.bed"
export AgtV5UTR="$EXOMRES/CaptureKitBeds/SureSelect_Human_All_Exon_V5_UTRs_Covered.ordered.bed"
export NbgV2="$EXOMRES/CaptureKitBeds/SeqCap_EZ_Exome_v2.b37.targets.bed"
export NbgV3="$EXOMRES/CaptureKitBeds/SeqCap_EZ_Exome_v3.b37.targets.bed"
export IllTS="$EXOMRES/CaptureKitBeds/truseq_exome_targeted_regions.b37.ordered.bed"
export BigTgt="$EXOMRES/CaptureKitBeds/custom_intervals/FullIntervals.bed"
export RefSeq="$EXOMRES/CaptureKitBeds/Target_RefSeq_Exon_UCSC.bed"
export VCRv2="$EXOMRES/CaptureKitBeds/Nimblegen_VCRome_v2.hg19.bed"
export TGTCODES="AgtV2:AgtV4:AgtV5:AgtV5UTR:NbgV2:NbgV3:IllTS:BigTgt:RefSeq:VCRv2"

