## Resource Directories
export EXOMPPLN="/ifs/scratch/c2b2/af_lab/ads2202/Exome_Seq/scripts/Exome_pipeline_scripts_GATKv3" # Directory containing pipeline shell scripts
export EXOMRES="/ifs/scratch/c2b2/af_lab/ads2202/Exome_Seq/resources" # Directory containing resources/references for pipeline

#jar files
GATKJAR="/ifs/scratch/c2b2/af_lab/ads2202/src/GenomeAnalysisTK_Current/GenomeAnalysisTK.jar" #Current GATK jar file
PICARD="/ifs/scratch/c2b2/af_lab/ads2202/src/picard-tools-1.101/" #directory containing Picard jar files

## References
export BUILD="b37" # shorthand for build
export DBSNP="$EXOMRES/hg19/dbsnp_137.hg19.vcf" # dbSNP vcf from GATK
export INDEL="$EXOMRES/hg19/Mills_and_1000G_gold_standard.indels.hg19.vcf" # Gold standard INDEL reference from GATK
export INDEL1KG="$EXOMRES/hg19/1000G_phase1.indels.hg19.vcf" # INDEL reference from 1000 genomes
export REF="$EXOMRES/hg19/ucsc.hg19.fasta" # UCSC genome assembly from GATK
export HAPMAP="$EXOMRES/hg19/hapmap_3.3.hg19.vcf" # hapmap vcf from GATK
export TGVCF="$EXOMRES/hg19/1000G_omni2.5.hg19.vcf" 
export ONEKG="$EXOMRES/hg19/1000G_phase1.snps.high_confidence.hg19.vcf" # 1000 genome SNPs vcf
export ANNHDB="/ifs/home/c2b2/af_lab/ads2202/scratch/src/annovar/humandb/" #Location of annovar databases

#GATK no-phone-home key
export ETKEY="$EXOMRES/ads2202_c2b2.columbia.edu.key"
