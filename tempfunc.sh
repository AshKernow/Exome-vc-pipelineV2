funcRunStep (){ #run each command
funcLogStepStart
echo $StepCmd >> $TmpLog
eval $StepCmd
funcLogStepFinit
}

#get ReadGroupHeader from input BAM
RgHeader=$(samtools view -H 1303362.bam | grep ^@RG | cut -f 2-)
echo "ReadGroup header: $RgHeader" >> $TmpLog
if [[ $RgHeader == "" ]]||[[ $(echo "$RgHeader" | wc -l) -gt 1 ]]; then #check that we have a  RG header and if not write a warning to the log file
	echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" >> $TmpLog
	echo "     Problem with ReadGroup header" >> $TmpLog
	echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" >> $TmpLog
fi

###Align using BWA mem algorithm
StepName="Align with BWA mem"
StepCmd='samtools bam2fq $BamFil | bwa mem -M -R \"$RgHeader\" -t 6 -p $REF - | samtools view -bS - > $AlnFil'
funcRunStep

#Mark the duplicates
StepName="Mark PCR Duplicates" >> $TmpLog
StepCmd="$JAVA7BIN -Xmx4G -Djava.io.tmpdir=$TmpDir -jar $PICARD/MarkDuplicates.jar INPUT=$AlnFil OUTPUT=$AlnFil METRICS_FILE=$BamNam.Dup.metrics.txt CREATE_INDEX=TRUE"
funcRunStep

#flag stats
StepName="Output flag stats"
funcLogStepStart
cmd="$SAMTOOLS flagstat $BamFil.bam > $BamFil.flagstat"
echo "    "$cmd >> $TmpLog
eval $cmd
funcLogStepFinit

#Index stats
StepName="Output idx stats"
funcLogStepStart
cmd="$SAMTOOLS idxstats $BamFil.bam > $BamFil.idxstats"
echo "    "$cmd >> $TmpLog
eval $cmd
funcLogStepFinit
echo "----------------------------------------------------------------" >> $TmpLog
