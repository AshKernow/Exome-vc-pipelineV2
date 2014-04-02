#Library of functions used throughout Exome analysis scriptd

#-------------------------------------------------------------------------------------------------------
#Function to get input file name from a list of files in an array job
funcFilfromList() {
ChecList=${InpFil##*.}
if [[ $ChecList == "list" ]];then
	echo $ChecList
	InpFil=$(head -n $SGE_TASK_ID $InpFil | tail -n 1)
fi
}
#-------------------------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------------------------
#Function to enter information about the script initiation into the log
funcWriteStartLog () {
uname -a >> $TmpLog
echo "Start "$ProcessName" - $0:`date`" >> $TmpLog
echo " Job name: "$JOB_NAME >> $TmpLog
echo " Job ID: "$JOB_ID >> $TmpLog
if [[ ! -z $SGE_TASK_ID ]]; then echo " Array task ID: "$SGE_TASK_ID >> $TmpLog; fi
echo " Input File: "$InpFil >> $TmpLog
if [[ ! -z $BamFil ]]; then echo " Bam File: "$BamFil >> $TmpLog; fi
if [[ ! -z $BamNam ]]; then echo " Base name for outputs: $BamNam" >> $TmpLog; fi
if [[ ! -z $Chr ]]; then echo " Chromosome: "$Chr >> $TmpLog; fi
echo "----------------------------------------------------------------" >> $TmpLog
}
#-------------------------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------------------------
#function to log the start of each step wtihin a script
funcLogStepStart () { echo "- Start $StepName `date`...">> $TmpLog ; } 
#-------------------------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------------------------
#func to trim GATK output log and write it to the temp log
funcTrimGATKlog (){
if [[ $StepCmd == *$GATKJAR* ]]; then
	echo "  --- GATK output log for $StepName ----------------" >> $TmpLog
	grep -vE "ProgressMeter - *[dc0-9XY]|Copyright|INITIALIZATION COMPLETE|----|For support and documentation|Done preparing for traversal" $GatkLog | awk '{ print "\t\t"$0 }' >> $TmpLog
	echo "  --- --- --- --- --- --- ---" >> $TmpLog
	rm $GatkLog
fi
}
#-------------------------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------------------------
#function checks that the step has completed successfully, if not it writes and error message to the log and exits the script, otherwise it logs the completion of the step
funcLogStepFinit () { 
if [[ $? == 1 ]]; then #check exit status and if error then...
	echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" >> $TmpLog
	echo "     $StepName failed `date`" >> $TmpLog
	qstat -j $JOB_ID | grep -E "usage *$SGE_TASK_ID:" >> $TmpLog #get cluster usage stats
	echo "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" >> $TmpLog
	echo "=================================================================" >> $TmpLog
	funcTrimGATKlog
	grep "ERROR MESSAGE" $SGE_STDERR_PATH | awk '{ print "\t\t"$0 }' >> $TmpLog
	cat $TmpLog >> $LogFil
	rm $TmpLog
	exit 1
fi
funcTrimGATKlog
echo "- End $StepName `date`...">> $TmpLog # if no error log the completion of the step
echo "-----------------------------------------------------------------------" >> $TmpLog
}
#-------------------------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------------------------
#function to run and log the initiation/completion/failure of each step in the script
funcRunStep (){
funcLogStepStart
echo $StepCmd >> $TmpLog
eval $StepCmd
funcLogStepFinit
}
#-------------------------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------------------------
#function to log the end of each script and transfer the contents of temporary log file to the main log file
funcWriteEndLog () {
echo "End "$ProcessName" $0:`date`" >> $TmpLog
qstat -j $JOB_ID | grep -E "usage *$SGE_TASK_ID:" >> $TmpLog #get cluster usage stats
echo "===========================================================================================" >> $TmpLog
echo "" >> $TmpLog
cat $TmpLog >> $LogFil
rm -r $TmpLog $TmpDir
}
#-------------------------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------------------------
# function for common additional arguments to GATK
funcGatkAddArguments (){
if [[ $AllowMisencoded == "true" ]]; then StepCmd=$StepCmd" -allowPotentiallyMisencodedQuals"; fi
if [[ $FixMisencoded == "true" ]]; then StepCmd=$StepCmd" -fixMisencodedQuals"; fi
if [[ $BadET == "true" ]]; then StepCmd=$StepCmd" -et NO_ET -K $ETKEY"; fi
}

#-------------------------------------------------------------------------------------------------------

#-------------------------------------------------------------------------------------------------------
# function for calling next step in pipeline
funcPipeLine (){
if [[ $PipeLine == "true" ]]; then
	mkdir -p stdostde
	echo "- Call $NextJob `date`:" >> $TmpLog
	echo "    "$QsubCmd  >> $TmpLog
	eval $QsubCmd >> $TmpLog
	echo "----------------------------------------------------------------" >> $TmpLog
fi
}
#-------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------
#hpc workarounds
if [[ /bin/hostname==*.hpc ]]; then 
source /etc/profile.d/sge.sh  # SGE commands from within node
source /ifs/home/c2b2/af_lab/ads2202/.bash_profile
fi
