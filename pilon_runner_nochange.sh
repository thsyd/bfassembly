#!/bin/bash
#Pilon polishing until no changes are made or 6 rounds have been run. 
#after https://github.com/nataliering/Resolving-the-complex-Bordetella-pertussis-genome-using-barcoded-nanopore-sequencing/blob/master/pilon_runner
#this is written to run on a HPC that manages software through Environment Modules (http://modules.sourceforge.net/).
#Other users may have to comment out the lines loading and unloading modules.

if [ "$1" == "-h" ] ; then
    echo "Usage: pilon_runner <in.fasta> <illumina_S1> <illumina_S2> <out_prefix>"
    exit 0
fi
if [ "$1" == "--help" ] ; then
    echo "Usage: pilon_runner <in_fasta> <illumina_S1> <illumina_S2> <out_prefix> "
    exit 0
fi

#The -z switch will test if the expansion of "$1" is a null string or not. 
if [ -z "$1" ] ; then
	echo "pilon_runner <in.fasta> <illumina_S1> <illumina_S2> <out_prefix>"
	echo "no argument supplied"
	echo "remember to load modules tools anaconda3/4.4.0 jre/1.8.0-openjdk pilon/1.22 bwa/0.7.15"
	exit 0
fi
echo "pilon_runner_nochange.sh with input:"
echo "1: $1"
echo "2: $2"
echo "3: $3"
echo "4: $4"

#set scratchdir
sp="$PBS_JOBID"

module purge
module load tools jre/1.8.0 pilon/1.22 bwa/0.7.15 samtools/1.9
module list

#set counts
base_1=$(basename $1)
ls $base_1.pilonchanges.txt
rm $base_1.pilonchanges.txt
touch $base_1.pilonchanges.txt
PILON_CHANG=1 
PILON_ROUND=0 

# set target fasta file
target=$1

#loop until 6 pilon rounds have passed or pilon changes = 0  
#to prevent continious piloning.
until [[ PILON_ROUND -eq "6" || $PILON_CHANG -eq "0" ]]; do

	PILON_ROUND=$((PILON_ROUND+1)) # increase by 1 
	dt=`date '+%d/%m/%Y %H:%M:%S'`; echo "$dt"
	echo "Pilon polishing round":$PILON_ROUND
	echo "target: $target"
	
	echo "bwa indexing"
	bwa index $target
	echo "bwa memming"
	bwa mem -t $NPROCS $target $2 $3 \
		| samtools sort -@ $NPROCS > $sp/pilon.aln.bam
	
	echo "samtools indexing"
	samtools index $sp/pilon.aln.bam
	echo "samtools faidxing"
	samtools faidx $target
	
	echo "piloning"
	java -jar /services/tools/pilon/1.22/pilon-1.22.jar \
		--genome $target \
		--frags $sp/pilon.aln.bam \
		--output $4.pilonround.$PILON_ROUND \
		--mindepth 0.5 --changes --threads $NPROCS

	#log changes in one file	
	echo "Pilon Polishing round":$PILON_ROUND >> \
		$base_1.pilonchanges.txt
	cat $4.pilonround.$PILON_ROUND.changes >> \
		$base_1.pilonchanges.txt	
	PILON_CHANG=$(cat $4.pilonround.$PILON_ROUND.changes | wc -l)
	echo "Pilon Changes:" $PILON_CHANG
	
	#make sure the next target is the product of this rounds pilon polish
	target=$4.pilonround.$PILON_ROUND.fasta
	echo "$PILON_ROUND" > pilon_round.txt #save round for later
done

#clean up
echo "cleaning after $PILON_ROUND round(s) of Pilon"
rm $sp/pilon.aln.bam
rm $sp/pilon.aln.bam.bai
rm *.sa
rm *.changes
rm *.amb
rm *.bwt
rm *.ann
rm *.pac
rm *.fai

echo "end"



