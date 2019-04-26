#!/bin/bash
# racon with illumina reads
#after racon_runner https://github.com/nataliering/Resolving-the-complex-Bordetella-pertussis-genome-using-barcoded-nanopore-sequencing/blob/master/racon_runner
#this is written to run on a HPC that manages software through Environment Modules (http://modules.sourceforge.net/).
#Other users may have to comment out the lines loading and unloading modules.

#Check if arguments have been stated.
if [ "$1" == "-h" ] ; then
    echo "Usage: racon_illumina_runner <assembly.fasta> <illuminareads1.fastq.gz> <reads2.fastq.gz> <out.fasta>"
    exit 0
fi

if [ "$1" == "--help" ] ; then
    echo "Usage: racon_illumina_runner <assembly.fasta> <illuminareads1.fastq> <reads2.fastq> <out.fasta>"
    exit 0
fi

#The -z switch will test if the expansion of "$1" is a null string or not. 
if [ -z "$1" ] ; then
	echo "Usage: racon_illumina_runner <assembly.fasta> <illuminareads1.fastq> <reads2.fastq> <out.fasta>"
	echo "no argument supplied"
	echo "remember to load modules tools anaconda3/4.4.0 minimap2/2.14r883 miniasm/0.3r179 racon/1.3.1 bwa/0.7.15"
	exit 0
fi

NPROCS=`wc -l < $PBS_NODEFILE` # number of processes in environment
dt=`date '+%d/%m/%Y %H:%M:%S'`; echo "$dt"
echo racon_runner.sh with input files:
echo 1: $1
echo 2: $2
echo 3: $3
echo 4: $4
echo "workning directory is $PBS_O_WORKDIR"

#load required modules
module purge
module load tools anaconda3/4.4.0 minimap2/2.6 racon/1.3.1
module list

#set path for prepilluminaracon.py
scriptpath=""

#set path for scratch / temp files (not backed up).
sp="$PBS_JOBID" 

# racon requires illumina fastq as one file
# Read data (as fastq.gz) into the named pipe and put the process in the background (&)
rm $sp/ill.temp.fastq
zcat $2 $3 > $sp/ill.temp.fastq

#fix reads names to R1 and R2 are not the same up to whitespace with prepilluminaracon.py
python $scriptpath/prepilluminaracon.py $sp/ill.temp.fastq \
	> $sp/ill.prep.temp.fastq

#map with minimap2 (instead of BWA requires same read names)
minimap2 -$NPROCS -ax sr $1 $sp/ill.prep.temp.fastq > $sp/$1.temp.sam

echo "Running racon"
racon --threads "$NPROCS" \
	--include-unpolished \
	 $sp/ill.prep.temp.fastq "$sp/$1.temp.sam" "$1" > "$4"

echo "cleaning"
rm "$sp/$1.temp.sam"
rm $sp/ill.sequence.temp.fastq

echo "Racon_illumina_runner.sh finished"
