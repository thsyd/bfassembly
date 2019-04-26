#!/bin/bash
#after racon_runner https://github.com/nataliering/Resolving-the-complex-Bordetella-pertussis-genome-using-barcoded-nanopore-sequencing/blob/master/racon_runner
#for debugging
#set -x # echo commands and output, expands variables. 
#shopt -s expand_aliases #keep aliases from calling environment
#exec 2> ./racon_runner_errorlog.txt # save sterr to logfile.

if [ "$1" == "-h" ] ; then
    echo "Usage: racon_runner <assembly.fasta> <reads.fastq> <out.fasta>"
    exit 0
fi

if [ "$1" == "--help" ] ; then
    echo "Usage: racon_runner <assembly.fasta> <ont_reads.fastq> <out.fasta>"
    exit 0
fi

#The -z switch will test if the expansion of "$1" is a null string or not. 
if [ -z "$1" ] ; then
	echo "Usage: racon_runner <assembly.fasta> <reads.fastq> <out.fasta>"
	echo "no argument supplied"
	echo "remember to load modules tools anaconda3/4.4.0 minimap2/2.14r883 miniasm/0.3r179 racon/1.3.1"
	exit 0
fi

module purge
module load tools anaconda3/4.4.0 minimap2/2.14r883 miniasm/0.3r179 racon/1.3.1
module list

#set scratchdir
sp="/home/projects/cu_10128/scratch/$PBS_JOBID"

echo racon_runner.sh with input files:
echo 1: $1
echo 2: $2
echo 3: $3
echo "" 
module list
dt=`date '+%d/%m/%Y %H:%M:%S'`; echo "$dt"
echo "Running minimap"
#minimap2 uses -t INT+1 thread therefore -t 7 (as 8 are used).
#$NPROCS = number of processes allocated on computerome
NPROCS=`wc -l < $PBS_NODEFILE` # number of processes in environment
mmthread=$(( NPROCS -1 ))
minimap2 -x ava-ont -t $mmthread $1 $2 > $sp/$1.temp.paf

echo "Running racon"
racon --threads $NPROCS --include-unpolished $2 $sp/$1.temp.paf $1 > $3

echo "cleaning"
rm $sp/$1.temp.paf

echo "end"
