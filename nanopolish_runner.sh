#! /bin/bash
#heavily inspired from https://github.com/nataliering/Resolving-the-complex-Bordetella-pertussis-genome-using-barcoded-nanopore-sequencing/blob/master/nanopolish_runner

set -x # echo commands and output, expands variables. 
#shopt -s expand_aliases #keep aliases from calling environment

logfile="./nanopolish_runner_errorlog.txt" # save sterr to logfile.
exec > $logfile 2>&1

if [ "$1" == "-h" ] ; then
    echo "Usage: run_nanopolish <directory_of_fast5s> <path/to/sequencing_summary.txt> <reads_from_albacore.fastq> <draft_assembly.fasta> <output_name.fasta>"
    exit 0
fi

if [ "$1" == "--help" ] ; then
    echo "Usage: run_nanopolish <directory_of_fast5s> <path/to/sequencing_summary.txt> <reads_from_albacore.fastq> <draft_assembly.fasta> <output_name.fasta>"
    exit 0
fi

#The -z switch will test if the expansion of "$1" is a null string or not. 
if [ -z "$1" ] ; then
	echo "Usage: run_nanopolish <directory_of_fast5s> <path/to/sequencing_summary.txt"
	echo "<reads_from_albacore.fastq> <draft_assembly.fasta> <output_name.fasta>"
	echo "No argument 1 supplied"
	echo "remember to load modules tools anaconda3/4.4.0 nanopolish/0.10.2 bwa/0.7.15"
	exit 0
fi

# load modules
module purge
module load tools anaconda3/4.4.0 nanopolish/0.10.2 \
		minimap2/2.6 samtools/1.9 parallel/20181222
module list

# set variables
NPROCS=`wc -l < $PBS_NODEFILE` # number of processes in environment

#set scratchdir
sp="/home/projects/cu_10128/scratch/$PBS_JOBID"
mkdir $sp
echo scratchdir is "$sp"

echo "using $NPROCS processes"
echo "Scratch temporary directory: $sp"
echo "Arguments:"
echo "<directory of fast5s>: "$1""
echo "<path/to/sequencing_summary.txt>: "$2""
echo "<reads_from_albacore.fastq>: "$3""
echo "<draft_assembly.fasta>: "$4""
echo "<output_name.fasta>: "$5"" 

dt=`date '+%d/%m/%Y %H:%M:%S'`; echo "$dt"
echo "Running nanopolish index"
#index has been done on the fasta.gz files so many times. 
#$NANOPOLISH/nanopolish index -d $1 -s "$2" "$3"

dt=`date '+%d/%m/%Y %H:%M:%S'`; echo "$dt"
echo "Running minimap2"
#removed -a before pipe
minimap2 -ax map-ont -t $NPROCS $4 $3 | \
	samtools sort -@ $NPROCS -o $sp/reads.sorted.bam -T $sp/reads.tmp

dt=`date '+%d/%m/%Y %H:%M:%S'`; echo "$dt"
echo "Running samtools index"
samtools index -@ $NPROCS $sp/reads.sorted.bam
#rm *.fai
dt=`date '+%d/%m/%Y %H:%M:%S'`
echo "$dt"

echo "Running parallel nanopolishing"
para_nr="4"
para_threads=$(( NPROCS / para_nr ))
python $NANOPOLISH/scripts/nanopolish_makerange.py $4 | \
	parallel --results $sp/nanopolish.results -P $para_nr \
	$NANOPOLISH/nanopolish variants --consensus \
			-o $sp/polished.{1}.vcf -w {1} \
			--reads $3 \
			--bam $sp/reads.sorted.bam \
			--genome $4 \
			--threads $para_threads \
			-q dcm,dam


dt=`date '+%d/%m/%Y %H:%M:%S'`; echo "$dt"
echo "Merging parallel fasta files"
nanopolish vcf2fasta -g "$4" $sp/polished.*.vcf > $5


dt=`date '+%d/%m/%Y %H:%M:%S'`; echo "$dt"
echo "nanopolish_runner.sh done"
echo "cleaning"
	rm $sp/*.vcf
	rm $sp/reads.sorted.bam
	rm $sp/reads.sorted.bam.bai
	rm -r $sp/nanopolish.results

echo "Output file is: $5"
