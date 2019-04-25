# Complete genome assembly of multidrug resistant *Bacteroides fragilis* reveals diverse plasmid and chromosomal context

This repository supports our paper, which is currenctly available as a preprint from BioRXiv: 

Content includes:
- code and commands for calling community-built and published tools
- links to data repositories, including raw and processed reads (fast5 and fastq) as well as genome assemblies
- supplementary results and links to supplementary files

## Abstract
*Bacteroides fragilis* is an opportunistic pathogen that resides in the human gut as part of the human commensal flora. Short DNA sequencing alone is unable to resolve the structure of the genome. VSeven complete genome assemblies have been deposited to public databases to date (2019-02-2018). To add to the public databases and to explore the location of antimicrobial resistance genes in the context of the complex genomes of B. fragilis, we aimed to complete genome assembly using a hybrid assembly approach. Using the Sanger sequenced *Bacteroides fragilis* NCTC 9343(NCBI RefSeq accession GCF_000025985.1) as a comparative reference, we sequenced *B. fragilis* CCUG4856T (=NCTC9343=ATCC25285) with Illumina 150 bp PE and Oxford Nanopore (ONT) MinION Rapid barcoding kit and assembled using various assemblers and polishing tools. The best assembly was with hybrid Unicycler with FilterByTile and Cutadapt filteres Illumina reads and FiltLong filtered and Canu corrected ONT reads. This approach was then applied to six multidrug resistant isolates. Resistance genes and IS elements were identified. Plasmids where charaterised. 

## Commands for tools mentioned in the manuscript

### Preparing Illumina reads

### Preparing Nanopore reads

#### Filtlong
Filtlong 0.2.0
```
filtlong --min_length 100 --keep_percent 90 --target_bases 500000000 input.fastq.gz | gzip > output.fastq.gz
```

### Running assemblies

#### Canu correction
Canu 1.8
run with standard options or
corOutCoverage=999 or
corMinCoverage=0
```
canuprefix=$(basename output.fastq.gz .fastq.gz)
canu -correct \
	-p $canuprefix \
	-d <outdir> \
	genomeSize=5.4m \
	-nanopore-raw output.fastq.gz \
	useGrid=false
```

### Unicycler hybrid assembly
Unicycler 0.4.7
Dependencies: jre/1.7.0 pilon/1.22 racon/1.3.1 spades/3.13.0 samtools/1.9
```
unicycler \
--short1 <path to short reads 1> \
--short2 <path to short reads 2> \
--long <path to long reads> \
--out <output path> \
--threads=<number of threads> \
--verbosity=2 \
--min_fasta_length=100 \
--mode normal
```

### Flye ONT read assembly
Flye 2.3.7
```
flye --nano-raw <long reads> \
	--genome-size 5.3m \
	--threads <number of threads> \
	--out-dir <outdir>
```

### Flye ONT read assembly with Flyepolish

Flye 2.3.7 run with three iterations of Flyes internal polishing engine
```
flye --nano-raw <long reads> \
	--genome-size 5.3m \
	--threads <number of threads> \
	--out-dir <outdir> \
	-i 3
```

### Minimap2 - miniasm assembly
minimap2/2.6 miniasm/0.3r179
```
minimap2 -x ava-ont -t $NPROCS ont_reads.fastq ont_reads.fastq | gzip -1 \
	> output.reads.paf.gz
miniasm -f <ONT reads> output.reads.paf.gz > output.miniasm.gfa

# convert to fasta after lh3 https://www.biostars.org/p/169516/
awk '/^S/{print">"$2"\n"$3}' output.miniasm.gfa | \
	fold > output.miniasm.fasta
```

### Skesa assembly
skesa/2.3.0
```
skesa 	--fastq <short reads 1> \
		--fastq <short reads 2> \
		--cores <number of cores> \
		--memory <allocated memory> \
		--contigs_out output.skesa.fasta
```

### SPAdes assembly
spades/3.13.0
```
spades.py 	--threads <number of cores> \
		--memory <allocated memory> \
		-k 31,51,71 \
		--pe1-1 <short reads 1> \
		--pe1-2 <short reads 2> \
		-o output.spades \
		--careful \
		--tmp-dir <scratch directory>
```
### HybridSPAdes assembly
spades/3.13.0
```
spades.py 	--threads <number of cores> \
		--memory <allocated memory> \
		-k 31,51,71 \
		--pe1-1 <short reads 1> \
		--pe1-2 <short reads 2> \
		--nanopore <nanopore reads> \
		-o output.hybridspades \
		--careful \
		--tmp-dir <scratch directory>
```
