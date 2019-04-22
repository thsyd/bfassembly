# Complete genome assembly of multidrug resistant *Bacteroides fragilis* reveals diverse plasmid and chromosomal context

This repository supports our paper, which is currenctly available as a preprint from BioRXiv: 

Content includes:
- code and commands for calling community-built and published tools
- links to data repositories, including raw and processed reads (fast5 and fastq) as well as genome assemblies
- supplementary results and links to supplementary files

## Abstract
*Bacteroides fragilis* is an opportunistic pathogen that resides in the human gut as part of the human commensal flora. Short DNA sequencing alone is unable to resolve the structure of the genome. VSeven complete genome assemblies have been deposited to public databases to date (2019-02-2018). To add to the public databases and to explore the location of antimicrobial resistance genes in the context of the complex genomes of B. fragilis, we aimed to complete genome assembly using a hybrid assembly approach. Using the Sanger sequenced *Bacteroides fragilis* NCTC 9343(NCBI RefSeq accession GCF_000025985.1) as a comparative reference, we sequenced *B. fragilis* CCUG4856T (=NCTC9343=ATCC25285) with Illumina 150 bp PE and Oxford Nanopore (ONT) MinION Rapid barcoding kit and assembled using various assemblers and polishing tools. The best assembly was with hybrid Unicycler with FilterByTile and Cutadapt filteres Illumina reads and FiltLong filtered and Canu corrected ONT reads. This approach was then applied to six multidrug resistant isolates. Resistance genes and IS elements were identified. Plasmids where charaterised. 

## Commands for tools mentioned in the manuscript.

### Preparing Illumina reads

### Filtlong
```
filtlong --min_length 100 --keep_percent 90 --target_bases 500000000 input.fastq.gz | gzip > output.fastq.gz
```
### Canu correction


canuprefix=$(basename output.fastq.gz .fastq.gz)
canu -correct \
	-p $canuprefix \
	-d <outdir> \
	genomeSize=5.2m \
	-nanopore-raw output.fastq.gz \
	useGrid=false
  
### Unicycler hybrid assembly



```

```
