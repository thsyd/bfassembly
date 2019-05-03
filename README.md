# Complete genome assembly of clinical multidrug resistant Bacteroides fragilis isolates enables comprehensive identification of antimicrobial resistance genes and plasmids.

This repository supports our paper, which is currenctly available as a preprint from BioRXiv: 

Content includes:
- Code and commands for calling community-built and published tools
- Links to data repositories, including raw and processed reads (fast5 and fastq) as well as genome assemblies
- Supplementary results and links to supplementary files

## Abstract
*Bacteroides fragilis* is an opportunistic pathogen that resides in the human gut as part of the human commensal flora. Short DNA sequencing alone is unable to resolve the structure of the genome. VSeven complete genome assemblies have been deposited to public databases to date (2019-02-2018). To add to the public databases and to explore the location of antimicrobial resistance genes in the context of the complex genomes of B. fragilis, we aimed to complete genome assembly using a hybrid assembly approach. Using the Sanger sequenced *Bacteroides fragilis* NCTC 9343(NCBI RefSeq accession GCF_000025985.1) as a comparative reference, we sequenced *B. fragilis* CCUG4856T (=NCTC9343=ATCC25285) with Illumina 150 bp PE and Oxford Nanopore (ONT) MinION Rapid barcoding kit and assembled using various assemblers and polishing tools. The best assembly was with hybrid Unicycler with FilterByTile and Cutadapt filteres Illumina reads and FiltLong filtered and Canu corrected ONT reads. This approach was then applied to six multidrug resistant isolates. Resistance genes and IS elements were identified. Plasmids where charaterised. 

## Commands for tools mentioned in the manuscript

### Preparing Illumina reads

#### filterbytile from [BBMap](https://sourceforge.net/projects/bbmap/)
bbmap/36.49
```
filterbytile.sh in1=<short reads 1> in2=<short reads 2> out1=filtered1.fq.gz out2=filtered2.fq.gz
```

#### [TrimGalore](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)
trim_galore/0.4.4
```
trim_galore --paired --quality 10 <input 1> <input 2>
```

#### Subsample to <100x read depth
seqtk/1.0-r82-dirty pigz/2.3.4
```
#set estimated genome length
GLEN=5.3
#set target max depth
DEPTH=100
BASES=$(zcat val_1.fq.gz | seqtk seq -A | grep -v "^>" | tr -dc "ACTGNactgn" | wc -m)
BASES=$(expr $BASES \* 2)	#total bases is double R1. need to escape the asterix

ORI_DEPTH=$(expr $BASES / $GLEN) #calculate original (crude) depth
echo DEPTH: $DEPTH
calc() { awk "BEGIN{print $*}"; } #need this short awk funktion as bash expr only calculates integer
FACTOR=$(calc $DEPTH / $ORI_DEPTH)

if [ "$ORI_DEPTH" -gt "$DEPTH" ]; then 		#subsample
	seqtk sample -s100 val_1.fq.gz $FACTOR | pigz --fast -c -p $NPROCS > R1.sub.fq.gz 
	seqtk sample -s100 val_2.fq.gz $FACTOR | pigz --fast -c -p $NPROCS > R2.sub.fq.gz
else
fi

```
### Preparing Nanopore reads
Nanopore reads were demultiplexed following the [notes by Ryan Wick](https://github.com/rrwick/Deepbinner/wiki/Using-Deepbinner-with-Albacore).
Raw reads were demultiplexed with [Deepbinner](https://github.com/rrwick/Deepbinner) and then basescalled with Albacore, retaining only reads where the two agreed on the barcode. 
[PoreChop](https://github.com/rrwick/Porechop) was then used to trim adapters and barcodes from reads.

Running PoreChop
```
porechop_out="$b"_trimmed.fastq.gz
porechop_log="$b"_porechop.log
porechop 	--threads $NPROCS \
		--check_reads 100 \
		--discard_middle \
		-i demultiplexed_fastqs/"$b"_untrimmed.fastq.gz \
			2> "$porechop_log" | pigz -p $NPROCS > "$porechop_out"
```
#### Concatenating fastq files from two runs
```
inp1=input1.fastq.gz
inp2=input2.fastq.gz
NPROCS=<number of processes/threads>
output=output.fastq.gz
zcat $inp1 $inp2 | pigz --fast -c -p $NPROCS > $output
```
#### [Filtlong](https://github.com/rrwick/Filtlong)
Filtlong 0.2.0
```
filtlong --min_length 100 --keep_percent 90 --target_bases 500000000 input.fastq.gz | gzip > output.fastq.gz
```

### Running assemblies

#### [Canu](https://canu.readthedocs.io/en/latest/index.html) correction of ONT reads
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
Canu assembly
```
canu \
	-p canu.ont \
	-d canu_ont \
	genomeSize=5.2m \
	useGrid=false \
	-nanopore-raw <path to nanopore reads>
 	maxMemory=<set memory> \
	maxThreads=<set threads>
```
#### [Unicycler](https://github.com/rrwick/Unicycler) hybrid assembly
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

#### [Flye](https://github.com/fenderglass/Flye) ONT read assembly
Flye 2.3.7
```
flye --nano-raw <long reads> \
	--genome-size 5.3m \
	--threads <number of threads> \
	--out-dir <outdir>
```

Flye ONT read assembly with Flyepolish

Flye 2.3.7 run with three iterations of Flyes internal polishing engine
```
flye --nano-raw <long reads> \
	--genome-size 5.3m \
	--threads <number of threads> \
	--out-dir <outdir> \
	-i 3
```

#### [miniasm](https://github.com/lh3/miniasm) assembly
minimap2/2.6 miniasm/0.3r179
```
minimap2 -x ava-ont -t $NPROCS ont_reads.fastq ont_reads.fastq | gzip -1 \
	> output.reads.paf.gz
miniasm -f <ONT reads> output.reads.paf.gz > output.miniasm.gfa

# convert to fasta after lh3 https://www.biostars.org/p/169516/
awk '/^S/{print">"$2"\n"$3}' output.miniasm.gfa | \
	fold > output.miniasm.fasta
```

#### [Skesa](https://github.com/ncbi/SKESA) assembly
skesa/2.3.0
```
skesa 	--fastq <short reads 1> \
		--fastq <short reads 2> \
		--cores <number of cores> \
		--memory <allocated memory> \
		--contigs_out output.skesa.fasta
```

#### [SPAdes](https://github.com/ablab/spades) assembly
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
#### [HybridSPAdes](https://github.com/ablab/spades) assembly
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

### Assembly polishing
Assemblies were polished with [Nanopolish](https://nanopolish.readthedocs.io/en/latest/), [Racon](https://github.com/isovic/racon), and/or [Pilon](https://github.com/broadinstitute/pilon).
[nanopolish_runner.sh](../master/nanopolish_runner.sh), [../master/pilon_runner_nochange.sh](../master/pilon_runner_nochange.sh), [racon_runner.sh](../master/racon_runner.sh), and [../master/racon_illumina_runner.sh](racon_illumina_runner.sh) were heavily inspired by @nataliering 's [scripts](https://github.com/nataliering/Resolving-the-complex-Bordetella-pertussis-genome-using-barcoded-nanopore-sequencing/).

Usage:
Nanopolish_runner.sh
```
run_nanopolish <directory_of_fast5s> <path/to/sequencing_summary.txt> <reads_from_albacore.fastq> <draft_assembly.fasta> <output_name.fasta>
```
pilon_runner_nochange.sh
```
pilon_runner <in.fasta> <illumina_S1> <illumina_S2> <out_prefix>
```
racon_runner.sh
```
racon_runner <assembly.fasta> <reads.fastq> <out.fasta>
```
racon_illumina_runner.sh
```
racon_illumina_runner <assembly.fasta> <illuminareads1.fastq.gz> <reads2.fastq.gz> <out.fasta>
```

### Evaluating assemblies
Assemblies is .fasta format were collected in the same folder, and the following scripts run in that folder.
This produced a number of output directories and .tsv files, with concatenated results for all .fasta files.
These .tsv files could then be merged, creating the supplementary file XXX with evaluation results.

#### [ALE](https://github.com/sc932/ALE)
The following runs ALE for all *fasta files in a folder, and extracts the score from the individual ALE results files and writes the scores to a seperate .tsv file.
```
mkdir ale
#create file for aggregation of results
printf "Assembly\total_ALE_score\n" > ale/totalALE.tsv
#analyse all fasta
for F in *.fasta; do
	
	N=$(basename $F .fasta )
	
	bwa index $F
	bwa mem -t $NPROCS $F $S1 $S2 \
		| samtools sort -@ $NPROCS > $sp/$N.ale.aln.bam

	ALE --nout $sp/$N.ale.aln.bam $F ale/$N.output.ale.meta.txt

	#cleanup
	rm $sp/$N.ale.aln.bam
	rm $F.amb $F.ann $F.bwt $F.pac $F.sa
	rm ale/$N.output.ale.meta.txt.param

	#extract info and print in results file
	totalALE=$(head -n 1 ale/$N.output.ale.meta.txt | tr ' ' \\t | cut -f3)
	printf "$F\t$totalALE\n" >> ale/totalALE.tsv
done
```

#### Average nucleotide identity was calculated using @chjp's [ANI.pl](https://github.com/chjp/ANI/blob/master/ANI.pl)
Reference genome was RefSeq GCF_000025985.1.
```
ref_genome=GCF_000025985.1_ASM2598v1_genomic.fna
allfasta=$(ls *.fasta)

printf "Assembly\tANI_pl\n" > ani.pl.results.txt
for F in $allfasta; do
	
	N=$(basename $F .fasta | sed s/\[.]/_/g)
	assembly=$(basename $F .fasta)
	perl ANI.pl \
		-bl <path to blastall> \
		-fd <path to formatdb> \
		-qr $F \
		-sb $ref_genome \
		-od "$N"_"ani" \
		>> $assembly.ani.pl.txt
		
		#concatenate ani.pl results
		printf "$assembly\t" >> ani.pl.results.txt
		grep -w -F -e ANI $assembly.ani.pl.txt | sed 's/.*\ //' >> "ani.pl.results.txt"
		#printf "\n" >> "ani.pl.results.txt" #need newline after

		#cleaning
		rm -r "$N"_"ani"
		rm $assembly.ani.pl.txt
done
```

#### [BUSCO](https://busco.ezlab.org/)
```
allfasta=$(ls *.fasta)
mkdir busco_out
cd busco_out
for F in $allfasta; do
	
	echo "BUSCO file $F"
	N=$(basename $F .fasta | sed s/\[.]/_/g)
	python run_BUSCO.py \
	--in ../$F \
	--out $N \
	--lineage_path <path to bacteroidetes_odb9> \
	--mode genome \
	--cpu $NPROCS \
	--tmp_path $sp/busco_tmp \
	--force

	#cleaning
	rm -r run_$N/augustus_output
	rm -r run_$N/blast_output
	rm -r run_$N/checkpoint.tmp
	rm -r run_$N/hmmer_output

done
```

#### [CheckM](https://ecogenomics.github.io/CheckM/)
```
binfolder="./"
outputfolder=checkm_out
markerfile=checkmmarkerfile

checkm tree $binfolder $outputfolder -x fasta --threads $NPROCS
checkm tree_qa $outputfolder 
checkm lineage_set $outputfolder $markerfile
checkm analyze $markerfile $binfolder $outputfolder -x fasta --threads $NPROCS
checkm qa $markerfile $outputfolder --out_format 2 \
	--threads $NPROCS --file $outputfolder/checkm_qa_out.tsv --tab_table

#cleaning
rm -r checkm_out/bins
rm -r checkm_out/storage
```

#### dnadiff from [MUMmer](http://mummer.sourceforge.net/)

```
for F in *.fasta; do
	N=$(basename $F .fasta | sed s/\[.]/_/g)
	dnadiff --prefix dnadiff/$N.dnadiff $ref_genome $F
	#get 1-to-1 AvgIdentity for QRY and print to file
	avgid_out=dnadiff/mummerdnadiff_results.txt
	printf "$F""\t">> $avgid_out
	grep -w -F -e AvgIdentity dnadiff/$N.dnadiff.report | \
	head -n 1 | tail -c 7 | tr -d '[:space:]' >> "$avgid_out"
	printf "\n" >> "$avgid_out" #need newline after
	
	#cleaning
	rm dnadiff/$N.dnadiff*
done
```

#### [Prokka](https://github.com/tseemann/prokka) annotation
Prokka is run for all assemblies in .fasta format.
The reference proteins of GCF_000025985.1 is used as protein database for prokka.
Counts for CDS, gene, RNA ect is then collected into one .tsv file.
```
ref_features=GCF_000025985.1_ASM2598v1_genomic.gbff
allfasta=$(ls *.fasta)

for F in $allfasta; do
	echo "PROKKA file $F"
	N=$(basename $F .fasta | sed s/\[.]/_/g)
	prokprefix=$(basename $F .fasta)
	prokka --compliant \
		--outdir prokka/$N --force \
		--prefix $prokprefix \
		--proteins $ref_genomic_gbff \
		--cpus $NPROCS \
		$F

	#cleaning
	rm prokka/$N/*.err
	rm prokka/$N/*.ffn
	rm prokka/$N/*.sqn
	rm prokka/$N/*.fsa

done

#Gather results in prokkaresults.tsv

#prepare column headers for prokkaresults.tsv
printf 	"Prokka bases\tProkka CDS\tProkka contigs\tProkka CRISPR\tProkka genes\tProkka isolate\tProkka rRNA\tProkka tmRNA\tProkka tRNA\n" > prokkaresults.tsv

for F in */*.txt; do

	F_base=$(basename $F .txt)
	# "isolate" and the basename to first line in tempfile.tsv
	printf "isolate\t$F_base\n" > tempfile.tsv
	#remove the first line of the prokka result.txt, convert space to tab, remove the colon
	tail -n +2 $F | tr ' ' \\t | tr -d : >> tempfile.tsv
	# we have to sort alfabetically as prokka does not output results to the txt in the same order for 
	# each annotation run!
	# sort the tempfile cut the second column and paste as a row (transpose), and add this to prokkaresults.tsv
	sort tempfile.tsv | cut -f2 | paste -s >> prokkaresults.tsv
done

#Reorder columns for more logical output.
awk 'BEGIN {FS="\t"}{OFS="\t"}{print $6,$2,$5,$1,$3,$7,$8,$9,$4}' prokkaresults.tsv > prokkaresults_reorder.tsv
```

#### Quast
Several statistics were collected from Quast
GCF_000025985.1 were used as reference for Quast.

```
allfasta=$(ls *.fasta)
ref_genome=GCF_000025985.1_ASM2598v1_genomic.fna
ref_features=GCF_000025985.1_ASM2598v1_genomic.gff

quast.py $allfasta \
        -r $ref_genome \
        --features $ref_features \
        -o quast_output \
	--threads $NPROCS \
	--min-contig 500 \
	--k-mer-stats \
	--plots-format png \
	--no-icarus \
	--no-html

#cleaning
cd quast_output

rm -r aligned_stats
rm -r k_mer_stats
rm -r basic_stats
rm -r contigs_reports
rm -r genome_stats
```

### Identification of resistance genes and IS elements with [ABRicate](https://github.com/tseemann/abricate)
ABRicate does not output the orientation of a hit in the results. To ascertain if an IS elements
is oriented upstream and in the opposite strand of the resistance gene of interest, a small change
was made to the ABRicate script prior to use. See ABRicate [issue 83](https://github.com/tseemann/abricate/issues/83).
Resistance genes were identified using the Resfinder, CARD and NCBI databases that come with ABRicate.
Gene sequences of the putative multidrug efflux pumps BexA (GenBank: AB067769.1:3564..4895) and BexB (GenBank: AY375536.1:4599..5963) were downloaded from the NCBI Nucleotide database.
The [ISFinder](https://isfinder.biotoul.fr/howto.php) Insertion Sequence database was downloaded from [here](https://github.com/thanhleviet/ISfinder-sequences) and duplicates removed (see [this issue](https://github.com/thanhleviet/ISfinder-sequences/issues/1).
The above databases and sequences were collated into one ABRicate database (abricate_tvs) following the guide on ABRicates [github site](https://github.com/tseemann/abricate).
ABRicate was run with very low --minid and --mincov settings, as we previously had experienced missing hits (from split genes or low homology sequences).
```
abricate --db abricate_tvs --minid 40 --mincov 25 --threads 3 <contigs.fasta> > out.tab
```

### Identification of plasmid replicon domain families.
We wanted further information and evidence supporting whether circularised contigs were indeed plasmid sequences.
Increased relative coverage compared to the main replicon/chromosome is indicative of plasmid sequence, as is circularisation.
Sequences of the plasmid replication domain families, listed in Table 1 in [Jørgensen et al 2014](https://doi.org/10.1371/journal.pone.0087924) were downloaded from the [Pfam database](https://pfam.xfam.org/).
Specifically the full-length sequences for all sequences in the full alignments were downloaded for:
[PF01051](https://pfam.xfam.org/family/PF01051#tabview=tab3)
[PF01402](https://pfam.xfam.org/family/PF01402#tabview=tab3)
[PF01446](https://pfam.xfam.org/family/PF01446#tabview=tab3)
[PF01719](https://pfam.xfam.org/family/PF01719#tabview=tab3)
[PF01815](https://pfam.xfam.org/family/PF01815#tabview=tab3)
[PF02486](https://pfam.xfam.org/family/PF02486#tabview=tab3)
[PF03090](https://pfam.xfam.org/family/PF03090#tabview=tab3)
[PF03428](https://pfam.xfam.org/family/PF03428#tabview=tab3)
[PF04796](https://pfam.xfam.org/family/PF04796#tabview=tab3)
[PF05732](https://pfam.xfam.org/family/PF05732#tabview=tab3)
[PF06504](https://pfam.xfam.org/family/PF06504#tabview=tab3)
[PF06970](https://pfam.xfam.org/family/PF06970#tabview=tab3)
[PF07042](https://pfam.xfam.org/family/PF07042#tabview=tab3)
[PF10134](https://pfam.xfam.org/family/PF10134#tabview=tab3)

The following awk one liner was used to improve the sequence names for use with ABRicata.
The example below is for PF03428_full_length_sequences.fasta
```
awk -F '>' '/^>/ { $1 = ">pfam~~~RP-C~~~" } { gsub (/[()]/,"",$2) } { gsub (" ","\t",$2) } { gsub ("\t"," PF03428.8 RP-C Replication protein C N-terminal domain 0 ",$2) }{ print $1 $2}' PF03428_full_length_sequences.fasta > PF03428_RP-C_abr.fasta
```
The individual *_abr.fasta files were concatenated, but the concatenated file contained duplicates.
[fasta_unique](https://github.com/b-brankovics/fasta_tools/blob/master/bin/fasta_unique) by @b-brankovics was used to remove duplicates.
The removed duplicates were writen to a .tab file for later review, and [Pfam_replicationdomain_unique.fasta](../master/Pfam_replicationdomain_unique.fasta) used to build a protein sequence database with ABRicate.
```
cat *_abr.fasta > Pfam_replicondomains.fasta
perl fasta_unique.perl Pfam_replicationdomain.fasta> Pfam_replicationdomain_unique.fasta 2>Pfam_replicationdomain_unique.tab
```


###

#### BLAST settings for ONT read mapping in CLC Genomics Workbench
```
# Local BLAST of np reads to assembly
'Filter low Complexity'
Expect	0,00000000000000000000000000000001
Word size		11
Match/mismatch	'Match 1, Mismatch -2'
Gap costs		Existence 1, Extension 1
Max number of hit sequences		500
```

## Data repository links
The NCBI/ENA Bioprojects contain sequence data and final assemblies:

| Strain | Alternative desigatores | Bioproject accession |
|--- | --- | --- |
| CCUG4856T | ATCC25285, NCTC9343| [PRJNA525024](http://www.ncbi.nlm.nih.gov/bioproject/PRJNA525024) |
| BF017 | O17, BF17, DCMOUH0017B| [PRJNA244943](http://www.ncbi.nlm.nih.gov/bioproject/PRJNA244943) |
| BFO18 | O18, BF18, DCMOUH0018B| [PRJNA244944](http://www.ncbi.nlm.nih.gov/bioproject/PRJNA244944) |
| S01 | DCMSKEJBY001B |[PRJNA244942](http://www.ncbi.nlm.nih.gov/bioproject/PRJNA244942) |
| BFO42| O42, BF42, DCMOUH0042B| [PRJNA253771](http://www.ncbi.nlm.nih.gov/bioproject/PRJNA253771) |
| BFO67| O67, BF67, DCMOUH0067B| [PRJNA254401](http://www.ncbi.nlm.nih.gov/bioproject/PRJNA254401) |
| BFO85| O85, BF85, DCMOUH0085B| [PRJNA254455](http://www.ncbi.nlm.nih.gov/bioproject/PRJNA254455) |

The 141 variatrions of assemblies and polished assembleis of CCUG4856T are available at [![DOI](https://www.zenodo.org/badge/DOI/10.5281/zenodo.2648546.svg)](https://doi.org/10.5281/zenodo.2648546)

The stages of genome assemblies for the complete circular assemblies of the seven isolates are available at [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.2661704.svg)](https://doi.org/10.5281/zenodo.2661704)

