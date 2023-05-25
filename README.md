# Genome assembly, annotation, and RNA-Seq scripts/code

This is all the scripts and code used in the paper __________

## Package versions

[Versions of all command-line software used are tabulated here](https://github.com/kd-lab/Genome_Assembly_Annotation/blob/main/Package_versions.xlsx)

## Check quality of reads

```bash
conda activate fastqc
fastqc -t 16 /data/hifi_reads.fastq.gz -o /step1_fastqc
```

## Generate initial stats

```bash
meryl count k=21 memory=50GB threads=20 /data/hifi_reads.fastq.gz output read-db.meryl
meryl histogram read-db.meryl > meryl_histogram.txt

genomescope2 --input meryl_histogram.txt --output genomescope2 --kmer_length 21 --testing
```

## Assemble genome with hifiasm

```bash
conda activate hifiasm
cd /step3_hifiasm
hifiasm -t 20 -o output -f 37 -l 0 -s 0.75 -O 1 --purge-max 114 --primary /data/hifi_reads.fastq
```

## Convert gfa to fasta

```bash
cd /step3_hifiasm
awk '/^S/{print ">"$2;print $3}' output.p_ctg.gfa > primary_assembly.fa
awk '/^S/{print ">"$2;print $3}' output.a_ctg.gfa > alternate_assembly.fa
```

## Check quality of initial assembly

```bash
cd /step4_quality
quast --pacbio /data/hifi_reads.fastq --labels Primary_assembly,Alternate_assembly -o quast_out --est-ref-size 1184000000 --eukaryote --min-contig 500 --min-alignment 65 --min-identity 95.0 --ambiguity-usage 'one' --ambiguity-score '0.99' --contig-thresholds '0,1000' --extensive-mis-size 1000 --scaffold-gap-max-size 1000 --unaligned-part-size 500 /step3_hifiasm/primary_assembly.fa /step3_hifiasm/alternate_assembly.fa --threads 16

busco --in /step3_hifiasm/primary_assembly.fa --update-data --mode 'geno' --out busco_out --cpu 8 --evalue 0.001 --limit 3 --auto-lineage-euk

ln -s $MERQURY/merqury.sh
cd /step4_quality/merqury_out
merqury.sh step2_initial_stats/read-db.meryl /step3_hifiasm/primary_assembly.fa /step3_hifiasm/primary_assembly.fa output_merqury
```

## Scaffold

No improvement observed so this step (step 5) was skipped

# Genome Annotation

## Identify repeats, mask

Create blast database which is used for later steps:

```bash
cd /step6_repeat_mask
BuildDatabase -name pokeweed2022 /step3_hifiasm/primary_assembly.fa
```

Run repeatmodeler

```bash
cd /step6_repeat_mask
nohup RepeatModeler -database pokeweed2022 -pa 14 -LTRStruct >& run.out &
```

Mask repeats in assembly

```bash
cd /step6_repeat_mask
nohup RepeatMasker -pa 3 -xsmall -a -s -gff -no_is -lib pokeweed2022-families.fa /step3_hifiasm/primary_assembly.fa &> RM.run.out &
```

## Gather and align RNA-Seq reads for gene prediction

Index genome assembly

```bash
cd /step7_genome_annotation/mRNA
ln /step3_hifiasm/primary_assembly.fa
mkdir star_indexes
STAR --runMode genomeGenerate --runThreadN 15 --genomeFastaFiles primary_assembly.fa
```

Save as file called new_data/accessions.txt

```
SRR11497117
SRR11497118
...
```

Download raw data, convert to fastq

```bash
cd /step7_genome_annotation/mRNA/new_data
conda activate sra-tools
prefetch --option-file accessions.txt &> prefetch_log.txt
cat accessions.txt | xargs fasterq-dump --outdir fasta &> fasterq_log.txt
```

## Check quality of downloaded data
Note: since commands were all the same for each fasta file subset what is shown below is the template used for each, square brackets would be where code changes for each bioproject/collection of files

```bash
fastqc -f fastq -t 6 [path to fasta file]
Compress files
pigz *.fastq
```

## Filter/trim
Shell script file makes it easier to trim so many fasta files:

```
#!/bin/bash
# arg1: number of threads
# to run: 
# chmod +x trim.sh
# <path>/trim.sh <number of threads>
# Example: ./trim.sh 40

for f in *_1.fastq.gz # for each sample

do
    n=${f%%_1.fastq.gz} # strip part of file name
    trimmomatic PE -threads $1 ${n}_1.fastq.gz  ${n}_2.fastq.gz \
    ${n}_1_trimmed.fastq.gz ${n}_1_unpaired.fastq.gz ${n}_2_trimmed.fastq.gz \
    ${n}_2_unpaired.fastq.gz \
    LEADING:3 TRAILING:3 SLIDINGWINDOW:4:25 MINLEN:36

done
```

Run the script to trim each file (note: substitute square brackets for path to fasta files)

```bash
conda activate trimmomatic
chmod +x trim.sh
cd [path/to/fasta/files]
./trim.sh 12
```

## Sort the files
```bash
mv *_trimmed.fastq.gz ./trimmed/
mv *_unpaired.fastq.gz ./unpaired/
```

## Map downloaded reads to genome, sort

Note: since commands were all the same for each fasta file subset what is shown below is the template used for each, square brackets would be where code changes for each bioproject/collection of files

```bash
cd /step7_genome_annotation/mRNA/[unique identifier]/trimmed
conda activate star
STAR --genomeDir /step7_genome_annotation/mRNA/GenomeDir --runThreadN 15 --readFilesIn [forward_trimmed.fastq.gz files, separates by comma] [reverse_trimmed.fastq.gz files] --readFilesCommand gunzip -c --readFilesType Fastx --outFileNamePrefix /step7_genome_annotation/mRNA/star_aligned/[unique identifier] --outSAMtype BAM Unsorted

conda activate samtools
cd /step7_genome_annotation/mRNA/star_aligned
samtools sort -@ 10 [unique identifier]Aligned.out.bam -o [unique identifier]_RNAseq_sorted.bam
samtools index [unique identifier]_RNAseq_sorted.bam
```

## Gather protein data for gene prediction

Get protein data (Orthodb)

```bash
#download from orthodb
cd /step7_genome_annotation/protein
wget https://v100.orthodb.org/download/odb10_plants_fasta.tar.gz
#unzip (output in file: /step7_genome_annotation/protein/plants/Rawdata)
tar xvf odb10_plants_fasta.tar.gz
```

Get protein data (Uniprot)

Download from: https://www.uniprot.org/taxonomy/3524

Then transfer to server in /step7_genome_annotation/protein/plants/Rawdata

```bash
cd /step7_genome_annotation/protein/plants/Rawdata
gunzip uniprot-compressed_true_download_true_format_fasta_includeIsoform_tr-2022.12.01-22.47.06.41.fasta.gz
cd /step7_genome_annotation/protein
#combine into one file
cat protein/plants/Rawdata/* > proteins.fasta
```

## Predict hints

```bash
cd /step7_genome_annotation/protein
prothint.py /step6_repeat_mask/primary_assembly_masked.fa /step7_genome_annotation/protein/proteins.fasta
BRAKER2 annotation for RNA
```

Note: square brackets denote where all the sorted RNA-seq BAM file paths would go, separated by commas

```bash
cd /step7_genome_annotation/mRNA
conda activate braker2
nohup braker.pl --cores=14 --species=Pokeweed --genome=/step6_repeat_mask/primary_assembly_masked.fa --bam=[path/to/_RNAseq_sorted.bam] --softmasking &> BRAKER_out.txt &
BRAKER2 annotation for protein
conda activate braker2

cd /step7_genome_annotation/protein
nohup braker.pl --species=Pokeweed_prot --cores=8  --genome=/step6_repeat_mask/primary_assembly_masked.fa --hints=prothint_augustus.gff --softmasking --AUGUSTUS_CONFIG_PATH=/step7_genome_annotation/mRNA/config --useexisting &> BRAKER_out.txt &
Combine the two
cd /step7_genome_annotation/combined
git clone https://github.com/Gaius-Augustus/TSEBRA

conda activate braker2
./TSEBRA/bin/tsebra.py -g /step7_genome_annotation/mRNA/braker/augustus.hints.gtf,/step7_genome_annotation/protein/braker/augustus.hints.gtf -c ./TSEBRA/config/default.cfg -e /step7_genome_annotation/mRNA/braker/hintsfile.gff,/step7_genome_annotation/protein/braker/hintsfile.gff -o annotations_combined.gtf
```

To convert it to GFF3, upload it to galaxy and convert it there
    1. upload combined gtf to galaxy (usegalaxy.org) 
    2. make sure file type is set to gtf
    3. use gffread to convert to gff

Input BED, GFF3 or GTF feature file 

```
filters 	Nothing selected.
Filter by genome region 	none
Filter out transcipts with large introns 	Not available.
Replace reference sequence names 	
Transcript merging 	none
Reference Genome 	none
Feature File Output 	gff
Ensembl GTF to GFF3 conversion 	False
Trackname to use in the second column of each GFF output line 	braker_merged
full GFF attribute preservation (all attributes are shown) 	False
decode url encoded characters within attributes 	False
warn about duplicate transcript IDs and other potential problems with the given GFF/GTF records 	False
```

    4. download the output

## Genome annotation quality check

make fasta file from annotations

```bash
cd step7_genome_annotation/combined
gffread -w annotations_combined.fa -g primary_assembly.fa annotations_combined.gff3
BUSCO
conda activate busco
cd /step9_annotation_quality

busco --in /step7_genome_annotation/mRNA/merged/annotations_combined.fa --update-data --mode 'transcriptome' --out busco_out_braker --cpu 8 --evalue 0.001 --limit 3 --auto-lineage-euk
```

## Try other annotation programs

Multiple attempts from MAKER produced lower quality results (step 8), so were not used.

# Protein domain and GO identification/annotation

## Interpro 

```bash
cd /step7_genome_annotation/combined
transeq -sequence annotations_combined.fa -frame 1 -clean Y -trim Y -outseq annotations_combined_pep.fa
```

Predict

```bash
cd /step10_interpro
ln -s /step7_genome_annotation/combined/annotations_combined_pep.fa

nohup /home/kyrad/my_interproscan/interproscan-5.59-91.0/interproscan.sh -i annotations_combined_pep.fa -goterms -iprlookup -pa -cpu 6 >& interpro.out &
```

## Blast against swissprot database

```bash
conda activate blast
cd /step11_blast_swisspro

ln -s /step7_genome_annotation/combined/annotations_combined_pep.fa

/home/miniconda3/envs/blast/bin/update_blastdb.pl --decompress swissprot

nohup blastp -db swissprot -query annotations_combined_pep.fa -evalue 1e-50 -outfmt 10 -out blast_swiss_annotations_combined.csv &> blast_out.txt &
```

## eggNOG

    1. go to http://eggnog-mapper.embl.de/
    2. upload translated transcript sequences
    3. adjust options:
        a. minimum hit e-value: 0.00001
        b. minimum hit bit-score: 60
        c. percentage identity: 40
        d. minimum% of query coverage: 60
        e. minimum % of subject coverage: 60
        f. taxonomic scope: viridiplantae
        g. orthology restrictions: transfer annotations from any ortholog
        h. gene ontology evidence: transfer non-electronic annotations
        i. pfam refinement: report pfam domains from orthologs
        j. smart annotation: skip smart annotation
        k. database: eggNOG5
    4. submit, download results

## Combine GO terms in R

[See this R markdown file for all R code](https://github.com/kd-lab/Genome_Assembly_Annotation/blob/main/R_code.rmd)

---

# RNA-Seq analysis

## Filter/trim (already done for genome assembly)
```bash
cd /RNAseq/Step1_filter_trim
# make symbolic link of trimmed reads
ln -s -t . /data/JA_timecourse/trimmed/*.fastq.gz
```

## Align reads with STAR
```bash
cd /RNAseq/Step2_align
conda activate star
mkdir star_indexes
STAR --runMode genomeGenerate --runThreadN 10 --genomeDir ./star_indexes --genomeFastaFiles primary_assembly.fa --sjdbGTFtagExonParentTranscript Parent --sjdbGTFfile /step7_genome_annotation/combined/annotations_combined.gff3 --sjdbOverhang 100
mkdir star_aligned
ln -s -t . /data/JA_timecourse/trimmed/*.fastq.gz
chmod +x align.sh
nohup ./align.sh &> star_out.txt &
```

code in align.sh
```
#!/bin/bash
# to run: 
# chmod +x align.sh
# <path>/align.sh
# Example: ./align.sh
for i in *1_trimmed.fastq.gz; do STAR --genomeDir ./star_indexes --runThreadN 10 --readFilesIn ${i} ${i%1_trimmed.fastq.gz}2_trimmed.fastq.gz --readFilesType Fastx --readFilesCommand zcat --readMapNumber -1 --clip3pNbases 0 --clip5pNbases 0 --outFileNamePrefix ./star_aligned/${i%1_trimmed.fastq.gz} --outSAMtype BAM Unsorted; done
```

Look at the files ending in Log.final.out to see how well the reads aligned (you want to make sure there wasn't a large amount of contamination in your reads)

## Sort with samtools

Sort each file by read name and genome position

```bash
cd /RNAseq/Step2_align/star_aligned
conda activate samtools
for i in *Aligned.out.bam; do samtools sort -@10 -n -o - ${i} | samtools sort -@10 -o /RNAseq/Step3_sort/${i%Aligned.out.bam}Aligned_sorted.bam -; done
```

## Count reads with htseq

```bash
cd /RNAseq/Step3_sort
conda activate htseq
chmod +x count.sh
nohup ./count.sh &> htseq_out.txt &
```

in count.sh

```
#!/bin/bash
# to run: 
# chmod +x count.sh
# <path>/count.sh
# Example: ./count.sh
for i in *A*Aligned_sorted.bam; do htseq-count -m union -s yes -r pos -i Parent -a 10 -f bam ./${i} /step7_genome_annotation/combined/annotations_combined.gff3 > /RNAseq/Step4_count/${i%_L003_Aligned_sorted.bam}.tsv; done
```

## Differential expression analysis in R

[See this R markdown file for all R code](https://github.com/kd-lab/Genome_Assembly_Annotation/blob/main/R_code.rmd)
