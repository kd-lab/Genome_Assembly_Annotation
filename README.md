# Genome assembly and annotation scripts readme

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

# Predict hints

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

Combine interproscan, eggnog, and blast results

```r
interpro <- read_tsv("interpro-blast/annotations_combined_pep.fa.tsv", col_names = c("Protein_accession", "Sequence_MD5_digest", "Sequence_length", "Analysis", "Signature_accession", "Signature_description", "Start_location", "Stop_location", "Evalue", "Status_of_match", "Date", "InterPro_annotations", "InterPro_annotation_description", "GO_annotations", "Pathways_annotations")) %>% 
  mutate(Protein_accession = str_remove(Protein_accession, "_[0-9]")) %>% 
  dplyr::select(!c(Sequence_MD5_digest,Pathways_annotations))
swi <- read_csv("blast_swiss_annotations_combined.csv", col_names = NULL) %>% 
  # select only top hit
  # group by transcript id
  group_by(X1) %>% 
  # select hit with lowest evalue
  dplyr::slice(which.min(X11)) %>% 
  # select the transcrtipt id and associated swissprot accession
  dplyr::select(c(X1, X2)) %>% 
  dplyr::rename("Protein_accession" = X1, "Entry" = X2) %>% 
  mutate(Protein_accession = str_remove(Protein_accession, "_[0-9]"))
write(paste(unique(swi$Entry), collapse = " "), "data_out/swissprot_hits.txt")
```

1. go to https://www.uniprot.org/id-mapping and load the text file
2. select from: UniProtKB_AC-ID     to: UniProtKB-Swiss-Prot
3. select the results when loaded
4. click 'customize columns' deselect all defaults, and select 'gene ontology IDs', "gene names", and "descriptions"
5. download results as tsv, uncompressed
6. rename to "get_swissprot_go.tsv"

```r
temp <- read_tsv("data_in/get_swissprot_go.tsv") %>% dplyr::select(c(From, Gene_Ontology_IDs)) %>% dplyr::rename("Entry" = From, "GO_sw" = "Gene_Ontology_IDs")

head(temp)
#check for improper concatenation
temp %>% filter(str_detect(GO_sw, "[0-9]GO"))
swi <- full_join(swi, temp, by = "Entry") %>% mutate(GO_sw = str_replace_all(GO_sw, "; ", ",")) %>% dplyr::select(!Entry) %>% mutate(GO_sw = ifelse(is.na(GO_sw), "", GO_sw))
#check for improper concatenation
swi %>% filter(str_detect(GO_sw, "[0-9]GO"))
```

add eggNOG

```r
eggnog <- read_tsv("C:/Users/Kyra/Desktop/Masters/improve_genome/interpro-blast/eggnog.tsv", comment = "#") %>% 
  dplyr::select(c(Protein_accession, GOs)) %>% 
  mutate(Protein_accession = str_remove(Protein_accession, "_[0-9]")) %>% 
  mutate(GOs = str_remove(GOs, "-"))
#check for improper concatenation
eggnog %>% filter(str_detect(GOs, "[0-9]GO"))
```

get GO terms per transcript

```r
GO_clean <- interpro[,c(1,13)] %>% 
  # group by transcriptID (as each has multiple lines)
  group_by(Protein_accession) %>% 
  # add a new column concatenating all that are in the same group, and remove the old GO column
  mutate(GO = paste0(GO_annotations, collapse = ","), .keep = "unused") %>% 
  # the previous step means that each duplicate will contain all the same data in the GO column, so now can delete duplicates
  distinct(Protein_accession, .keep_all = T) %>% 
  # clean up GO column to delete "NA", "-", and "|"
  mutate(GO = str_remove_all(GO, "NA[,]*")) %>% 
  mutate(GO = str_remove_all(GO, "-[,]*")) %>% 
  mutate(GO = str_remove_all(GO, ",$")) %>% 
  mutate(GO = str_replace_all(GO, "\\|", ","))
#check for improper concatenation
GO_clean %>% filter(str_detect(GO, "[0-9]GO"))
```

Combine

```r
GO_clean <- left_join(GO_clean, eggnog, by = "Protein_accession")

GO_clean <- full_join(GO_clean, swi, by = "Protein_accession") %>% 
  # convert added NAs to blanks
  replace(is.na(.), "")
GO_clean <- GO_clean %>% 
  #mutate(GO = paste0(GO, GOs, GO_sw, collapse = ",")) %>% 
  #mutate(GO = paste0(GO, GOs, collapse = ",")) %>%
  unite("GO", c(GO, GOs, GO_sw), sep = ",", remove = T) %>% 
  # remove extra commas
  # doubles
  mutate(GO = str_replace(GO, ",,", ",")) %>% 
  # at start
  mutate(GO = str_remove(GO, "^,")) %>% 
  # at end
  mutate(GO = str_remove(GO, ",$")) %>% 
  # remove go_sw and GOs now that part of GO
  #dplyr::select(c(Protein_accession, GO)) %>% 
  # split cleaned GO terms for removal of duplicates (next step)
  mutate(GO_str = strsplit(GO, ","))
#check for improper concatenation
GO_clean %>% filter(str_detect(GO, "[0-9]GO"))
Get rid of duplicate GO terms, refine
for (i in 1:nrow(GO_clean)){
  GO_clean$GO[i] <- paste0(unique(unlist(GO_clean$GO_str[i])), collapse = ",")
}

GO_clean <- GO_clean %>% 
  # convert blanks to NAs
  mutate_all(na_if,"") %>% 
  # remove GO_str
  dplyr::select(!GO_str)
GO_names <- GO_clean$GO[!is.na(GO_clean$GO)]
GO_names <- paste0(GO_names, collapse = ",")
GO_names <- strsplit(GO_names, ",")
GO_names <- as.data.frame(GO_names, col.names = "GO")
GO_names <- as.data.frame(Term(GOTERM)) %>% rename("Term(GOTERM)" = "GO_term") %>% rownames_to_column(var = "GO") %>% left_join(GO_names, ., by = "GO") %>% distinct(GO, .keep_all = T) 
bad_GO <- GO_names %>% filter(!is.na(GO_term)) %>% filter(str_detect(GO_term, "aging|larval|nervous system|neur[oa]|behavior|mating|animal|glial|imaginal disc|wing|Schwann|cardiac|atrial|kidney|intestine|stomach|eye|retina|spleen|axon |myoblast|bone|osteo|male|MHC"))

GO_clean <- GO_clean %>% 
  # remove incorrect GO terms
  mutate(GO = str_remove_all(GO, paste0(bad_GO$GO, collapse = "|"))) %>% 
  # remove multiple commas
  mutate(GO = str_replace_all(GO, ",[,]+", ","))
write_tsv(GO_clean, "data_out/GO.tsv", na = "", col_names = F)
```
