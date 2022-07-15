#!/bin/zsh

SECONDS=0

####################################################################################################
## RNA seq analysis pipeline for fastq files
####################################################################################################

# change working directory
cd /Users/ryan/test_analysis/bioinformatics_101/


### STEP 1: Quality Control with fastqc

# help file: fastqc -h
fastqc data/fastq/demo.fastq -o data/fastq/
echo "fastqc finished"

# use trimmomatic to trim reads with poor quality, if applicable
java -jar /usr/local/bin/Trimmomatic-0.39/trimmomatic-0.39.jar SE -threads 4 data/fastq/demo.fastq data/fastq/demo_trimmed.fastq TRAILING:10 -phred33
echo "trimmomatic finished"

# re-run fastqc with trimmed data file
fastqc data/fastq/demo_trimmed.fastq -o data/fastq/
echo "fastqc trimmed finished"


### STEP 2: Alignment with HISAT2

# mkdir HISAT2

# download the genome indices and extract files
# wget https://genome-idx.s3.amazonaws.com/hisat/grch38_genome.tar.gz
# tar -xvf grch38_genome.tar.gz

# run hisat2 alignment
hisat2 -q --rna-strandness R -x grch38/genome -U data/fastq/demo_trimmed.fastq | samtools sort -o HISAT2/demo_trimmed.bam
echo "HISAT2 finished"


### STEP 3: Quantification with featurecounts

# get gtf
# wget http://ftp.ensembl.org/pub/release-106/gtf/homo_sapiens/Homo_sapiens.GRCh38.106.gtf.gz
# gunzip Homo_sapiens.GRCh38.106.gtf.gz

# mkdir quants

featureCounts -S 2 -a Homo_sapiens.GRCh38.106.gtf -o quants/demo_featurecounts.txt HISAT2/demo_trimmed.bam
echo "featureCounts finished"


duration=$SECONDS
echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed."

# preview counts file
cat quants/demo_featurecounts.txt | cut -f1,7 | head
