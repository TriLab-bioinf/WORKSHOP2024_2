#!/bin/bash
## This is for Paired End data
#SBATCH --cpus-per-task=16
#SBATCH --mem=32g

sample=example
adapter=/usr/local/apps/trimmomatic/0.39/adapters/TruSeq3-PE.fa
module load trimmomatic
java -Djava.io.tmpdir=. -jar $TRIMMOJAR PE -phred33 -threads 12 \
    ${sample}_R1.fastq.gz ${sample}_R2.fastq.gz \
    ${sample}_forward_paired.fq.gz ${sample}_forward_unpaired.fq.gz \
    ${sample}_reverse_paired.fq.gz ${sample}_reverse_unpaired.fq.gz \
    ILLUMINACLIP:$adapter:2:30:10 LEADING:3 TRAILING:3 \
    SLIDINGWINDOW:4:15 MINLEN:36
