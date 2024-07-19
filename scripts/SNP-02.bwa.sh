#!/bin/bash
#SBATCH --cpus-per-task=16
#SBATCH --mem=32g

module load bwa
module load samtools

genome=hg38_chr17.fa
id=example
sm=example
bwa mem -t 32 \
        -R "@RG\tID:$id\tPL:ILLUMINA\tSM:$sm" \
        $genome \
        ${sm}_forward_paired.fq.gz \
        ${sm}_reverse_paired.fq.gz \
        | samtools sort \
        -@ 12 \
        -O BAM \
        -o ${sm}_bwa_sorted.bam
