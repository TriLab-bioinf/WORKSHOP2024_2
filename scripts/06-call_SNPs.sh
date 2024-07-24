#!/bin/bash
#SBATCH --cpus-per-task=16
#SBATCH --mem=32g

GENOME=./WGS_data/hg38_chr17.fa
bamfile=$1

OUTDIR=${WORKSHOPDIR}/06-SNPcalling
mkdir -p $OUTDIR

READ_NAME=$(basename $bamfile)
PREFIX=${READ_NAME%.dedup.bam}

module load bcftools
bcftools mpileup -f $GENOME $bamfile | bcftools call -mv -Oz -o $OUTDIR/$PREFIX.vcf.gz
