#!/bin/bash
#SBATCH --cpus-per-task=8 --mem=16g

GENOME=$1
GTF=$2
OUTDIR=Step3-star_mapping
READLEN=48

module load STAR

# Create genome index
time STAR --runMode genomeGenerate \
    --runThreadN 16 \
    --genomeDir ${OUTDIR} \
    --sjdbGTFfile ${GTF} \
    --sjdbGTFfeatureExon exon \
    --sjdbGTFtagExonParentGeneType protein_coding \
    --sjdbOverhang ${READLEN} \
    --genomeFastaFiles ${GENOME}

