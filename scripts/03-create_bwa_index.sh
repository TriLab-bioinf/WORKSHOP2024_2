#!/bin/bash
#SBATCH --cpus-per-task=16 --mem=32g

genome=$1
OUTDIR=03-reference_index
mkdir -p ${OUTDIR}

module load bwa
bwa index $genome -p ${OUTDIR}/hg38_chr17