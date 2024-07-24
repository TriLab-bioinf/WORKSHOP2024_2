#!/bin/bash
#SBATCH --cpus-per-task=16
#SBATCH --mem=32g

# Load read files from the command line
READ1=$1
READ2=$2

# variables setting
THREADS=16
REFERENCE=${WORKSHOPDIR}/03-reference_index/hg38_chr17
OUTDIR=${WORKSHOPDIR}/04-mapping_bwa

mkdir -p ${OUTDIR}

READ_NAME=$(basename $READ1)
PREFIX=${READ_NAME%_forward_paired.fq.gz}


id=$PREFIX
sm=$PREFIX

module load bwa
module load samtools

bwa mem -t ${THREADS} \
        -R "@RG\tID:$id\tPL:ILLUMINA\tSM:$sm" \
        $REFERENCE \
        $READ1 \
        $READ2 \
        | samtools sort \
        -@ ${THREADS} \
        -O BAM \
        -o ${OUTDIR}/${PREFIX}_bwa_sorted.bam
