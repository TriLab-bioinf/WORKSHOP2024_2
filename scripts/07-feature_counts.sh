#!/bin/bash
#SBATCH --cpus-per-task=8

BAM=$1
GENOME=${WORKSHOPDIR}/data/GRCh38.chr17.fa
THREADS=8
GTF=${WORKSHOPDIR}/data/gencode.v45.annotation.chr17.gtf
OUTPUT_FILE=${WORKSHOPDIR}/Step7-read_counts/read_counts.txt

# Make output direcory
mkdir -p ${WORKSHOPDIR}/Step7-read_counts

module load subread

featureCounts -G ${GENOME} -T ${THREADS}\
  -a ${GTF} \
  -t exon \
  -g gene_id \
  -O \
  -s 2 \
  -p --countReadPairs -C \
  --ignoreDup \
  -M --fraction \
  -o ${OUTPUT_FILE} ${BAM}
