#!/bin/bash
#SBATCH --cpus-per-task=16 --mem=32g

# Load read files from the command line
READ1=$1
READ2=$2

# STAR-specific variables
THREADS=16
REFERENCE=${WORKSHOPDIR}/Step3-reference_index
OUTDIR=${WORKSHOPDIR}/Step4-mapping_star
READ_NAME=$(basename $READ1)
PREFIX=${READ_NAME%.paired.R1.fastq.gz}

module load STAR

STAR --runMode alignReads \
  --runThreadN ${THREADS} \
  --genomeDir ${REFERENCE} \
  --outFilterMismatchNmax 5 \
  --alignEndsType EndToEnd \
  --readFilesIn ${READ1} ${READ2} \
  --readFilesCommand zcat \
  --outFileNamePrefix ${OUTDIR}/${PREFIX}.sorted. \
  --quantMode GeneCounts \
  --outSAMtype BAM SortedByCoordinate \
  --outSAMattributes All \
  --twopassMode Basic
