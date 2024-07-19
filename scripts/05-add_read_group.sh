#!/bin/bash
#SBATCH

# Get bam file from the command line
INPUT_BAM=$1
DIRNAME=$(dirname $INPUT_BAM)
PREFIX=example.sorted

module load samtools

samtools addreplacerg -w -r "@RG\tID:RG1\tSM:SampleName\tPL:Illumina\tLB:Library.fa" -o ${DIRNAME}/${PREFIX}.bam ${INPUT_BAM}
