#!/bin/bash
#SBATCH --cpus-per-task=16

READ1=$1
READ2=$2
OUTDIR=Step1-fastqc

# Create output directory
mkdir -p $OUTDIR

module load fastqc
fastqc -o $OUTDIR $READ1 $READ2
