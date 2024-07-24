#!/bin/bash
#SBATCH

# Enter path to read files from STDIN 
READ1=$1
READ2=$2

# Set name of output directory
OUTDIR=01-fastqc

# Create output directory
mkdir -p $OUTDIR

module load fastqc
fastqc -o $OUTDIR $READ1 $READ2
