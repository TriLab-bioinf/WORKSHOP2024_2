#!/bin/bash
#SBATCH --cpus-per-task=16
#SBATCH --mem=32g

## This is for Paired End data

# Load read files from the command line
READ1=$1
READ2=$2

# Load some other variables
adapter=/usr/local/apps/trimmomatic/0.39/adapters/TruSeq3-PE.fa
THREADS=16
OUTDIR_TRIMMOMATIC=${WORKSHOPDIR}/02-trimming_trimmomatic
READNAME=$(basename ${READ1})
PREFIX=${READNAME%_R1.fastq.gz}

# Create output directory
mkdir -p ${OUTDIR_TRIMMOMATIC}

# run trimmomatic
module load trimmomatic
java -Djava.io.tmpdir=. -jar $TRIMMOJAR PE -phred33 -threads $THREADS \
    $READ1 $READ2 \
    $OUTDIR_TRIMMOMATIC/${PREFIX}_forward_paired.fq.gz $OUTDIR_TRIMMOMATIC/${PREFIX}_forward_unpaired.fq.gz \
    $OUTDIR_TRIMMOMATIC/${PREFIX}_reverse_paired.fq.gz $OUTDIR_TRIMMOMATIC/${PREFIX}_reverse_unpaired.fq.gz \
    ILLUMINACLIP:$adapter:2:30:10 LEADING:3 TRAILING:3 \
    SLIDINGWINDOW:4:15 MINLEN:36
