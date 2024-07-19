#!/bin/bash
#SBATCH --cpus-per-task=16

echo; echo RUNNING FASTP; echo

# Load read files from the command line
READ1=$1
READ2=$2

# Load some other variables used by fastp
ADAPTERS=${WORKSHOPDIR}/data/TruSeq3-PE.fa
THREADS=16
OUTDIR_FASTP=${WORKSHOPDIR}/Step2-trimming_fastp
READNAME=$(basename ${READ1})
PREFIX=${READNAME%.R1.fastq.gz}

# Create output directory
mkdir -p ${OUTDIR_FASTP}

# Run fastp
module load fastp
fastp -i ${READ1} -I ${READ2} \
  -o ${OUTDIR_FASTP}/${PREFIX}.paired.R1.fastq.gz -O ${OUTDIR_FASTP}/${PREFIX}.paired.R2.fastq.gz \
  --unpaired1 ${OUTDIR_FASTP}/${PREFIX}.unpaired.R1.fastq.gz --unpaired2 ${OUTDIR_FASTP}/${PREFIX}.unpaired.R2.fastq.gz \
  --qualified_quality_phred 20 \
  --length_required 25 \
  --adapter_fasta ${ADAPTERS} \
  --trim_front1 1 \
  --trim_front2 1 \
  --cut_right --cut_mean_quality 20 --cut_window_size 5 \
  --thread ${THREADS} &> ${OUTDIR_FASTP}/fastp.log
#  --umi --umi_loc=read1 --umi_len=8  # Required if reads contain UMIs
