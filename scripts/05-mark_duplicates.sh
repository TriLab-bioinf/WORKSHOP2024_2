#!/bin/bash
#SBATCH --cpus-per-task=5 --mem=32g --gres=lscratch:40

# Flag duplicated reads with Picard
module load picard/3.2.0

# Enter bam file from the command line
BAM=$1
PREFIX=example
OUTDIR=${WORKSHOPDIR}/05-markduplicates

# Make output directory
mkdir -p $OUTDIR

java -Xmx32g -XX:ParallelGCThreads=5 -jar $PICARDJARPATH/picard.jar MarkDuplicates \
  --INPUT ${BAM} \
  --OUTPUT ${OUTDIR}/${PREFIX}.dedup.bam \
  --METRICS_FILE ${OUTDIR}/${PREFIX}.metrix.txt \
  --CREATE_INDEX true \
  --REMOVE_DUPLICATES false \
  --TMP_DIR /lscratch/$SLURM_JOBID
