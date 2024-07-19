#!/bin/bash
#SBATCH --cpus-per-task=16
#SBATCH --mem=32g

module load picard
sample=example
time java -Xmx8g -XX:ParallelGCThreads=5 -jar $PICARDJARPATH/picard.jar MarkDuplicates \
CREATE_INDEX=true \
METRICS_FILE=$sample.metrix.txt \
INPUT=${sample}_bwa_sorted.bam \
REMOVE_DUPLICATES=false \
OUTPUT=${sample}_bwa_sorted_dedup.bam
