## Workshop guide

### A- File formats

1. Fasta
2. Fastq
3. SAM/BAM
4. GTF/GFF
5. VCF

## B- Proprocessing of sequencing reads
### B.1 Quality control of sequencing data 
```
READ1=my_reads.R1.fastq.gz
READ2=my_reads.R2.fastq.gz
OUTDIR=FASTQC

module load fastqc
fastqc -o $OUTDIR $READ1 $READ2
```

Going through a fastqc report


### B.2 Trimming reads

Single-end, paired ends
Preparing adaptor fasta file for triiming?


### B.3 Mapping reads 
