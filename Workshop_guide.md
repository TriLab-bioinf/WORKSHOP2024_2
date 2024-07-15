## Workshop guide

### A- File formats

1. Fasta
2. Fastq
  Forward and reverse
  Stranded: how to check if seq data is stranded or not?
  
4. SAM/BAM
5. GTF/GFF
6. VCF

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

[Good Illumina data](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/good_sequence_short_fastqc.html)

[Bad Illumina data](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/bad_sequence_fastqc.html)

[Others examples fromt he fastqc website](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

### B.2 Trimming reads

Single-end, paired ends
Preparing adaptor fasta file for triming?
```
```

### B.3 Dealing with UMIs

5

[UMI-tools](https://umi-tools.readthedocs.io/en/latest/index.html) can handle any UMI tagged sequencing data where deduplication happens after mapping.

The process is to extract the UMIs from the read sequence and add it to the read names. There are two ways to do this, and between them provide the flexibility to handle any read configuration I can think of (see https://umi-tools.readthedocs.io/en/latest/regex.html)

You then map your reads with your favourite mapper.

The next step depends on whether your technique fragments the cDNA before or after PCR. If fragmentation happens after PCR, then the next step is to assign reads to features (e.g. genes) using featureCounts. If PCR happened after fragmentation, then you do the read assignment/quantification after deduping.

Then you group/dedup/count (depending on your downstream application). If fragmentation happened after PCR then you need to do this on a per-gene basis.


### B.4 Mapping reads 

```
module load STAR

```
