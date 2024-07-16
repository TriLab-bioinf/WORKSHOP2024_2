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
READ1=example.R1.fastq.gz
READ2=example.R2.fastq.gz
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

THREADS=16
READ_PREFIX=example
OUTPUT_PREFIX=example
ADAPTERS=data/00adapters/truseq.fa.gz
LOG=bbduk.log

module load bbtools/39.06
bbduk.sh -Xmx1g threads= \
  in1=${READ_PREFIX}.R1.fastq.gz in2=${READ_PREFIX}.R2.fastq.gz \
  out1=${READ_PREFIX}.paired.R1.fastq.gz out2=${READ_PREFIX}.paired.R2.fastq.gz outs=$${READ_PREFIX}.unpaired.fastq.gz \
  ref=${ADAPTERS} \
  ktrim=r k=23 mink=11 hdist=1 tpe tbo qtrim=rl trimq=20 overwrite=t \
  stats=${LOG}

# Trimmomatic
MIN_READ_LEN=50
LEADING_BASES=0
TRALING_BASES=0

module load trimmomatic
java -jar $TRIMMOMATIC_JAR PE ${READ_PREFIX}.R1.fastq.gz in2=${READ_PREFIX}.R2.fastq.gz \
  ${READ_PREFIX}.paired.R1.fastq.gz ${READ_PREFIX}.unpaired.R1.fastq.gz \
  ${READ_PREFIX}.paired.R2.fastq.gz ${READ_PREFIX}.unpaired.R2.fastq.gz \
  ILLUMINACLIP:${ADAPTERS}:2:30:10:2:True LEADING:${LEADING_BASES} TRAILING:${TRALING_BASES} MINLEN:${MIN_READ_LEN}

# Fastp

module load fastp
fastp -i ${READ_PREFIX}.R1.fastq.gz -I ${READ_PREFIX}.R2.fastq.gz \
  -o ${READ_PREFIX}.paired.R1.fastq.gz -O ${READ_PREFIX}.paired.R2.fastq.gz \
  --unpaired1 ${READ_PREFIX}.unpaired.R1.fastq.gz --unpaired2 ${READ_PREFIX}.unpaired.R2.fastq.gz \
  --qualified_quality_phred 20 \
  --length_required 50 \
  --adapter_fasta ${ADAPTERS} \
  --cut_right --cut_mean_quality 20 --cut_window_size 5 \
  --umi --umi_loc=read1 --umi_len=0
```

### B.3 Dealing with UMIs


[UMI-tools](https://umi-tools.readthedocs.io/en/latest/index.html) can handle any UMI tagged sequencing data where deduplication happens after mapping.

The process is to extract the UMIs from the read sequence and add it to the read names. There are two ways to do this, and between them provide the flexibility to handle any read configuration I can think of (see https://umi-tools.readthedocs.io/en/latest/regex.html)

You then map your reads with your favourite mapper.

The next step depends on whether your technique fragments the cDNA before or after PCR. If fragmentation happens after PCR, then the next step is to assign reads to features (e.g. genes) using featureCounts. If PCR happened after fragmentation, then you do the read assignment/quantification after deduping.

Then you group/dedup/count (depending on your downstream application). If fragmentation happened after PCR then you need to do this on a per-gene basis.

```
# Extract UMI info from reads

module load umitools
# For Single-end reads
# umi_tools extract --stdin=${READ_PREFIX}.paired.R1.fastq.gz --bc-pattern=NNNNNNNNN --log=processed.log --stdout ${READ_PREFIX}.umi.fastq.gz 
# For Paired-end reads
# umi_tools extract [OPTIONS] -p PATTERN [-I IN_FASTQ[.gz]] [-S OUT_FASTQ[.gz]] --read2-in=IN2_FASTQ[.gz] --read2-out=OUT2_FASTQ[.gz]

READ_PREFIX=example
OUTPUT_PREFIX=example
umi_tools extract -I ${READ_PREFIX}.paired.R1.fastq.gz --bc-pattern=NNNXXXXNN --bc-pattern2=NNNXXXXNN \ 
  --read2-in=${READ_PREFIX}.paired.R2.fastq.gz --stdout=${OUTPUT_PREFIX}.umi.paired.R1.fastq.gz \
  --read2-out=${OUTPUT_PREFIX}.umi.paired.R2.fastq.gz
```

### B.4 Mapping reads 
sjdbGTFfeatureExon                      exon
    string: feature type in GTF file to be used as exons for building transcripts

sjdbGTFtagExonParentTranscript          transcript_id
    string: GTF attribute name for parent transcript ID (default "transcript_id" works for GTF files)

sjdbGTFtagExonParentGene                gene_id
    string: GTF attribute name for parent gene ID (default "gene_id" works for GTF files)

sjdbGTFtagExonParentGeneName            gene_name
    string(s): GTF attribute name for parent gene name

sjdbGTFtagExonParentGeneType            gene_type gene_biotype
    string(s): GTF attribute name for parent gene type

sjdbOverhang                            100
    int>0: length of the donor/acceptor sequence on each side of the junctions, ideally = (mate_length - 1)



In the folder "scripts" create a script named "create_star_index.sh" with the following commands. STAR is pretty memmory intensive and therefore, we will runit with 64G of RAM and 16 CPUs.
```
#!/bin/bash
#SBATCH --cpus-per-task=16 --mem=32g

GENOME=$1
GTF=$2
OUTDIR=$3
READLEN=$4

module load STAR

# Create genome index
time STAR --runMode genomeGenerate \
    --runThreadN 16 \
    --genomeDir ${OUTDIR} \
    --sjdbGTFfile ${GTF} \
    --sjdbGTFfeatureExon exon \
    --sjdbGTFtagExonParentGeneType protein_coding \
    --sjdbOverhang ${READLEN} \
    --genomeFastaFiles ${GENOME}
```
And then run the script on the Biowulf grid like this:
```
# Make script executable
chmod a+x ./create_star_index.sh

# Run script with sbatch
sbatch ./create_star_index.sh ../data/GRCh38.chr1.fa ../data/gencode.v45.annotation.gtf ../reference 50
sbatch ./create_star_index.sh
```
### B.5 Deduplicate reads
```
# Sort bam file

# Index bam file

# If you used UMIs
umi_tools dedup -I mapped.bam --paired -S deduplicated.bam

```
