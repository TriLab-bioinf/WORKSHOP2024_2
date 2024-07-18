## Workshop guide

### A- File formats

1. Fasta (examples: file.fa ; file.fasta ; file.fsa)
```
>Seq1
GCTACGCTAGTGTTTTATGCTGATTATAATTATTTTTTAGCTAGCT
GGTCTAGGG
>Seq2
GGGGGTCGATTTTGTCTTT
>Seq3
GCGCGGGGGGGCTAGTGATGTATTGTGTGTGTTTTTGGGACGCGGT
GCTTAGTGTTTCGGGGTGCTTTGGTCGTTTTTTGCAGAAGATTCGT
AGATGT
```
2. Fastq (examples: file.fastq ; files.fq)
```
@A00941:835:H2WLKDRX2:1:2101:1226:1016 1:N:0:GACACCCTGT+TACACGTTGA
CGTTAAAGGTCAGCTTCCCGCAGGCTGGCCTCAGGCGGAGTCTGGGTCAG
+
:FFFF,FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFF
@A00941:835:H2WLKDRX2:1:2101:1389:1016 1:N:0:GACACCATGT+TACACGTTGA
GGGACCATGTCTTCTTAGTCTGCCCTCTCCAGGGTGCTCGGCTACATGCC
+
FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF,FFFF
@A00941:835:H2WLKDRX2:1:2101:1714:1016 1:N:0:GACACCATGT+TACACGTTGA
CTTAAGATGTCTTATTCCTCAAGCCAGAGCCAATACAAGGATAAGAACTG
+
FFFFFFFFFFFFFFFFFFFFFFFFF:FFFFF:FFFFFFFFFFFFFFFFFF
```
Fastq files can be for single-end or paired-end.

For paired-end data, read-1 is usually `forward` and read-2 is `reverse`.

For RNA-seq experiements, fastq files can be strand-specific (stranded) or not. To check if the sequencing data is stranded or not, you can use the following command using the `rseqc` Biowulf module on a bam file:


```
module load rseqc
infer_experiment.py -i <bam file> -r <gene annotations in bed format>

# Example Output
Tris is a PairEnd Data
Fraction of reads failed to determine: 0.057
Fraction of reads explained by "1++, 1--, 2+-, 2-+": 0.0123
Fraction of reads explained by "1+-, 1-+, 2++, 2--": 0.9371
```

## B- Preparing the working directory in Biowulf

**Step 1:**

Login into Biowulf and start a new interactive session having 8 CPUs and 16GB of memory like so:
```
sinteractive --cpus-per-task=8 --mem=16g --time=3:00:00
```

**Step 2:**

Download the `WORKSHOP2024_2` data from the Workshop's GitHub repository into your Biowulf data directory like so:
```
cd /data/$USER/

git clone https://github.com/TriLab-bioinf/WORKSHOP2024_2.git
```
You should see now a new directory named `WORKSHOP2024_2`. Go to that directory `cd WORKSHOP2024_2` and check the content of it by typing `tree`.

**Step 3:**

Create an environmental variable named `WORKSHOPDIR` with the path of your WORKSHOP2024_2 directory. We will use this variable during the processing of sequencing reads:
```
export WORKSHOPDIR=/data/$USER/WORKSHOP2024_2
```


## C- Proprocessing of sequencing reads for RNAseq analysis

### C.1 Workflow:
![](https://github.com/TriLab-bioinf/WORKSHOP2024_2/blob/main/figures/RNAseq.png)

### C.2 Quality control of sequencing data 
Before you start processing your sequencing files it is important to first check the quality of the sequencing data. Among the things we want to check for are:
1. Total number of reads per file 
2. Read base quality along reads
3. Presence of sequencing adaptors
4. GC-content
5. Distribution of bases [ACGT] along reads
6. Distribution of overrepresented sequences

There are several tools that you can use for this. In our case we will you a program called `fastqc`.

**Step 1:**

Create a script named `01-fastqc.sh` with your favorite editor containing the following code:
```
#!/bin/bash
#SBATCH

# Enter path to read files from STDIN 
READ1=$1
READ2=$2

# Set name of output directory
OUTDIR=Step1-fastqc

# Create output directory
mkdir -p $OUTDIR

module load fastqc
fastqc -o $OUTDIR $READ1 $READ2
```
**Note:** Remember to make the script execulatble with `chmod +x 01-fastqc.sh`.


Now run the script localy like this:
```
./01-fastqc.sh ${WORKSHOPDIR}/data/example.R1.fastq.gz ${WORKSHOPDIR}/data/example.R2.fastq.gz
```
You can also have tried to run the same script remotely on a Biowulf node like this:
```
sbatch ./01-fastqc.sh ${WORKSHOPDIR}/data/example.R1.fastq.gz ${WORKSHOPDIR}/data/example.R2.fastq.gz
```

Now you should see the new directory `Step1-fastqc` containing the fastq output. In particular we are interested in the html report files named `example.R1_fastqc.html` and `example.R2_fastqc.html`. You can open them in a browser. To do so, first you need to link your Biowulf directories to the ones in your computer like so:

1. For Macs, go to Finder > Go > Connect to Server
2. Enter `smb://hpcdrive.nih.gov/data`
3. press `Connect`

A new Finder Windows should appear showing your Biowulf data directory. Navigate to the fastqc report files into `WORKSHOP2024_2/Step1-fastqc/` and open one of the html reports. How the reads look like?


Other examples of fastqc reports:

[Good Illumina data](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/good_sequence_short_fastqc.html)

[Bad Illumina data](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/bad_sequence_fastqc.html)

[Others examples from the fastqc website](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

### C.3 Trimming reads

During this step we will try to accomplish the following:

1. To trim reads based on their sequencing quality by chopping off the portions of the read that falls below a particular quality cutoff (usually a value between 20-30).
2. Cut bases at the ends of the read that might be consistently low quality.
3. Remove contaminating sequencing adapters, usually from 3'end of reads.
4. Move UMI sequences, is present, from the read sequence to the fastq header.

**Trimming reads with fastp:**

Create the script `02-trim_reads_fastp.sh` with the following code:
```
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
```
Make the script executable with `chmod +x 02-trim_reads_fastp.sh` and run it locally like this:
```
./02-trim_reads_fastp.sh ./data/example.R1.fastq.gz ./data/example.R2.fastq.gz
```

As before, if you want to use sbatch you can run the script by adding the `sbatch` command at the front:
```
sbatch ./02-trim_reads_fastp.sh ./data/example.R1.fastq.gz ./data/example.R2.fastq.gz
```
In the later case, you can follow the status of the sbatch command using the following:
```
showq -u $USER
```

A nice feature about `fastp` is that it will also generate a QC report about the reads before and after trimming, similar to the one generated by `fastqc`.

You can check the report in html format with your browser, as you did with `fastqc`. The report should be located in your workinf directory and named `fastp.html`.

### C.4 Dealing with UMIs

`fastp` or [`umi_tools`](https://umi-tools.readthedocs.io/en/latest/index.html) can handle any UMI tagged sequencing data where deduplication happens after mapping.

For bulk RNAseq, the process is to extract the UMIs from the read sequence and add it to the read names. Then you map your reads with your favourite mapper, remove duplicated reads from the resulting bam file with, for example, `umi_tools dedup`, and finally you count the number of reads per feature using the deduplicated bam file.


### C.5 Mapping reads 

Next step is to map the trimmed (filtered) reads to the reference genome. In our case we will use the human chromosome 17 as reference. In general, any mapping program require the referece sequence to be indexed before it can start mapping the reads. Genome indexing might take several hours for large genomes as the Human or Mouse genomes. Fortunately, index files are already available for most common model organisms, which can be found in Biowulf under the following path:

```
/fdb/STAR_indices/2.7.11b/GENCODE/Gencode_human/release_39/ref/Genome
```

However, many times there is no available index files for our reference of interest and therefore, we will need to create is from scratch by using as input the fasta file of the referennce genome. In our case, we will create the index files for our reference, the human chr17.

**Step 1:**

Create a script named "03-create_star_index.sh" with the following commands. STAR indexing and mapping is pretty memmory intensive and therefore, we usually run it with 64G of RAM and 16 CPUs for large genomes. In our case, the reference is not that large and therefore, we will use just 8 CPUs and 16Gb of memory.

```
#!/bin/bash
#SBATCH --cpus-per-task=8 --mem=16g

GENOME=$1
GTF=$2
OUTDIR=Step3-reference_index
READLEN=48

module load STAR

# Create genome index
STAR --runMode genomeGenerate \
    --runThreadN 16 \
    --genomeDir ${OUTDIR} \
    --sjdbGTFfile ${GTF} \
    --sjdbGTFfeatureExon exon \
    --sjdbGTFtagExonParentGeneType protein_coding \
    --sjdbOverhang ${READLEN} \
    --genomeFastaFiles ${GENOME}
```
And then run the script, after making it executable, like this:
```
./03-create_star_index.sh ./data/GRCh38.chr17.fa ./data/gencode.v45.annotation.chr17.gtf
```
Once done, take a look at the output file `Step3-reference_index`. You will see a bunch of files that corrspond to the reference indexes used for read mapping.

**Step 2:**

Now we are ready for mapping the trimmed reads to the reference genome with STAR. Create the script `04-mapping_reads_star.sh` with the following code and make it executable, as before:
```
#!/bin/bash
#SBATCH --cpus-per-task=16 --mem=32g

# Load read files from the command line
READ1=$1
READ2=$2

# STAR-specific variables
THREADS=16
REFERENCE=${WORKSHOPDIR}/Step3-reference_index
OUTDIR=${WORKSHOPDIR}/Step4-mapping_star
READ_NAME=$(basename $READ1)
PREFIX=${READ_NAME%.paired.R1.fastq.gz}

module load STAR

STAR --runMode alignReads \
  --runThreadN ${THREADS} \
  --genomeDir ${REFERENCE} \
  --outFilterMismatchNmax 5 \
  --alignEndsType EndToEnd \
  --readFilesIn ${READ1} ${READ2} \
  --readFilesCommand zcat \
  --outFileNamePrefix ${OUTDIR}/${PREFIX}.sorted. \
  --quantMode GeneCounts \
  --outSAMtype BAM SortedByCoordinate \
  --outSAMattrRGline ID:RG1 , SM:SampleName , PL:Illumina , LB:Library.fa \
  --outSAMattributes All &> ${OUTDIR}/star.log
```

Then run the script `04-mapping_reads_star.sh` in a Biowulf node like this:
```
sbatch ./04-mapping_reads_star.sh ./data/example.R1.fastq.gz ./data/example.R2.fastq.gz
```

STAR wil output the bam file `example.bam` containing read mapping information. 


**Note:** Another popular mapper for RNAseq analysis is [HISAT2](https://daehwankimlab.github.io/hisat2/).

### B.5 Deduplicate reads

The next step requires the bam file to be sorted by read coordinates and indexed. In our case, STAR sorted the reads during the mapping step, but if that wasn't the case, then you can sort the bam file by read coordinate using samtools and the following command:
```
# Sort bam file by coordinate

GENOME=${WORSHOPDIR}/data/GRCh38.chr17.fa
PREFIX=example
UNSORTED_BAM=$(WORKSHOPDIR}/Step4-mapping_star/${PREFIX}.bam

module load samtools

samtools sort --threads 8 \
  -O BAM \
  --reference ${GENOME} \
  -T tmp_file \
  -o $(WORKSHOPDIR}/Step4-mapping_star/${PREFIX}.sorted.bam \
  ${UNSORTED_BAM}
```

**Step 3:**

Once the bam file is sorted, it is necessary to create an index of it. FOr this we will use the followig `samtools` command:
```
module load samtools

samtools index $(WORKSHOPDIR}/Step4-mapping_star/example.sorted.bam
```


# If you used UMIs
umi_tools dedup -I mapped.bam --paired -S deduplicated.bam

```
#!/bin/bash
#SBATCH --cpus-per-task=5 --mem=32g --gres=lscratch:40

# Flag duplicated reads with Picard
module load picard/3.2.0

# Enter bam file from the command line
BAM=$1

OUTDIR=${WORKSHOPDIR}/Step5-markduplicates

java -Xmx32g -XX:ParallelGCThreads=5 -jar $PICARDJARPATH/picard.jar MarkDuplicates \
  --INPUT ${BAM} \
  --OUTPUT ${OUTDIR} \
  --METRICS_FILE <File> \
  --CREATE_INDEX true \
  --REMOVE_DUPLICATES false \
  --CREATE_INDEX true \
  --TMP_DIR /lscratch/$SLURM_JOBID

```
Run the 05-mark_duplicates.sh script remotely with sbatch:
```
sbatch lscratch:40
```


samtools index {output}

TMP_DIR=/lscratch/$SLURM_JOBID

```
resources:
        cpus_per_task = 4,
        mem_mb = 128000,
        partition = "quick",
        time = "4:00:00",
        gres = "lscratch:40"
    shell:
        """
        picard -Xmx32g MarkDuplicates \
         I={input} \
         O={output} \
         M={log} \
         {params}
        samtools index {output}
```

### B.6 Count reads per feature
Count reads as fragments for PE reads
```
featureCounts
```
