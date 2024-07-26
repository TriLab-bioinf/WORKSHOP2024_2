## Workshop guide for RNAseq read-processing workflow

## A- File formats

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

For RNA-seq experiements, fastq files can be strand-specific (stranded) or not. To check if the sequencing data is stranded or not, you can map a small fraction of reads to reference genome and then use the following command using the `rseqc` Biowulf module on the resulting bam file:

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

Now we need to copy the data files into the Login into `WORKSHOP2024_2` directory. Data files are stored in the folder `/scratch/lorenziha/data`, that is not accessible from the interactive session. Therefore, open another terminal, and this time login into the `Helix` like so:
```
ssh $USER@helix.nih.gov 
```
Then, run the following command:
```
cp -r /scratch/lorenziha/data /data/$USER/WORKSHOP2024_2/

# Check that the "data" directory has been copied just fine
ls -lrt /data/$USER/WORKSHOP2024_2/data/
```
Now you can close the Helix terminal.

**Step 4:**

In your interactive session, create an environmental variable named `WORKSHOPDIR` with the path of your WORKSHOP2024_2 directory. We will use this variable during the processing of sequencing reads:
```
export WORKSHOPDIR=/data/$USER/WORKSHOP2024_2
```


## C- Proprocessing of sequencing reads for RNAseq analysis

### C.1 Workflows:
There are many different RNAseq processing workflows that you can use to analyze you data. The main factor to consider when choosing the read processing approach is what the goal is of your project or what is the question you are trying to answer. The figure below depicts a number of potential workflows that can be used, but many more combinations are possible!
![](https://github.com/TriLab-bioinf/WORKSHOP2024_2/blob/main/figures/alternative_RNAseq_pipelines.png)
*From: Corchete, L.A., Rojas, E.A., Alonso-López, D. et al. Systematic comparison and assessment of RNA-seq procedures for gene expression quantitative analysis. Sci Rep 10, 19737 (2020). https://doi.org/10.1038/s41598-020-76881-x*

The workflow we will be using during the workshop for processing bulk RNAseq data is summarized below:
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

Next step is to map the trimmed (filtered) reads to the reference genome. In our case, we will use the human chromosome 17 as reference and the read-mapper `STAR`. In general, any read-mapper program requires the referece sequence to be indexed before it can start mapping the reads. Genome indexing might take several hours for large genomes as the Human or Mouse genomes. Fortunately, index files are already available for most common model organisms, which can be found in Biowulf under the following path:

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
  --outSAMattributes All \
  --twopassMode Basic # Optional, it imprives support for novel splice junctions
```

Then run the script `04-mapping_reads_star.sh` in a Biowulf node like this:
```
sbatch ./04-mapping_reads_star.sh ./Step2-trimming_fastp/example.paired.R1.fastq.gz ./Step2-trimming_fastp/example.paired.R2.fastq.gz
```

This command will produce a bam file, which contains read mapping information. In this case, STAR wil output a bam file name `example.sorted.Aligned.sortedByCoord.out.bam` within the directory `Step4-mapping_star`.

SAM/BAM files can be visualized with the command `samtools view <BAM_FILE>`. These files are composed of two parts:
1. A header at the top, composed of several lines starting with the `@` symbol, which contain general information such as whether the file is sorted or not, sizes of the reference chromosomes/scaffolds, commands used to generate the SAM/BAM file, and names of the sample, sequencing library, etc. For instance:
```
@HD	VN:1.4	SO:coordinate
@SQ	SN:chr17	LN:83257441
@PG	ID:STAR	PN:STAR	VN:2.7.11b	CL:STAR   --runMode alignReads      --runThreadN 16   --genomeDir /home/lorenziha/data/WORKSHOP2024_2/Step3-reference_index   --readFilesIn ./data/example.R1.fastq.gz   ./data/example.R2.fastq.gz      --readFilesCommand zcat      --outFileNamePrefix /home/lorenziha/data/WORKSHOP2024_2/Step4-mapping_star/example.R1.fastq.gz.sorted.   --outSAMtype BAM   SortedByCoordinate      --outSAMattributes All      --outSAMattrRGline "ID:RG1	SM:SampleName	PL:Illumina	LB:Library.fa"      --outFilterMismatchNmax 5   --alignEndsType EndToEnd   --quantMode GeneCounts
@RG	ID:RG1	SM:SampleName	PL:Illumina	LB:Library.fa
```

2. After the header comes the mapping information about the reads organized in a tab-delimited manner, where each line stores mapping data for one single read. Mapping data looks like this:
```
A00941:835:H2WLKDRX2:1:2271:9200:35524	99	chr17	59989	255	50M	=	60092	153	TGTGCCGGCCCTGATCATGCAGCTCTTCCAGGCCCACTGCTTCTTCCTGT	,FFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFF	NH:i:1	HI:i:1	AS:i:84	nM:i:1	NM:i:12	MD:Z:0N0N0N0N0N0N0N0N0N0N0N0N38	jM:B:c,-1	jI:B:i,-1	MC:Z:50M	RG:Z:RG1
A00941:835:H2WLKDRX2:1:2115:12048:34601	99	chr17	59992	255	50M	=	60037	95	ACCGGCCCTGATCATGCAGCTCTTCCAGGCCCACTGCTTCTTCCTGTCCA	,FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF,FFF	NH:i:1	HI:i:1	AS:i:87	nM:i:1	NM:i:9	MD:Z:0N0N0N0N0N0N0N0N0N41	jM:B:c,-1	jI:B:i,-1	MC:Z:50M	RG:Z:RG1
```
The first 11 columns of the SAM/BAM file contain the following information:
![](https://github.com/TriLab-bioinf/WORKSHOP2024_2/blob/main/figures/SAM_fileds.png)
*From: SAM (file format). Wikipedia. https://en.wikipedia.org/wiki/SAM_(file_format)*

A more in-depth description of what information is stored in each column can be found [here](https://en.wikipedia.org/wiki/SAM_(file_format)).
 
For RNAseq data, reads spanning a splice junction are represented by Ns in the CIGAR string of the sam/bam file (column 6) with a length equal to the span of the intronic region. For instance:
```
A00941:835:H2WLKDRX2:1:2106:14886:24987	99	chr17	213795	255	49M	=	219698	33326	CGGGCAGGGTCTGGCAGGAATCCTCCACAGGGAAGTCTGTTCCAGGCAC	FFFFFFFFFFFFFFFFFFFFFFFFFFFF,FFFFFFFFFFFFFFFFFFFF	NH:i:1	HI:i:1	AS:i:96	nM:i:0	NM:i:0	MD:Z:49	jM:B:c,-1	jI:B:i,-1	MC:Z:39M27374N10M
A00941:835:H2WLKDRX2:1:2106:14886:24987	147	chr17	219698	255	39M27374N10M	=	213795	-33326	AGCTAAGATCCGAGTCACTGTCACTGTCACTGGAAACCACTCTTCCTCG	F:FFFF::FFFFFFFFF:FF:FFFFFFF:FFFFFFFFFFFFFF,FFFFF	NH:i:1	HI:i:1	AS:i:96	nM:i:0	NM:i:0	MD:Z:49	jM:B:c,22	jI:B:i,219737,247110	MC:Z:49M
```

**Example of spliced read aligment:**

![](https://github.com/TriLab-bioinf/WORKSHOP2024_2/blob/main/figures/Slide1.png)

To find out the significance of specific FLAGS in column 2 of the sam/bam files you can use the following tool: [Decoding SAM flags](https://broadinstitute.github.io/picard/explain-flags.html)

**Note:** Another popular mapper for RNAseq analysis is [HISAT2](https://daehwankimlab.github.io/hisat2/).

### C.6 Remove or flag deduplicated reads

The next step requires the bam file to be sorted by read coordinates. In our case, STAR sorted the reads during the mapping step, but if that wasn't the case, then you can sort the bam file by read coordinate using samtools and the following command:
```
# Sort bam file by coordinate

module load samtools
samtools sort --threads 8 \
  -O BAM \
  --reference ${WORSHOPDIR}/data/GRCh38.chr17.fa \
  -T tmp_file \
  -o $(WORKSHOPDIR}/Step4-mapping_star/${PREFIX}.sorted.bam \
  $(WORKSHOPDIR}/Step4-mapping_star/${PREFIX}.Aligned.sortedByCoord.out.bam
```

**Step 3:**

Now we need to flag the reads that are duplicated. PCR duplicated reads are those whose 5'end coordinate are the same (for paired-end reads, both 5'ends have to match), also share the same UMI (if they have one). The second type of duplication are optical duplicates, that share the same coordinate on the flow cell. 

**Case 1: UMIs were added during library preparation**

If your reads have UMIs, then you can deduplicate the bam file from STAR with `umi_tools dedup` like so:
```

# If you used UMIs
module load umitools
umi_tools dedup -I $(WORKSHOPDIR}/Step4-mapping_star/${PREFIX}.Aligned.sortedByCoord.out.bam \
  --paired -S example.dedup.bam
```

**Case 2: No UMIs were added during library preparation**

In this case, one option is to mark (an optionally remove) duplicated reads with `picard MarkDuplicates`. One of the issues of this tool, is that requires the bam file to contain a SM field in the @RG (Read Group) tag of the bam header, which is not done by STAR. Therefore, we will need to add it with `samtools addreplacerg`. To do this, let's create the script `05-add_read_group.sh` with the following code:
```
#!/bin/bash
#SBATCH

# Get bam file from the command line
INPUT_BAM=$1
DIRNAME=$(dirname $INPUT_BAM)
PREFIX=example.sorted

module load samtools

samtools addreplacerg -w -r "@RG\tID:RG1\tSM:SampleName\tPL:Illumina\tLB:Library.fa" -o ${DIRNAME}/${PREFIX}.bam ${INPUT_BAM}
```
This will create a new bam file named `example.sorted.bam` in the same directory as the input bam file, containing a new @RG line that looks like this:
```
@RG	ID:RG1	SM:SampleName	PL:Illumina	LB:Library.fa
```

Now we are ready for flagging duplicated reads with picard. To do so, prepare a script named `06-mark_duplicates.sh` with the following lines:
```diff
#!/bin/bash
#SBATCH --cpus-per-task=5 --mem=32g --gres=lscratch:40

# Flag duplicated reads with Picard
module load picard/3.2.0

# Enter bam file from the command line
BAM=$1
PREFIX=example
OUTDIR=${WORKSHOPDIR}/Step6-markduplicates

# Make output directory
mkdir -p ${WORKSHOPDIR}/Step6-markduplicates

java -Xmx32g -XX:ParallelGCThreads=5 -jar $PICARDJARPATH/picard.jar MarkDuplicates \
  --INPUT ${BAM} \
  --OUTPUT ${OUTDIR}/${PREFIX}.dedup.bam \
  --METRICS_FILE ${OUTDIR}/${PREFIX}.metrix.txt \
  --CREATE_INDEX true \
  --REMOVE_DUPLICATES false \
  --TMP_DIR /lscratch/$SLURM_JOBID

```
Next, run 06-mark_duplicates.sh remotely with sbatch (remember to make it executable first!!):
```
sbatch ./06-mark_duplicates.sh ./Step4-mapping_star/example.sorted.bam 
```
Picard MarkDuplicates requires a lot of memmory and disk space to run. Therefore, you need to make sure that you request enough memory if you are working with a big genome (human, mouse). This tool though will only require 5 CPUs, as configured above. You can monitor how much resources your program is using while runnign with the command `jobload -j job_#`. For a human genome, I usually request the following to the sbatch job:
```
sbatch --cpus-per-task=5 --mem=128g --gres=lscratch:40 <my picard MarkDuplicates scrip here>
```


### C.7 Count reads per genomic feature

The last step os to generate a tabulated reads containing reads counts per feature from the `dedup.bam` file. . For bulk RNAseq, features are usually either genes or transcripts. For this task we will use a tool from the `subread` module called `featureCounts`. Let's create the script `07-feature_counts.sh` with the following code:

```diff
#!/bin/bash
#SBATCH --cpus-per-task=8

BAM=$1
GENOME=${WORKSHOPDIR}/data/GRCh38.chr17.fa
THREADS=8
GTF=${WORKSHOPDIR}/data/gencode.v45.annotation.chr17.gtf
OUTPUT_FILE=${WORKSHOPDIR}/Step7-read_counts/read_counts.txt

# Make output direcory
mkdir -p ${WORKSHOPDIR}/Step7-read_counts

module load subread

featureCounts -G ${GENOME} -T ${THREADS}\
  -a ${GTF} \
  -t exon \
  -g gene_id \
  -O \
  -s 2 \
  -p --countReadPairs -C \
  --ignoreDup \
  -M --fraction \
  -o ${OUTPUT_FILE} ${BAM}

```
Now run the script using the following command:
```
sbatch ./07-feature_counts.sh ./Step6-markduplicates/example.dedup.bam
```
Within the `Step7-read_counts` output directory, you will find a file named `read_counts.txt`  containing the read counts per feature, that can be imported into R for differential gene expression analysis. Also, you will find a second file named `read_counts.txt.summary` that contains a summary of the reads that were mapped or not to thhe genome features of interest.

You can check the content of the `read_counts.txt` file with the following linux command:
```
cut -f 1,6- Step7-read_counts/read_counts.txt | head
```

### C.8 Generate a report summarizing the results from all steps above

To build a graphical report we will use the utility `multiqc`. Run the following commands on the `WORKSHOP2024_2` directory:
```
module load multiqc
multiqc ${WORKSHOPDIR}
```
This command will generate a html file named `multiqc_report.html` that can be inspected with a browser. Analysis of this report can help to identify potential issues in the sequencing data that are not that obvious from the raw (text) log outputs for the programs above. For instance, here is an example of how the log file output by STAR helped to identify a RNAseq dataset contaminated with rRNA reads, due to problems in the rRNA-depletion step during sequencing library preparation:

![](https://github.com/TriLab-bioinf/WORKSHOP2024_2/blob/main/figures/rRNA_contamination_example.png)

This rRNA contamination was further confirmed by visualization of mapped reads in IGV.
