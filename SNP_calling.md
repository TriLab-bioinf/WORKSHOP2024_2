## SNP identification
![pipeline](https://github.com/TriLab-bioinf/WORKSHOP2024_2/blob/main/figures/pipeline.jpg)


### get example data

```
ssh $USER@helix.nih.gov
mkdir -p /data/$USER/WORKSHOP2024_2-2
cp -r /scratch/wangy80/WGS_data /data/$USER/WORKSHOP2024_2-2/WGS_data

# Check that the "data" directory has been copied just fine
ls -lrt /data/$USER/WORKSHOP2024_2-2/WGS_data/

```

### 1) Quality Control
#### 1.1 Fastqc

Create a script named "01-fastqc.sh" with the following commands.
```
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
```
Make the script executable with chmod +x 01-fastqc.sh and run it locally like this:
```
export WORKSHOPDIR=/gpfs/gsfs12/users/$USER/WORKSHOP2024_2-2
./01-fastqc.sh ${WORKSHOPDIR}/WGS_data/example_R1.fastq.gz ${WORKSHOPDIR}/WGS_data/example_R2.fastq.gz
```
Then for Macs, you can go to **Finder > Go > Connect** to Server

Enter **smb://hpcdrive.nih.gov/data**

press **Connect** to see the fastqc outputs.

#### 1.2 Trimming reads with trimmomatic

Create the script 02-trim_reads_trimmomatic.sh with the following code:

```
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
```

This will perform the following:

1. Remove adapters (ILLUMINACLIP:TruSeq3-PE.fa:2:30:10)

   ILLUMINACLIP: fastaWithAdaptersEtc: seed mismatches : palindrome clip threshold : simple clip threshold
   
    ● **fastaWithAdaptersEtc**: specifies the path to a fasta file containing all the adapters, PCR sequences etc. The naming of the various sequences within this file determines how they are used.

    ● **seedMismatches**: specifies the maximum mismatch count which will still allow a full match to be performed

    ● **palindromeClipThreshold**: specifies how accurate the match between the two 'adapter ligated' reads must be for PE palindrome read alignment.

    ● **simpleClipThreshold**: specifies how accurate the match between any adapter etc. sequence must be against a read.

2. Remove leading low quality or N bases (below quality 3) (LEADING:3)
3. Remove trailing low quality or N bases (below quality 3) (TRAILING:3)
4. Scan the read with a 4-base wide sliding window, cutting when the average quality per base drops below 15 (SLIDINGWINDOW:4:15)
5. Drop reads below the 36 bases long (MINLEN:36)

Make the script executable with chmod +x 02-trim_reads_trimmomatic.sh and run it locally like this:

```
./02-trim_reads_trimmomatic.sh ./WGS_data/example_R1.fastq.gz ./WGS_data/example_R2.fastq.gz
```
As before, if you want to use sbatch you can run the script by adding the sbatch command at the front:
```
sbatch ./02-trim_reads_trimmomatic.sh ./WGS_data/example_R1.fastq.gz ./WGS_data/example_R2.fastq.gz
```
For Single End:
```
java -jar trimmomatic-0.35.jar SE -phred33 input.fq.gz output.fq.gz ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
```

### 2)	Mapping (using BWA)

BWA is a software package for mapping low-divergent sequences against a large reference genome, such as the human genome. It consists of three algorithms: BWA-backtrack, BWA-SW and BWA-MEM. The first algorithm is designed for Illumina sequence reads up to 100bp, while the rest two for longer sequences ranged from 70bp to 1Mbp. BWA-MEM and BWA-SW share similar features such as long-read support and split alignment, but BWA-MEM, which is the latest, is generally recommended for high-quality queries as it is faster and more accurate. BWA-MEM also has better performance than BWA-backtrack for 70-100bp Illumina reads.

2.1.	Index the reference (genome) sequence

Create a script named "03-create_bwa_index.sh" with the following commands.
```
#!/bin/bash
#SBATCH --cpus-per-task=16 --mem=32g

genome=$1
OUTDIR=03-reference_index
mkdir -p ${OUTDIR}

module load bwa
bwa index $genome -p ${OUTDIR}/hg38_chr17
```

And then run the script, after making it executable, like this:
```
./03-create_bwa_index.sh ./WGS_data/hg38_chr17.fa
```

2.2.	Perform the alignment

Create the script 04-mapping_reads_bwa.sh with the following code:
```
#!/bin/bash
#SBATCH --cpus-per-task=16
#SBATCH --mem=32g

# Load read files from the command line
READ1=$1
READ2=$2

# variables setting
THREADS=16
REFERENCE=${WORKSHOPDIR}/03-reference_index/hg38_chr17
OUTDIR=${WORKSHOPDIR}/04-mapping_bwa
mkdir -p ${OUTDIR}

READ_NAME=$(basename $READ1)
PREFIX=${READ_NAME%_forward_paired.fq.gz}

id=$PREFIX
sm=$PREFIX

module load bwa
module load samtools

bwa mem -t ${THREADS} \
        -R "@RG\tID:$id\tPL:ILLUMINA\tSM:$sm" \
        $REFERENCE \
        $READ1 \
        $READ2 \
        | samtools sort \
        -@ ${THREADS} \
        -O BAM \
        -o ${OUTDIR}/${PREFIX}_bwa_sorted.bam
```
Then run the script 04-mapping_reads_bwa.sh in a Biowulf node like this:
```
sbatch ./04-mapping_reads_bwa.sh ./02-trimming_trimmomatic/example_forward_paired.fq.gz ./02-trimming_trimmomatic/example_reverse_paired.fq.gz
```

### 3) Mark Duplicates with Picard

Create the script 05-mark_duplicates.sh with the following code:
```
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
```
Next, run 05-mark_duplicates.sh remotely with sbatch:
```
sbatch ./05-mark_duplicates.sh ./04-mapping_bwa/example_bwa_sorted.bam
```

### 4) Call SNPs (using bcftools) 

Create the script 06-call_SNPs.sh with the following code:

```
#!/bin/bash
#SBATCH --cpus-per-task=16
#SBATCH --mem=32g

GENOME=./WGS_data/hg38_chr17.fa
bamfile=$1

OUTDIR=${WORKSHOPDIR}/06-SNPcalling
mkdir -p $OUTDIR

READ_NAME=$(basename $bamfile)
PREFIX=${READ_NAME%.dedup.bam}

module load bcftools
bcftools mpileup -f $GENOME $bamfile | bcftools call -mv -Oz -o $OUTDIR/$PREFIX.vcf.gz
```

Next, run 06-call_SNPs.sh remotely with sbatch:
```
sbatch ./06-call_SNPs.sh ./05-markduplicates/example.dedup.bam
```

● bcftools mpileup

 Collects summary information in the input BAMs, computes the likelihood of data given each possible genotype and stores the likelihoods in the BCF format.
 
● bcftools call

 Applies the prior and does the actual calling. The -m switch tells the program to use the default calling method, the -v option asks to output only variant sites, finally the -O option selects the output format.

**Note:** Another popular software for SNP calling is GATK, please refer to https://hpc.nih.gov/training/gatk_tutorial/workflow-overview.html.

## VCF file formats
```
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##bcftoolsVersion=1.19+htslib-1.19
##bcftoolsCommand=mpileup -f hg38_chr17.fa example_bwa_sorted_dedup.bam
##reference=file://hg38_chr17.fa
##contig=<ID=chr17,length=83257441>
##ALT=<ID=*,Description="Represents allele(s) other than observed.">
##INFO=<ID=INDEL,Number=0,Type=Flag,Description="Indicates that the variant is an INDEL.">
##INFO=<ID=IDV,Number=1,Type=Integer,Description="Maximum number of raw reads supporting an indel">
##INFO=<ID=IMF,Number=1,Type=Float,Description="Maximum fraction of raw reads supporting an indel">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Raw read depth">
##INFO=<ID=VDB,Number=1,Type=Float,Description="Variant Distance Bias for filtering splice-site artefacts in RNA-seq data (bigger is better)",Version="3">
##INFO=<ID=RPBZ,Number=1,Type=Float,Description="Mann-Whitney U-z test of Read Position Bias (closer to 0 is better)">
##INFO=<ID=MQBZ,Number=1,Type=Float,Description="Mann-Whitney U-z test of Mapping Quality Bias (closer to 0 is better)">
##INFO=<ID=BQBZ,Number=1,Type=Float,Description="Mann-Whitney U-z test of Base Quality Bias (closer to 0 is better)">
##INFO=<ID=MQSBZ,Number=1,Type=Float,Description="Mann-Whitney U-z test of Mapping Quality vs Strand Bias (closer to 0 is better)">
##INFO=<ID=SCBZ,Number=1,Type=Float,Description="Mann-Whitney U-z test of Soft-Clip Length Bias (closer to 0 is better)">
##INFO=<ID=SGB,Number=1,Type=Float,Description="Segregation based metric, http://samtools.github.io/bcftools/rd-SegBias.pdf">
##INFO=<ID=MQ0F,Number=1,Type=Float,Description="Fraction of MQ0 reads (smaller is better)">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="List of Phred-scaled genotype likelihoods">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes for each ALT allele, in the same order as listed">
##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
##INFO=<ID=DP4,Number=4,Type=Integer,Description="Number of high-quality ref-forward , ref-reverse, alt-forward and alt-reverse bases">
##INFO=<ID=MQ,Number=1,Type=Integer,Description="Average mapping quality">
##bcftools_callVersion=1.19+htslib-1.19
##bcftools_callCommand=call -mv -Oz -o calls.vcf.gz; Date=Fri Jul 19 09:50:42 2024
##bcftools_filterVersion=1.19+htslib-1.19
##bcftools_filterCommand=filter -e QUAL<20 -O z -o my.var-final.vcf.gz calls.vcf.gz; Date=Fri Jul 19 09:54:01 2024
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  example
chr17   62608   .       G       A       21.9555 PASS    DP=3;VDB=0.28;SGB=-0.453602;RPBZ=-1.22474;MQBZ=0;BQBZ=0;SCBZ=0;MQ0F=0;AC=1;AN=2;DP4=0,1,0,2;MQ=60       GT:PL   0/1:55,0,26
chr17   76752   .       C       A       22.4195 PASS    DP=4;VDB=0.52;SGB=-0.453602;RPBZ=0;MQBZ=0;MQSBZ=0;BQBZ=-0.408248;SCBZ=0;MQ0F=0;AC=1;AN=2;DP4=0,2,2,0;MQ=60      GT:PL   0/1:55,0,61
chr17   114276  .       A       G       61.4147 PASS    DP=3;VDB=0.0618664;SGB=-0.511536;MQ0F=0;AC=2;AN=2;DP4=0,0,3,0;MQ=60     GT:PL   1/1:91,9,0
chr17   114551  .       G       C       40.4148 PASS    DP=2;VDB=0.5;SGB=-0.453602;MQ0F=0;AC=2;AN=2;DP4=0,0,0,2;MQ=60   GT:PL   1/1:70,6,0
chr17   298712  .       aa      aACa    189.416 PASS    INDEL;IDV=4;IMF=0.8;DP=5;VDB=0.709369;SGB=-0.590765;RPBZ=1.41421;MQBZ=0;MQSBZ=0;SCBZ=0;MQ0F=0;AC=2;AN=2;DP4=0,0,2,3;MQ=60       GT:PL   1/1:219,15,0
chr17   355539  .       gt      gTt     37.7211 PASS    INDEL;IDV=2;IMF=0.666667;DP=3;VDB=0.28;SGB=-0.453602;RPBZ=0;MQBZ=0;MQSBZ=0;BQBZ=-1.22474;SCBZ=0;MQ0F=0;AC=1;AN=2;DP4=0,1,1,1;MQ=60      GT:PL   0/1:71,0,26
chr17   671239  .       agaaag  aGAAAGgaaag     116.327 PASS    INDEL;IDV=3;IMF=0.75;DP=4;VDB=0.470313;SGB=-0.511536;RPBZ=-1.34164;MQBZ=0;MQSBZ=0;SCBZ=0;MQ0F=0;AC=2;AN=2;DP4=0,1,0,3;MQ=60     GT:PL   1/1:143,4,0
```

In a nutshell, VCF format is tab-separated text file having the following columns:

Chromosome name

Position

Variant's ID

Reference genome

Alternative (i.e. variant)

Quality score

Filter (whether or not the variant passed quality filters)

INFO : Generic information about this variant. SnpEff adds annotation information in this column.

Information about the following columns - The GT in the FORMAT column tells us to expect genotypes in the following columns.

Individual identifier: The previous column told us to expect to see genotypes here. The genotype is in the form 0|1, where 0 indicates the reference allele and 1 indicates the alternative allele, i.e it is heterozygous.

### 5) Filter SNPs 
Variant filtering is not easy. The variant callers provide a quality score (the QUAL) column, which gives an estimate of how likely it is to observe a call purely by chance. An easy way to filter low quality calls is
```
module load bcftools
bcftools filter -e 'QUAL<20' 06-SNPcalling/example.vcf.gz -O z -o 06-SNPcalling/my.var-final.vcf.gz
```

Other useful metrics are:

● sequencing depth (DP bigger than twice the average depth indicates problematic regions and is often enriched for artefacts)

● the minimum number of high-quality non-reference reads

● proximity to indels (bcftools filter -g)

● etc.

### 6) SNP annotations
SnpEff: Genetic variant annotation, and functional effect prediction toolbox. It annotates and predicts the effects of genetic variants on genes and proteins (such as amino acid changes).

Create the script 07-annotate_SNPs.sh with the following code:

```
#!/bin/bash
# -- this file is snpEff.sh --
#SBATCH --cpus-per-task=16
#SBATCH --mem=32g

module load snpEff
VCFfile=$1

OUTDIR=${WORKSHOPDIR}/07-SNPannotation
mkdir -p $OUTDIR

java -Xmx${SLURM_MEM_PER_NODE}m -jar $SNPEFF_JAR -v hg38 $VCFfile > $OUTDIR/example.eff.vcf
cat $OUTDIR/example.eff.vcf | java -jar $SNPSIFT_JAR filter "( EFF[*].IMPACT = 'HIGH' )" > $OUTDIR/example.filtered.vcf
#java -jar $SNPSIFT_JAR dbnsfp -v -db /fdb/dbNSFP2/dbNSFP3.2a.txt.gz file.eff.vcf > example.annotated.vcf
```

Next, run 07-annotate_SNPs.sh remotely with sbatch:
```
sbatch ./07-annotate_SNPs.sh ./06-SNPcalling/my.var-final.vcf.gz
```


After SnpEff annotation, the information is added into INFO column
```
##SnpEffVersion="5.2a (build 2023-10-24 14:24), by Pablo Cingolani"
##SnpEffCmd="SnpEff  hg38 my.var-final.vcf.gz "
##INFO=<ID=ANN,Number=.,Type=String,Description="Functional annotations: 'Allele | Annotation | Annotation_Impact | Gene_Name | Gene_ID | Feature_Type | Feature_ID | Transcript_BioType | Rank | HGVS.c | HGVS.p | cDNA.pos / cDNA.length | CDS.pos / CDS.length | AA.pos / AA.length | Distance | ERRORS / WARNINGS / INFO' ">
##INFO=<ID=LOF,Number=.,Type=String,Description="Predicted loss of function effects for this variant. Format: 'Gene_Name | Gene_ID | Number_of_transcripts_in_gene | Percent_of_transcripts_affected'">
##INFO=<ID=NMD,Number=.,Type=String,Description="Predicted nonsense mediated decay effects for this variant. Format: 'Gene_Name | Gene_ID | Number_of_transcripts_in_gene | Percent_of_transcripts_affected'">
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  example
chr17   62608   .       G       A       21.9555 PASS    DP=3;VDB=0.28;SGB=-0.453602;RPBZ=-1.22474;MQBZ=0;BQBZ=0;SCBZ=0;MQ0F=0;AC=1;AN=2;DP4=0,1,0,2;MQ=60;ANN=A|intergenic_region|MODIFIER|CHR_START-SCGB1C2|CHR_START-SCGB1C2|intergenic_region|CHR_START-SCGB1C2|||n.62608G>A||||||   GT:PL   0/1:55,0,26
chr17   76752   .       C       A       22.4195 PASS    DP=4;VDB=0.52;SGB=-0.453602;RPBZ=0;MQBZ=0;MQSBZ=0;BQBZ=-0.408248;SCBZ=0;MQ0F=0;AC=1;AN=2;DP4=0,2,2,0;MQ=60;ANN=A|intergenic_region|MODIFIER|CHR_START-SCGB1C2|CHR_START-SCGB1C2|intergenic_region|CHR_START-SCGB1C2|||n.76752C>A||||||  GT:PL   0/1:55,0,61
chr17   148758  .       C       T       160.416 PASS    DP=7;VDB=0.933109;SGB=-0.636426;MQSBZ=0;MQ0F=0;AC=2;AN=2;DP4=0,0,2,5;MQ=60;ANN=T|intron_variant|MODIFIER|DOC2B|DOC2B|transcript|NM_003585.5|protein_coding|7/8|c.1005+353G>A||||||      GT:PL   1/1:190,21,0
chr17   150509  .       T       TA      43.4147 PASS    INDEL;IDV=2;IMF=1;DP=2;VDB=0.28;SGB=-0.453602;MQSBZ=0;BQBZ=1;MQ0F=0;AC=2;AN=2;DP4=0,0,2,0;MQ=60;ANN=TA|intron_variant|MODIFIER|DOC2B|DOC2B|transcript|NM_003585.5|protein_coding|6/8|c.924-1318_924-1317insT||||||      GT:PL   1/1:73,6,0
chr17   156324  .       G       C       42.4147 PASS    DP=2;VDB=0.06;SGB=-0.453602;MQSBZ=0;MQ0F=0;AC=2;AN=2;DP4=0,0,1,1;MQ=60;ANN=C|synonymous_variant|LOW|DOC2B|DOC2B|transcript|NM_003585.5|protein_coding|6/9|c.819C>G|p.Ser273Ser|990/6062|819/1239|273/412||      GT:PL   1/1:72,6,0
chr17   17042759        .       CG      CGG     26.4242 PASS    INDEL;IDV=2;IMF=1;DP=2;VDB=0.42;SGB=-0.453602;BQBZ=0;MQ0F=0;AC=2;AN=2;DP4=0,0,0,2;MQ=60;ANN=CGG|frameshift_variant|HIGH|MPRIP|MPRIP|transcript|NM_001364716.2|protein_coding|1/24|c.8dupG|p.Gly4fs|305/15419|9/7374|3/2457||,CGG|5_prime_UTR_variant|MODIFIER|MPRIP|MPRIP|transcript|NM_015134.4|protein_coding|1/23|c.-89dupG|||||88|,CGG|5_prime_UTR_variant|MODIFIER|MPRIP|MPRIP|transcript|NM_201274.4|protein_coding|1/24|c.-89dupG|||||88|        GT:PL   1/1:56,6,0
chr17   75589275        .       G       A       25.4267 PASS    DP=2;VDB=0.5;SGB=-0.453602;MQ0F=0;AC=2;AN=2;DP4=0,0,0,2;MQ=60;ANN=A|stop_gained|HIGH|MYO15B|MYO15B|transcript|NM_001309242.1|protein_coding|1/63|c.1218G>A|p.Trp406*|1218/9540|1218/9177|406/3058||;LOF=(MYO15B|MYO15B|1|1.00);NMD=(MYO15B|MYO15B|1|1.00)     GT:PL   1/1:55,6,0
chr17   21300898        .       C       T       74.4149 PASS    DP=3;VDB=0.0785113;SGB=-0.511536;MQSBZ=0;MQ0F=0;AC=2;AN=2;DP4=0,0,1,2;MQ=60;ANN=T|stop_gained|HIGH|MAP2K3|MAP2K3|transcript|NM_145109.3|protein_coding|5/12|c.304C>T|p.Gln102*|514/2239|304/1044|102/347||,T|stop_gained|HIGH|MAP2K3|MAP2K3|transcript|NM_001316332.2|protein_coding|6/13|c.217C>T|p.Gln73*|639/2364|217/957|73/318||,T|stop_gained|HIGH|MAP2K3|MAP2K3|transcript|NM_002756.4|protein_coding|5/12|c.217C>T|p.Gln73*|311/2060|217/957|73/318||;LOF=(MAP2K3|MAP2K3|3|1.00);NMD=(MAP2K3|MAP2K3|3|1.00) GT:PL   1/1:104,9,0
chr17   75589275        .       G       A       25.4267 PASS    DP=2;VDB=0.5;SGB=-0.453602;MQ0F=0;AC=2;AN=2;DP4=0,0,0,2;MQ=60;ANN=A|stop_gained|HIGH|MYO15B|MYO15B|transcript|NM_001309242.1|protein_coding|1/63|c.1218G>A|p.Trp406*|1218/9540|1218/9177|406/3058||;LOF=(MYO15B|MYO15B|1|1.00);NMD=(MYO15B|MYO15B|1|1.00)     GT:PL   1/1:55,6,0
```

Here is a description of the meaning of each sub-field:
1. Allele (or ALT): In case of multiple ALT fields, this helps to identify which ALT we are referring to. E.g.:
```
# CHROM  POS     ID  REF  ALT    QUAL  FILTER  INFO
chr1    123456  .   C    A      .     .       ANN=A|...
chr1    234567  .   A    G,T    .     .       ANN=G|... , T|...
```
2. Annotation (a.k.a. effect): Annotated using Sequence Ontology terms. Multiple effects can be concatenated using '&'.
```
#CHROM  POS     ID  REF  ALT  QUAL  FILTER  INFO
chr1    123456  .   C    A    .     .      ANN=A|intron_variant&nc_transcript_variant|...
```
3. Putative_impact: A simple estimation of putative impact / deleteriousness : {HIGH, MODERATE, LOW, MODIFIER}
   
| **Impact**  | **Meaning** | **Example** |
| ------------- | ------------- | ------------- |
| HIGH  | The variant is assumed to have high (disruptive) impact in the protein, probably causing protein truncation, loss of function or triggering nonsense mediated decay. | stop_gained, frameshift_variant |
| MODERATE  | A non-disruptive variant that might change protein effectiveness. | missense_variant, inframe_deletion |
| LOW | Assumed to be mostly harmless or unlikely to change protein behavior. | 	synonymous_variant |
| MODIFIER | Usually non-coding variants or variants affecting non-coding genes, where predictions are difficult or there is no evidence of impact. | exon_variant, downstream_gene_variant |

4. Gene Name: Common gene name (HGNC). Optional: use closest gene when the variant is "intergenic".
5. Gene ID: Gene ID
6. Feature type: Which type of feature is in the next field (e.g. transcript, motif, miRNA, etc.). It is preferred to use Sequence Ontology (SO) terms, but 'custom' (user defined) are allowed.
7. Feature ID: Depending on the annotation, this may be: Transcript ID (preferably using version number), Motif ID, miRNA, ChipSeq peak, Histone mark, etc. Note: Some features may not have ID (e.g. histone marks from custom Chip-Seq experiments may not have a unique ID).
8. Transcript biotype: The bare minimum is at least a description on whether the transcript is {"Coding", "Noncoding"}. Whenever possible, use ENSEMBL biotypes.
9. Rank / total: Exon or Intron rank / total number of exons or introns.
10. HGVS.c: Variant using HGVS notation (DNA level)
11. HGVS.p: If variant is coding, this field describes the variant using HGVS notation (Protein level). Since transcript ID is already mentioned in 'feature ID', it may be omitted here.
12. cDNA_position / cDNA_len: Position in cDNA and trancript's cDNA length (one based).
13. CDS_position / CDS_len: Position and number of coding bases (one based includes START and STOP codons).
14. Protein_position / Protein_len: Position and number of AA (one based, including START, but not STOP).
15. Distance to feature: All items in this field are options, so the field could be empty.
    Up/Downstream: Distance to first / last codon
    Intergenic: Distance to closest gene
    Distance to closest Intron boundary in exon (+/- up/downstream). If same, use positive number.
    Distance to closest exon boundary in Intron (+/- up/downstream)
    Distance to first base in MOTIF
    Distance to first base in miRNA
    Distance to exon-intron boundary in splice_site or splice _region
    ChipSeq peak: Distance to summit (or peak center)
    Histone mark / Histone state: Distance to summit (or peak center)
16. Errors, Warnings or Information messages: Add errors, warnings or informative message that can affect annotation accuracy.

 
## Bam files conversion to wig, bigwig and tdf
### wig: computes average alignment or feature density for over a specified window size across the genome
```
module load IGVTools
igvtools --memory 20g count 05-markduplicates/example.dedup.bam example.dedup.wig WGS_data/hg38_chr17.fa
```

wig example
```
track type=wiggle_0
variableStep chrom=chr17 span=25
60051   1.0
60076   1.76
60101   2.0
60126   2.0
60151   1.12
60176   1.28
```

### convert an wig file to tiled data format (tdf)
```
igvtools toTDF example.dedup.wig example.dedup.tdf WGS_data/hg38_chr17.fa
```

### wig to bigwig
create chrom.sizes
```
chr1    248956422
chr2    242193529
chr3    198295559
chr4    190214555
chr5    181538259
chr6    170805979
chr7    159345973
chr8    145138636
chr9    138394717
chr10   133797422
chr11   135086622
chr12   133275309
chr13   114364328
chr14   107043718
chr15   101991189
chr16   90338345
chr17   83257441
chr18   80373285
chr19   58617616
chr20   64444167
chr21   46709983
chr22   50818468
chrX    156040895
chrY    57227415
chrM    16569
```
We can use samtools faidx reference genome and get fai file, then cut the first two columns into chrom.sizes.

```
module load samtools
module load ucsc

## create fasta index file
samtools faidx WGS_data/hg38_chr17.fa
cut -f1,2 WGS_data/hg38_chr17.fa.fai >chrom.sizes

wigToBigWig example.dedup.wig chrom.sizes example.dedup.bw
```



## gff/gtf files conversion to bed
```
module load bedops
convert2bed --input=gtf --output=bed < data/gencode.v45.annotation.chr17.gtf > gencode.v45.annotation.chr17.bed
```

## IGV visualization
```
module load IGV
igv --memory 20g
```


