## SNP identification
![pipeline](https://github.com/TriLab-bioinf/WORKSHOP2024_2/blob/main/figures/pipeline.jpg)


### 1) Quality Control
1.1 Fastqc

1.2 trimmomatic
```
#!/bin/bash
## This is for Paired End data
module load trimmomatic || exit 1
java -Djava.io.tmpdir=. -jar $TRIMMOJAR PE -phred33 -threads $SLURM_CPUS_PER_TASK \
    ERR194160_1.fastq.gz ERR194160_2.fastq.gz \
    output_forward_paired.fq.gz output_forward_unpaired.fq.gz \
    output_reverse_paired.fq.gz output_reverse_unpaired.fq.gz \
    ILLUMINACLIP:/usr/local/apps/trimmomatic/0.39/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 \
    SLIDINGWINDOW:4:15 MINLEN:36
```

This will perform the following:

Remove adapters (ILLUMINACLIP:TruSeq3-PE.fa:2:30:10)

Remove leading low quality or N bases (below quality 3) (LEADING:3)

Remove trailing low quality or N bases (below quality 3) (TRAILING:3)

Scan the read with a 4-base wide sliding window, cutting when the average quality per base drops below 15 (SLIDINGWINDOW:4:15)

Drop reads below the 36 bases long (MINLEN:36)

```
# Single End:
java -jar trimmomatic-0.35.jar SE -phred33 input.fq.gz output.fq.gz ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
```

### 2)	Mapping (using BWA)

BWA is a software package for mapping low-divergent sequences against a large reference genome, such as the human genome. It consists of three algorithms: BWA-backtrack, BWA-SW and BWA-MEM. The first algorithm is designed for Illumina sequence reads up to 100bp, while the rest two for longer sequences ranged from 70bp to 1Mbp. BWA-MEM and BWA-SW share similar features such as long-read support and split alignment, but BWA-MEM, which is the latest, is generally recommended for high-quality queries as it is faster and more accurate. BWA-MEM also has better performance than BWA-backtrack for 70-100bp Illumina reads.

2.1.	Index the reference (genome) sequence 
```
module load bwa
bwa index hg38_chr17.fa
```

2.2.	Perform the alignment 
```
genome=hg38_chr17.fa
bwa mem -t 32 \
        -R "@RG\tID:$id\tPL:ILLUMINA\tLB:$lb\tSM:$sm" \
        $genome \
        output_forward_paired.fq.gz \
        output_reverse_paired.fq.gz \
        | samtools sort \
        -@ 12 \
        -O BAM \
        -o ERR194160_bwa_sorted.bam \
        2> mapping.log
```

### 3) Mark Duplicates
```
module load picard
java -Xmx4g -XX:ParallelGCThreads=5 -jar $PICARDJARPATH/picard.jar MarkDuplicates \
CREATE_INDEX=true \
INPUT=<input.bam> \
VALIDATION_STRINGENCY=STRICT
```

### 4) Call SNPs (using bcftools) 
```
bcftools mpileup -f reference.fa alignments.bam | bcftools call -mv -Oz -o calls.vcf.gz
```

● bcftools mpileup

 Collects summary information in the input BAMs, computes the likelihood of data given each possible genotype and stores the likelihoods in the BCF format.
 
● bcftools call

 Applies the prior and does the actual calling. The -m switch tells the program to use the default calling method, the -v option asks to output only variant sites, finally the -O option selects the output format.

### 5) Filter SNPs 
Variant filtering is not easy. The variant callers provide a quality score (the QUAL) column, which gives an estimate of how likely it is to observe a call purely by chance. An easy way to filter low quality calls is
```
bcftools filter -i'%QUAL>20' calls.vcf.gz -O z -o my.var-final.vcf.gz
```
Other useful metrics are:

● sequencing depth (DP bigger than twice the average depth indicates problematic regions and is often enriched for artefacts)

● the minimum number of high-quality non-reference reads

● proximity to indels (bcftools filter -g)

● etc.

To give a concrete example, the following filter seemed to work quite well for one particular dataset (human data, exomes):

```
bcftools filter -sLowQual -g3 -G10 \
    -e'%QUAL<10 || (RPB<0.1 && %QUAL<15) || (AC<2 && %QUAL<15) || %MAX(DV)<=3 || %MAX(DV)/%MAX(DP)<=0.3' \
    calls.vcf.gz
```

### 6) SNP annotations
SnpEff: Genetic variant annotation, and functional effect prediction toolbox. It annotates and predicts the effects of genetic variants on genes and proteins (such as amino acid changes).

```
#!/bin/bash
# -- this file is snpEff.sh --

module load snpEff
ln -s $SNPEFF_HOME/example/file.vcf .
java -Xmx${SLURM_MEM_PER_NODE}m -jar $SNPEFF_JAR -v hg19 file.vcf > file.eff.vcf
cat file.eff.vcf | java -jar $SNPSIFT_JAR filter "( EFF[*].IMPACT = 'HIGH' )" > file.filtered.vcf
java -jar $SNPSIFT_JAR dbnsfp -v -db /fdb/dbNSFP2/dbNSFP3.2a.txt.gz file.eff.vcf > file.annotated.vcf
```

## VCF file formats
```
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##bcftoolsVersion=1.9+htslib-1.9
##bcftoolsCommand=mpileup -f hg38.fa NA12878.bam
##reference=file://hg38.fa
##contig=<ID=chr1,length=248956422>
...
...
##ALT=<ID=*,Description="Represents allele(s) other than observed.">
##INFO=<ID=INDEL,Number=0,Type=Flag,Description="Indicates that the variant is an INDEL.">
##INFO=<ID=IDV,Number=1,Type=Integer,Description="Maximum number of reads supporting an indel">
##INFO=<ID=IMF,Number=1,Type=Float,Description="Maximum fraction of reads supporting an indel">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Raw read depth">
##INFO=<ID=VDB,Number=1,Type=Float,Description="Variant Distance Bias for filtering splice-site artefacts in RNA-seq data (bigger is better)",Version="3">
##INFO=<ID=RPB,Number=1,Type=Float,Description="Mann-Whitney U test of Read Position Bias (bigger is better)">
##INFO=<ID=MQB,Number=1,Type=Float,Description="Mann-Whitney U test of Mapping Quality Bias (bigger is better)">
##INFO=<ID=BQB,Number=1,Type=Float,Description="Mann-Whitney U test of Base Quality Bias (bigger is better)">
##INFO=<ID=MQSB,Number=1,Type=Float,Description="Mann-Whitney U test of Mapping Quality vs Strand Bias (bigger is better)">
##INFO=<ID=SGB,Number=1,Type=Float,Description="Segregation based metric.">
##INFO=<ID=MQ0F,Number=1,Type=Float,Description="Fraction of MQ0 reads (smaller is better)">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="List of Phred-scaled genotype likelihoods">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##INFO=<ID=ICB,Number=1,Type=Float,Description="Inbreeding Coefficient Binomial test (bigger is better)">
##INFO=<ID=HOB,Number=1,Type=Float,Description="Bias in the number of HOMs number (smaller is better)">
##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes for each ALT allele, in the same order as listed">
##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
##INFO=<ID=DP4,Number=4,Type=Integer,Description="Number of high-quality ref-forward , ref-reverse, alt-forward and alt-reverse bases">
##INFO=<ID=MQ,Number=1,Type=Integer,Description="Average mapping quality">
##bcftools_callVersion=1.9+htslib-1.9
##bcftools_callCommand=call -mv -Oz -o calls.vcf.gz; Date=Wed Jul 17 15:46:50 2024
##bcftools_filterVersion=1.9+htslib-1.9
##bcftools_filterCommand=filter -i%QUAL>20 -O z -o my.var-final.vcf.gz calls.vcf.gz; Date=Wed Jul 17 15:50:54 2024
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NA12878
chr1	5501186	.	C	G	30.4183	PASS	DP=2;SGB=-0.379885;MQ0F=0;AC=2;AN=2;DP4=0,0,1,0;MQ=60	GT:PL	1/1:60,3,0
chr1	6706325	.	Catatatatatatata	CATATatatatatatatata	30.4183	PASS	
chr1	8123500	.	G	A	30.4183	PASS	DP=2;SGB=-0.379885;MQ0F=0;AC=2;AN=2;DP4=0,0,1,0;MQ=60	GT:PL	1/1:60,3,0
chr1	10564183	.	ctt	cTTtt	30.4183	PASS	INDEL;IDV=1;IMF=1;DP=1;SGB=-0.379885;MQ0F=0;AC=2;AN=2;DP4=0,0,0,1;MQ=60	GT:PL	1/1:60,3,0
chr1	17091463	.	C	T	30.4183	PASS	DP=2;SGB=-0.379885;MQ0F=0;AC=2;AN=2;DP4=0,0,1,0;MQ=60	GT:PL	1/1:60,3,0
chr1	17304501	.	C	G	21.4353	PASS	DP=4;SGB=-0.379885;RPB=1;MQB=1;BQB=1;MQ0F=0;ICB=1;HOB=0.5;AC=1;AN=2;DP4=1,0,1,0;MQ=60	GT:PL	0/1:54,0,54
chr1	17304527	.	C	A	80	PASS	DP=4;VDB=0.7;SGB=-0.453602;MQ0F=0;AC=2;AN=2;DP4=0,0,2,0;MQ=60	GT:PL	1/1:110,6,0
chr1	24866346	.	G	T	30.4183	PASS	DP=2;SGB=-0.379885;MQ0F=0;AC=2;AN=2;DP4=0,0,0,1;MQ=60	GT:PL	1/1:60,3,0
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

 
## Bam files conversion to wig, bigwig and tdf
### wig: computes average alignment or feature density for over a specified window size across the genome
```
module load IGVTools
igvtools count NA12878.bam NA12878.wig hg38.fa
```

### convert an wig file to tiled data format (tdf)
```
igvtools toTDF NA12878.wig NA12878.tdf hg38.fa
```

### wig to bigwig
```
module load ucsc 
wigToBigWig in.wig chrom.sizes out.bw
```

## gff/gtf files conversion to bed
```
module load bedops
convert2bed --input=fmt [--output=fmt] [options] < input > output
```

## IGV visualization
```
module load IGV
igv --memory 20g
```


