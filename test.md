## SNP identification
 
### 1)	Align reads to reference (using BWA)
Bwa-mem2 is the next version of the bwa-mem algorithm in bwa. It produces alignment identical to bwa and is ~1.3-3.1x faster depending on the use-case, dataset and the running machine.

1.1.	Index the reference (genome) sequence 
```
bwa-mem2 index [-p prefix] <in.fasta>
##on biowulf, we can find pre-index genome: /fdb/bwa-mem2/hg38/genome.fa
```

1.2.	Perform the alignment 
```
bwa-mem2 mem -t 32 \
        -R "@RG\tID:$id\tPL:ILLUMINA\tLB:$lb\tSM:$sm" \
        $genome \
        output_forward_paired.fq.gz \
        output_reverse_paired.fq.gz \
        | samtools sort \
        -m 1706M \
        -@ 12 \
        -O BAM \
        -o ERR194160_bwa_sorted.bam \
        2> mapping.log
```
  	
### 2) Call SNPs (using bcftools) 
```
bcftools mpileup -f reference.fa alignments.bam | bcftools call -mv -Oz -o calls.vcf.gz
```

● bcftools mpileup 
 Collects summary information in the input BAMs, computes the likelihood of data given each possible genotype and stores the likelihoods in the BCF format.
 
● bcftools call 
 Applies the prior and does the actual calling.

### 3) Filter SNPs 
```
bcftools filter -i'%QUAL>20' calls.vcf.gz -O z -o my.var-final.vcf.gz
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
 
## Bam files conversion to bed, wig, bigwig and tdf
### bed
```
bedtools bamtobed [OPTIONS] -i <bam>
```

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


## IGV visualization
```
module load IGV IGVTools
igv --memory 20g
```


