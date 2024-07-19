## SNP identification
![pipeline](https://github.com/TriLab-bioinf/WORKSHOP2024_2/blob/main/figures/pipeline.jpg)


### get example data

```
ssh $USER@helix.nih.gov
cp -r /scratch/wangy80/WGS_data ./
cd WGS_data
```

### 1) Quality Control
1.1 Fastqc
```
module load fastqc
fastqc example_R1.fastq.gz
fastqc example_R2.fastq.gz
```

1.2 trimmomatic
```
#!/bin/bash
## This is for Paired End data
sample=example
adapter=/usr/local/apps/trimmomatic/0.39/adapters/TruSeq3-PE.fa
module load trimmomatic || exit 1
java -Djava.io.tmpdir=. -jar $TRIMMOJAR PE -phred33 -threads 12 \
    ${sample}_R1.fastq.gz ${sample}_R2.fastq.gz \
    ${sample}_forward_paired.fq.gz ${sample}_forward_unpaired.fq.gz \
    ${sample}_reverse_paired.fq.gz ${sample}_reverse_unpaired.fq.gz \
    ILLUMINACLIP:$adapter:2:30:10 LEADING:3 TRAILING:3 \
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
#!/bin/bash
#SBATCH --cpus-per-task=16
#SBATCH --mem=32g

module load bwa
module load samtools

genome=hg38_chr17.fa
id=example
sm=example
bwa mem -t 32 \
        -R "@RG\tID:$id\tPL:ILLUMINA\tSM:$sm" \
        $genome \
        ${sm}_forward_paired.fq.gz \
        ${sm}_reverse_paired.fq.gz \
        | samtools sort \
        -@ 16 \
        -O BAM \
        -o ${sm}_bwa_sorted.bam
```

### 3) Mark Duplicates
```
#!/bin/bash
module load picard
sample=example
time java -Xmx8g -XX:ParallelGCThreads=5 -jar $PICARDJARPATH/picard.jar MarkDuplicates \
CREATE_INDEX=true \
METRICS_FILE=$sample.metrix.txt \
INPUT=${sample}_bwa_sorted.bam \
REMOVE_DUPLICATES=false \
OUTPUT=${sample}_bwa_sorted_dedup.bam
```

### 4) Call SNPs (using bcftools) 
```
module load bcftools
bcftools mpileup -f hg38_chr17.fa example_bwa_sorted_dedup.bam | bcftools call -mv -Oz -o calls.vcf.gz
```

● bcftools mpileup

 Collects summary information in the input BAMs, computes the likelihood of data given each possible genotype and stores the likelihoods in the BCF format.
 
● bcftools call

 Applies the prior and does the actual calling. The -m switch tells the program to use the default calling method, the -v option asks to output only variant sites, finally the -O option selects the output format.

### 5) Filter SNPs 
Variant filtering is not easy. The variant callers provide a quality score (the QUAL) column, which gives an estimate of how likely it is to observe a call purely by chance. An easy way to filter low quality calls is
```
bcftools filter -e 'QUAL<20' calls.vcf.gz -O z -o my.var-final.vcf.gz
```
Other useful metrics are:

● sequencing depth (DP bigger than twice the average depth indicates problematic regions and is often enriched for artefacts)

● the minimum number of high-quality non-reference reads

● proximity to indels (bcftools filter -g)

● etc.

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


### 6) SNP annotations
SnpEff: Genetic variant annotation, and functional effect prediction toolbox. It annotates and predicts the effects of genetic variants on genes and proteins (such as amino acid changes).

```
#!/bin/bash
# -- this file is snpEff.sh --
#SBATCH --cpus-per-task=12
#SBATCH --mem=24g

module load snpEff
#ln -s $SNPEFF_HOME/example/file.vcf .
java -Xmx${SLURM_MEM_PER_NODE}m -jar $SNPEFF_JAR -v hg38 my.var-final.vcf.gz > example.eff.vcf
cat example.eff.vcf | java -jar $SNPSIFT_JAR filter "( EFF[*].IMPACT = 'HIGH' )" > example.filtered.vcf
#java -jar $SNPSIFT_JAR dbnsfp -v -db /fdb/dbNSFP2/dbNSFP3.2a.txt.gz file.eff.vcf > example.annotated.vcf
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
chr17   156366  .       A       G       81.415  PASS    DP=3;VDB=0.301053;SGB=-0.511536;MQSBZ=0;MQ0F=0;AC=2;AN=2;DP4=0,0,2,1;MQ=60;ANN=G|synonymous_variant|LOW|DOC2B|DOC2B|transcript|NM_003585.5|protein_coding|6/9|c.777T>C|p.Thr259Thr|948/6062|777/1239|259/412||  GT:PL   1/1:111,9,0
chr17   17042759        .       CG      CGG     26.4242 PASS    INDEL;IDV=2;IMF=1;DP=2;VDB=0.42;SGB=-0.453602;BQBZ=0;MQ0F=0;AC=2;AN=2;DP4=0,0,0,2;MQ=60;ANN=CGG|frameshift_variant|HIGH|MPRIP|MPRIP|transcript|NM_001364716.2|protein_coding|1/24|c.8dupG|p.Gly4fs|305/15419|9/7374|3/2457||,CGG|5_prime_UTR_variant|MODIFIER|MPRIP|MPRIP|transcript|NM_015134.4|protein_coding|1/23|c.-89dupG|||||88|,CGG|5_prime_UTR_variant|MODIFIER|MPRIP|MPRIP|transcript|NM_201274.4|protein_coding|1/24|c.-89dupG|||||88|        GT:PL   1/1:56,6,0
```

 
## Bam files conversion to wig, bigwig and tdf
### wig: computes average alignment or feature density for over a specified window size across the genome
```
module load IGVTools
igvtools count example_bwa_sorted_dedup.bam example_bwa_sorted_dedup.wig hg38_chr17.fa
```

### convert an wig file to tiled data format (tdf)
```
igvtools toTDF example_bwa_sorted_dedup.wig example_bwa_sorted_dedup.tdf hg38_chr17.fa
```

### wig to bigwig
```
module load ucsc 
wigToBigWig example_bwa_sorted_dedup.wig chrom.sizes example_bwa_sorted_dedup.bw
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


