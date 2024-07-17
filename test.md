### SNP identification
 
1)	Align reads to reference (using BWA)
Bwa-mem2 is the next version of the bwa-mem algorithm in bwa. It produces alignment identical to bwa and is ~1.3-3.1x faster depending on the use-case, dataset and the running machine.

1.	Index the reference (genome) sequence 
```
bwa-mem2 index [-p prefix] <in.fasta> ##on biowulf, we can find pre-index genome: /fdb/bwa-mem2/hg38/genome.fa
```

2.	Perform the alignment 
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
  	
2) Call SNPs (using bcftools) 
```
bcftools mpileup -f reference.fa alignments.bam | bcftools call -mv -Oz -o calls.vcf.gz
```

● bcftools mpileup 
 Collects summary information in the input BAMs, computes the likelihood of data given each possible genotype and stores the likelihoods in the BCF format. 
● bcftools call 
 Applies the prior and does the actual calling.

3) Filter SNPs 
```
bcftools filter -i'%QUAL>20' calls.vcf.gz -O z -o my.var-final.vcf.gz
```


 
### Bam files conversion to bed, wig, bigwig and tdf
# bed
```
bedtools bamtobed [OPTIONS] -i <bam>
```

# wig: computes average alignment or feature density for over a specified window size across the genome
```
igvtools count NA12878.bam NA12878.wig hg38.fa
```

# convert an wig file to tiled data format (tdf)
```
igvtools toTDF NA12878.wig NA12878.tdf hg38.fa
```

# wig to bigwig
```
module load ucsc 
wigToBigWig in.wig chrom.sizes out.bw
```


### IGV visualization


