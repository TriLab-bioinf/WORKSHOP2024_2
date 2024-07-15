### WORKSHOP2024_2

#### Workshop outline:


1. File formats (fastq, sam/bam/, gff/gtf/, VCF)
2. Read quality check (fastqc: example of rRNA contamination from yeast RNAseq from Orna -maybe should be in multiQC), example of low quality reads, example of reads with adapter contamination)

3. Read trimming and filtering (trimmomatic? Try others as well?)
4. Mapping reads to a reference genome (STAR, bwa/bowtie2?) Should we mention mapping of long reads?
5. Remove duplicated reads (picard, have you tried samtools for this?)
6. Quantification of reads per gene. (featureCounts)
  . Counting reads or fragments
  . Dealing with multimappers   
8. SNP/INDEL calls and annotation.
9. MultiQC to summarize program outputs.
