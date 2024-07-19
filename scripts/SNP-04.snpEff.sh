#!/bin/bash
#SBATCH --cpus-per-task=16
#SBATCH --mem=32g

# -- this file is snpEff.sh --

module load snpEff
#ln -s $SNPEFF_HOME/example/file.vcf .
java -Xmx${SLURM_MEM_PER_NODE}m -jar $SNPEFF_JAR -v hg38 my.var-final.vcf.gz > example.eff.vcf
cat example.eff.vcf | java -jar $SNPSIFT_JAR filter "( EFF[*].IMPACT = 'HIGH' )" > example.filtered.vcf
#java -jar $SNPSIFT_JAR dbnsfp -v -db /fdb/dbNSFP2/dbNSFP3.2a.txt.gz file.eff.vcf > example.annotated.vcf
