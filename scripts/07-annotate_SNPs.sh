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
