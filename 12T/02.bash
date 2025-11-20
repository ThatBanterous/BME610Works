#!/bin/bash
#SBATCH --chdir=.
#SBATCH --cpus-per-task=10
#SBATCH --error='%x-%A.txt'
#SBATCH --output='%x-%A.txt'
#SBATCH --job-name=BME610-02
#SBATCH --mem=20G
#SBATCH --export=ALL
# DO NOT MODIFY the above lines
/BiO/Share/Tools/samtools-1.21/samtools sort -l 9 --threads 10 --reference /BiO/Share/Tools/gatk-bundle/hg38/Homo_sapiens_assembly38.fasta --write-index -o 12T.Sort.bam 12T.bam
/BiO/Share/Tools/gatk-4.6.1.0/gatk MarkDuplicatesSpark --input 12T.Sort.bam --output 12T.Sort.MarkDuplicates.bam --reference /BiO/Share/Tools/gatk-bundle/hg38/Homo_sapiens_assembly38.fasta --metrics-file 12T.Sort.MarkDuplicates.metrics --duplicate-tagging-policy 'OpticalOnly' -- --spark-master 'local[10]' --spark-verbosity 'INFO'
/BiO/Share/Tools/gatk-4.6.1.0/gatk BaseRecalibrator --input 12T.Sort.MarkDuplicates.bam --reference /BiO/Share/Tools/gatk-bundle/hg38/Homo_sapiens_assembly38.fasta --output 12T.Sort.MarkDuplicates.BQSR.table --create-output-bam-index true --known-sites /BiO/Share/Tools/gatk-bundle/hg38/1000G_phase1.snps.high_confidence.hg38.vcf --known-sites /BiO/Share/Tools/gatk-bundle/hg38/1000G_omni2.5.hg38.vcf.gz --known-sites /BiO/Share/Tools/gatk-bundle/hg38/Homo_sapiens_assembly38.dbsnp138.vcf --known-sites /BiO/Share/Tools/gatk-bundle/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf
/BiO/Share/Tools/gatk-4.6.1.0/gatk ApplyBQSR --bqsr-recal-file 12T.Sort.MarkDuplicates.BQSR.table --input 12T.Sort.MarkDuplicates.bam --output 12T.Sort.MarkDuplicates.BQSR.bam --reference /BiO/Share/Tools/gatk-bundle/hg38/Homo_sapiens_assembly38.fasta --create-output-bam-index true

/BiO/Share/Tools/samtools-1.21/samtools sort -l 9 --threads 10 --reference /BiO/Share/Tools/gatk-bundle/hg38/Homo_sapiens_assembly38.fasta --write-index -o TUMOR.Sort.bam TUMOR.bam
/BiO/Share/Tools/gatk-4.6.1.0/gatk MarkDuplicatesSpark --input TUMOR.Sort.bam --output TUMOR.Sort.MarkDuplicates.bam --reference /BiO/Share/Tools/gatk-bundle/hg38/Homo_sapiens_assembly38.fasta --metrics-file TUMOR.Sort.MarkDuplicates.metrics --duplicate-tagging-policy 'OpticalOnly' -- --spark-master 'local[10]' --spark-verbosity 'INFO'
/BiO/Share/Tools/gatk-4.6.1.0/gatk BaseRecalibrator --input TUMOR.Sort.MarkDuplicates.bam --reference /BiO/Share/Tools/gatk-bundle/hg38/Homo_sapiens_assembly38.fasta --output TUMOR.Sort.MarkDuplicates.BQSR.table --create-output-bam-index true --known-sites /BiO/Share/Tools/gatk-bundle/hg38/1000G_phase1.snps.high_confidence.hg38.vcf --known-sites /BiO/Share/Tools/gatk-bundle/hg38/1000G_omni2.5.hg38.vcf.gz --known-sites /BiO/Share/Tools/gatk-bundle/hg38/Homo_sapiens_assembly38.dbsnp138.vcf --known-sites /BiO/Share/Tools/gatk-bundle/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf
/BiO/Share/Tools/gatk-4.6.1.0/gatk ApplyBQSR --bqsr-recal-file TUMOR.Sort.MarkDuplicates.BQSR.table --input TUMOR.Sort.MarkDuplicates.bam --output TUMOR.Sort.MarkDuplicates.BQSR.bam --reference /BiO/Share/Tools/gatk-bundle/hg38/Homo_sapiens_assembly38.fasta --create-output-bam-index true
